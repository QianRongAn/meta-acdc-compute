#!/usr/bin/env python3
"""da234 watchdog: resume the 234-positive 3-seed DA retrain on the AutoDL
instance when it powers back on, pull the checkpoints, and (optionally) run
the local stability gate on them.

The 2026-10-15 03:57 instance run was interrupted by a shutdown before any
checkpoint was saved (train_domain_adapt.py saves only after each seed's 150
epochs), so the whole 3-seed run must be redone. This watchdog does:

  WAIT    poll connect.nmb2.seetacloud.com:PORT until reachable
  ENSURE  if /root/autodl-tmp/meta-acdc-compute or its base dataset is
          missing (empty data disk after a power cycle), re-upload
          /tmp/opencode/repo.bundle + /tmp/opencode/data_pack.tar.gz and
          unpack (git clone + tar xzf)
  TRAIN   launch the 3-seed training detached if it is not running and the
          checkpoints da_seed{0,1,2}.model.pt are missing (max 3 relaunches)
  PULL    once all 3 checkpoints exist, SFTP them + train.log back into
          data/processed/da_v2_crystals_instance/ locally
  DONE    stop (local gate runs via scripts/run_da234_gate.sh <dir>)

Long-lived background process, same pattern as scripts/fetch_colabfold_run2.py:
  nohup .venv/bin/python scripts/da234_watchdog.py &
Status: /tmp/opencode/da234.status   Log: logs/da234_watchdog.log
"""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
from pathlib import Path

HOST = "connect.nmb2.seetacloud.com"
PORT = 39776
USER = "root"
PW = os.environ.get("AUTODL_PW", "j6Wtd36duL2q")
REMOTE_ROOT = "/root/autodl-tmp"
REMOTE_REPO = f"{REMOTE_ROOT}/meta-acdc-compute"
REMOTE_OUT = f"{REMOTE_REPO}/data/processed/da_v2_crystals"
CKPTS = ["da_seed0.model.pt", "da_seed1.model.pt", "da_seed2.model.pt"]

REPO = Path("/home/administrator/docs/meta-acdc")
BUNDLE = Path("/tmp/opencode/repo.bundle")
DATA_PACK = Path("/tmp/opencode/data_pack.tar.gz")
LOCAL_OUT = REPO / "data/processed/da_v2_crystals_instance"
STATUS = Path("/tmp/opencode/da234.status")
LOG = REPO / "logs/da234_watchdog.log"

INTERVAL = 90                 # seconds between checks
CONNECT_TIMEOUT = 15
MAX_RELAUNCH = 3

import paramiko  # venv only


def log(msg: str) -> None:
    line = f"[{time.strftime('%F %T')}] {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def set_status(state: str) -> None:
    STATUS.write_text(f"{time.strftime('%F %T')} {state}\n",
                      encoding="utf-8")


def port_open() -> bool:
    try:
        with socket.create_connection((HOST, PORT), timeout=CONNECT_TIMEOUT):
            return True
    except OSError:
        return False


def connect():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PW,
                timeout=CONNECT_TIMEOUT, look_for_keys=False,
                allow_agent=False)
    return ssh


def remote_sh(ssh, cmd: str, timeout: int = 120) -> str:
    _, out, err = ssh.exec_command(cmd, timeout=timeout)
    o = out.read().decode("utf-8", "replace")
    e = err.read().decode("utf-8", "replace")
    if e.strip():
        o += "\n[stderr] " + e[-1500:]
    return o


def ensure_repo(ssh) -> bool:
    """True if repo+data present after this call."""
    check = (f"test -d {REMOTE_REPO}/data/processed "
             f"&& test -f {REMOTE_REPO}/data/processed/egnn_dataset.pt "
             f"&& ls {REMOTE_REPO}/data/processed/da_positives_v2/*_model_0.cif "
             f"| wc -l")
    out = remote_sh(ssh, check).strip()
    out = out.splitlines()[-1].strip() if out else ""
    if out.isdigit() and int(out) >= 200:
        return True
    log("instance repo/data missing (empty disk?) -> re-upload bundles")
    if not (BUNDLE.exists() and DATA_PACK.exists()):
        log(f"local bundles missing: {BUNDLE} {DATA_PACK} -> need human")
        set_status("needs-human: local bundles missing")
        sys.exit(1)
    sftp = ssh.open_sftp()
    for local, remote in [(BUNDLE, f"{REMOTE_ROOT}/repo.bundle"),
                          (DATA_PACK, f"{REMOTE_ROOT}/data_pack.tar.gz")]:
        log(f"uploading {local.name} ({local.stat().st_size/1048576:.0f} MB)")
        t0 = time.time()
        sftp.put(str(local), remote)
        log(f"  done in {(time.time()-t0)/60:.1f} min")
    sftp.close()
    out = remote_sh(
        ssh,
        f"set -e; cd {REMOTE_ROOT}; rm -rf meta-acdc-compute; "
        f"git clone {REMOTE_ROOT}/repo.bundle meta-acdc-compute; "
        f"cd meta-acdc-compute; git checkout main; "
        f"tar xzf {REMOTE_ROOT}/data_pack.tar.gz; "
        f"ls data/processed/da_positives_v2/*_model_0.cif | wc -l",
        timeout=600)
    log(f"unpack result: {out[-400:]}")
    last = [ln.strip() for ln in out.strip().splitlines() if ln.strip()]
    return bool(last) and last[-1].isdigit() and int(last[-1]) >= 200


def train_state(ssh) -> tuple[bool, str]:
    """(done, tail) — done if all 3 checkpoints exist."""
    out = remote_sh(ssh, f"ls {REMOTE_OUT}/ 2>/dev/null; "
                         f"pgrep -f train_domain_adapt || true")
    done = all(c in out for c in CKPTS)
    return done, out


def launch_train(ssh) -> None:
    cmd = (
        f"cd {REMOTE_REPO} && mkdir -p data/processed/da_v2_crystals && "
        f"nohup setsid env PYTHONPATH=src /root/miniconda3/bin/python "
        f"src/meta_acdc/models/train_domain_adapt.py "
        f"--af3 data/processed/da_positives_v2 --epochs 150 --seeds 0 1 2 "
        f"--out-dir data/processed/da_v2_crystals --device cuda "
        f"> data/processed/da_v2_crystals/train.log 2>&1 < /dev/null &"
    )
    out = remote_sh(ssh, cmd)
    log(f"launched: {out[-200:]}")


def pull_checkpoints(ssh) -> None:
    LOCAL_OUT.mkdir(parents=True, exist_ok=True)
    sftp = ssh.open_sftp()
    for name in CKPTS + ["train.log"]:
        remote = f"{REMOTE_OUT}/{name}"
        local = LOCAL_OUT / name
        sftp.get(remote, str(local))
        log(f"pulled {name} ({local.stat().st_size/1024:.0f} KB)")
    sftp.close()


def main() -> int:
    log(f"watchdog start: host={HOST}:{PORT} interval={INTERVAL}s "
        f"max_relaunch={MAX_RELAUNCH}")
    relaunches = 0
    was_up = False
    while True:
        if not port_open():
            if was_up:
                log("instance went offline")
                was_up = False
            set_status("waiting-instance")
            time.sleep(INTERVAL)
            continue
        was_up = True
        try:
            ssh = connect()
        except Exception as e:
            log(f"connect failed: {e}")
            time.sleep(INTERVAL)
            continue
        try:
            if not ensure_repo(ssh):
                log("repo ensure failed, retry next cycle")
                time.sleep(INTERVAL)
                continue
            done, tail = train_state(ssh)
            if done:
                log("all 3 checkpoints present -> pull")
                pull_checkpoints(ssh)
                set_status("done: checkpoints pulled")
                log(f"FINISHED. local gate: scripts/run_da234_gate.sh "
                    f"{LOCAL_OUT.relative_to(REPO)}")
                return 0
            running = "train_domain_adapt" in tail
            if running:
                set_status("training-running")
                log(f"training running: {tail[-200:]}")
            else:
                relaunches += 1
                if relaunches > MAX_RELAUNCH:
                    set_status("needs-human: training keeps dying")
                    log("too many relaunches -> needs human")
                    return 1
                log(f"training not running (relaunch {relaunches}/"
                    f"{MAX_RELAUNCH}); tail: {tail[-200:]}")
                launch_train(ssh)
        except Exception as e:
            log(f"cycle error: {e}")
        finally:
            ssh.close()
        time.sleep(INTERVAL)


if __name__ == "__main__":
    sys.exit(main())
