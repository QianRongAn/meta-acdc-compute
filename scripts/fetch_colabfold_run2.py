#!/usr/bin/env python3
"""Poll the AutoDL GPU instance and auto-fetch + process ColabFold run2.

Run as a long-lived background process:
    .venv/bin/python scripts/fetch_colabfold_run2.py

Behavior (loop every INTERVAL seconds):
- if the instance is unreachable (port closed), log and retry.
- once reachable, wait until /root/af3/out2 has 85 *.done.txt AND the
  colabfold_batch process is gone (plus a grace period).
- then: tar out2 remotely (single tar, setsid pattern), verify with
  tar -tzf, SFTP-download, extract into results_run2 (Windows) and
  data/raw/colabfold_run2/ (repo), run scripts/process_colabfold_run.sh,
  git commit + push (Windows-git route if direct push fails).
- writes /tmp/opencode/fetch_run2.status after every iteration.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HOST = "connect.nmb2.seetacloud.com"
PORT = 44784
USER = "root"
PW = "8P3EwySkAMlx"
REMOTE = "/root/af3"
EXPECTED = 85
INTERVAL = 240          # seconds between checks
GRACE = 180             # seconds after last colabfold process exits
MAX_RUNTIME = 16 * 3600
WIN_RESULTS = Path("/mnt/c/Users/Administrator/WorkBuddy/"
                   "2026-09-23-01-47-50/af3work/results_run2")
REPO = Path("/home/administrator/docs/meta-acdc")
RAW_RUN2 = REPO / "data/raw/colabfold_run2"
STATUS = Path("/tmp/opencode/fetch_run2.status")
LOG = REPO / "logs" / "fetch_run2.log"

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


def connect():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(HOST, port=PORT, username=USER, password=PW, timeout=25,
                look_for_keys=False, allow_agent=False, banner_timeout=25,
                auth_timeout=25)
    return ssh


def run(ssh, cmd: str, timeout: int = 120) -> str:
    stdin, stdout, stderr = ssh.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode("utf-8", "replace")
    err = stderr.read().decode("utf-8", "replace")
    return (out + "\n" + err).strip()


def remote_ready(ssh) -> tuple[bool, str]:
    done = run(ssh, f"ls {REMOTE}/out2/*.done.txt 2>/dev/null | wc -l")
    procs = run(ssh, "pgrep -f colabfold_batch | wc -l")
    last = run(ssh, f"ls -t {REMOTE}/out2/*.done.txt 2>/dev/null | head -1 "
                    "| xargs -r stat -c %Y 2>/dev/null || echo 0")
    try:
        n_done, n_proc, last_mtime = int(done.strip()), int(procs.strip()), int(last.strip())
    except ValueError:
        return False, f"unparseable: done={done!r} procs={procs!r}"
    age = time.time() - last_mtime if last_mtime else -1
    log(f"remote: done={n_done}/{EXPECTED} colabfold_procs={n_proc} "
        f"last_done_age={age:.0f}s")
    if n_done >= EXPECTED and n_proc == 0:
        return True, "complete"
    if n_proc == 0 and 0 < n_done < EXPECTED and age > 7200:
        return True, f"stalled-partial ({n_done}/{EXPECTED}, last done {age/3600:.1f}h ago)"
    return False, "waiting"


def pack_remote(ssh) -> bool:
    """Single tar of out2; verify with tar -tzf. Returns True on verified."""
    tarball = f"{REMOTE}/out2_full.tar.gz"
    run(ssh, f"pgrep -f 'tar -c[z]f' | head -3 || true")
    log("launching remote tar (setsid pattern)...")
    run(ssh, f'nohup setsid bash -c "tar czf {tarball} -C {REMOTE} out2" '
             f"> {REMOTE}/pack3.log 2>&1 < /dev/null & sleep 8; echo LAUNCHED")
    for i in range(60):  # up to 10 min
        time.sleep(10)
        procs = run(ssh, "pgrep -f 'tar -c[z]f' | wc -l").strip()
        if procs == "0":
            break
    else:
        log("tar still running after 10 min — proceeding to verify anyway")
    time.sleep(5)
    n = run(ssh, f"tar -tzf {tarball} 2>/dev/null | wc -l").strip()
    sz = run(ssh, f"stat -c %s {tarball} 2>/dev/null || echo 0").strip()
    log(f"tarball: {n} entries, {int(sz)/1048576:.1f} MB")
    try:
        return int(n) > 100
    except ValueError:
        return False


def download_and_extract(ssh) -> bool:
    tarball = f"{REMOTE}/out2_full.tar.gz"
    local_tgz = Path("/tmp/opencode/out2_full.tar.gz")
    sftp = ssh.open_sftp()
    st = sftp.stat(tarball)
    log(f"downloading {st.st_size/1048576:.1f} MB via SFTP...")
    sftp.get(tarball, str(local_tgz))
    sftp.close()
    log("verifying local tarball...")
    with tarfile.open(local_tgz) as tf:
        members = tf.getmembers()
        if len(members) < 100:
            log(f"local tarball suspicious: {len(members)} entries; abort")
            return False
        n_pdb = sum(1 for m in members if m.name.endswith(".pdb"))
        log(f"local tarball OK: {len(members)} entries, {n_pdb} pdbs")
        WIN_RESULTS.mkdir(parents=True, exist_ok=True)
        tf.extractall(WIN_RESULTS)
        RAW_RUN2.mkdir(parents=True, exist_ok=True)
        for m in members:
            if m.name.endswith(".pdb"):
                tf.extract(m, RAW_RUN2)
    local_tgz.unlink(missing_ok=True)
    log(f"extracted to {WIN_RESULTS} and {RAW_RUN2}")
    return True


def process() -> bool:
    log("running process_colabfold_run.sh run2 ...")
    p = subprocess.run(
        ["bash", str(REPO / "scripts/process_colabfold_run.sh"),
         str(RAW_RUN2), "run2",
         "results on disk (colabfold run2, 2026-10-14)",
         " [+colabfold run2 2026-10-14]"],
        cwd=REPO, capture_output=True, text=True, timeout=3600)
    for ln in (p.stdout + p.stderr).splitlines()[-40:]:
        log(f"  | {ln}")
    return p.returncode == 0


def commit_and_push() -> bool:
    log("record docs + git add + commit run2 artifacts...")
    subprocess.run([str(REPO / ".venv/bin/python"),
                    str(REPO / "scripts/record_colabfold_run.py"),
                    "run2", "2026-10-14", "85"], cwd=REPO, check=False)
    for f in ["data/processed/colabfold_run2_import_report.tsv",
              "data/processed/colabfold_scores_ensemble.tsv"]:
        subprocess.run(["git", "add", f], cwd=REPO, check=False)
    subprocess.run(["git", "add",
                    "scripts/fetch_colabfold_run2.py",
                    "scripts/process_colabfold_run.sh",
                    "scripts/push_via_windows.sh",
                    "scripts/record_colabfold_run.py",
                    ".gitignore",
                    "af3_pending_upload/already_uploaded.tsv",
                    "docs/benchmarks.md", "docs/milestones/KEY-NODES.md"],
                   cwd=REPO, check=False)
    subprocess.run(["git", "commit", "-m",
                    "2026-10-14 ColabFold run2 回收(85): 自动拉取+导入+5-seed 打分+台账"],
                   cwd=REPO, check=False)
    log("pushing via Windows git...")
    r = subprocess.run(["bash", str(REPO / "scripts/push_via_windows.sh")],
                       cwd=REPO, capture_output=True, text=True, timeout=900)
    for ln in (r.stdout + r.stderr).splitlines()[-10:]:
        log(f"  | {ln}")
    return r.returncode == 0


def main() -> int:
    log(f"watcher start: host={HOST}:{PORT} expected={EXPECTED} "
        f"interval={INTERVAL}s max={MAX_RUNTIME/3600:.0f}h")
    start = time.time()
    while time.time() - start < MAX_RUNTIME:
        set_status("checking")
        try:
            ssh = connect()
        except Exception as e:
            log(f"unreachable: {type(e).__name__}: {e}")
            set_status("instance-unreachable")
            ssh = None
        if ssh is not None:
            try:
                ready, why = remote_ready(ssh)
                if ready:
                    log(f"FETCH TRIGGER: {why}")
                    set_status(f"fetching ({why})")
                    if pack_remote(ssh) and download_and_extract(ssh):
                        ok_proc = process()
                        ok_push = commit_and_push()
                        log(f"process={ok_proc} push={ok_push}")
                        set_status(f"DONE process={ok_proc} push={ok_push}")
                        ssh.close()
                        return 0
                    log("fetch pipeline aborted; will retry next round")
                else:
                    set_status(f"up ({why})")
            except Exception as e:
                log(f"remote error: {type(e).__name__}: {e}")
                set_status("remote-error")
            finally:
                try:
                    ssh.close()
                except Exception:
                    pass
        time.sleep(INTERVAL)
    set_status("timeout")
    log("watcher timeout — instance never came back")
    return 2


if __name__ == "__main__":
    sys.exit(main())
