#!/usr/bin/env bash
# Push current HEAD to GitHub via Windows git (github.com git-smart-HTTP is
# blocked from WSL; api.github.com works, but big blobs get corrupted by the
# middlebox). Proven route: bundle -> Windows clone -> git push with token.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
TOKEN="$(gh auth token 2>/dev/null || true)"
if [ -z "$TOKEN" ]; then echo "push_via_windows: no gh token"; exit 1; fi

BUNDLE_WIN='C:\Users\Administrator\AppData\Local\Temp\meta-acdc-push.bundle'
DEST_WIN='C:\Users\Administrator\AppData\Local\Temp\meta-acdc-push-repo'
BUNDLE_WSL='/mnt/c/Users/Administrator/AppData/Local/Temp/meta-acdc-push.bundle'

echo "== bundle $(date +%F\ %T)"
git bundle create "$BUNDLE_WSL" main || exit 1

echo "== clone + push via Windows git"
cmd.exe /c "cd /d C:\Users\Administrator\AppData\Local\Temp && rmdir /s /q meta-acdc-push-repo >nul 2>&1 & D:\Git\cmd\git.exe clone -q C:\Users\Administrator\AppData\Local\Temp\meta-acdc-push.bundle meta-acdc-push-repo && cd /d meta-acdc-push-repo && D:\Git\cmd\git.exe checkout main >nul 2>&1 && D:\Git\cmd\git.exe push https://x-access-token:${TOKEN}@github.com/QianRongAn/meta-acdc-compute.git main:main" 2>&1 | grep -av "UNC\|不支持\|当前目录\|卷 的" | tail -6
RC=${PIPESTATUS[0]}

echo "== cleanup"
cmd.exe /c "cd /d C:\Users\Administrator\AppData\Local\Temp && rmdir /s /q meta-acdc-push-repo >nul 2>&1 & del /q meta-acdc-push.bundle >nul 2>&1 & echo ok" >/dev/null 2>&1
exit $RC
