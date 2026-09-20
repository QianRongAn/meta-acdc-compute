#!/usr/bin/env bash
# Meta-ACDC 计算环境安装脚本(KN-1)
# 前提: 系统 Python 3.14(无需降级)、NVIDIA 驱动 >= 525(CUDA 12.6 兼容)。
# 用法: bash scripts/setup_env.sh
set -euo pipefail

cd "$(dirname "$0")/.."
PYPI_MIRROR="${PYPI_MIRROR:-https://pypi.tuna.tsinghua.edu.cn/simple}"
PT_MIRROR="${PT_MIRROR:-https://mirror.sjtu.edu.cn/pytorch-wheels/cu126}"

echo "==> 1/4 创建 venv(无 ensurepip 系统用 --without-pip + get-pip.py 引导)"
if [ ! -x .venv/bin/python ]; then
  python3 -m venv --without-pip .venv
fi
if [ ! -x .venv/bin/pip ]; then
  curl -sL https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py
  .venv/bin/python /tmp/get-pip.py -q
fi
.venv/bin/pip config set global.index-url "$PYPI_MIRROR" >/dev/null 2>&1 || true

echo "==> 2/4 下载 torch/torchvision 轮子(本地安装,断点续传)"
mkdir -p .cache
TORCH_WHL="torch-2.9.1%2Bcu126-cp314-cp314-manylinux_2_28_x86_64.whl"
TV_WHL="torchvision-0.24.1%2Bcu126-cp314-cp314-manylinux_2_28_x86_64.whl"
[ -f .cache/torch.whl ] || curl -sSL --retry 5 -C - -o .cache/torch.whl "$PT_MIRROR/$TORCH_WHL"
[ -f .cache/torchvision.whl ] || curl -sSL --retry 5 -C - -o .cache/torchvision.whl "$PT_MIRROR/$TV_WHL"

echo "==> 3/4 安装 torch + torchvision(本地轮子)"
.venv/bin/pip install --no-index .cache/torch.whl .cache/torchvision.whl

echo "==> 4/4 安装科学计算栈"
.venv/bin/pip install \
  numpy pandas scipy scikit-learn matplotlib tqdm \
  torch_geometric pyg_lib 2>/dev/null || \
.venv/bin/pip install numpy pandas scipy scikit-learn matplotlib tqdm torch_geometric

echo "==> 验证"
.venv/bin/python -c "import torch; print('torch', torch.__version__, '| cuda available:', torch.cuda.is_available(), '| device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')"
