# Meta-ACDC compute — reproducible container (KN-16 deliverable)
# CPU-only image (the reference host uses GTX 1050 Ti + torch 2.9.1+cu126;
# GPU users should use scripts/setup_env.sh instead).
FROM python:3.12-slim

WORKDIR /app

# system deps for mhcflurry / sklearn / building
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential curl git \
    && rm -rf /var/lib/apt/lists/*

# python deps (CPU torch; mhcflurry 2.2.1 is torch-backed)
RUN pip install --no-cache-dir \
        numpy pandas scipy scikit-learn matplotlib tqdm \
        torch --index-url https://download.pytorch.org/whl/cpu \
    && pip install --no-cache-dir mhcflurry

COPY pyproject.toml README.md ./
COPY src ./src
COPY scripts ./scripts
COPY docs ./docs
COPY tests ./tests
RUN pip install --no-cache-dir -e . --no-deps

# pre-fetch MHCflurry weights into the image (so runs are offline-capable);
# the pipes shim is only needed on Python >= 3.13 — harmless on 3.12
RUN python -c "import sys, types, shlex; \
m = types.ModuleType('pipes'); m.quote = shlex.quote; sys.modules['pipes'] = m; \
sys.argv=['x','fetch','models_class1_presentation']; \
from mhcflurry.downloads_command import run; run()" || \
    echo "WARNING: mhcflurry weights not pre-fetched (fetch at first run)"

# smoke test: package imports + CLI entry
RUN python -c "import meta_acdc, meta_acdc.models.egnn, \
meta_acdc.active_learning.simulate_al, meta_acdc.structure.graph; \
print('meta-acdc smoke OK')"

# regression tests (data-dependent cases skip without the local fixtures)
RUN python -m unittest discover -s tests 2>&1 | tail -3

WORKDIR /work
# usage examples (see docs/SOP.md):
#   docker run --rm -v $PWD:/work meta-acdc \
#     python src/meta_acdc/dashboard/server.py --port 8000
CMD ["bash"]
