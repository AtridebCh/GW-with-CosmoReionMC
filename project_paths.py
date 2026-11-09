"""
project_paths.py
================
Single source of truth for every filesystem path used in this repository.
Lives at the repository root, next to `.env`.

Resolution order for each configurable path:
    1. An environment variable already set in the shell / job script
    2. A value in the `.env` file at the repo root (machine-specific, NOT committed)
    3. A default inside the repository

Usage from any subfolder (mcmc/, GravWave/, ...):
    from project_paths import CAMB_DIR, PLANCK_DATA_DIR, CHAINS_DIR

For the import to work, the repo root must be on PYTHONPATH, e.g. in a job script:
    export PYTHONPATH=/path/to/repo:$PYTHONPATH

Check what got resolved on the current machine:
    python project_paths.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Folder containing this file = repository root.
PROJECT_ROOT = Path(__file__).resolve().parent


# --- .env loading --------------------------------------------------------
# Variables already set in the environment take precedence, so a job
# script can always override what's in .env. No external package needed.
def _load_env_file(env_path: Path) -> None:
    if not env_path.is_file():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key.startswith("export "):
            key = key[len("export "):].strip()
        value = value.strip().strip('"').strip("'")
        os.environ.setdefault(key, value)


_load_env_file(PROJECT_ROOT / ".env")


def _path_from_env(var: str, default: Path) -> Path:
    """Return the path in env var `var`, or `default` if it isn't set.

    Supports ~ and $VARS inside the value, e.g. PLANCK_DATA_DIR=$SCRATCH/planck
    """
    value = os.environ.get(var)
    if value:
        return Path(os.path.expandvars(value)).expanduser().resolve()
    return default.resolve()


# --- Code folders (inside the repository) --------------------------------
# These are part of the repo, so they normally need no .env entry.
# They can still be overridden, e.g. to use a CAMB build elsewhere.

CAMB_DIR = _path_from_env("CAMB_DIR", PROJECT_ROOT / "CAMB-1.1.0_modified")
PLANCK_CODE_DIR = _path_from_env("PLANCK_CODE_DIR", PROJECT_ROOT / "COM_Likelihood_Code-v3.0_R3.10")
GLOBAL21CM_DIR = PROJECT_ROOT / "global_21cm_signal"
GRAVWAVE_DIR = PROJECT_ROOT / "GravWave"
MCMC_DIR = PROJECT_ROOT / "mcmc"
REIONIZATION_DIR = PROJECT_ROOT / "ReionizationWithF2Py"
REIONIZATION_PY_DIR = REIONIZATION_DIR / "python_dir"   # contains reion_f (F2Py)


# --- Data (large, usually outside the repository) ------------------------

DATA_DIR = _path_from_env("DATA_DIR", PROJECT_ROOT / "data")
CATALOGUE_DIR = _path_from_env("CATALOGUE_DIR", GRAVWAVE_DIR / "catalogues")
PLANCK_DATA_DIR = _path_from_env(
    "PLANCK_DATA_DIR", PROJECT_ROOT / "COM_Likelihood_Data-baseline_R3.00"
)

# Individual Planck likelihood files (.clik), built from PLANCK_DATA_DIR.
# clik expects plain strings, so pass e.g. str(CLIK_PLIK_LITE_TTTEEE).
PLANCK_PLC_DIR = PLANCK_DATA_DIR / "baseline" / "plc_3.0"
CLIK_COMMANDER_LOWL_TT = PLANCK_PLC_DIR / "low_l" / "commander" / "commander_dx12_v3_2_29.clik"
CLIK_SIMALL_LOWL_EE = PLANCK_PLC_DIR / "low_l" / "simall" / "simall_100x143_offlike5_EE_Aplanck_B.clik"
CLIK_PLIK_LITE_TTTEEE = PLANCK_PLC_DIR / "hi_l" / "plik_lite" / "plik_lite_v22_TTTEEE.clik"

# Observational data files (small, kept inside the repository).
REIONIZATION_OBS_DIR = REIONIZATION_DIR / "ObsData"
REIONIZATION_SPEC_DIR = REIONIZATION_DIR / "SpecData"   # stellar SEDs (dndm_PopII_*, dndm_PopIII_*)
LYMAN_LIMIT_DATA_FILE = REIONIZATION_OBS_DIR / "Lyman_limit.dat"
GAMMA_DATA_FILE = REIONIZATION_OBS_DIR / "gamma_data_all_combined.dat"
TB_DATA_FILE = GLOBAL21CM_DIR / "brightness_temp.txt"


# --- Outputs (generated, never committed) --------------------------------

OUTPUT_DIR = _path_from_env("OUTPUT_DIR", PROJECT_ROOT / "outputs")
CHAINS_DIR = _path_from_env("CHAINS_DIR", OUTPUT_DIR / "chains")
FIGURES_DIR = _path_from_env("FIGURES_DIR", OUTPUT_DIR / "figures")
LOG_DIR = _path_from_env("LOG_DIR", OUTPUT_DIR / "logs")

# Directories the code is allowed to create automatically.
_WRITABLE_DIRS = [OUTPUT_DIR, CHAINS_DIR, FIGURES_DIR, LOG_DIR]

# Folders whose Python / F2Py modules should be importable by name.
_CODE_DIRS = [CAMB_DIR, GLOBAL21CM_DIR, GRAVWAVE_DIR, MCMC_DIR, REIONIZATION_PY_DIR]


# --- Helpers -------------------------------------------------------------

def ensure_dirs() -> None:
    """Create the output directories if they don't exist. Call once at startup."""
    for d in _WRITABLE_DIRS:
        d.mkdir(parents=True, exist_ok=True)


def add_code_dirs_to_syspath() -> None:
    """Optional: make modules in the code folders importable by name.

    Replaces the hardcoded sys.path.append(...) lines, so that
        from gw_likelihood import GWLikelihood      # GravWave/
        from gw_model import GWModel                # GravWave/
        import reion_f                              # ReionizationWithF2Py/python_dir/
        from brightnessTemp import Generate21cmSignal   # global_21cm_signal/
    work on any machine.
    """
    for d in _CODE_DIRS:
        s = str(d)
        if d.is_dir() and s not in sys.path:
            sys.path.insert(0, s)


_ALL_NAMES = [
    "PROJECT_ROOT",
    "CAMB_DIR", "PLANCK_CODE_DIR", "GLOBAL21CM_DIR", "GRAVWAVE_DIR",
    "MCMC_DIR", "REIONIZATION_DIR", "REIONIZATION_PY_DIR",
    "DATA_DIR", "CATALOGUE_DIR" "PLANCK_DATA_DIR",
    "CLIK_COMMANDER_LOWL_TT", "CLIK_SIMALL_LOWL_EE", "CLIK_PLIK_LITE_TTTEEE",
    "REIONIZATION_OBS_DIR", "REIONIZATION_SPEC_DIR", "LYMAN_LIMIT_DATA_FILE", "GAMMA_DATA_FILE", "TB_DATA_FILE",
    "OUTPUT_DIR", "CHAINS_DIR", "FIGURES_DIR", "LOG_DIR",
]


def describe() -> str:
    """Human-readable summary of resolved paths (handy at the top of logs)."""
    g = globals()
    width = max(len(n) for n in _ALL_NAMES)
    lines = []
    for n in _ALL_NAMES:
        p = g[n]
        flag = "" if p.exists() else "   (missing)"
        lines.append(f"{n:<{width}} = {p}{flag}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
