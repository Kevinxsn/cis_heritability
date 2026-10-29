"""Pick a fast scratch directory for the transient LD / GRM traffic.

Every region hands R a p x p float64 LD matrix (620 MB at p = 8.8k) and GCTA a
binary GRM.  With dozens of concurrent workers that traffic must land on
node-local storage: the project filesystem here is NFS at ~117 MB/s, the local
NVMe at ~1.5 GB/s, and 48 workers writing LD matrices to NFS turn a
CPU-bound job into an I/O-bound one.

Preference order: an allocation-provided scratch, then node-local NVMe, then
/tmp, and only as a last resort a directory beside the results.
"""

from __future__ import annotations

import getpass
import os
import shutil
from typing import Optional

# candidate roots, best first; formatted with the current user
_CANDIDATES = (
    "$SLURM_TMPDIR",
    "/scratch/$USER/$SLURM_JOB_ID",
    "/scratch/$USER",
    "/data/tmp",
    "/tmp",
)

MIN_FREE_GB = 50.0


def _expand(template: str) -> Optional[str]:
    path = os.path.expandvars(template)
    if "$" in path:  # an undefined variable, e.g. no SLURM_JOB_ID
        return None
    return path


def _usable(path: str, min_free_gb: float) -> bool:
    if not os.path.isdir(path):
        return False
    if not os.access(path, os.W_OK):
        return False
    try:
        free_gb = shutil.disk_usage(path).free / 1e9
    except OSError:
        return False
    return free_gb >= min_free_gb


def resolve_scratch(
    explicit: Optional[str] = None,
    fallback: Optional[str] = None,
    tag: str = "cish2",
    min_free_gb: float = MIN_FREE_GB,
) -> str:
    """Return a writable scratch directory, creating it if needed.

    :param explicit: honour this path if the caller configured one.
    :param fallback: used only when no fast candidate qualifies.
    """
    if explicit:
        os.makedirs(explicit, exist_ok=True)
        return explicit

    user = os.environ.get("USER") or getpass.getuser()
    for template in _CANDIDATES:
        root = _expand(template)
        if root is None or not _usable(root, min_free_gb):
            continue
        path = root if root.rstrip("/").endswith(user) else os.path.join(root, f"{tag}-{user}")
        try:
            os.makedirs(path, exist_ok=True)
        except OSError:
            continue
        if os.access(path, os.W_OK):
            return path

    path = fallback or os.path.join(os.getcwd(), "scratch")
    os.makedirs(path, exist_ok=True)
    return path


def describe(path: str) -> str:
    try:
        usage = shutil.disk_usage(path)
        return f"{path} ({usage.free / 1e9:.0f} GB free)"
    except OSError:
        return path
