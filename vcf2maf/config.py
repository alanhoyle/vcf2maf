"""
Load user/site configuration for vcf2maf tools.

The config file uses TOML format with a [defaults] table.  Keys are
argparse dest names (underscores).  String values starting with ~ are
expanded.  Numeric values are passed through as-is.

Search order (first found wins):
  1. Path passed to load_config()
  2. ~/.vcf2maf.toml
  3. .vcf2maf.toml  (current working directory)

Requires Python ≥ 3.11 (tomllib stdlib) or the ``tomli`` back-port for
earlier versions.
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import tomllib  # Python ≥ 3.11
except ModuleNotFoundError:
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ModuleNotFoundError:
        tomllib = None  # type: ignore[assignment]

_SEARCH_PATHS = [
    Path("~/.vcf2maf.toml").expanduser(),
    Path(".vcf2maf.toml"),
]


def load_config(path: Optional[str] = None) -> Dict[str, Any]:
    """
    Return a dict of defaults from the config file.

    Keys match argparse dest names (underscores).  String values starting
    with ~ are expanded.  Returns an empty dict when no config file is found
    or when TOML support is unavailable.
    """
    if tomllib is None:
        return {}

    candidates = [Path(path)] if path else _SEARCH_PATHS
    for candidate in candidates:
        if candidate.exists():
            with candidate.open("rb") as fh:
                data = tomllib.load(fh)
            defaults = data.get("defaults", {})
            return {
                k: (os.path.expanduser(v) if isinstance(v, str) and v.startswith("~") else v)
                for k, v in defaults.items()
            }
    return {}
