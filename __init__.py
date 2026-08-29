"""Shim so this git checkout is importable as ``molecular_qm_util``.

The installable package lives in ``./molecular_qm_util/``. When pytest (or
PyCharm) puts the parent repo on ``sys.path``, this directory would otherwise
become an empty namespace package and shadow that inner package — including
stale leftover folders such as ``pubchempy_scripts/`` that contain only
``__pycache__``.
"""
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
_INNER = _ROOT / "molecular_qm_util"
# Inner package first so real submodules win; checkout root second so pytest can
# still import ``molecular_qm_util.tests`` from this tree.
__path__ = [str(_INNER), str(_ROOT)]

_init = _INNER / "__init__.py"
exec(compile(_init.read_text(encoding="utf-8"), str(_init), "exec"))
