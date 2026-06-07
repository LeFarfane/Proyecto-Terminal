"""Compatibility shims for upsetplot 0.9.0 under modern pandas / matplotlib.

upsetplot 0.9.0 is the latest release and predates both pandas 3.0 and
matplotlib 3.10. Two independent incompatibilities break ``UpSet.plot()``:

1. **Marker colours (pandas >= 3.0).** ``plot_matrix`` fills default dot colours
   with chained inplace fillna, e.g. ``styles["facecolor"].fillna(c, inplace=True)``.
   Under pandas 3.0 Copy-on-Write those operate on a temporary copy and silently
   do nothing, leaving NaN colours that matplotlib rejects with
   ``ValueError: Invalid RGBA argument: nan``. (Also emits noisy deprecation
   warnings under pandas 2.x.)

2. **Count labels (matplotlib >= 3.10).** ``_label_sizes`` computes a text margin
   as ``0.01 * abs(np.diff(ax.get_xlim()))``. ``np.diff`` returns a 1-element
   array, so the label x/y position is an array. matplotlib < 3.10 tolerated
   that; 3.10 raises ``TypeError: only 0-dimensional arrays can be converted to
   Python scalars`` at draw time, breaking ``show_counts`` / ``show_percentages``.

There is no fixed upstream release to upgrade to, so we re-derive each method's
source and apply the minimal, behaviour-preserving rewrite upstream would.
Deriving from source (rather than vendoring whole methods) keeps this working if
the surrounding logic changes in a future release. Import this module once; the
patch is applied as an import side effect and is idempotent.
"""
from __future__ import annotations

import inspect
import re
import textwrap

from upsetplot import plotting as _up

_FLAG = "_upset_compat_patched"


def _patch_method(method_name: str, subs: list[tuple[str, str]]) -> bool:
    """Re-derive `UpSet.<method_name>` source, apply regex `subs`, and rebind.

    Returns True if the method is patched (or was already), False if it could
    not be sourced or no substitution matched (left untouched).
    """
    method = getattr(_up.UpSet, method_name, None)
    if method is None:
        return False
    if getattr(method, _FLAG, False):
        return True
    try:
        src = textwrap.dedent(inspect.getsource(method))
    except (OSError, TypeError):
        return False

    total = 0
    for pattern, repl in subs:
        src, n = re.subn(pattern, repl, src)
        total += n
    if total == 0:
        return False  # nothing matched — upstream may have already fixed it

    namespace: dict = {}
    exec(compile(src, f"<upset_compat:{method_name}>", "exec"), _up.__dict__, namespace)
    new_method = namespace[method_name]
    new_method.__dict__[_FLAG] = True
    setattr(_up.UpSet, method_name, new_method)
    return True


def apply() -> dict[str, bool]:
    """Apply all shims. Returns {method_name: patched?} for diagnostics."""
    results = {}

    # (1) pandas COW: styles["col"].fillna(VALUE, inplace=True)
    #                  -> styles["col"] = styles["col"].fillna(VALUE)
    results["plot_matrix"] = _patch_method(
        "plot_matrix",
        [(
            r'(styles\["[^"]+"\])\.fillna\((.+?),\s*inplace=True\)',
            r"\1 = \1.fillna(\2)",
        )],
    )

    # (2) matplotlib >=3.10: np.diff(...) 1-element array used as a text position.
    #     Coerce the margin basis to a scalar with .item().
    results["_label_sizes"] = _patch_method(
        "_label_sizes",
        [(
            r"np\.diff\((ax\.get_[xy]lim\(\))\)",
            r"np.diff(\1).item()",
        )],
    )
    return results


apply()
