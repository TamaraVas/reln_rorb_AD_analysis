"""
figure_style — minimal matplotlib styling helpers for the RELN-RORB analysis suite.

apply_figure_style() sets publication rcParams; panel_letter() stamps bold
capital panel labels top-left outside the axes. Appearance only — not required
to reproduce any statistic.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.family"] = font
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.12, dy=1.03, case="upper", fontsize=13):
    l = letter.upper() if case == "upper" else letter.lower()
    ax.text(dx, dy, l, transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", va="bottom", ha="right")


def panel_letter(ax, letter, dx=-0.12, dy=1.03, case="upper", fontsize=13):
    l = letter.upper() if case == "upper" else letter.lower()
    ax.text(dx, dy, l, transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", va="bottom", ha="right")

