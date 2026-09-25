"""
    Contact sheet of what each method keeps at the MIDOG box, for eyeballing: RGB crop of the 73 px patch with production's
    template footprint (yellow), the saliency kept mask (cyan fill) and its template (cyan), and the 1 um template
    (magenta dotted). Picks cases from out/geometry_midog_box.csv by category.
    Usage: contact_sheet.py  -> out/contact_sheet.png
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import common as C  # noqa: E402

PER_CATEGORY = 6


def pick(geo: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    """
        A few annotations per category, deterministically.

        geo (pd.DataFrame): geometry_midog_box.csv.

        Returns list[tuple[str, pd.Series]]: (category, row).
    """
    rng = np.random.default_rng(7)
    cats = {
        "same object, typical": geo[(geo.best_in_gray >= 0.5) & (geo.sal_n_kept == 1)],
        "saliency keeps >= 2 masses": geo[geo.sal_n_kept >= 2],
        "gray gate refused": geo[geo.prod_status != "refined"],
        "different object (best < 50% in gray)": geo[geo.best_in_gray < 0.5],
    }
    out = []
    for name, g in cats.items():
        idx = rng.permutation(len(g))[:PER_CATEGORY]
        out += [(name, g.iloc[i]) for i in idx]
    return out


def square(ax, base_size, tx, ty, ox, oy, **kw):
    """
        Draw a template footprint on patch axes.

        ax (matplotlib Axes): target axes.
        base_size (int): template side.
        tx (float): template centre x, image coordinates.
        ty (float): template centre y.
        ox (int): patch origin x in image coordinates.
        oy (int): patch origin y.
        kw: Rectangle styling.

        Returns None.
    """
    h = int(base_size) // 2
    ax.add_patch(Rectangle((int(round(tx)) - h - ox - 0.5, int(round(ty)) - h - oy - 0.5), int(base_size), int(base_size), fill=False, **kw))


def main():
    """
        Render the sheet.

        Returns None.
    """
    geo = pd.read_csv(C.OUT_DIR / "geometry_midog_box.csv")
    picks = pick(geo)
    rois = C.roi_table().set_index("file_name")
    ncol = PER_CATEGORY
    nrow = int(np.ceil(len(picks) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.9 * nrow), squeeze=False)
    loaded = {}
    for i, (cat, r) in enumerate(picks):
        ax = axes[i // ncol][i % ncol]
        if r.file_name not in loaded:
            loaded.clear()
            loaded[r.file_name] = C.load_roi(rois.loc[r.file_name, "path"])
        roi = loaded[r.file_name]
        ix, iy = int(round(r.cx)), int(round(r.cy))
        ox, oy = ix - C.HALF_PATCH, iy - C.HALF_PATCH
        ax.imshow(roi["rgb"][oy:oy + C.PATCH, ox:ox + C.PATCH])
        ann = roi["gt"][roi["gt"].ann_id == r.ann_id].iloc[0]
        _, sal = C.saliency_spec(roi["hem"], *C.midog_box(ann), roi["mpp"], 51, dilate=False)
        if sal["mask"] is not None:
            ov = np.zeros((C.PATCH, C.PATCH, 4))
            ov[sal["mask"]] = (0, 0.9, 1, 0.35)
            ax.imshow(ov)
        ax.add_patch(Rectangle((C.HALF_PATCH - 25 - 0.5, C.HALF_PATCH - 25 - 0.5), 51, 51, fill=False, edgecolor="white", lw=0.8, ls="--"))
        square(ax, r.prod_base, r.prod_tx, r.prod_ty, ox, oy, edgecolor="yellow", lw=1.4)
        square(ax, r.sal_base, r.sal_tx, r.sal_ty, ox, oy, edgecolor="cyan", lw=1.4)
        square(ax, r.sal1_base, r.sal1_tx, r.sal1_ty, ox, oy, edgecolor="magenta", lw=1.0, ls=":")
        ax.plot(C.HALF_PATCH, C.HALF_PATCH, "w+", ms=6)
        ax.set_title(f"{cat}\n{r.file_name} {int(r.ann_id)}  prod {int(r.prod_base)}{'' if r.prod_status == 'refined' else ' (fallback)'} / sal {int(r.sal_base)} / kept {int(r.sal_n_kept)}", fontsize=6.5)
        ax.axis("off")
    for j in range(len(picks), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("yellow = production template (gray Otsu gate), cyan = saliency template (cyan fill = kept chromatin), magenta dotted = saliency + 1 um, white dashed = MIDOG 51 px box", fontsize=8)
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "contact_sheet.png", dpi=140)
    print("->", C.OUT_DIR / "contact_sheet.png")


if __name__ == "__main__":
    main()
