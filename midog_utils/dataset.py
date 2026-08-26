"""Loading MIDOG++ annotations and ROI images.

Two dataset conventions matter here and are both easy to get wrong:

1. ``bbox`` in ``databases/MIDOG++.json`` is ``[x1, y1, x2, y2]`` in absolute pixels --
   *not* the COCO ``[x, y, width, height]`` the file otherwise imitates. Values are
   floats; 3114 annotations carry ``.5`` coordinates.
2. Every one of the 26286 boxes is exactly 50x50 and synthetic. The real ground truth
   is a single SlideRunner point click and the box is ``point +- 25``. Everything
   downstream of this module is therefore point-based, never IoU-based.

ROI images are read with ``tifffile`` rather than ``openslide``. These files are flat
~2 mm^2 regions, not pyramidal whole-slide images, so there is nothing to gain from a
slide reader -- and ``openslide`` is not installed in the environment that has a
working ``cv2``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

BOX_SIZE = 50  # every MIDOG++ bbox, mitotic and look-alike alike

MITOTIC = 1
LOOKALIKE = 2  # "not mitotic figure" -- pathologist-marked hard negative
CATEGORY_NAMES = {MITOTIC: "mitotic figure", LOOKALIKE: "not mitotic figure"}

# Image ids 151-200 are in the JSON but carry zero annotations, are absent from
# datasets_xvalidation.csv, and are not downloadable via Setup.ipynb. 553 - 50 = 503,
# which is the specimen count the README quotes.
UNANNOTATED_IMAGE_IDS = range(151, 201)

# The JSON and datasets_xvalidation.csv disagree on two spellings. Normalise to the
# JSON's, since that is what this package joins on.
TUMOR_ALIASES = {"canine lymphoma": "canine lymphosarcoma"}
SCANNER_ALIASES = {"Hamammatsu XR": "Hamamatsu XR"}

# TIFF ResolutionUnit tag values -> microns per unit.
_RESOLUTION_UNIT_UM = {2: 25400.0, 3: 10000.0}  # 2 = inch, 3 = centimetre


def canonical_tumor(name: str) -> str:
    return TUMOR_ALIASES.get(name, name)


def canonical_scanner(name: str) -> str:
    return SCANNER_ALIASES.get(name, name)


def load_annotations(json_path="databases/MIDOG++.json", drop_unannotated=True):
    """Read the COCO-ish JSON into two tidy frames.

    Returns ``(images, annotations)`` where

    * ``images``      -- image_id, file_name, width, height, tumor_type
    * ``annotations`` -- ann_id, image_id, file_name, cx, cy, category_id, n_votes,
                         n_mitotic_votes, unanimous

    ``cx``/``cy`` are the recovered point clicks (box centres), as floats.
    ``unanimous`` marks annotations every expert scored the same way -- useful for
    splitting recall into "missed an obvious mitosis" vs "missed one the pathologists
    themselves argued about".
    """
    raw = json.loads(Path(json_path).read_text())

    images = pd.DataFrame(
        [
            {
                "image_id": im["id"],
                "file_name": im["file_name"],
                "width": im["width"],
                "height": im["height"],
                "tumor_type": canonical_tumor(im["tumor_type"]),
            }
            for im in raw["images"]
        ]
    )

    id_to_name = dict(zip(images["image_id"], images["file_name"]))

    rows = []
    for a in raw["annotations"]:
        x1, y1, x2, y2 = a["bbox"]  # TLBR, not xywh
        votes = a.get("labels", [])
        rows.append(
            {
                "ann_id": a["id"],
                "image_id": a["image_id"],
                "file_name": id_to_name[a["image_id"]],
                "cx": (x1 + x2) / 2.0,
                "cy": (y1 + y2) / 2.0,
                "w": x2 - x1,
                "h": y2 - y1,
                "category_id": a["category_id"],
                "n_votes": len(votes),
                "n_mitotic_votes": sum(1 for v in votes if v == MITOTIC),
                "unanimous": len(set(votes)) == 1 if votes else False,
            }
        )
    annotations = pd.DataFrame(rows)

    if drop_unannotated:
        keep = ~images["image_id"].isin(UNANNOTATED_IMAGE_IDS)
        images = images[keep].reset_index(drop=True)

    return images, annotations


def load_slide_metadata(csv_path="datasets_xvalidation.csv"):
    """Scanner / origin / species, which live only in the CSV -- not in the JSON."""
    df = pd.read_csv(csv_path, delimiter=";")
    df = df.rename(columns={"Slide": "image_id"})
    df["Tumor"] = df["Tumor"].map(canonical_tumor)
    df["Scanner"] = df["Scanner"].map(canonical_scanner)
    return df


def check_invariants(annotations: pd.DataFrame) -> None:
    """Fail loudly if the two conventions this package relies on ever stop holding."""
    bad_size = annotations[(annotations["w"] != BOX_SIZE) | (annotations["h"] != BOX_SIZE)]
    if len(bad_size):
        raise AssertionError(f"{len(bad_size)} annotations are not {BOX_SIZE}x{BOX_SIZE}")

    cats = set(annotations["category_id"].unique())
    if not cats <= {MITOTIC, LOOKALIKE}:
        raise AssertionError(f"unexpected category ids: {cats - {MITOTIC, LOOKALIKE}}")


def image_annotations(annotations: pd.DataFrame, file_name: str, category_id=None):
    sel = annotations[annotations["file_name"] == file_name]
    if category_id is not None:
        sel = sel[sel["category_id"] == category_id]
    return sel.reset_index(drop=True)


def points(df: pd.DataFrame) -> np.ndarray:
    """(N, 2) float array of (x, y) centres."""
    if len(df) == 0:
        return np.zeros((0, 2), dtype=np.float64)
    return df[["cx", "cy"]].to_numpy(dtype=np.float64)


def load_roi(path) -> np.ndarray:
    """Full-resolution RGB uint8 array for one ROI.

    Handles both TIFF layouts present in this dataset: flat uncompressed single-page
    files (001-406) and 7-level LZW pyramids (505/506). Alpha is stripped -- the files
    are RGBA but the alpha channel is uniform.
    """
    with tifffile.TiffFile(str(path)) as tf:
        series = tf.series[0]
        levels = getattr(series, "levels", None)
        arr = levels[0].asarray() if levels else series.asarray()
    if arr.ndim != 3:
        raise ValueError(f"expected an HxWxC image, got shape {arr.shape}")
    return np.ascontiguousarray(arr[:, :, :3])


def roi_mpp(path) -> float:
    """Microns per pixel from the TIFF resolution tags.

    Honouring ``ResolutionUnit`` is not optional. 505.tiff and 506.tiff store unit 3
    (centimetre) while every other image stores unit 2 (inch). Assuming inches gives
    0.578 um/px for those two -- implying a 13 mm^2 ROI -- instead of the correct
    0.227 um/px and 2.0 mm^2.
    """
    with tifffile.TiffFile(str(path)) as tf:
        page = tf.pages[0]
        num, den = page.tags["XResolution"].value
        unit = int(page.tags["ResolutionUnit"].value)
    if unit not in _RESOLUTION_UNIT_UM:
        raise ValueError(f"unsupported ResolutionUnit {unit} in {path}")
    pixels_per_unit = num / den
    return _RESOLUTION_UNIT_UM[unit] / pixels_per_unit


def roi_area_mm2(path, shape=None) -> float:
    mpp = roi_mpp(path)
    if shape is None:
        with tifffile.TiffFile(str(path)) as tf:
            s = tf.series[0]
            levels = getattr(s, "levels", None)
            shape = (levels[0].shape if levels else s.shape)[:2]
    h, w = shape[0], shape[1]
    return h * w * mpp * mpp / 1e6


def check_roi_scale(path, shape=None, expected=(1.9, 2.1)) -> float:
    """MIDOG++ ROIs are all ~2 mm^2. This is what catches a misread resolution tag."""
    area = roi_area_mm2(path, shape)
    if not (expected[0] <= area <= expected[1]):
        raise AssertionError(
            f"{path}: computed ROI area {area:.2f} mm^2 outside {expected} -- "
            "check the ResolutionUnit tag"
        )
    return area
