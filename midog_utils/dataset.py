"""Loading MIDOG++ annotations and ROI images."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

BOX_SIZE = 50  # every MIDOG++ bbox is a synthetic point +- 25, mitotic and look-alike alike

MITOTIC = 1
LOOKALIKE = 2  # "not mitotic figure" -- pathologist-marked hard negative
CATEGORY_NAMES = {MITOTIC: "mitotic figure", LOOKALIKE: "not mitotic figure"}

UNANNOTATED_IMAGE_IDS = range(151, 201)  # in the JSON but carry zero annotations

TUMOR_ALIASES = {"canine lymphoma": "canine lymphosarcoma"}
SCANNER_ALIASES = {"Hamammatsu XR": "Hamamatsu XR"}

_RESOLUTION_UNIT_UM = {2: 25400.0, 3: 10000.0}  # TIFF ResolutionUnit tag: 2 = inch, 3 = centimetre


def canonical_tumor(name: str) -> str:
    return TUMOR_ALIASES.get(name, name)


def canonical_scanner(name: str) -> str:
    return SCANNER_ALIASES.get(name, name)


def load_annotations(json_path="databases/MIDOG++.json", drop_unannotated=True):
    """
        Read the COCO-ish JSON into two tidy frames.

        json_path (str): path to MIDOG++.json.
        drop_unannotated (bool): drop image ids with zero annotations.

        images (pd.DataFrame): image_id, file_name, width, height, tumor_type.
        annotations (pd.DataFrame): ann_id, image_id, file_name, cx, cy, category_id,
            n_votes, n_mitotic_votes, unanimous.
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
    """
        Scanner / origin / species, which live only in the CSV -- not in the JSON.

        csv_path (str): path to the crossvalidation CSV.

        Returns pd.DataFrame: the CSV, with image_id/Tumor/Scanner canonicalised.
    """
    df = pd.read_csv(csv_path, delimiter=";")
    df = df.rename(columns={"Slide": "image_id"})
    df["Tumor"] = df["Tumor"].map(canonical_tumor)
    df["Scanner"] = df["Scanner"].map(canonical_scanner)
    return df


def check_invariants(annotations: pd.DataFrame) -> None:
    """
        Fail loudly if the two conventions this package relies on ever stop holding.

        annotations (pd.DataFrame): must have w, h, category_id columns.

        Returns None. Raises AssertionError naming the violated invariant.
    """
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
    """
        (N, 2) float array of (x, y) centres.

        df (pd.DataFrame): must have cx, cy columns.

        Returns np.ndarray: shape (N, 2), float64; empty when ``df`` is empty.
    """
    if len(df) == 0:
        return np.zeros((0, 2), dtype=np.float64)
    return df[["cx", "cy"]].to_numpy(dtype=np.float64)


def load_roi(path) -> np.ndarray:
    """
        Full-resolution RGB uint8 array for one ROI. Handles both flat single-page
        TIFFs and pyramidal ones; the alpha channel is stripped.

        path (str or Path): path to the ROI TIFF.

        Returns np.ndarray: HxWx3 uint8 array.
    """
    with tifffile.TiffFile(str(path)) as tf:
        series = tf.series[0]
        levels = getattr(series, "levels", None)
        arr = levels[0].asarray() if levels else series.asarray()
    if arr.ndim != 3:
        raise ValueError(f"expected an HxWxC image, got shape {arr.shape}")
    return np.ascontiguousarray(arr[:, :, :3])


def roi_mpp(path) -> float:
    """
        Microns per pixel from the TIFF resolution tags. Honours ResolutionUnit --
        some scanners in this dataset store centimetres rather than inches.

        path (str or Path): path to the ROI TIFF.

        Returns float: microns per pixel.
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
    """
        MIDOG++ ROIs are all ~2 mm^2. This is what catches a misread resolution tag.

        path (str or Path): path to the ROI TIFF.
        shape (tuple or None): (height, width) to use instead of reading the TIFF.
        expected (tuple[float, float]): accepted (min, max) area range, in mm^2.

        Returns float: the computed ROI area, in mm^2.
        Raises AssertionError when outside ``expected``.
    """
    area = roi_area_mm2(path, shape)
    if not (expected[0] <= area <= expected[1]):
        raise AssertionError(
            f"{path}: computed ROI area {area:.2f} mm^2 outside {expected} -- "
            "check the ResolutionUnit tag"
        )
    return area
