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

UNANNOTATED_IMAGE_IDS = range(151, 201)  # in the JSON but carry zero annotations

TUMOR_ALIASES = {"canine lymphoma": "canine lymphosarcoma"}

_RESOLUTION_UNIT_UM = {2: 25400.0, 3: 10000.0}  # TIFF ResolutionUnit tag: 2 = inch, 3 = centimetre


def canonical_tumor(name: str) -> str:
    return TUMOR_ALIASES.get(name, name)


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


def image_annotations(annotations: pd.DataFrame, file_name: str, category_id=None):
    sel = annotations[annotations["file_name"] == file_name]
    if category_id is not None:
        sel = sel[sel["category_id"] == category_id]
    return sel.reset_index(drop=True)


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
