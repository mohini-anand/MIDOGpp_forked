"""Step 2: embed every candidate crop with each encoder under test.

Runs under the **system** python (torch 2.2.2), not the anaconda one -- see the module
docstring of `click_rank_stage1.py` for why the experiment is split across interpreters.
`torch.from_numpy` is unusable here (torch built for numpy 1.x, system numpy is 2.x), so
arrays cross into torch through `torch.frombuffer`.

Embeddings are a function of the ROI alone -- they do not depend on which cell was clicked
-- so they are computed once per (ROI, encoder) and reused by all five seeds. That is also
the product architecture: this is ingest-time work, and only the cosine is click-time.

The 16 um arm is the exact centre crop of the stored 32 um crop, so "does the answer depend
on how much context the encoder sees" is a control that costs one extra forward pass and no
extra extraction.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

OUT = os.environ.get("CLICK_RANK_DIR", os.path.expanduser("~/.cache/annotatedx_click_rank"))
ROIS = ("301.tiff", "246.tiff")
BATCH = 64

# (arm name, timm spec, input side, crop micron field, extra create_model kwargs).
# Ordered most-informative-first: the two pathology-domain encoders are the arms the
# conclusion turns on, so they finish first if the run has to be cut short.
#
# **Every encoder sees the identical stored crop at its stored resolution** (128 px over a
# 32 um field; 126 for the patch-14 model, which needs a multiple of 14). Nothing is
# resampled, so no arm can win or lose on interpolation.
#
# Why 128 rather than the checkpoints' native 224: what a ViT sees is set by
# microns-per-token, not by image size. These encoders were trained near 0.5 um/px, so a
# patch-16 token covered ~8 um of tissue. Feeding this 0.25 um/px crop at 224 px would make
# a token cover 2.3 um -- 3.5x too zoomed -- while feeding it at 128 px gives 4 um, 2x too
# zoomed, and costs a quarter as much. Both DINO and DINOv2 pretraining used 96-98 px local
# crops, so 128 px is inside these models' training experience; timm interpolates the
# position embeddings at `create_model` time.
ENCODERS = [
    ("lunit_dino_vits16", "hf-hub:1aurent/vit_small_patch16_224.lunit_dino",  128, 32, {"img_size": 128}),
    ("kaiko_vits16", "hf-hub:1aurent/vit_small_patch16_224.kaiko_ai_towards_large_pathology_fms", 128, 32, {"img_size": 128}),
    # 16 um variants: the same encoder on the cell alone rather than cell-plus-context,
    # so a negative result cannot be blamed on the embedding being dominated by
    # surrounding tissue. The 16 um crop is the exact centre of the stored 32 um one.
    ("lunit_dino_vits16_16um", "hf-hub:1aurent/vit_small_patch16_224.lunit_dino", 128, 16, {"img_size": 128}),
    ("resnet18_in1k",    "resnet18.a1_in1k",                                  128, 32, {}),
    ("resnet18_in1k_16um", "resnet18.a1_in1k",                                128, 16, {}),
    # A mitosis-trained backbone, as a **leaky upper bound**: the ResNet-18 body of the
    # MIT-licensed FCOS detector from `jonas-amme/FCOS_Inference_CLI` (F1 0.737-0.753 on
    # the MIDOG 2022 test set), which was trained on MIDOG++ -- including 301.tiff and
    # 246.tiff. Its features have seen these ROIs' mitotic figures with labels, so it
    # cannot support a product claim. It answers a different question: if even a
    # representation built for this exact task does not make the click pay, the negative
    # result is about the click rather than about encoder choice.
    ("fcos_midogpp_r18", "local:FCOS_18.ckpt",                                 128, 32, {}),
]


def centre_crop(a: np.ndarray, frac: float) -> np.ndarray:
    if frac >= 1.0:
        return np.asarray(a)
    p = a.shape[1]
    k = int(round(p * frac))
    o = (p - k) // 2
    return np.asarray(a[:, o:o + k, o:o + k, :])


def from_tensor(t: torch.Tensor) -> np.ndarray:
    """torch -> numpy without the numpy bridge (unavailable: torch built for numpy 1.x)."""
    try:
        return np.from_dlpack(t.contiguous()).copy()
    except Exception:
        return np.array(t.contiguous().tolist(), dtype=np.float32)


def to_tensor(batch_u8: np.ndarray) -> tuple:
    """uint8 NHWC -> float32 NCHW tensor, through the buffer protocol.

    Returns the backing array as well: `torch.frombuffer` does not copy, so the numpy
    array must outlive the tensor. Two copies were removed here (`tobytes()` and
    `clone()`); on this CPU they cost about as much as the convolutions.
    """
    b = np.ascontiguousarray(batch_u8.transpose(0, 3, 1, 2), dtype=np.float32)
    t = torch.frombuffer(memoryview(b), dtype=torch.float32).view(b.shape)
    return t / 255.0, b


def embed(model, crops: np.ndarray, side: int, mean, std, frac: float) -> np.ndarray:
    mean_t = torch.tensor(mean).view(1, 3, 1, 1)
    std_t = torch.tensor(std).view(1, 3, 1, 1)
    outs = []
    with torch.no_grad():
        for i in range(0, len(crops), BATCH):
            x, _keepalive = to_tensor(centre_crop(crops[i:i + BATCH], frac))
            if x.shape[-1] != side:
                x = F.interpolate(x, size=(side, side), mode="bicubic",
                                  align_corners=False, antialias=x.shape[-1] > side)
            x = (x - mean_t) / std_t
            f = model(x)
            f = F.normalize(f.float(), dim=1)
            outs.append(from_tensor(f))
    return np.concatenate(outs).astype(np.float32)


def load_local(fname: str):
    """The FCOS detector's ResNet-18 body as a pooled-feature encoder.

    The checkpoint is a Lightning bundle whose `model.backbone.body.*` keys are exactly a
    torchvision ResNet-18 (standard BN, not FrozenBatchNorm -- `num_batches_tracked` is
    present), so they load into `torchvision.models.resnet18` directly; `fc` is the only
    unmatched module and is replaced by identity. Normalisation is ImageNet, which is what
    torchvision's detection transform uses.
    """
    import torchvision as tv
    path = os.path.join(os.environ.get("CLICK_RANK_CKPT", OUT), fname)
    ck = torch.load(path, map_location="cpu")
    sd = ck.get("state_dict", ck)
    pre = "model.backbone.body."
    body = {k[len(pre):]: v for k, v in sd.items() if k.startswith(pre)}
    if not body:
        raise KeyError(f"no {pre}* keys in {path}")
    m = tv.models.resnet18(weights=None)
    missing, unexpected = m.load_state_dict(body, strict=False)
    if unexpected or [k for k in missing if not k.startswith("fc.")]:
        raise ValueError(f"backbone load mismatch: missing={missing} unexpected={unexpected}")
    m.fc = torch.nn.Identity()
    return m.eval(), (0.485, 0.456, 0.406), (0.229, 0.224, 0.225)


def main():
    import timm
    torch.set_grad_enabled(False)
    torch.set_num_threads(6)
    only = sys.argv[1:] or None

    for name, spec, side, crop_um, kw in ENCODERS:
        if only and name not in only:
            continue
        if spec.startswith("local:"):
            model, mean, std = load_local(spec[len("local:"):])
        else:
            model = timm.create_model(spec, pretrained=True, num_classes=0, **kw).eval()
            cfg = timm.data.resolve_model_data_config(model)
            mean, std = cfg["mean"], cfg["std"]
        frac = crop_um / 32.0
        for fn in ROIS:
            dst = f"{OUT}/emb_{fn}_{name}.npy"
            if os.path.exists(dst):
                print(f"skip {dst}", flush=True)
                continue
            t0 = time.time()
            # mmap: a resident 1 GB copy per job halves this CPU's throughput
            crops = np.load(f"{OUT}/crops_{fn}.npy", mmap_mode="r")
            e = embed(model, crops, side, mean, std, frac)
            np.save(dst, e)
            seeds = np.load(f"{OUT}/seedcrops_{fn}.npy")
            np.save(f"{OUT}/embseed_{fn}_{name}.npy",
                    embed(model, seeds, side, mean, std, frac))
            json.dump({"encoder": name, "timm_spec": spec, "input_side": side,
                       "crop_um": crop_um, "dim": int(e.shape[1]),
                       "mean": list(map(float, mean)), "std": list(map(float, std)),
                       "n": int(e.shape[0]), "seconds": round(time.time() - t0, 1)},
                      open(f"{OUT}/embmeta_{fn}_{name}.json", "w"), indent=1)
            print(f"[{name}] {fn}: {e.shape} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
