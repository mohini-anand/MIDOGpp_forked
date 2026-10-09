"""v2_variants.py with production hem computed on 1024-row blocks (bit-identical, verified), to measure the full worker's peak RSS."""
import sys
import v2lib
from midog_utils import channels as _ch
_orig = _ch.to_channel
def _blocked(rgb, name, rows=1024):
    import numpy as np
    out = np.empty(rgb.shape[:2], np.float32)
    for r0 in range(0, rgb.shape[0], rows):
        out[r0:r0 + rows] = _orig(rgb[r0:r0 + rows], name)
    return out
_ch.to_channel = _blocked
exec(open('v2_variants.py').read())
