#!/usr/bin/env python
"""Retag SAW annotated BAM with scAPAtrap/umi_tools-compatible tags. v2.

v1 (pilot): CB / UB / GX / GN. v2 adds RX = same fixed-length UMI as UB —
umi_tools dedup (called by scAPAtrap dedupByPos) reads UMI from tag RX by
default, and SAW's native RX can be heterogeneous (len 21-22 on the rat
thymus run) which trips umi_tools' same-length assertion.

  CB = "{Cx}_{Cy}"           (spatial barcode string)
  UB = RX = (UR + 'NNNNN')[:5]  (UMI padded to fixed length 5)
  GX = GI  (gene ID)
  GN = GS  (gene symbol)
"""
import pysam, sys, time

IN, OUT = sys.argv[1], sys.argv[2]

t0 = time.time()
bin_ = pysam.AlignmentFile(IN)
bout = pysam.AlignmentFile(OUT, "wb", template=bin_)
n = n_gi = 0
for r in bin_.fetch(until_eof=True):
    if r.has_tag('Cx') and r.has_tag('Cy'):
        r.set_tag('CB', "{}_{}".format(r.get_tag('Cx'), r.get_tag('Cy')), 'Z')
    if r.has_tag('UR'):
        ub = (r.get_tag('UR') + 'NNNNN')[:5]
        r.set_tag('UB', ub, 'Z')
        r.set_tag('RX', ub, 'Z')
    if r.has_tag('GI'):
        r.set_tag('GX', r.get_tag('GI'), 'Z'); n_gi += 1
    if r.has_tag('GS'):
        r.set_tag('GN', r.get_tag('GS'), 'Z')
    bout.write(r)
    n += 1
    if n % 10000000 == 0:
        print("  {} reads... ({:.1f} min, {:.0f}% gene-tagged)".format(
            n, (time.time() - t0) / 60.0, 100.0 * n_gi / n), flush=True)
bout.close()
print("retag done: {} reads ({} gene-tagged, {:.1f}%) in {:.1f} min".format(
    n, n_gi, 100.0 * n_gi / n, (time.time() - t0) / 60.0), flush=True)
