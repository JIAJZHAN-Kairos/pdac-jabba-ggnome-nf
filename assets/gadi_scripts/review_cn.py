import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

root, sample = Path(sys.argv[1]), sys.argv[2]
out = root / "graphs" / sample
nodes = pd.read_csv(out / "nodes.tsv", sep="\t")
purple = pd.read_csv(root / "inputs" / f"{sample}.purple.cnv.somatic.tsv", sep="\t")
nodes["chromosome"] = nodes.seqnames.astype(str).str.replace(r"^chr", "", regex=True)
purple["chromosome"] = purple.chromosome.astype(str).str.replace(r"^chr", "", regex=True)
chroms = [str(x) for x in range(1, 23)]
pairs = []
for chrom in chroms:
    ns = nodes.loc[nodes.chromosome.eq(chrom)].sort_values("start")
    ps = purple.loc[purple.chromosome.eq(chrom)].sort_values("start")
    assert not ns.empty and not ps.empty
    assert np.all(ns.start.to_numpy()[1:] > ns.end.to_numpy()[:-1])
    j = 0
    nr = list(ns.itertuples())
    for p in ps.itertuples():
        while j < len(nr) and nr[j].end < p.start:
            j += 1
        k = j
        while k < len(nr) and nr[k].start <= p.end:
            n = nr[k]
            left, right = max(n.start, p.start), min(n.end, p.end)
            if left <= right:
                pairs.append((chrom, left, right, p.copyNumber, n.cn))
            k += 1
comparison = pd.DataFrame(pairs, columns=["chromosome", "start", "end", "PURPLE_CN", "JaBbA_CN"])
comparison["bp"] = comparison.end - comparison.start + 1
comparison["residual"] = comparison.JaBbA_CN - comparison.PURPLE_CN
comparison.to_csv(out / "CN_comparison.tsv.gz", sep="\t", index=False)
known = comparison.loc[np.isfinite(comparison.residual)].copy()
assert not known.empty
weights = known.bp.to_numpy()
residual = known.residual.to_numpy()
metrics = {
    "sample": sample,
    "comparison_scope": "Autosomes; overlap with original PURPLE segments; consistency only, not independent truth",
    "compared_bp": int(weights.sum()),
    "PURPLE_autosomal_bp": int((purple.loc[purple.chromosome.isin(chroms), "end"] -
                                purple.loc[purple.chromosome.isin(chroms), "start"] + 1).sum()),
    "overlap_with_unknown_graph_CN_bp": int(comparison.loc[~np.isfinite(comparison.JaBbA_CN), "bp"].sum()),
    "length_weighted_CN_MAE": float(np.average(np.abs(residual), weights=weights)),
    "length_weighted_CN_bias": float(np.average(residual, weights=weights)),
    "bp_fraction_abs_residual_le_1": float(np.average(np.abs(residual) <= 1, weights=weights)),
    "bp_fraction_abs_residual_gt_2": float(np.average(np.abs(residual) > 2, weights=weights)),
    "length_weighted_JaBbA_CN": float(np.average(known.JaBbA_CN, weights=weights)),
    "length_weighted_PURPLE_CN": float(np.average(known.PURPLE_CN, weights=weights)),
}
(out / "CN_fit_QC.json").write_text(json.dumps(metrics, indent=2) + "\n")

fig = plt.figure(figsize=(15, 8), layout="constrained")
axes = fig.subplot_mosaic([["profile", "profile"], ["scatter", "residual"]])
offset, ticks, labels = 0, [], []
for chrom in chroms:
    ns = nodes.loc[nodes.chromosome.eq(chrom)]
    ps = purple.loc[purple.chromosome.eq(chrom)]
    length = max(ns.end.max(), ps.end.max())
    for table, field, color, lw in [(ps, "copyNumber", "#777777", 1.5), (ns, "cn", "#0072B2", .8)]:
        lines = [((offset + r.start, getattr(r, field)), (offset + r.end, getattr(r, field)))
                 for r in table.itertuples() if np.isfinite(getattr(r, field))]
        axes["profile"].add_collection(LineCollection(lines, colors=color, linewidths=lw))
    ticks.append(offset + length / 2)
    labels.append(chrom)
    offset += length
    axes["profile"].axvline(offset, lw=.3, color="#cccccc")
axes["profile"].set_xlim(0, offset)
axes["profile"].set_ylim(-.1, max(known.PURPLE_CN.max(), known.JaBbA_CN.max()) * 1.1 + 1)
axes["profile"].set_yscale("symlog", linthresh=3)
axes["profile"].set_xticks(ticks, labels)
axes["profile"].set(xlabel="Chromosome", ylabel="Total CN (symlog)",
                    title=f"{sample}: JaBbA (blue), PURPLE (grey); autosomes")
axes["scatter"].scatter(known.PURPLE_CN, known.JaBbA_CN,
                         s=3 + 35 * np.sqrt(weights / weights.max()), alpha=.35, color="#0072B2")
high = max(known.PURPLE_CN.max(), known.JaBbA_CN.max()) * 1.05 + 1
axes["scatter"].plot([0, high], [0, high], "--", lw=.8, color="#555555")
axes["scatter"].set(xlabel="PURPLE CN", ylabel="JaBbA integer CN")
axes["scatter"].set_xscale("symlog", linthresh=3)
axes["scatter"].set_yscale("symlog", linthresh=3)
axes["residual"].hist(residual, bins=80, weights=weights / weights.sum(), color="#0072B2")
axes["residual"].set(xlabel="JaBbA minus PURPLE CN", ylabel="Fraction of compared bases",
                     title=f"Weighted MAE {metrics['length_weighted_CN_MAE']:.3f}; bias {metrics['length_weighted_CN_bias']:.3f}")
fig.savefig(out / "CN_fit_QC.png", dpi=150)
fig.savefig(out / "CN_fit_QC.pdf")
plt.close(fig)
print(json.dumps(metrics), flush=True)
