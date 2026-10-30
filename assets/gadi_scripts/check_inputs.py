"""Apply the pilot input checks to one authorized cohort sample."""
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1]).resolve()
out = root / "prepared"
out.mkdir(exist_ok=True)
samples = list(csv.DictReader((root / "pilot_samples.tsv").open(), delimiter="\t"))
selected = sys.argv[2]
samples = [s for s in samples if s["sample"] == selected]
assert len(samples) == 1
autosomes = [f"chr{i}" for i in range(1, 23)]

for r in csv.DictReader((root / "input_sha256.tsv").open(), delimiter="\t"):
    if not r["file"].startswith(selected + "."):
        continue
    p = root / "inputs" / r["file"]
    assert p.stat().st_size == int(r["bytes"]), p
    with p.open("rb") as handle:
        assert hashlib.file_digest(handle, "sha256").hexdigest() == r["sha256"], p

summaries = []
for sample in samples:
    name = sample["sample"]
    dest = out / name
    dest.mkdir(exist_ok=True)
    prefix = root / "inputs" / name
    headers, records, lengths = [], {}, {}
    with gzip.open(str(prefix) + ".esvee.somatic.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                match = re.match(r"##contig=<ID=([^,]+),length=(\d+)", line)
                if match:
                    lengths[match[1]] = int(match[2])
            else:
                row = line.rstrip().split("\t")
                assert row[2] not in records, (name, "duplicate VCF ID", row[2])
                info = dict(x.split("=", 1) for x in row[7].split(";") if "=" in x)
                records[row[2]] = (row, info, line)
    assert all(chrom in lengths for chrom in autosomes), name
    eligible, exclusions = [], {}
    for ident, (row, info, line) in records.items():
        reason = None
        mate = records.get(info.get("MATEID", ""))
        alt = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", row[4])
        if row[6] != "PASS":
            reason = "non_PASS"
        elif not mate or not alt:
            reason = "unpaired_or_non_BND"
        elif mate[0][6] != "PASS":
            reason = "mate_non_PASS"
        elif mate[1].get("MATEID") != ident:
            raise ValueError((name, ident, "nonreciprocal MATEID"))
        elif row[0] not in autosomes or alt[1] not in autosomes:
            reason = "outside_autosome_QC_scope"
        else:
            assert alt[1] == mate[0][0] and int(alt[2]) == int(mate[0][1])
            assert 1 <= int(row[1]) <= lengths[row[0]]
        if reason:
            exclusions[reason] = exclusions.get(reason, 0) + 1
        else:
            eligible.append(line)
    assert len(eligible) % 2 == 0
    with gzip.open(dest / "paired_PASS.autosomes.vcf.gz", "wt") as handle:
        handle.writelines(headers + eligible)

    coverage = pd.read_csv(str(prefix) + ".cobalt.ratio.tsv.gz", sep="\t")
    input_rows = len(coverage)
    coverage = coverage.loc[coverage.chromosome.isin(autosomes)].copy()
    assert set(coverage.chromosome) == set(autosomes)
    assert not coverage.duplicated(["chromosome", "position"]).any()
    assert coverage.position.mod(1000).eq(1).all()
    assert coverage.groupby("chromosome", sort=False).position.diff().dropna().eq(1000).all()
    assert coverage.position.le(coverage.chromosome.map(lengths)).all()
    valid = (np.isfinite(coverage.tumorGCRatio) & (coverage.tumorGCRatio >= 0)
             & np.isfinite(coverage.referenceGCDiploidRatio)
             & (coverage.referenceGCDiploidRatio > 0)
             & (coverage.referenceReadDepth > 0))
    coverage["ratio"] = np.where(valid,
        coverage.tumorGCRatio / coverage.referenceGCDiploidRatio, np.nan)
    assert np.isfinite(coverage.loc[valid, "ratio"]).all()
    assert coverage.loc[valid].groupby("chromosome").size().reindex(autosomes).gt(0).all()
    coverage["end"] = np.minimum(coverage.position + 999, coverage.chromosome.map(lengths))
    covout = coverage.rename(columns={"chromosome": "seqnames", "position": "start"})
    covout[["seqnames", "start", "end", "ratio"]].to_csv(
        dest / "coverage.autosome_candidate.tsv.gz", sep="\t", index=False, na_rep="NA")

    cn = pd.read_csv(str(prefix) + ".purple.cnv.somatic.tsv", sep="\t")
    cn = cn.loc[cn.chromosome.isin(autosomes)].copy()
    assert set(cn.chromosome) == set(autosomes)
    assert cn.start.ge(1).all() and cn.end.ge(cn.start).all()
    assert cn.end.le(cn.chromosome.map(lengths)).all()
    purity, ploidy = float(sample["purity"]), float(sample["ploidy"])
    assert 0 < purity <= 1 and ploidy > 0
    purity_source = pd.read_csv(str(prefix) + ".purple.purity.tsv", sep="\t").iloc[0]
    assert np.isclose(float(purity_source.purity), purity)
    assert np.isclose(float(purity_source.ploidy), ploidy)
    global_mean = coverage.loc[valid, "ratio"].mean()
    segment_tables = []
    for chrom in autosomes:
        seg = cn.loc[cn.chromosome == chrom].sort_values("start").reset_index(drop=True)
        assert (seg.start.to_numpy()[1:] > seg.end.to_numpy()[:-1]).all()
        cov = coverage.loc[coverage.chromosome == chrom]
        midpoint = (cov.position.to_numpy() + cov.end.to_numpy()) // 2
        ix = np.searchsorted(seg.start.to_numpy(), midpoint, side="right") - 1
        inside = (ix >= 0) & (midpoint <= seg.end.to_numpy()[np.maximum(ix, 0)])
        values = cov.ratio.to_numpy()
        usable = inside & np.isfinite(values)
        count = np.bincount(ix[usable], minlength=len(seg))
        total = np.bincount(ix[usable], weights=values[usable], minlength=len(seg))
        seg["valid_bin_count"] = count
        seg["candidate_ratio_mean"] = np.divide(total, count,
            out=np.full(len(seg), np.nan), where=count > 0)
        seg["candidate_CN"] = ((purity * ploidy + 2 * (1 - purity)) / purity
            * seg.candidate_ratio_mean / global_mean - 2 * (1 - purity) / purity)
        seg["CN_residual"] = seg.candidate_CN - seg.copyNumber
        segment_tables.append(seg)
    comparison = pd.concat(segment_tables, ignore_index=True)
    comparison.to_csv(dest / "PURPLE_segment_coverage_comparison.tsv", sep="\t", index=False)
    fit = comparison.loc[comparison.valid_bin_count >= 5]
    weights = fit.end - fit.start + 1
    summary = dict(sample=name, cohort=sample["cohort"], purity=purity, ploidy=ploidy,
        cobalt_total_rows=input_rows, autosomal_bins=len(coverage), valid_bins=int(valid.sum()),
        valid_fraction=float(valid.mean()), mean_paired_ratio=float(global_mean),
        PASS_paired_autosomal_junctions=len(eligible) // 2, exclusions=exclusions,
        compared_segments=len(fit),
        length_weighted_CN_MAE=float(np.average(np.abs(fit.CN_residual), weights=weights)),
        segment_median_abs_CN_residual=float(np.median(np.abs(fit.CN_residual))),
        state="input_QC_only; not a fitted or validated JaBbA graph")
    (dest / "input_QC.json").write_text(json.dumps(summary, indent=2) + "\n")
    summaries.append(summary)
    print(json.dumps(summary), flush=True)

(out / selected / "INPUT_QC_COMPLETE").write_text("One sample input checks passed; graph QC remains required.\n")
