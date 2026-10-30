args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
sample <- args[2]
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(VariantAnnotation)
})
source(file.path(root, "scripts", "cbs_coverage.R"))
dest <- file.path(root, "graph_inputs", sample)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
vcf <- readVcf(file.path(root, "inputs", paste0(sample, ".esvee.somatic.vcf.gz")), genome = "GRCh38")
sl <- seqlengths(rowRanges(vcf))
canonical <- intersect(c(paste0("chr", 1:22), "chrX", "chrY"), names(sl))
cov <- fread(cmd = paste("gzip -dc", shQuote(file.path(root, "inputs",
                       paste0(sample, ".cobalt.ratio.tsv.gz")))))
cov <- cov[chromosome %in% canonical]
cov[, valid.normal := is.finite(referenceGCDiploidRatio) &
        referenceGCDiploidRatio > 0 & referenceReadDepth > 0]
cov[, block := floor((position - 1) / 5e6)]
support <- cov[, .(valid.fraction = mean(valid.normal)), by = .(chromosome, block)]
supported <- unique(support[valid.fraction > 0.90, chromosome])
baseline <- cov[valid.normal == TRUE,
               .(median.normal.ratio = median(referenceGCDiploidRatio)), by = chromosome]
baseline[, cn := round(2 * median.normal.ratio)]
baseline[, included := chromosome %in% supported & cn > 0]
fwrite(baseline, file.path(dest, "normal_chromosome_baselines.tsv"), sep = "\t")
keep.chr <- canonical[canonical %in% baseline[included == TRUE, chromosome]]
stopifnot(all(paste0("chr", 1:22) %in% keep.chr), "chrX" %in% keep.chr)
purity <- fread(file.path(root, "inputs", paste0(sample, ".purple.purity.tsv")))
stopifnot(nrow(purity) == 1)
normal.x <- baseline[chromosome == "chrX", cn]
expected.x <- switch(purity$gender, MALE = 1, FEMALE = 2, NA_real_)
stopifnot(is.finite(expected.x), identical(as.numeric(normal.x), expected.x))
cov <- cov[chromosome %in% keep.chr]
stopifnot(!anyDuplicated(cov[, .(chromosome, position)]))
end <- pmin(cov$position + 999L, unname(sl[cov$chromosome]))
si <- Seqinfo(keep.chr, seqlengths = sl[keep.chr], genome = "GRCh38")
tumor <- GRanges(cov$chromosome, IRanges(cov$position, end), seqinfo = si)
normal <- tumor
tumor$cobalt_gcratio <- cov$tumorGCRatio
tumor$cobalt_gcratio[!is.finite(tumor$cobalt_gcratio) | tumor$cobalt_gcratio < 0] <- NA_real_
normal$cobalt_gcratio <- cov$referenceGCDiploidRatio
normal$cobalt_gcratio[!cov$valid.normal] <- NA_real_
ncn <- baseline$cn[match(cov$chromosome, baseline$chromosome)]
paired <- prepare_cbs_coverage(tumor, normal, ncn, input_field = "cobalt_gcratio")
nseg <- GRanges(keep.chr, IRanges(1, sl[keep.chr]), seqinfo = si)
nseg$cn <- baseline$cn[match(keep.chr, baseline$chromosome)]
nseg$mean <- baseline$median.normal.ratio[match(keep.chr, baseline$chromosome)]
purple <- fread(file.path(root, "inputs", paste0(sample, ".purple.cnv.somatic.tsv")))
purple <- purple[chromosome %in% keep.chr]
seg <- GRanges(purple$chromosome, IRanges(purple$start, purple$end), seqinfo = si)
seg$purple.copyNumber <- purple$copyNumber
seg$ncn <- cbs_segment_ncn(seg, paired$cov)
bundle <- finalize_cbs_bundle(paired$cov, seg, nseg, paired$mode,
                             paired$input_field, paired$foreground_scale)
metadata(bundle$cov)$input_source <- list(tool = "COBALT", sample = sample,
  normal_cn = "rounded 2x calibrated chromosome median; chrX verified against PURPLE gender",
  normal_support_rule = ">90% valid normal bins in at least one 5Mb block, native nf-gos CBS rule",
  helper_revision = "89524ac3aacd2fc321121e2898a91bda58d5b691")
saveRDS(bundle$cov, file.path(dest, "cov.rds"), compress = FALSE)
saveRDS(bundle$seg, file.path(dest, "seg.rds"))
saveRDS(bundle$nseg, file.path(dest, "nseg.rds"))

ids <- names(rowRanges(vcf))
mates <- as.character(info(vcf)$MATEID)
mate.ix <- match(mates, ids)
alt <- as.character(unlist(alt(vcf)))
pass <- as.character(fixed(vcf)$FILTER) == "PASS"
in.scope <- as.character(seqnames(rowRanges(vcf))) %in% keep.chr
eligible <- which(pass & grepl("[\\[\\]]", alt, perl = TRUE) & !is.na(mate.ix))
eligible <- eligible[pass[mate.ix[eligible]] & in.scope[eligible] & in.scope[mate.ix[eligible]]]
stopifnot(length(eligible) > 0, all(mates[mate.ix[eligible]] == ids[eligible]),
          length(eligible) %% 2 == 0)
paired.vcf <- vcf[eligible]
seqlevels(paired.vcf, pruning.mode = "coarse") <- keep.chr
stopifnot(identical(seqlevels(paired.vcf), keep.chr))
paired.header <- header(paired.vcf)
meta(paired.header)$contig <- meta(paired.header)$contig[keep.chr, , drop = FALSE]
header(paired.vcf) <- paired.header
# ESVEE uses BND ALT syntax even when its event annotation is DEL/DUP/INV.
# Dispatch the derived file through gGnome's BND reader; preserve original types below.
info(paired.vcf)$SVTYPE <- rep("BND", length(eligible))
writeVcf(paired.vcf, file.path(dest, "paired_PASS.vcf"))
stopifnot(setequal(rownames(meta(scanVcfHeader(file.path(dest, "paired_PASS.vcf")))$contig),
                   keep.chr))
fwrite(data.table(id = ids, mate = mates, filter = fixed(vcf)$FILTER,
                 original_svtype = as.character(info(vcf)$SVTYPE),
                 chromosome = as.character(seqnames(rowRanges(vcf))),
                 included = seq_along(ids) %in% eligible),
       file.path(dest, "junction_inclusion_audit.tsv"), sep = "\t")
summary <- list(sample = sample, chromosomes = keep.chr, purity = purity$purity,
  ploidy = purity$ploidy, normal_cn = baseline, paired_junctions = length(eligible) / 2,
  coverage_info = cbs_coverage_info(bundle$cov))
saveRDS(summary, file.path(dest, "input_summary.rds"))
writeLines(capture.output(str(summary)), file.path(dest, "input_summary.txt"))
