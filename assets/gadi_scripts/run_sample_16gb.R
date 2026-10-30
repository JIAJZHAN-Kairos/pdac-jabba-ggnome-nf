args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
sample <- args[2]
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
  library(JaBbA)
  library(gGnome)
  library(data.table)
})
setDTthreads(4)
set.seed(42)
inputs <- file.path(root, "graph_inputs", sample)
out <- file.path(root, "graphs", sample)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
info <- readRDS(file.path(inputs, "input_summary.rds"))
parsed.junctions <- gGnome::jJ(file.path(inputs, "paired_PASS.vcf"), chr.convert = FALSE)
stopifnot(length(parsed.junctions$grl) == info$paired_junctions,
          setequal(seqlevels(parsed.junctions$grl), info$chromosomes))
saveRDS(parsed.junctions$grl, file.path(out, "native_parsed_junctions.rds"))
model.nseg <- getFromNamespace("seqlevels_to_ncbi", "JaBbA")(readRDS(file.path(inputs, "nseg.rds")))
saveRDS(model.nseg, file.path(out, "normal_segments_NCBI.rds"))

# The pinned JaBbA source refers to the CLI opt object in its final gap message.
# Supply that one option without changing the upstream inference implementation.
opt <- list(epgap = 1e-4)
gg <- JaBbA::JaBbA(junctions = file.path(inputs, "paired_PASS.vcf"),
  coverage = file.path(inputs, "cov.rds"), seg = file.path(inputs, "seg.rds"),
  nseg = model.nseg, field = "foreground",
  purity = info$purity, ploidy = info$ploidy, name = sample, outdir = out,
  max.threads = 4, max.mem = 12, mc.cores = 1, tilim = 6000,
  epgap = opt$epgap, slack.penalty = 100, drop.chr = TRUE)
saveRDS(gg, file.path(out, "fitted.gg.rds"))
