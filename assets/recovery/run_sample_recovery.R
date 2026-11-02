args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
sample <- args[2]
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({library(JaBbA); library(gGnome); library(data.table)})
spec <- jsonlite::fromJSON(file.path(root, "recovery_spec.json"))
stopifnot(spec$sample == sample, spec$mode %in% c("LONG_FRESH", "NATIVE_WARMSTART"))
setDTthreads(spec$cpus)
set.seed(42)
inputs <- file.path(root, "graph_inputs", sample)
out <- file.path(root, "graphs", sample)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
info <- readRDS(file.path(inputs, "input_summary.rds"))
opt <- list(epgap = 1e-4)

if (spec$mode == "NATIVE_WARMSTART") {
  old <- file.path(root, "previous_native", "graphs", sample)
  stopifnot(file.copy(file.path(old, "karyograph.rds"), file.path(out, "karyograph.rds")))
  kag <- readRDS(file.path(out, "karyograph.rds"))
  incumbent <- readRDS(file.path(old, "jabba.raw.rds"))
  stopifnot(nrow(incumbent$adj) == length(incumbent$segstats))
  # Same native recovery entrypoint and settings as the existing Gadi recovery.
  getFromNamespace("ramip_stub", "JaBbA")(
    kag.file = file.path(out, "karyograph.rds"), out.file = file.path(out, "jabba.raw.rds"),
    mc.cores = 1, max.threads = spec$cpus, mem = spec$solver_memory_GB, tilim = spec$tilim,
    slack.prior = 1 / 100, customparams = TRUE,
    purity.min = info$purity, purity.max = info$purity,
    ploidy.min = info$ploidy, ploidy.max = info$ploidy,
    mipstart = incumbent, epgap = opt$epgap, verbose = TRUE,
    edge.nudge = readRDS(file.path(old, "edge.nudge.rds")),
    ab.force = readRDS(file.path(old, "ab.force.rds")),
    ab.exclude = readRDS(file.path(old, "ab.exclude.rds")),
    loose.penalty.mode = "boolean", dyn.tuning = TRUE,
    debug.ix = seq_len(2L * length(kag$segstats)))
}

parsed <- gGnome::jJ(file.path(inputs, "paired_PASS.vcf"), chr.convert = FALSE)
stopifnot(length(parsed$grl) == info$paired_junctions,
          setequal(seqlevels(parsed$grl), info$chromosomes))
saveRDS(parsed$grl, file.path(out, "native_parsed_junctions.rds"))
model.nseg <- getFromNamespace("seqlevels_to_ncbi", "JaBbA")(readRDS(file.path(inputs, "nseg.rds")))
saveRDS(model.nseg, file.path(out, "normal_segments_NCBI.rds"))
gg <- JaBbA::JaBbA(junctions = file.path(inputs, "paired_PASS.vcf"),
  coverage = file.path(inputs, "cov.rds"), seg = file.path(inputs, "seg.rds"),
  nseg = model.nseg, field = "foreground", purity = info$purity, ploidy = info$ploidy,
  name = sample, outdir = out, max.threads = spec$cpus,
  max.mem = spec$solver_memory_GB, mc.cores = 1, tilim = spec$tilim,
  epgap = opt$epgap, slack.penalty = 100, drop.chr = TRUE,
  overwrite = spec$mode != "NATIVE_WARMSTART")
saveRDS(gg, file.path(out, "fitted.gg.rds"))
