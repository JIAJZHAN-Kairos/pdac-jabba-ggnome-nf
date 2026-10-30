args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
sample <- args[2]
.libPaths(c(Sys.getenv("GUTILS_COMPAT_LIB"), Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({library(gGnome); library(gTrack); library(gUtils); library(data.table)})
out <- file.path(root, "graphs", sample)
gg <- readRDS(file.path(out, "events.gg.rds"))
nd <- gg$nodes$dt
ed <- gg$edges$dt
loose <- c(nd$loose.cn.left, nd$loose.cn.right)
report <- fread(file.path(out, "optimization_report.tsv"))
review <- list(sample = sample,
  unknown_node_bp = sum(nd$end[!is.finite(nd$cn)] - nd$start[!is.finite(nd$cn)] + 1),
  total_node_bp = sum(nd$end - nd$start + 1),
  nonzero_loose_end_sides = sum(loose > 0, na.rm = TRUE),
  total_loose_end_CN = sum(loose, na.rm = TRUE),
  ALT_edges = sum(ed$type == "ALT"),
  positive_CN_ALT_edges = sum(ed$type == "ALT" & ed$cn > 0, na.rm = TRUE),
  zero_CN_ALT_edges = sum(ed$type == "ALT" & ed$cn == 0, na.rm = TRUE),
  unknown_CN_ALT_edges = sum(ed$type == "ALT" & !is.finite(ed$cn)),
  optimization_status_counts = as.list(table(report$status)),
  convergence_category_counts = as.list(table(report$converge)))
jsonlite::write_json(review, file.path(out, "graph_review_QC.json"), auto_unbox = TRUE, pretty = TRUE)
ev <- as.data.table(readRDS(file.path(out, "native_events.rds")))
selected <- ev[type %in% c("bfb", "tyfonas", "rigma", "pyrgo", "chromoplexy", "dm", "cpxdm", "chromothripsis")]
windows <- list()
if (nrow(selected)) for (i in seq_len(nrow(selected))) {
  fp <- gUtils::parse.gr(unlist(strsplit(selected$footprint[i], "[;,]")),
                        seqlengths = GenomeInfoDb::seqlengths(gg$nodes$gr))
  windows[[paste0(selected$type[i], "_", selected$ev.id[i])]] <- fp
}
top <- gg$nodes$gr[which.max(gg$nodes$gr$cn)]
windows[["highest_CN_region"]] <- top
dir.create(file.path(out, "native_plots"), showWarnings = FALSE)
plot.windows <- list()
for (label in names(windows)) {
  fp <- windows[[label]]
  reference.length <- GenomeInfoDb::seqlengths(gg$nodes$gr)[as.character(seqnames(fp))]
  window <- fp
  ranges(window) <- IRanges::IRanges(start = pmax(1L, start(fp) - 100000L),
                                    end = pmin(reference.length, end(fp) + 100000L))
  stopifnot(all(is.finite(reference.length)), all(start(window) >= 1),
            all(end(window) <= reference.length),
            identical(unname(GenomeInfoDb::seqlengths(window)[as.character(seqnames(window))]),
                      unname(reference.length)))
  plot.windows[[label]] <- data.table(plot = label, seqnames = as.character(seqnames(window)),
                                     start = start(window), end = end(window),
                                     reference_seqlength = unname(reference.length))
  track <- gg$gtrack(y.field = "cn", name = paste(sample, label))
  png(file.path(out, "native_plots", paste0(label, ".png")), width = 1800, height = 900, res = 150)
  plot(track, window)
  dev.off()
  pdf(file.path(out, "native_plots", paste0(label, ".pdf")), width = 12, height = 6)
  plot(track, window)
  dev.off()
}
fwrite(rbindlist(plot.windows), file.path(out, "native_plots", "plotted_windows.tsv"), sep = "\t")
writeLines(names(windows), file.path(out, "native_plots", "plotted_events.txt"))
