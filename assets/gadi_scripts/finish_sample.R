args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
sample <- args[2]
.libPaths(c(Sys.getenv("GUTILS_COMPAT_LIB"), Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({library(gGnome); library(data.table)})
setDTthreads(4)
set.seed(42)
inputs <- file.path(root, "graph_inputs", sample)
out <- file.path(root, "graphs", sample)
info <- readRDS(file.path(inputs, "input_summary.rds"))
gg <- readRDS(file.path(out, "fitted.gg.rds"))
opt <- list(epgap = 1e-4)
# Native RDS retains NULL metadata; normalize only missing TSV list cells.
write_tsv <- function(value, path, ...) {
  value <- copy(as.data.table(value))
  for (column in names(value)) if (is.list(value[[column]])) {
    value[[column]] <- lapply(value[[column]], function(cell)
      if (is.null(cell)) NA_character_ else cell)
  }
  data.table::fwrite(value, path, ...)
}
nodes <- gg$nodes$dt
edges <- gg$edges$dt
write_tsv(nodes, file.path(out, "nodes.tsv"), sep = "\t")
write_tsv(edges, file.path(out, "edges.tsv"), sep = "\t")
balance <- getFromNamespace("bcheck", "gGnome")(gg)
node.index <- match(balance$node.id, nodes$node.id)
chromosome.length <- seqlengths(gg$nodes$gr)[as.character(nodes$seqnames[node.index])]
stopifnot(all(is.finite(chromosome.length)))
# JaBbA leaves physical chromosome ends unconstrained (see native telo.in/out).
# Preserve raw bcheck residuals; enforce zero residual on every finite internal side.
balance[, reference_boundary_left := nodes$start[node.index] == 1]
balance[, reference_boundary_right := nodes$end[node.index] == chromosome.length]
write_tsv(balance, file.path(out, "junction_balance.tsv"), sep = "\t")
known.cn <- nodes$cn[is.finite(nodes$cn)]
known.edges <- edges$cn[is.finite(edges$cn)]
stopifnot(length(known.cn) > 0, all(known.cn >= -1e-6),
          all(abs(known.cn - round(known.cn)) <= 1e-6),
          all(known.edges >= -1e-6), all(abs(known.edges - round(known.edges)) <= 1e-6))
finite.balance <- is.finite(balance$left) & is.finite(balance$right)
internal.left <- is.finite(balance$left) & !balance$reference_boundary_left
internal.right <- is.finite(balance$right) & !balance$reference_boundary_right
stopifnot(any(internal.left), any(internal.right),
          all(abs(balance$left[internal.left]) <= 1e-6),
          all(abs(balance$right[internal.right]) <= 1e-6))
report.files <- unique(normalizePath(list.files(out, "^opt.report.rds$", recursive = TRUE,
                                               full.names = TRUE)))
stopifnot(length(report.files) > 0)
report <- rbindlist(lapply(report.files, function(f) {
  value <- as.data.table(readRDS(f)); value[, source := f]; value
}), fill = TRUE)
write_tsv(report, file.path(out, "optimization_report.tsv"), sep = "\t")
stopifnot(all(is.finite(report$epgap)), all(report$status %in% c(101, 102)),
          all(report$epgap <= opt$epgap * 1.001))

event.gg <- gGnome::events(readRDS(file.path(out, "fitted.gg.rds")))
saveRDS(event.gg, file.path(out, "events.gg.rds"))
events <- event.gg$meta$events
if (is.null(events) || nrow(events) == 0)
  events <- data.table(type = character(), footprint = character())
saveRDS(events, file.path(out, "native_events.rds"))
write_tsv(events, file.path(out, "native_events.tsv"), sep = "\t")
write_tsv(events[, .N, by = type], file.path(out, "event_type_counts.tsv"), sep = "\t")
stats <- list(sample = sample, node_count = nrow(nodes), edge_count = nrow(edges),
  unknown_node_CN = sum(!is.finite(nodes$cn)), unknown_edge_CN = sum(!is.finite(edges$cn)),
  checked_balance_nodes = sum(finite.balance), unknown_balance_nodes = sum(!finite.balance),
  max_abs_balance_residual = max(abs(c(balance$left[finite.balance], balance$right[finite.balance]))),
  checked_internal_balance_sides = sum(internal.left) + sum(internal.right),
  reference_chromosome_boundary_sides = sum(balance$reference_boundary_left) + sum(balance$reference_boundary_right),
  max_abs_internal_balance_residual = max(abs(c(balance$left[internal.left], balance$right[internal.right]))),
  optimization_rows = nrow(report),
  optimization_unknown_gaps = sum(!is.finite(report$epgap)),
  max_known_optimization_gap = if (any(is.finite(report$epgap))) max(report$epgap[is.finite(report$epgap)]) else NA_real_,
  events = nrow(events), chromosomes = info$chromosomes,
  status = "fit_and_native_classification_finished; requires pilot review before cohort expansion")
jsonlite::write_json(stats, file.path(out, "sample_QC.json"), auto_unbox = TRUE,
                     pretty = TRUE, na = "null", digits = NA)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(out, "SAMPLE_COMPLETE"))
