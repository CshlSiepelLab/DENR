library(tibble)
library(magrittr)
library(dplyr)
library(plyranges)
library(GenomicFeatures)
library(Gviz)

current_dir <- rstudioapi::getSourceEditorContext()$path
inst_dir <- dirname(dirname(current_dir))
extdata_dir <- file.path(inst_dir, "extdata")

################################################################################
# generate test txdb
g1 <- tibble(
    start = c(1000, 1000, 1020, 1400, 2500),
    width = c(5000, 5000, 4980, 4600, 3500),
    seqnames = "1",
    strand = "+",
    type = c("gene", rep("transcript", 4)),
    ID = c("g1", paste0(rep("t1.", 4), seq(1:4))),
    Parent = c(NA, rep("g1", 4))
)

g2 <- tibble(
    start = c(7000, 7000, 8000),
    width = c(5000, 5000, 3000),
    seqnames = "1",
    strand = "+",
    type = c("gene", rep("transcript", 2)),
    ID = c("g2", paste0(rep("t2.", 2), seq(1:2))),
    Parent = c(NA, rep("g2", 2))
)

grng <- dplyr::bind_rows(g1, g2) %>% as_granges()
txdb <- makeTxDbFromGRanges(grng)

saveDb(txdb, file.path(extdata_dir, "txdb_test_single_strand"))

################################################################################
# generate test bigwig file
set.seed(2019)

tx <- transcripts(txdb)
tx_subset_plus <- tx[tx$tx_name %in% c("t1.1", "t1.4", "t2.1")]

wig <- tibble(seqnames = 1,
    start = 1:13000,
    width = 1,
    strand = "+",
    base = rbinom(n = 13000, size = 1, prob = 0.3),
    pausing_peak = 0,
    tx_sum = 0)

# generate read count near TSS region
add_peak_count <- function(wig, start_pos, end_pos) {
    wig[(wig$start >= start_pos & wig$start < end_pos), ] %<>%
        mutate(pausing_peak = abs(rnorm(200, mean = 200, sd = 20)))
    wig
}

# add pausing peak read counts 
peak_se <- list(peak_starts = start(tx_subset_plus) - 100,
              peak_ends = start(tx_subset_plus) + 100)

for (i in seq_len(length(peak_se$peak_starts))) {
    wig <- add_peak_count(wig, peak_se$peak_starts[i], peak_se$peak_ends[i])
}

# add transcript expression
tx_se <- list(tx_starts = start(tx_subset_plus),
              tx_ends = end(tx_subset_plus) + 1,
              bin_height = c(50, 150, 200))

for (i in seq_len(length(tx_se$tx_starts))) {
    tx_start <- tx_se$tx_starts[i]
    tx_end <- tx_se$tx_ends[i]
    wig[(wig$start >= tx_start & wig$start < tx_end), ] %<>%
        mutate(tx_sum = abs(rnorm(tx_end - tx_start,
                              mean = tx_se$bin_height[i],
                              sd = tx_se$bin_height[i] / 2)))
    }

wig <- wig %>%
    mutate(score = Reduce("+", .[c("base", "pausing_peak", "tx_sum")])) %>%
    dplyr::select(seqnames, start, width, score) %>%
    as_granges()
    
plot(start(wig), wig$score)

rtracklayer::export.wig(wig, file.path(extdata_dir, "test_single_strand.wig"))
rtracklayer::wigToBigWig(file.path(extdata_dir, "test_single_strand.wig"),
                         seqinfo = Seqinfo(seqnames = c("1"),
                                           seqlengths = 1e5,
                                           isCircular = NA, genome = NA))
file.remove(file.path(extdata_dir, "test_single_strand.wig"))

################################################################################
# Use Gviz to visualize
txdb <- AnnotationDbi::loadDb(file.path(extdata_dir, "txdb_test_single_strand"))

options(ucscChromosomeNames = FALSE)

axistrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(txdb, chromosome = 1, shape = "arrow")
dtrack <- DataTrack(range = file.path(extdata_dir, "test_single_strand.bw"),
                    type = "h", name = "PROseq_plus", stream = T)

plotTracks(c(axistrack, grtrack, dtrack), from = 1, to = 13000,
           showId = T, transcriptAnnotation = "symbol")
