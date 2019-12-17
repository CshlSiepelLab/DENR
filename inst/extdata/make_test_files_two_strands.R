library(tibble)
library(magrittr)
library(dplyr)
library(plyranges)
library(GenomicFeatures)
library(Gviz)

current_file <- rstudioapi::getSourceEditorContext()$path
inst_dir <- dirname(dirname(current_file))
extdata_dir <- file.path(inst_dir, "extdata")

################################################################################
# generate test txdb
g1 <- tibble(
    start = c(5200, 5200, 7000),
    width = c(5000, 5000, 3200),
    seqnames = "1",
    strand = "+",
    type = c("gene", rep("transcript", 2)),
    ID = c("g1", paste0(rep("t1.", 2), seq(1:2))),
    Parent = c(NA, rep("g1", 2))
)

g2 <- tibble(
    end = c(6300, 6300, 5100),
    width = c(5000, 5000, 3000),
    seqnames = "1",
    strand = "-",
    type = c("gene", rep("transcript", 2)),
    ID = c("g2", paste0(rep("t2.", 2), seq(1:2))),
    Parent = c(NA, rep("g2", 2))
)

grng <- dplyr::bind_rows(g1, g2) %>% as_granges()
txdb <- makeTxDbFromGRanges(grng)

saveDb(txdb, file.path(extdata_dir, "test_double_strand.txdb"))

################################################################################
# generate test bigwig file
set.seed(2019)

tx <- transcripts(txdb)
tx_subset_plus <- tx[tx$tx_name %in% c("t1.1")]
tx_subset_minus <- tx[tx$tx_name %in% c("t2.2")]

get_wig_template <- function(strand) {
    wig <- tibble(seqnames = 1,
                  start = 1:13000,
                  width = 1,
                  strand = strand,
                  # set the base expression
                  base = rbinom(n = 13000, size = 1, prob = 0.3),
                  pausing_peak = 0,
                  tx_sum = 0)
    return(wig)
}

wig_plus <- get_wig_template("+")
wig_minus <- get_wig_template("-")

# generate read count near TSS region
add_peak_count <- function(wig, start_pos, end_pos) {
    wig[(wig$start >= start_pos & wig$start < end_pos), ] %<>%
        # set the height of the peak
        mutate(pausing_peak = abs(rnorm(200, mean = 200, sd = 20)))
    return(wig)
}

# add pausing peak read counts
get_peak_wig <- function(tx_subset, wig, strand) {
    if (strand == "+") {
        peak_se <- list(peak_starts = start(tx_subset) - 100,
                        peak_ends = start(tx_subset) + 100)
    } else if (strand == "-") {
        peak_se <- list(peak_starts = end(tx_subset) - 100,
                        peak_ends = end(tx_subset) + 100)
    }
    for (i in seq_len(length(peak_se$peak_starts))) {
        wig <- add_peak_count(wig, peak_se$peak_starts[i], peak_se$peak_ends[i])
    }
    return(wig)
}

wig_plus <- get_peak_wig(tx_subset_plus, wig_plus, "+")
wig_minus <- get_peak_wig(tx_subset_minus, wig_minus, "-")

# add transcript expression
get_tx_wig <- function(tx_subset, wig, bin_height_vec) {
    tx_se <- list(tx_starts = start(tx_subset),
                  tx_ends = end(tx_subset) + 1,
                  bin_height = bin_height_vec)
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
    return(wig)
}

#
wig_plus <- get_tx_wig(tx_subset_plus, wig_plus, c(150))
wig_minus <- get_tx_wig(tx_subset_minus, wig_minus, c(200))

plot(start(wig_plus), wig_plus$score)
plot(start(wig_minus), wig_minus$score)

get_bigwig <- function(wig, name) {
    rtracklayer::export.wig(wig, file.path(extdata_dir, name))
    rtracklayer::wigToBigWig(file.path(extdata_dir, name),
                             seqinfo = Seqinfo(seqnames = c("1"),
                                               seqlengths = 1e5,
                                               isCircular = NA, genome = NA))
    file.remove(file.path(extdata_dir, name))
}

get_bigwig(wig_plus, "test_double_strand_plus.wig")
get_bigwig(wig_minus, "test_double_strand_minus.wig")

################################################################################
# Use Gviz to visualize
txdb <- AnnotationDbi::loadDb(file.path(extdata_dir, "test_double_strand.txdb"))

options(ucscChromosomeNames = FALSE)

proseq_ylim <- c(0, 600)

axistrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(txdb, chromosome = 1, shape = "arrow")
dtrack_plus <- DataTrack(range = file.path(extdata_dir,
                                           "test_double_strand_plus.bw"),
                    type = "h", name = "PROseq_plus",
                    stream = T, ylim = proseq_ylim)
dtrack_minus <- DataTrack(range = file.path(extdata_dir,
                                            "test_double_strand_minus.bw"),
                         type = "h", name = "PROseq_minus",
                         stream = T, col = "red", ylim = rev(proseq_ylim))

plotTracks(c(axistrack, grtrack, dtrack_plus, dtrack_minus),
           from = 1, to = 13000,
           showId = T, transcriptAnnotation = "symbol")
