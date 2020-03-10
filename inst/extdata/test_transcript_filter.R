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
    start = c(1, 1, 301),
    width = c(1000, 1000, 200),
    seqnames = "1",
    strand = "+",
    type = c("gene", rep("transcript", 2)),
    ID = c("g1", paste0(rep("t1.", 2), seq(1:2))),
    Parent = c(NA, rep("g1", 2))
)

grng <- g1 %>% as_granges()
txdb <- makeTxDbFromGRanges(grng)

saveDb(txdb, file.path(extdata_dir, "test_transcript_filter.txdb"))

################################################################################
# generate test bigwig file
wig_1_1000 <- tibble(seqnames = 1,
              start = 1:1000,
              width = 1,
              score = 1) %>% as_granges()

wig_301_500 <- tibble(seqnames = 1,
                      start = 1:1000,
                      width = 1,
                      score = c(rep(0, 300), rep(1, 200), rep(0, 500))) %>%
    as_granges()

plot(start(wig_1_1000), wig_1_1000$score)
plot(start(wig_301_500), wig_301_500$score)

generate_bigwig <- function(wig, wig_name) {
    rtracklayer::export.wig(wig, file.path(extdata_dir, wig_name))
    rtracklayer::wigToBigWig(file.path(extdata_dir, wig_name),
                             seqinfo = Seqinfo(seqnames = c("1"),
                                               seqlengths = 1e3,
                                               isCircular = NA, genome = NA))
    file.remove(file.path(extdata_dir, wig_name))
}

generate_bigwig(wig = wig_1_1000, "test_transcript_filter_1_1000.wig")
generate_bigwig(wig = wig_301_500, "test_transcript_filter_301_500.wig")

################################################################################
# Use Gviz to visualize
txdb <- AnnotationDbi::loadDb(file.path(extdata_dir,
                                        "test_transcript_filter.txdb"))

options(ucscChromosomeNames = FALSE)

axistrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(txdb, chromosome = 1, shape = "arrow")
dtrack_1_1000 <- DataTrack(
    range = file.path(extdata_dir, "test_transcript_filter_1_1000.bw"),
    type = "h", name = "1 - 1000", stream = T
)
dtrack_301_500 <- DataTrack(
    range = file.path(extdata_dir, "test_transcript_filter_301_500.bw"),
    type = "h", name = "301 - 500", stream = T)
