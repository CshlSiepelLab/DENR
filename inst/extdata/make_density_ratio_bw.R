library(tibble)
library(magrittr)
library(dplyr)
library(plyranges)
library(GenomicFeatures)
library(Gviz)

as_granges <- function(x) {
  return(GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
}

current_dir <- rstudioapi::getSourceEditorContext()$path
inst_dir <- dirname(dirname(current_dir))
extdata_dir <- file.path(inst_dir, "extdata")

################################################################################
# generate test txdb
g1 <- tibble(
  start = c(1e4, 1.5e4, 3e4),
  end = start + 1e4 - 1,
  seqnames = "1",
  strand = c("+", "+", "-"),
  type = c(rep("transcript", 3)),
  ID = paste0(rep("t1.", 3), seq_len(3)),
)

grng <- g1 %>% as_granges()
txdb <- makeTxDbFromGRanges(grng)

saveDb(txdb, file.path(extdata_dir, "test_upr.txdb"))

################################################################################
# generate test bigwig file
wig_plus <- tibble(seqnames = 1,
                   start = c(start(grng)[1], end(grng)[3]),
                   end = c(end(grng)[2], end(grng)[3] + 1e2),
                   strand = "+",
                   score = c(1, 1)) %>% as_granges()
seqlengths(wig_plus) <- c(1e5)

wig_minus <- tibble(seqnames = 1,
                   start = c(start(grng)[1] - 1e2, start(grng)[3]),
                   end = c(start(grng)[1], end(grng)[3]),
                   strand = "-",
                   score = c(-2, -2)) %>% as_granges()
seqlengths(wig_minus) <- c(1e5)

generate_bigwig <- function(wig, wig_name) {
  rtracklayer::export.bw(object = wig, con = file.path(extdata_dir, wig_name))
}

generate_bigwig(wig = wig_plus, "test_upr_plus.bw")
generate_bigwig(wig = wig_minus, "test_upr_minus.bw")
