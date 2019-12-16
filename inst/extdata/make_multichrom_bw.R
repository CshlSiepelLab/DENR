library(GenomicFeatures)
library(rtracklayer)

current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# generate test bigwig file
wig_multi <- GRanges(seqnames = rep(c(1, "Z"), each = 1e3),
              IRanges(start = rep(1:1e3, 2),
                      width = 1),
              strand = "+",
              score = c(rep(1, 1.5e3), rep(0, 500)))
rtracklayer::export.wig(wig_multi,
                        file.path(current_dir, "test_multichrom.wig"))
rtracklayer::wigToBigWig(file.path(current_dir, "test_multichrom.wig"),
                         seqinfo = Seqinfo(seqnames = c("1", "Z"),
                                           seqlengths = rep(1e5, 2),
                                           isCircular = NA, genome = NA))
unlink(file.path(current_dir, "test_multichrom.wig"),
       recursive = TRUE, force = TRUE)
