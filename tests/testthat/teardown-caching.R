db_string <- "homo_sapiens-grch38.p12-ensembl_genes_96_test_Fh4Jf93l1.txdb"
cache_loc <- cache_get_dir()

tmp_db <- file.path(cache_loc, db_string)

if (file.exists(tmp_db)) {
  file.remove(tmp_db)
}
