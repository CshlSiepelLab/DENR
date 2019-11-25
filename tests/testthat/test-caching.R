context("Test txdb caching")

# Load in test txdb
txdb_path = system.file("extdata", "hsapiens_gene_ensembl_96_txdb_test",
                        package = "tuSelecter2")
txdb = AnnotationDbi::loadDb(file = txdb_path)
db_string = 'homo_sapiens-grch38.p12-ensembl_genes_96_test_Fh4Jf93l1.txdb'

test_that("Cache location can be set and retrieved", {
    cache_set_dir()
    expect_silent({cache_loc <- cache_get_dir()})
    expect_true(file.exists(cache_loc))
})

test_that("TxDb can be cached and retrieved",{
    cache_loc <- cache_get_dir()
    file_path = paste0(file.path(cache_loc,db_string))
    expect_message(save_txdb(txdb = txdb,
                             annotation_version = 'ensembl_genes_96_test_Fh4Jf93l1'),
                   regexp = paste0("Saving db to ",file_path))
    expect_s4_class({txdb2 <<- load_txdb(db_string)},
                    class = "TxDb")
    AnnotationDbi::dbFileDisconnect(dbconn =  AnnotationDbi::dbconn(txdb2))
    # Force release of txdb file
    invisible(gc())
})

test_that("TxDb listed correctly",{
    expect_s3_class({txdb_list <- list_cached_txdb()},
                    "data.frame")
    expect_true(db_string %in% txdb_list$file_name)
})
