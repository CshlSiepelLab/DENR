library(biomaRt)
library(GenomicFeatures)
library(stringr)

extdata_dir <- '/home/yixin/NutstoreFiles/Nutstore/package_development/tuSelecter2/inst/extdata'

# generate txdb for some genes for testing 
dataset_name <- "hsapiens_gene_ensembl"
release <- 96

# test gene list for genes with more than one TSS
gl_moretss <- c("ENSG00000142632", "ENSG00000100029", "ENSG00000162909", 
                "ENSG00000130751", "ENSG00000250021")

gl_onetss <- c("ENSG00000136942")

gl_embedded <- c("ENSG00000197312", "ENSG00000215695")

gene_list <- c(gl_moretss, gl_onetss, gl_embedded) 
# gene_filter <- list(ensembl_gene_id = "ENSG00000142632")

get_txdb <- function(extdata_dir, dataset_name, release, gene_list){
    # get host name based on Ensembl release version
    df <- biomaRt::listEnsemblArchives()
    host_name <- df[df$version == release, "url"]
    
    txdb_name <-
        file.path(extdata_dir, str_c(dataset_name, release, "txdb_test", sep = "_"))
    
    if (file.exists(txdb_name)){
        txdb <- AnnotationDbi::loadDb(txdb_name)
    } else{
        txdb <- GenomicFeatures::makeTxDbFromBiomart(dataset = dataset_name, 
                                                 host = host_name,
                                                 filter = list(ensembl_gene_id = gene_list))
        AnnotationDbi::saveDb(txdb, txdb_name)
    }
    return(txdb)
}

get_txdb(extdata_dir, dataset_name, release, gene_list)
