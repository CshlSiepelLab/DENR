#' Get transcript model from Ensembl
#'
#' This function creates a txDb object if it doesn't exist locally.
#'
#' @param txDb_dir A directory to save and retrieve txDb object.
#' @param dataset_name A dataset name used by Ensembl, like
#'   "hsapiens_gene_ensembl".
#' @param release The Ensembl release version.
#' @param host_name The host name used to connect to Ensembl.
#' @return A txDB object.
#' @export

get_tx_model <- function(txDb_dir, dataset_name, release){
    # get host name based on Ensembl release version
    df <- biomaRt::listEnsemblArchives()
    host_name <- df[df$version == release, "url"]

    txDb_name <-
        file.path(txDb_dir, str_c(dataset_name, release, "txDb", sep = "_"))

    if (file.exists(txDb_name)){
        txDb <- AnnotationDbi::loadDb(txDb_name)
    } else{
        txDb <-
            GenomicFeatures::makeTxDbFromBiomart(dataset = dataset_name, host = host_name)
        AnnotationDbi::saveDb(txDb, txDb_name)
    }
    return(txDb)
}
