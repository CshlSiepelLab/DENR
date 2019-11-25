#' Retrieve ensembl biomart
#'
#' Convenience function for retrieving information about available ensembl marts
#' that can be turned into TxDbs.
#'
#' @param species the species to download the TxDb for. Use the first letter
#' of the genus and the full species as one string (e.g. homo sapiens ->
#' hsapiens). If not specified will return all species.
#' @param ensembl_version The Ensembl release version. Defaults to current
#' release.
#' @rdname view_available_txdb
#' @return A data.frame of available datasets
#' @export
view_available_txdb <- function(species = NULL, ensembl_version = NULL) {
  # List ensembl versions
  ensembl_version_df <- biomaRt::listEnsemblArchives()
  # Default url
  ensembl_url <- "www.ensembl.org"
  # Look for version specified
  if (!is.null(ensembl_version)) {
    ensembl_url <- ensembl_version_df[ensembl_version_df$version ==
                                       ensembl_version, ]$url
    if (length(ensembl_url) == 0) {
      stop("Invalid ensembl version specified.")
    }
    mart <- biomaRt::useMart(biomart = "ensembl", host = ensembl_url)
  } else {
    mart <- biomaRt::useMart(biomart = "ensembl")
  }

  species_df <- biomaRt::listDatasets(mart)

  # Filter returned mart by species if specified
  if (!is.null(species)) {
    query_string <- paste0(species, "_gene_ensembl")
    species_df <- species_df[species_df$dataset == query_string, ]
    if (nrow(species_df) == 0) {
      stop("Invalid species specified.")
    }
  }
  species_df$host <- ensembl_url
  return(species_df)
}

#' Download transcript model from biomart
#'
#' This function creates a TxDb object from the biomaRt ensembl database,
#' converts it to a TxDb, and returns it (WARNING: this may take a few
#' minutes.) To view available species for a given ensembl version run
#' view_available_txdb().
#'
#' @param species the species to download the TxDb for. Use the first letter
#' of the genus and the full species as one string (e.g. homo sapiens ->
#' hsapiens).
#' @param ensembl_version The Ensembl release version. Defaults to current
#' release.
#' @return TxDb object
#' @export
download_txdb <- function(species, ensembl_version = NULL, cache = TRUE) {
  if (length(species) > 1) {
    stop("Only one species may be specified.")
  }
  # Retrieve relevant mart
  tx_df <- view_available_txdb(species = species,
                              ensembl_version = ensembl_version)
  # Convert it to a txdb
  txdb <- GenomicFeatures::makeTxDbFromBiomart(dataset = tx_df$dataset,
                                              biomart = "ensembl",
                                              host = tx_df$host)
  # Cache the downloaded txdb
  if (isTRUE(cache)) {
    cache_set_dir()
    save_txdb(txdb = txdb)
  }
  return(txdb)
}
