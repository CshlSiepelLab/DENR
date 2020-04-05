#' @title Split path
#'
#' @description Splits path into vector of path elements
#'
#' @param path file path
#'
#' @return A vector of path elements.
#'
split_path <- function(path) {
  # https://stackoverflow.com/questions/29214932/
  # split-a-file-path-into-folder-names-vector
  if (dirname(path) %in% c(".", path)) {
    return(basename(path))
  }
  return(c(split_path(dirname(path)), basename(path)))
}

#' @title Create cache defaults
#'
#' @description Creates paths for default cache_dir and cache_path.txt
#' file
#' @return A list with two elements \code{dir} and \code{path_file}
#'
cache_defaults <- function() {
  # Set default cache dir location
  default_cache_dir <- rappdirs::user_data_dir(appname = "tuselecter")
  default_cache_dir <- do.call(
    file.path,
    as.list(split_path(path = default_cache_dir))
  )

  # File that stores location of cache path
  cache_path_file <- file.path(default_cache_dir, "cache_path.txt")
  return(list(dir = default_cache_dir, path_file = cache_path_file))
}

#' @title Create local cache
#'
#' @description
#' Creates local data file cache. If the directory does not exist, it
#' will be created recursively. If no custom path is set, the
#' default user data directory for the package will be used. See
#' \code{\link[rappdirs]{user_data_dir}} for details.
#'
#' @param cache_dir Character path to new local files path. If null,
#' path will be reset to default user data directory location.
#' @param force boolean that forces overwriting of cache location (default:
#' FALSE)
#'
#' @rdname cache_create_dir
#' @seealso \code{\link{cache_set_dir}} \link{cache_get_dir}
#'
#' @examples
#' \dontrun{
#' # set a different cache path
#' set_cache_path("z:/transcript_db/tuselecter")
#' }
#'
#' @export
cache_create_dir <- function(cache_dir = NULL, force = FALSE) {
  defaults <- cache_defaults()
  # Set default cache dir location
  default_cache_dir <- defaults$dir
  # File that stores location of cache path
  cache_path_file <- defaults$path_file

  # Create cache_dir to serve as true cache dir or place to store
  # cache path file to record cache location between sessions
  if (!dir.exists(default_cache_dir)) {
    dir.create(default_cache_dir, recursive = TRUE)
  }

  # If cache_dir is unset use the default directory as the cache_dir
  if (is.null(cache_dir)) {
    cache_dir <- default_cache_dir
  }

  # Reparse cache_dir
  cache_dir <- do.call(file.path, as.list(split_path(path = cache_dir)))

  # Check to see if cache path file already exists, if so, require
  # user input to overwrite/remove if set is different from recorded
  # location
  if (file.exists(cache_path_file)) {
    cache_con <- file(cache_path_file)
    cache_path_string <- readLines(cache_con,
      n = 1,
      warn = FALSE
    )
    close(cache_con)
    counter <- 0
    if (isFALSE(force) && cache_dir != cache_path_string) {
      usr_text <- paste0(
        "Specified cache directory does not match previous",
        " path directory:
",
        cache_path_string, " -> ", cache_dir
      )
      cat(usr_text)
      overwrite_cache_string <- tolower(
        readline(prompt = "Continue setting cache path [y/n]:")
      )
      while (counter < 3 && !overwrite_cache_string %in% c("y", "n")) {
        cat("Invalid input. Input must be \'y\' or \'n\'")
        overwrite_cache_string <- tolower(
          readline(prompt = "Continue setting cache path [y/n]:")
        )
        counter <- counter + 1
      }
      # If not overwriting cache location set cache_dir to previous location
      if (overwrite_cache_string == "n") {
        cache_dir <- cache_path_string
      }
    }
  }

  # Create cache_dir if does not exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Record cache dir path
  write(x = normalizePath(cache_dir, winslash = "/"), file = cache_path_file)
}


#' @title Set local cache path
#'
#' @description
#' Downloaded databases are cached locally. This function
#' sets the R environmental variable \code{CACHE_DIR_TUSEL}
#' with the location of the cache dir. If \code{cache_dir} is
#' specified, it sets the cache directory to that location,
#' creating the specified directory if one does not exit.
#' Otherwise it just uses the previously specified cache dir.
#' If no cache directory has been created  and \code{cache_dir}
#' is unspecified it creates one in the default location.
#'
#' @param cache_dir Character path to new local files path. If null,
#' path will be reset to default user data directory location.
#'
#' @seealso \code{\link{cache_get_dir}}
#' @rdname cache_set_dir
#' @export
cache_set_dir <- function(cache_dir = NULL) {
  defaults <- cache_defaults()
  # File that stores location of cache path
  cache_path_file <- defaults$path_file

  # Run create cache dir function which will do nothing if it already
  # exists
  cache_create_dir(cache_dir)

  # Get cache path
  cache_con <- file(cache_path_file)
  cache_path_string <- readLines(cache_con,
    n = 1,
    warn = FALSE
  )
  close(cache_con)

  # Set environmental variable CACHE_DIR_TUSEL to the cache file path
  Sys.setenv("CACHE_DIR_TUSEL" = cache_path_string)
}


#' @title Get local cache path
#'
#' @description
#' Downloaded databases are cached locally. This function
#' gets the local cache path from the  the R environmental variable
#' \code{CACHE_DIR_TUSEL} or returns an error if it does not exist.
#'
#' @rdname cache_get_dir
#' @seealso \code{\link{cache_set_dir}}
#'
#' @export
cache_get_dir <- function() {
  cache_path <- Sys.getenv(
    x = "CACHE_DIR_TUSEL",
    unset = "unset"
  )
  if (cache_path == "unset") {
    stop("Cache directory not set. Run cache_set_dir()")
  } else if (!dir.exists(cache_path)) {
    stop(paste(
      "Cache directory does not exist:", cache_path,
      "Run set_cache_path(cache_directory) to fix."
    ))
  }
  return(cache_path)
}


#' @title Cache a TxDb
#'
#' @description
#' Save a txdb in the local cache.
#'
#' @param txdb a TxDb object
#' @param species_name the name of the species (defaults to using
#' the organism name from the txdb if none specified).
#' @param genome_name the name of the genome (defaults to using
#' the BioMart dataset version from the txdb if none specified).
#' @param annotation_version a version id (defaults to using the BioMart
#' database version if none specified)
#' @param force boolean that forces writing of txdb even if matching one
#' already exists (default: FALSE)
#'
#' @rdname save_txdb
#' @seealso \code{\link{cache_get_dir}}
#'
#' @export
save_txdb <- function(txdb, species_name = NULL, genome_name = NULL,
                      annotation_version = NULL, force = FALSE) {
  # Check that txdb is an actual transcript database object
  if (!methods::is(txdb, "TxDb")) {
    stop("txdb is not a TxDb object")
  }
  # Get metadata from TxDb
  txdb_metadata <- S4Vectors::metadata(txdb)
  txdb_metadata$name <- tolower(gsub(" ", "_", txdb_metadata$name))
  # Set defaults if options not specified
  if (is.null(species_name)) {
    species_name <- tolower(GenomeInfoDb::organism(txdb))
  }
  if (is.null(genome_name)) {
    genome_name <- tolower(txdb_metadata[txdb_metadata$name ==
      "biomart_dataset_version", ]$value)
  }
  if (is.null(annotation_version)) {
    annotation_version <- tolower(txdb_metadata[txdb_metadata$name ==
      "biomart_database_version", ]$value)
  }

  # Remove all spaces and '-' from names
  species_name <- gsub("-|\\s", "_", species_name)
  genome_name <- gsub("-|\\s", "_", genome_name)
  annotation_version <- gsub("-|\\s", "_", annotation_version)

  # Create output name
  txdb_name <- paste0(
    paste(species_name, genome_name, annotation_version, sep = "-"),
    ".txdb"
  )
  cache_dir <- cache_get_dir()
  out_path <- do.call(file.path, as.list(split_path(path = cache_dir)))
  out_path <- file.path(out_path, txdb_name)
  write <- TRUE
  # Check if file already exists
  if (isFALSE(force) && file.exists(out_path)) {
    counter <- 0
    usr_text <- paste0(
      "File already exists:",
      txdb_name
    )
    cat(usr_text)
    overwrite <- tolower(
      readline(prompt = "Overwrite? [y/n]:")
    )
    while (counter < 3 && !overwrite %in% c("y", "n")) {
      cat("Invalid input. Input must be \'y\' or \'n\'")
      overwrite <- tolower(
        readline(prompt = "Overwrite? [y/n]:")
      )
      counter <- counter + 1
    }
    write <- (overwrite == "y")
  }
  # If still okay to write, then do so
  if (isTRUE(write)) {
    message(paste("Saving db to", out_path))
    invisible(AnnotationDbi::saveDb(txdb, file = out_path))
  }
}


#' @title View cache
#'
#' @description
#' View a list of txdbs that have been cached
#'
#' @rdname list_cached_txdb
#' @seealso \code{\link{save_txdb}}
#'
#' @export
list_cached_txdb <- function() {
  cache_dir <- cache_get_dir()
  # List all txdb
  txdb_list <- dir(path = cache_dir, pattern = "*.txdb")
  txdb_table <- stringr::str_match(
    string = txdb_list,
    pattern = "^(.*)-(.*)-(.*).txdb"
  )
  txdb_table <- as.data.frame(txdb_table[, c(2:4, 1), drop = FALSE])
  colnames(txdb_table) <- c("species", "genome", "annotation_version",
                            "file_name")
  return(txdb_table)
}

#' @title Load txdb from cache
#'
#' @description
#' View a list of txdbs that have been cached
#'
#' @param file_name the full file name of the txdb
#'
#' @rdname load_txdb
#' @seealso \code{\link{save_txdb}},\code{\link{list_cached_txdb}}
#'
#' @export
load_txdb <- function(file_name) {
  cache_dir <- cache_get_dir()
  in_path <- file.path(cache_dir, file_name)
  if (file.exists(in_path)) {
    return(invisible(AnnotationDbi::loadDb(file = in_path)))
  } else {
    stop("TxDb does not exist. Run list_cached_txdb() to view available db.")
  }
}

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
#' @param cache whether to automatically cache the downloaded TxDb (Default:
#' TRUE)
#'
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
