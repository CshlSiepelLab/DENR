#' @title Create cache defaults
#'
#' @description Creates paths for default cache_dir and cache_path.txt
#' file
#' @return A list with two elements \code{dir} and \code{path_file}
#'
cache_defaults <- function(){
    # Set default cache dir location
    default_cache_dir = rappdirs::user_data_dir(appname = "tuselecter")
    # File that stores location of cache path
    cache_path_file = file.path(default_cache_dir,"cache_path.txt")
    return(list(dir = default_cache_dir,path_file=cache_path_file))
}


#' @title Create local cache
#'
#' @description
#' Creates local data file cache. If the directory does not exist, it
#' will be created recursively. If no custom path is set, the
#' default user data directory for the package will be used. See
#' \code{\link[rappdirs]{user_data_dir}} for details.
#'
#' @param path
#' Character path to new local files path. If null, path will be
#' reset to default user data directory location.
#'
#' @rdname cache_create_dir
#' @seealso \code{\link{cache_set_dir}} \link{cache_get_dir}
#'
#' @examples
#' \dontrun{
#'   #set a different cache path
#'   set_cache_path('z:/transcript_db/tuselecter')
#' }
#'
#' @export
cache_create_dir <- function(cache_dir = NULL){
    defaults = cache_defaults()
    # Set default cache dir location
    default_cache_dir = defaults$dir
    # File that stores location of cache path
    cache_path_file = defaults$path_file

    # Create cache_dir to serve as true cache dir or place to store
    # cache path file to record cache location between sessions
    if(!dir.exists(default_cache_dir)){
        dir.create(default_cache_dir, recursive = TRUE)
    }

    # If cache_dir is unset use the default directory as the cache_dir
    if(is.null(cache_dir)){
        cache_dir = default_cache_dir
    }

    # Check to see if cache path file already exists, if so, require
    # user input to overwrite/remove if set is different from recorded
    # location
    if(file.exists(cache_path_file)){
        cache_con = file(cache_path_file)
        cache_path_string = readLines(cache_con, n = 1,
                                      warn = FALSE)
        close(cache_con)
        if(suppressWarnings(normalizePath(cache_dir)) != cache_path_string){
            usr_text = paste0('Specified cache directory does not match previous',
                              ' path directory:\n',
                              cache_path_string, ' -> ', cache_dir)
            cat(usr_text)
            overwrite_cache_string = tolower(
                readline(prompt = 'Continue setting cache path [y/n]:')
                )
            while(!overwrite_cache_string %in% c("y","n")){
                cat("Invalid input. Input must be \'y\' or \'n\'")
                overwrite_cache_string = tolower(
                    readline(prompt = 'Continue setting cache path [y/n]:')
                    )
            }
            # If not overwriting cache location set cache_dir to previous location
            if(overwrite_cache_string == "n") {
                cache_dir = cache_path_string
            }
        }
    }

    # Create cache_dir if does not exist
    if(!dir.exists(cache_dir)){
        dir.create(cache_dir, recursive = TRUE)
    }

    # Record cache dir path
    write(x = normalizePath(cache_dir), file = cache_path_file)
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
#' @seealso \code{\link{cache_set_dir}}
#' @rdname cache_set_dir
#' @examples
#'
#' \dontshow{cache_set_dir(temppath=TRUE)}
#'
#' print(cache_get_dir())
#'
#'
#' @export
cache_set_dir <- function(cache_dir = NULL){
    defaults = cache_defaults()
    # File that stores location of cache path
    cache_path_file = defaults$path_file

    # Run create cache dir function which will do nothing if it already
    # exists
    cache_create_dir(cache_dir)

    # Get cache path
    cache_con = file(cache_path_file)
    cache_path_string = readLines(cache_con, n = 1,
                                  warn = FALSE)
    close(cache_con)

    # Set environmental variable CACHE_DIR_TUSEL to the cache file path
    Sys.setenv("CACHE_DIR_TUSEL"= cache_path_string)
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
cache_get_dir <- function(){
    cache_path = Sys.getenv(x = "CACHE_DIR_TUSEL",
                          unset = "unset")
    if(cache_path == "unset"){
        stop("Cache directory not set. Run cache_set_dir()")
    } else if (!dir.exists(cache_path)){
        stop(paste("Cache directory does not exist:", cache_path,
                   "Run set_cache_path(cache_directory) to fix."))
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
#'
#' @rdname cache_txdb
#' @seealso \code{\link{cache_get_dir}}
#'
#' @export
cache_txdb <- function(txdb, species_name = NULL, genome_name = NULL,
                     annotation_version = NULL){
    # Check that txdb is an actual transcript database object
    if(!methods::is(txdb,"TxDb")){
        stop("txdb is not a TxDb object")
    }
    # Get metadata from TxDb
    txdb_metadata = S4Vectors::metadata(txdb)
    txdb_metadata$name = tolower(gsub(" ", "_", txdb_metadata$name))
    # Set defaults if options not specified
    if(is.null(species_name)){
        species_name=organism(txdb)
    }
    if(is.null(species_name)){
        genome_name=txdb_metadata[txdb_metadata$name ==
                                      "biomart_dataset_version",]$value
    }
    if(is.null(species_name)){
        annotation_version=txdb_metadata[txdb_metadata$name ==
                                             "biomart_database_version",]$value
    }

    # Remove all spaces and '-' from names
    species_name = tolower(gsub("-|\\s","_",species_name))
    genome_name = tolower(gsub("-|\\s","_",genome_name))
    annotation_version = tolower(gsub("-|\\s","_",annotation_version))

    # Create output name
    txdb_name = paste0(paste(species_name,genome_name,annotation_version,sep = "-")
                       ,".txdb")
    cache_dir = cache_get_dir()
    out_path =file.path(cache_dir,txdb_name)
    message(paste("Saving db to", out_path))
    invisible(AnnotationDbi::saveDb(txdb,file = out_path))
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
list_cached_txdb <- function(){
    cache_dir = cache_get_dir()
    # List all txdb
    txdb_list = dir(path = cache_dir, pattern = "*.txdb")
    txdb_table = stringr::str_match(string = txdb_list,
                                    pattern = "^(.*)-(.*)-(.*).txdb")
    txdb_table = as.data.frame(txdb_table[,c(2:4,1),drop=FALSE])
    colnames(txdb_table) = c("species","genome","annotation_version","file_name")
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
load_txdb <- function(file_name){
    cache_dir = cache_get_dir()
    in_path = file.path(cache_dir,file_name)
    if(file.exists(in_path)){
        return(invisible(AnnotationDbi::loadDb(file = in_path)))
    } else {
        stop("TxDb does not exist. Run list_cached_txdb() to view available db.")
    }
}


