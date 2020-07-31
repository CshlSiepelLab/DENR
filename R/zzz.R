.onLoad <- function(libname, pkgname){
  # Avoids any issues with non-standard chromosome names in plotting
  options(ucscChromosomeNames=FALSE)
}