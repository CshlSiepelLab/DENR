.onLoad <- function(libname, pkgname) { # nolint
  # Avoids any issues with non-standard chromosome names in plotting
  options(ucscChromosomeNames = FALSE)
}
