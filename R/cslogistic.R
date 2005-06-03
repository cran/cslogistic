
.onAttach <- function(...) {
     cat("##\n## Conditionally Specified Logistic Regression Model Package (cslogistic)\n")
     cat("## Copyright (C) 2005, Alejandro Jara and Maria Jose Garcia-Zattera \n")
     cat("##\n## Support provided by the Katholieke Universiteit Leuven \n")
     cat("## (Research Grant OT / 00 / 35) \n##\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("cslogistic", libpath)
}
