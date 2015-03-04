# Sources for needed classes
source("Class_paramList.R")

# Class backward
setClass(
	Class="Backward", 
	contains="ParamList",
    validity=function(object) {
    	cat("---------- Backward : verification ----------\n")
    	# add verification if needed
    }
)

# Constructor of backward
setMethod(
	f="initialize",
	signature="Backward",
	definition=function(.Object, niche, dispersion, mutation) {
		cat("---------- Backward : initiation ----------\n")
		if(!missing(niche) && !missing(dispersion) && !missing(mutation)) {
			.Object@niche = niche
			.Object@dispersion = dispersion
			.Object@mutation = mutation
			validObject(.Object)
		}
		else {
            stop("[Backward initiation] Missing argument(s)\n")
        }
        return(.Object)
	}
)

# User-friendly constructor of backward
backward = function(niche, dispersion, mutation) {
	cat("---------- Backward : construction ----------\n")
	new(Class="Backward", niche=niche, dispersion=dispersion, mutation=mutation)
}