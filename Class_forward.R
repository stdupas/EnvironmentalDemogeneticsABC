# Sources for needed classes
source("Class_composante.R")
source("Class_paramList.R")

# Class forward
setClass(
	Class="Forward",
	representation=representation(
		generation="Composante"
	),
	contains="ParamList",
    validity=function(object) {
    	cat("---------- Forward : verification ----------\n")
    	# add verification if needed
    }
)

# Constructor of forward
setMethod(
	f="initialize",
	signature="Forward",
	definition=function(.Object, niche, dispersion, mutation, generation) {
		cat("---------- Forward : initiation ----------\n")
		if(!missing(niche) && !missing(dispersion) && !missing(mutation) && !missing(generation)) {
			.Object@niche = niche
			.Object@dispersion = dispersion
			.Object@mutation = mutation
			.Object@generation = generation
			validObject(.Object)
		}
		else {
            stop("[Forward initiation] Missing argument(s)\n")
        }
        return(.Object)
	}
)

# User-friendly constructor of forward
forward = function(niche, dispersion, mutation, generation) {
	cat("---------- Forward : construction ----------\n")
	new(Class="Forward", niche=niche, dispersion=dispersion, mutation=mutation, generation=generation)
}