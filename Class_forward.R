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
    	# add verification if needed
    }
)

# Constructor of forward
setMethod(
	f="initialize",
	signature="Forward",
	definition=function(.Object) {
		.Object@niche = composante("niche")
		.Object@dispersion = composante("dispersion")
		.Object@mutation = composante("mutation")
		.Object@generation = composante("generation")
		validObject(.Object)
        return(.Object)
	}
)

# User-friendly constructor of forward
forward = function() {
	new(Class="Forward")
}