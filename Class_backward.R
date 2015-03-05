# Sources for needed classes
source("Class_paramList.R")

# Class backward
setClass(
    Class="Backward", 
    contains="ParamList",
    validity=function(object) {
        # add verification if needed
    }
)

# Constructor of backward
setMethod(
    f="initialize",
    signature="Backward",
    definition=function(.Object) {
        .Object@niche = composante("niche")
        .Object@dispersion = composante("dispersion")
        .Object@mutation = composante("mutation")
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of backward
backward = function() {
    new(Class="Backward")
}