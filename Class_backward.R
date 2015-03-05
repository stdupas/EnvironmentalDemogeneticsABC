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
    definition=function(.Object) {
        cat("---------- Backward : initiation ----------\n")
        niche = composante("niche")
        dispersion = composante("dispersion")
        mutation = composante("mutation")
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of backward
backward = function() {
    cat("---------- Backward : construction ----------\n")
    new(Class="Backward")
}