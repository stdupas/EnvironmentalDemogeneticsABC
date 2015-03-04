# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
	Class="ParamList", 
    representation=representation(
    	niche="Composante",
    	dispersion="Composante",
    	mutation="Composante"
    ),
    prototype=prototype(
    	niche=NULL,
    	dispersion=NULL,
    	mutation=NULL
    ),
    validity=function(object) {
    	cat("---------- ParamList : verification ----------\n")
    	# add verification if needed
    }
)
