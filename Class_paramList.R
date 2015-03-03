# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
	Class="paramList", 
    representation=representation(
    	niche="composante",
    	dispersion="composante",
    	mutation="composante"
    )
)