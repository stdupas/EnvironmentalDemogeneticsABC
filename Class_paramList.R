# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
	Class="ParamList", 
    representation=representation(
    	niche="Composante",
    	dispersion="Composante",
    	mutation="Composante"
    )
)