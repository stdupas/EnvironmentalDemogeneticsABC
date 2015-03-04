# Sources for needed classes
source("Class_composante.R")
source("Class_paramList.R")

# Class forward
setClass(
	Class="Forward",
	representation=representation(
		generation="Composante"
	),
	contains="ParamList"
)