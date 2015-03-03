# Sources for needed classes
source("Class_composante.R")
source("Class_paramList.R")

# Class forward
setClass(
	Class="forward",
	representation=representation(
		generation="composante"
	),
	contains="paramList"
)