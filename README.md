EnvironmentalDemogeneticsABC
========

Statistical models of coalescent on a graph

## Coalescence in an Approximate Bayesian Computation framework

### Model
Use *askListOfParameters* function to precise the model. The function will successively ask for a model, for the prior distribution of each parameter, and finally for values of hyperparameters. Then parameters values are drawn and stored with the model in a list object.

### Coalescent simulations
Use *simSpatialCoal* to simulate spatial coalescence.
#### Required objects
	The model list, a rasterStack of environment variables and a **spatial genetic dataset** :
| x | y | Locus 1 | Locus 2 | ... | Locus 3 |
|---|---|---------|---------|-----|---------|
|0.5|1|287|180|242|136|

#### Function
The function makes use of *mclapply* (package *parallel*) to perform parallel computing. This function relies on forking implementation, but **Windows does not support forking**. Furthermore, even on Linux, *mclapply* seems to behave quite weirdly. This will be settled later. 

#### Results
*simSpatialCoal* function creates a repertory in the current directory, named SimulResults, where are written the simulated genetic values for each locus and individuals.

### ABC analysis 
Makes use of *abc* package functions.
#### Summary statistics
Use *pca4abc* function to computes summary statistics of simulated genetic dataset. Firstly, PC are computed according to the observed genetic data, then transformation is applied to all simulated genetic data. This function returns a list of three objects that can be used in *abc* package functions : a vector of observed summary statistics, a matrix of simulated summary statistics, and a matrix of parameters values used for simulations.
 
#### Cross validation
*abc* package can perform cross validation, using *cv4abc*
#### ABC
*abc* function of *abc* package performs a classic analysis


## Forward models
