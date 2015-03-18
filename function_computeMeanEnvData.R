# Function to compute the mean of var ("pr", "tasmin", "tasmax")
# 	nbDays: number of days

computeMeanEnvData = function(data, var, nbDays) {
	dataOld = data
	for(d in 1:length(data[,1,var])) {
		for(i in (nbDays+1):length(data[d,,var])) {
			m = mean(dataOld[d,(i-nbDays):(i-1),var])
			data[d,i,var] = m
		}
	}
	return(data)
}

