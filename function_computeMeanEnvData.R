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


aggregateDays = function(data, nbDays, Dates) {
	nbDates = round(length(Dates)/nbDays)-1
	dataNew = array(0, dim=c(dim(data)[1],nbDates,dim(data)[3]))
	# pour chaque couche
	for(c in 1:dim(data)[3]) {
		# pour chaque deme
		for(d in 1:dim(data)[1]) {
			j = 1
			# pour toutes les nbDays dates
			for(i in seq(nbDays,length(Dates),by=nbDays)) {
				m = mean(data[d,(i-nbDays+1):i,c])
				dataNew[d,j,c] = m
				j = j+1
			}
		}
	}
	dimnames(dataNew)[1] = dimnames(data)[1]
	dimnames(dataNew)[2] = dimnames(data[,seq(nbDays,length(Dates),by=nbDays),])[2]
	dimnames(dataNew)[3] = dimnames(data)[3]

	return(dataNew)
}