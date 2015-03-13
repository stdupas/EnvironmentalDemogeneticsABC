res = NULL
lat = c(48.8566,49.8854,43.6046,43.2964,41.9027)
lon = c(2.3522,2.2875,1.4442,5.3697,12.4963)
mat = NULL
for(i in 1:5) {
	res = NULL
	for(j in 1:5) {
		res = c(res,distance(lat[i],lat[j],lon[i],lon[j]))
	}
	mat = rbind(mat,res)
}

