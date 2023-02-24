library(netcdf)

ncin <- nc_open("~/Downloads/sresa1b_ncar_ccsm3-example.nc")

lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
tmp_array <- ncvar_get(ncin, "tas")

image(lon,lat,tmp_array, col=rev(brewer.pal(10,"RdBu")))