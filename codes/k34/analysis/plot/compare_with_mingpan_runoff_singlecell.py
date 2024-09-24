import netCDF4 as nc
#use the mingpan runoff data to compare with the simulated total runoff

sFolder_mingpan_runoff = ''

#we will focus on the time period from 2000 to 2009
iYear_start = 2000
iYear_end = 2009

#use a for loop read the mingpan runoff


#use the first year data to retrieve the domain information
sYear = "{:04d}".format(iYear_start)
sFilename = sFolder_mingpan_runoff + 'runoff_' + sYear + '.nc'
pData_domain = nc.Dataset(sFilename)
#read the longitude and latitude
aLon = pData_domain.variables['lon'][:]
aLat = pData_domain.variables['lat'][:]
#get the dimension of the domain
ncolumn = len(aLon)
nrow = len(aLat)
#get the resolution of the data in x and y direction
dResolution_x = aLon[1] - aLon[0]
dResolution_y = aLat[1] - aLat[0]
#close the file
pData_domain.close()

for iYear in range(iYear_start, iYear_end+1):
    sYear = "{:04d}".format(iYear)
    sFilename = sFolder_mingpan_runoff + 'runoff_' + sYear + '.nc'

    pass