import numpy as np

import datetime
import netCDF4 as nc
from pyearth.system.define_global_variables import *
 
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase

from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file

from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data


sDate = '20230401'
iCase_index = 2

#site location
dLatitude = -2.6091
dLongitude = -60.2093

iYear_start = 1980
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'

sVariable = 'discharge'
sVariable = 'Main_Channel_Water_Depth_LIQ' #the river gage height
sUnit = r'Units: m3/s'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/case.xml'

aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,
                                                       iCase_index_in =  iCase_index ,
                                                       iYear_start_in = iYear_start, 
                                                       iYear_end_in = iYear_end,
                                                        iYear_subset_start_in = iYear_start, 
                                                         iYear_subset_end_in = iYear_end, 
                                                       sDate_in= sDate,
                                                       sModel_in = sModel, 
                                                           sRegion_in=sRegion,
                                                       sVariable_in = sVariable )

oCase = pycase(aParameter_case)
sFilename_domain = '/compyfs/liao313/00raw/drof/mosart_domain_amazon_2d_halfdegree.nc'
sWorkspace_drof ='/compyfs/liao313/00raw/drof'

#read the domain file
aDatasets = nc.Dataset(sFilename_domain)
netcdf_format = aDatasets.file_format    
print(netcdf_format)

for sKey, aValue in aDatasets.variables.items():
    if "mask" == sKey:
        aMask = (aValue[:]).data            
    if "xc" == sKey:
        aLon = (aValue[:]).data            
    if "yc" == sKey:
        aLat = (aValue[:]).data 

aMask=np.flip(aMask, 0)   
aLon=np.flip(aLon, 0) 
aLat=np.flip(aLat, 0) 

dLat_max = np.max(aLat)
dLat_min = np.min(aLat)
dLon_min = np.min(aLon)

nrow = len(aLat)
#get the index of the site
iIndex_lon = int((dLongitude - dLon_min) / 0.5)

iIndex_lat_flipped = int(( dLatitude- dLat_min) / 0.5) 
nflip = iIndex_lat_flipped + 1

iIndex_lat =  nrow - nflip 

print(iIndex_lon, iIndex_lat)
aDate=list()
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        for iDay in range(1, 32):
            try:
                dSimulation = datetime.datetime(iYear, iMonth, iDay)
                aDate.append( dSimulation )
            except ValueError:
                break
            pass
        
        pass




aTS=list()

for iYear in range(iYear_start, iYear_end+1):
    sYear = '{0:04d}'.format(iYear)
    #for iMonth in range(1, 12+1):
    sFilename = sWorkspace_drof + slash + 'drof_' + sYear + '.nc'
    aDatasets = nc.Dataset(sFilename)
    for sKey, aValue in aDatasets.variables.items():
        if "main_channel_water_depth_liq" == sKey:
            aData0= (aValue[:]).data 
            
            break
        
    aData0=np.flip(aData0, 1)
    #find if there is a value close to 0.8714
   
    aData = aData0[:, iIndex_lat, iIndex_lon]
    #attach each element to a list
    for k in range(0, len(aData)):
        aTS.append(aData[k])
#reshape the array as 1D
aTS = np.array(aTS)

#save the time series as a csv file
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/gage_height_daily'+'.csv'
np.savetxt(sFilename_out, aTS, delimiter=",")

#plot the time series
aTime = np.array([aDate])
aData = np.array([aTS])
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/gage_height_daily'+'.png'
print(np.max(aData))
plot_time_series_data( aTime , aData,
                           sFilename_out,
                           iFlag_scientific_notation_in=0,
                           iReverse_y_in = 0,
                           sTitle_in = 'River gage height',
                           sLabel_y_in= 'River gage height (m)',
                           dMax_y_in = 5,
                           dMin_y_in =0,
                           aColor_in = ['blue'] ,
                           aLabel_legend_in=['River gage height'],
                           iSize_x_in = 12,
                           iSize_y_in = 5)

print('finished')
