import numpy as np
from pathlib import Path
from os.path import realpath
from datetime import datetime
from netCDF4 import Dataset 
from pathlib import Path
from pyearth.system.define_global_variables import *
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from pyearth.gis.location.convert_lat_lon_range import convert_360_to_180
from pyearth.gis.location.find_index_by_latlon import find_index_by_latlon
from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pye3sm.elm.mesh.elm_retrieve_case_dimension_info import elm_retrieve_case_dimension_info
from pye3sm.elm.general.structured.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d

sPath_parent = str(Path(__file__).parents[4]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_figure = realpath( sPath_parent +  '/figures/' )

#set up latlon

dLon = 299.777709960938 
dLat = -2.62420988082886 

dLon1 = convert_360_to_180(dLon)

#set up time series

sFilename='/qfs/people/lili400/datashare/pf_simulation/par_35_date.csv'

date_all = text_reader_string(sFilename, iSkipline_in=1, cDelimiter_in=',')
sYear ='1905'
nts = date_all.shape[0]
aDate_parflow=list()
aDate_index=list()
for i in range(nts):
    dummy = date_all[i,0]
    if sYear in dummy:
        aDate_index.append(i)
        iMonth=int(dummy.split('-')[1])
        iDay=int(dummy.split('-')[2])
        dSimulation =datetime(2004, iMonth, iDay)
        aDate_parflow.append(dSimulation)
    pass

nts_subset = len(aDate_index)
aDate_index=np.array(aDate_index)
aDate_parflow=np.array(aDate_parflow)

#read in parflow ts
# 2.0) input the lat and lon
sFilename = '/qfs/people/lili400/datashare/pf_simulation/k34_dem1_adj.nc'
aDatasets = Dataset(sFilename)
    
for sKey, aValue in aDatasets.variables.items():
    if (sKey == 'lon'):                   
        aLongitude = (aValue[:]).data
        continue
    if (sKey == 'lat'):                    
        aLatitude = (aValue[:]).data
        continue


nrow=1
ncolumn=1

expnm = 'par_35'
ymdat = '1905'
sFilename = '/qfs/people/lili400/datashare/pf_simulation/par_35_wtd.npy'

wtd_mod = np.load(sFilename) # (nt, ny, nx)
print(wtd_mod.shape)  ### this is what you need ###


mod_day=np.full(nts_subset, -9999, dtype=float)
for i in range(nts_subset):
    a= wtd_mod[aDate_index[i],:,:]
    mod_day[i] =np.mean(a)

#read in hlgf using latlon


#using a case extract
iFlag_monthly=1
iFlag_log=0
sDate = '20230101'


#start loop
iCase_index=4




iYear_start = 2000
iYear_end = 2009
iYear_start_subset=2004
iYear_end_subset=2004
sModel = 'e3sm'
sRegion='k34'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)



sVariable='zwt'

aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                           iCase_index_in =  iCase_index ,\
                                                           iYear_start_in = iYear_start, \
                                                           iYear_end_in = iYear_end,\
                                                            iYear_subset_start_in = iYear_start_subset, \
                                                             iYear_subset_end_in = iYear_end_subset, \
                                                           sDate_in= sDate,\
                                                           sModel_in = sModel, \
                                                               sRegion_in=sRegion,\
                                                           sVariable_in = sVariable )

oCase_x = pycase(aParameter_case)

aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,
                                                           iCase_index_in =  59 ,
                                                           iYear_start_in = iYear_start, \
                                                           iYear_end_in = iYear_end,\
                                                            iYear_subset_start_in = iYear_start_subset, \
                                                             iYear_subset_end_in = iYear_end_subset, \
                                                           sDate_in= sDate,\
                                                           sModel_in = sModel, \
                                                               sRegion_in=sRegion,\
                                                           sVariable_in = sVariable )

oCase_y = pycase(aParameter_case)

aLon, aLat, aMask_ll =elm_retrieve_case_dimension_info(oCase_x)
#dimension
nrow = np.array(aMask_ll).shape[0]
ncolumn = np.array(aMask_ll).shape[1]
aMask_ul = np.flip(aMask_ll, 0)
nrow = np.array(aMask_ll).shape[0]
ncolumn = np.array(aMask_ll).shape[1]
aMask_ll_index = np.where(aMask_ll==1)
aMask_ul_index = np.where(aMask_ul==1)
nrow_extract, ncolumn_extract = aLon.shape

#resolution
dLon_min = np.min(aLon)
dLon_max = np.max(aLon)
dLat_min = np.min(aLat)
dLat_max = np.max(aLat)
dResolution_x = (dLon_max - dLon_min) / (ncolumn_extract-1)
dResolution_y = (dLat_max - dLat_min) / (nrow_extract-1)

row_index, column_index= find_index_by_latlon(dLon1, dLat, dLon_min,dLat_max,dResolution_x, dResolution_y, iFlag_center_in= 1)

iYear_start = oCase_x.iYear_start
iYear_end = oCase_x.iYear_end    
iYear_subset_start = oCase_x.iYear_subset_start
iYear_subset_end = oCase_x.iYear_subset_end
dates = list()
dates_year=list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    dSimulation0 = datetime(iYear, 6, 30)
    dates_year.append( dSimulation0 )
    for iMonth in range(1,13):
        dSimulation = datetime(iYear, iMonth, 1)
        dates.append( dSimulation )

nstress = nyear * nmonth
#take the subset
iMonth = 1
subset_index_start = (iYear_subset_start - iYear_start) * 12 + iMonth-1
subset_index_end = (iYear_subset_end + 1 - iYear_start) * 12 + iMonth-1
subset_index = np.arange( subset_index_start,subset_index_end, 1 )
dates=np.array(dates)
dates_subset = dates[subset_index]
nstress_subset= len(dates_subset)
if iFlag_monthly ==1:
    aData_monthly_x= np.full(nstress_subset, -9999, dtype=float)
    aData_monthly_y= np.full(nstress_subset, -9999, dtype=float)
    aData_ret_x = elm_retrieve_variable_2d( oCase_x, iFlag_monthly_in = 1)        
    aData_ret_y = elm_retrieve_variable_2d( oCase_y, iFlag_monthly_in = 1)    
    for i in np.arange(0, nstress_subset, 1):
        aImage = aData_ret_x[i]     
        dummy1 = np.reshape(aImage, (nrow, ncolumn))
        aImage = aData_ret_y[i]     
        dummy2 = np.reshape(aImage, (nrow, ncolumn))
        aData_monthly_x[i]=dummy1[row_index, column_index]
        aData_monthly_y[i]=dummy2[row_index, column_index]

aDate_all=[aDate_parflow,dates_subset, dates_subset]
aData_all=[mod_day, aData_monthly_x, aData_monthly_y]
sFilename_out=sWorkspace_figure + slash \
                + sVariable.lower() + '_tsplots_parflow' +'.png' 
aColor_in = ['black','red','blue']
sFormat_y='{:.1f}' 
sLabel_y=r'Water table depth (m)'
sTitle=''
aLegend = ['Parflow','ELM-Default','ELM-HLGF']
plot_time_series_data(aDate_all,
                          aData_all,\
                          sFilename_out,\
                          iReverse_y_in = 1, \
                          iFlag_log_in = 0,\
                          iFlag_scientific_notation_in=0,\
                          iFlag_miniplot_in = 0,\
                          ncolumn_in =int(5),\
                          dMax_y_in = 40,\
                          dMin_y_in = 0,\
                          dSpace_y_in = 5, \
                          sTitle_in = sTitle, \
                          sLabel_y_in= sLabel_y,\
                          sFormat_y_in= sFormat_y ,\
                          aLabel_legend_in = aLegend, \
                          aColor_in =aColor_in,\
                          aMarker_in = None,\
                          sLocation_legend_in = 'lower left' ,\
                          aLocation_legend_in = (0.0, 1.0),\
                          aLinestyle_in = None,\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)  