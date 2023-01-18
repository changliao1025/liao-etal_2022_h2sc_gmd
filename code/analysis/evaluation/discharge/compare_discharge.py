import os
from os.path import realpath
import numpy as np
from datetime import datetime
from pyearth.system.define_global_variables import *
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file, pye3sm_read_case_configuration_file
from pye3sm.mosart.evaluation.halfdegree.mosart_evaluate_stream_discharge_gsim import mosart_evaluate_stream_discharge_gsim
from pye3sm.mosart.mesh.mosart_retrieve_case_dimension_info import mosart_retrieve_case_dimension_info 
from pye3sm.mosart.general.mosart_prepare_outlet_coordinates_with_gsim_filenames import mosart_prepare_outlet_coordinates_with_gsim_filenames
from pye3sm.tools.mosart.gsim.read_gsim_data import read_gsim_data
from pye3sm.mosart.general.structured.twod.retrieve.mosart_retrieve_variable_2d import mosart_retrieve_variable_2d
from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data

sPath_parent = str(Path(__file__).parents[4]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_figure = realpath( sPath_parent +  '/figures/' )

sWorkspace_scratch_in = '/compyfs/liao313/'
sModel = 'e3sm'
sRegion = 'amazon'

sDate = '20220701'

sVariable = 'discharge'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'


#dLongitude = -55.5131  #make sure east and west
#dLatitude = -1.9192


aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration )
oE3SM = pye3sm(aParameter_e3sm)
iYear_start = 2000
iYear_end = 2009
iYear_subset_start = 1990
iYear_subset_end =1998
sLabel_y = r'River discharge ($m^{3} s^{-1}$)'

sFilename_mosart_gsim_info = '/qfs/people/liao313/data/h2sc/global/auxiliary/basin_ind.txt'

#the default case
aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                           iCase_index_in =  62 ,\
                                                           iYear_start_in = iYear_start, \
                                                           iYear_end_in = iYear_end,\
                                                           iYear_subset_start_in = iYear_subset_start, \
                                                           iYear_subset_end_in = iYear_subset_end, \
                                                           sDate_in= sDate,\
                                                           sLabel_y_in =  sLabel_y, \
                                                           sVariable_in = sVariable,\
                                                                 sModel_in = sModel, \
                                                           sRegion_in=sRegion,\
                                                           sWorkspace_scratch_in = sWorkspace_scratch_in )

oCase_x = pycase(aParameter_case)
#the hlgf case
aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                           iCase_index_in =  59 ,\
                                                           iYear_start_in = iYear_start, \
                                                           iYear_end_in = iYear_end,\
                                                           iYear_subset_start_in = iYear_subset_start, \
                                                           iYear_subset_end_in = iYear_subset_end, \
                                                           sDate_in= sDate,\
                                                           sLabel_y_in =  sLabel_y, \
                                                           sVariable_in = sVariable,\
                                                                 sModel_in = sModel, \
                                                           sRegion_in=sRegion,\
                                                           sWorkspace_scratch_in = sWorkspace_scratch_in )
oCase_y = pycase(aParameter_case)
#read data
iYear_start = oCase_x.iYear_start
iYear_end = oCase_x.iYear_end
iYear_subset_start = oCase_x.iYear_subset_start
iYear_subset_end = oCase_x.iYear_subset_end
sLabel_Y = oCase_x.sLabel_y
dConversion = oCase_x.dConversion
sVariable = oCase_x.sVariable

sWorkspace_analysis_case = oCase_x.sWorkspace_analysis_case
sWorkspace_variable_data = sWorkspace_analysis_case + slash + sVariable +  slash + 'tiff'

 #new approach
aLon, aLat , aMask_ll= mosart_retrieve_case_dimension_info(oCase_x)
#dimension
aMask_ul = np.flip(aMask_ll, 0)
nrow = np.array(aMask_ll).shape[0]
ncolumn = np.array(aMask_ll).shape[1]
aMask_ll_index = np.where(aMask_ll==0)
aMask_ul_index = np.where(aMask_ul==0)
#resolution
dLon_min = np.min(aLon)
dLon_max = np.max(aLon)
dLat_min = np.min(aLat)
dLat_max = np.max(aLat)
dResolution_x = (dLon_max - dLon_min) / (ncolumn-1)
dResolution_y = (dLat_max - dLat_min) / (nrow-1)
#the gsim file to be read

aBasin, aLat, aLon, aID, aFilename_gsim = mosart_prepare_outlet_coordinates_with_gsim_filenames(sFilename_mosart_gsim_info)

nsite = len(aFilename_gsim)
sWorkspace_gsim = '/compyfs/liao313/00raw/hydrology/gsim/GSIM_indices/TIMESERIES/monthly'

#read basin mask

dates = list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dSimulation = datetime(iYear, iMonth, 15)
        dates.append( dSimulation )
nstress = nyear * nmonth
#take the subset
iMonth = 1
subset_index_start = (iYear_subset_start - iYear_start) * 12 + iMonth-1
subset_index_end = (iYear_subset_end + 1 - iYear_start) * 12 + iMonth-1
subset_index1 = np.arange( subset_index_start,subset_index_end, 1 )
subset_index2 = np.arange( subset_index_start+1,subset_index_end+1, 1 )
dates=np.array(dates)
dates_subset1 = dates[subset_index1]
dates_subset2 = dates[subset_index2]
nstress_subset= len(dates_subset1)

aData_ret_x = mosart_retrieve_variable_2d( oCase_x, iFlag_monthly_in = 1)

aData_ret_y = mosart_retrieve_variable_2d( oCase_y, iFlag_monthly_in = 1)

for iSite in np.arange(1, nsite+1, 1):
    dLatitude = aLat[iSite-1]
    dLongitude = aLon[iSite-1]
    if dLatitude >=dLat_min and dLatitude<=dLat_max and dLongitude >=dLon_min and dLongitude<=dLon_max:
        pass
    else:
        continue
    sFilename_gsim = aFilename_gsim[iSite-1] + '.mon'
    sFilename_gsim = 'BR_0000244.mon'
    sFilename_gsim =  os.path.join(sWorkspace_gsim, sFilename_gsim) 
    print(sFilename_gsim)
    #continue
    aData = read_gsim_data(sFilename_gsim,iYear_start, iYear_end )
    row_index = int( (dLat_max-dLatitude)/0.5 )
    column_index =int ((dLongitude -dLon_min)/0.5)

    aVariable2 = np.full(nstress_subset, -9999, dtype=float)
    aVariable3 = np.full(nstress_subset, -9999, dtype=float)
    for i in np.arange(0, nstress_subset, 1):
        dummy1 = aData_ret_x[i]
        dummy2 = dummy1[ row_index, column_index ]
        #dummy3 =  area[row_index, column_index]
        aVariable2[i] = dummy2 #* dummy3 /1000.0

        dummy1 = aData_ret_y[i]
        dummy2 = dummy1[ row_index, column_index ]
        #dummy3 =  area[row_index, column_index]
        aVariable3[i] = dummy2 #* dummy3 /1000.0
        pass
    sTitle_in=''
    sLabel_Y=''
    #generate the time series plot using the pyes library
    aDate_all = np.array([dates_subset1, dates_subset1, dates_subset1])
    aData_all = np.array([aVariable2,aVariable3, aData[subset_index2] ])
    aLabel_legend=['ELM default - MOSART river discharge','ELM HLGF - MOSART river discharge','Observed river discharge']
    sFilename_out=sWorkspace_figure + slash \
                + sVariable.lower() + '_tsplots' +'.png' 
    sLabel_Y = r'River discharge (m3/s)'
    aColor = ['red', 'blue', 'black']
    plot_time_series_data(aDate_all,
                          aData_all ,\
                          sFilename_out,\
                          iReverse_y_in = 0, \
                            iFlag_scientific_notation_in=1,\
                            ncolumn_in= 3,\
                            aColor_in= aColor,\
                          sTitle_in = sTitle_in, \
                          sLabel_y_in= sLabel_Y,\
                          aLabel_legend_in = aLabel_legend, \
                          sLocation_legend_in = 'upper right' ,\
                          aLocation_legend_in = (1.0, 1.0),\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)


print('finished')