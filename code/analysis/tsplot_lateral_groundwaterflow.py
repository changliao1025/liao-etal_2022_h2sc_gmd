
import os, sys

import numpy as np
import datetime
from os.path import realpath

from pyearth.system.define_global_variables import *
 
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase

from pye3sm.elm.general.structured.twod.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.timeseries.plot_time_series_data_w_variation import plot_time_series_data_w_variation

sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = sPath_parent +  '/figures/'

sDate = '20220701'


iIndex_start=52
iIndex_end=58
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end

iFlag_scientific_notation_colorbar_in = 1

aVariable = ['QDRAI']
aFlag_scientific_notation_colorbar=[1]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'


aTitle= [ r'Subsurface runoff' ]
aUnit = [ r'Units: mm/s']
aData_min = [-1E-5]
aData_max = [1E-5]
aConversion = [1]

sExtend='both'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
ncase= len(aCase_index)
nvariable = len(aVariable)

#gather data

aData_all = list()
aData_upper_all = list()
aData_lower_all = list()

for i in range(ncase):

    iCase_index =  aCase_index[i]
    
    
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
        dData_min = aData_min[iVariable]
        dData_max = aData_max[iVariable]
        dConversion = aConversion[iVariable]
        iFlag_scientific_notation_colorbar = aFlag_scientific_notation_colorbar[iVariable]    

        aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                       iCase_index_in =  iCase_index ,\
                                                       iYear_start_in = iYear_start, \
                                                       iYear_end_in = iYear_end,\
                                                        iYear_subset_start_in = iYear_start, \
                                                         iYear_subset_end_in = iYear_end, \
                                                            dConversion_in= dConversion,\
                                                       sDate_in= sDate,\
                                                       sModel_in = sModel, \
                                                           sRegion_in=sRegion,\
                                                       sVariable_in = sVariable )

        oCase_in = pycase(aParameter_case)    
        aData_out_x = elm_retrieve_variable_2d(  oCase_in,\
                iFlag_monthly_in =1,\
                    iFlag_annual_mean_in=0,\
                iFlag_annual_total_in= 0     )

        iYear_start = oCase_in.iYear_start
        iYear_end = oCase_in.iYear_end
        iYear_subset_start = oCase_in.iYear_subset_start
        iYear_subset_end = oCase_in.iYear_subset_end
        dates = list()
        nyear = iYear_end - iYear_start + 1
        for iYear in range(iYear_start, iYear_end + 1):
            for iMonth in range(1,13):
                dSimulation = datetime.datetime(iYear, iMonth, 15)
                dates.append( dSimulation )

        dates=np.array(dates)
        nstress = nyear * nmonth
        iMonth = 1
        subset_index_start = (iYear_subset_start - iYear_start) * 12 + iMonth-1
        subset_index_end = (iYear_subset_end + 1 - iYear_start) * 12 + iMonth-1
        subset_index = np.arange( subset_index_start,subset_index_end, 1 )
        dates_subset = dates[subset_index]
        nstress_subset= len(dates_subset)
        aMean=np.full(nstress_subset, -9999, dtype=float)
        aUpper=np.full(nstress_subset, -9999, dtype=float)
        aLower=np.full(nstress_subset, -9999, dtype=float)

        aData_out_x=np.array(aData_out_x)
        for iStress in range(1, nstress_subset+1,1):
            dummy= aData_out_x[iStress-1,:,:]
            good_index = np.where(dummy !=  -9999)
            dummy4 = dummy[good_index]

            dummy1 = np.mean(dummy4)

            aPercentiles_in = np.arange(25, 76, 50)

            aInterval = cgpercentiles(dummy4, aPercentiles_in, missing_value_in = -9999) 

            dummy5 = np.log10(dummy1)  

            #dummy0 = np.std(dummy4)
            
            dummy2 = np.log10(aInterval[1])
            dummy3 = np.log10(aInterval[0])

            aMean[iStress -1 ]= dummy5
            aUpper[iStress -1 ]= dummy2
            aLower[iStress -1 ]= dummy3

        aData_all.append(aMean)
        aData_upper_all.append(aUpper)
        aData_lower_all.append(aLower)
        
        pass

#get 




aDate_all = np.tile([dates_subset], [ncase,1])

sFilename_out= sWorkspace_output + '/tsplot_qdrai.png'
sTitle_in =''
sLabel_Y=''
aLabel_legend=list()
for i in range(ncase):
    sCase = "{:0d}".format(i +1)
    sText = 'Case ' + sCase 
    aLabel_legend.append( sText )
aColor = create_diverge_rgb_color_hex(ncase)
plot_time_series_data_w_variation(aDate_all,
                          aData_all,\
                            aData_upper_all,\
                                aData_lower_all,\
                          sFilename_out,\
                            iFlag_scientific_notation_in=0,\
                          iReverse_y_in = 0, \
                          iFlag_log_in = 1,\
                          ncolumn_in = 4,\
                          #dMax_y_in = -5,\
                         # dMin_y_in = -7,\
                          #dSpace_y_in = dSpace_y_in, \
                          sTitle_in = sTitle_in, \
                          sLabel_y_in= sLabel_Y,\
                          sFormat_y_in= '{:.3f}' ,\
                          aLabel_legend_in = aLabel_legend, \
                          aColor_in = aColor,\
                          #aMarker_in = ['o','.','*','+'],\
                          sLocation_legend_in = 'lower right' ,\
                          aLocation_legend_in = (1.0, 0.0),\
                          #aLinestyle_in = ['-','--','-.' ,'solid'],\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)

print('finished')


