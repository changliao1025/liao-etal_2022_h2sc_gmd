

import os, sys
from re import I
import argparse
import subprocess
import numpy as np
import multiprocessing


from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.structured.twod.map.elm_map_variable_difference_dc_2d import elm_map_variable_difference_dc_2d
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pye3sm.elm.general.structured.twod.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file


sDate = '20220701'


iIndex_start=51
iIndex_end=58
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end



aVariable = ['ZWT']
aFlag_scientific_notation_colorbar=[0]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'


aTitle= [ 'Water table depth' ]
aUnit = [r'Unit: m']
aData_min = [-3]
aData_max = [3]
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

aData_all = list()

for i in range(ncase):

    iCase_index =  aCase_index[i]
    if iCase_index == 52:
        continue
    
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
        dData_min = aData_min[iVariable]
        dData_max = aData_max[iVariable]
        dConversion = aConversion[iVariable]
        iFlag_scientific_notation_colorbar = aFlag_scientific_notation_colorbar[iVariable]    

        aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                       iCase_index_in =  52 ,\
                                                       iYear_start_in = iYear_start, \
                                                       iYear_end_in = iYear_end,\
                                                        iYear_subset_start_in = iYear_start, \
                                                         iYear_subset_end_in = iYear_end, \
                                                            dConversion_in= dConversion,\
                                                       sDate_in= sDate,\
                                                       sModel_in = sModel, \
                                                           sRegion_in=sRegion,\
                                                       sVariable_in = sVariable )

        oCase_x = pycase(aParameter_case)    
        aData_out_x = elm_retrieve_variable_2d(  oCase_x,\
                iFlag_monthly_in =0,\
                    iFlag_annual_mean_in=1,\
                iFlag_annual_total_in= 0     )

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
        oCase_y = pycase(aParameter_case)

        aData_out_y = elm_retrieve_variable_2d(  oCase_y,\
                iFlag_monthly_in =0,\
                    iFlag_annual_mean_in=1,\
                iFlag_annual_total_in= 0     )

        dummy_x = aData_out_x.pop()
        dummy_y = aData_out_x.pop()

        nan_index = np.where(dummy_x ==  -9999)

        dummy=dummy_y-dummy_x
        dummy[nan_index] = -9999
        aData_all.append(dummy)
        
        pass

#get 

aPercentiles_in = np.arange(10, 90, 15)
aInterval = cgpercentiles(aData_all, aPercentiles_in, missing_value_in = -9999)  

aColor = create_diverge_rgb_color_hex(  len(aInterval) +1  )
aColor.reverse()

for i in range(0,ncase):

    iCase_index =  aCase_index[i]
    if iCase_index == 52:
        continue
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
        dData_min = aData_min[iVariable]
        dData_max = aData_max[iVariable]
        dConversion = aConversion[iVariable]
        iFlag_scientific_notation_colorbar = aFlag_scientific_notation_colorbar[iVariable]

        aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                       iCase_index_in =  52 ,\
                                                       iYear_start_in = iYear_start, \
                                                       iYear_end_in = iYear_end,\
                                                        iYear_subset_start_in = iYear_start, \
                                                         iYear_subset_end_in = iYear_end, \
                                                            dConversion_in= dConversion,\
                                                       sDate_in= sDate,\
                                                       sModel_in = sModel, \
                                                           sRegion_in=sRegion,\
                                                       sVariable_in = sVariable )

        oCase_x = pycase(aParameter_case)

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
        oCase_y = pycase(aParameter_case)

        aLegend = list()
        sCase = "{:0d}".format(iCase_index - 50)
        sText = 'Case ' + sCase + ' - ' + 'Case 2' 
        aLegend.append(sText)
        elm_map_variable_difference_dc_2d(oE3SM, oCase_x , oCase_y, dData_min_in=dData_min, dData_max_in=dData_max,\
            iFlag_scientific_notation_colorbar_in = iFlag_scientific_notation_colorbar, \
                   iFlag_monthly_in =0,\
                    iFlag_annual_mean_in=1,\
                iFlag_annual_total_in= 0,\
                     aInterval_in = aInterval,\
                sExtend_in = sExtend,\
            sUnit_in = sUnit,\
            sTitle_in=  sTitle,\
                aColor_in=aColor,\
            aLegend_in = aLegend )
        pass

print('finished')
