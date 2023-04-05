


import numpy as np



from pyearth.system.define_global_variables import *

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.unstructured.singlegrid.plot.elm_tsplot_variable_singlegrid import elm_tsplot_variable_singlegrid
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file


sDate = '20230101'
iYear_start = 1980
iYear_end = 2009
sModel = 'e3sm'
sRegion='k34'
for i in range(4, 5,1):
    iCase_index = i

    
    sVariable = 'zwt'    
    sLabel_y = r'Water table depth (m)'    
    sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
    sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'


    aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
    print(aParameter_e3sm)
    oE3SM = pye3sm(aParameter_e3sm)
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
    #print(aParameter_case)
    oCase = pycase(aParameter_case)  

    #subsurface runoff
    sVariable = 'zwt'  
    sLabel_y = r'Water table depth (m)'
    iReverse_y=1
    iFlag_scientific_notation_in=0  
    dMin_y=0
    dMax_y=80

    elm_tsplot_variable_singlegrid(oCase,
                                   iReverse_y_in = iReverse_y,
                                   iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                   dMax_y_in = dMax_y,                                   
                                   dMin_y_in = dMin_y,
                                   aLabel_legend_in=[sLabel_y],
                                   sLabel_y_in= sLabel_y,
                                   sVariable_in = sVariable )
    
    continue
    #subsurface runoff
    sVariable='qdrai'  
    sLabel_y=r'Subsurface runoff (mm/s)'
    iReverse_y=0
    iFlag_scientific_notation_in=1   
    dMin_y=None
    dMax_y=None
    elm_tsplot_variable_singlegrid(oCase,
                                   iReverse_y_in= iReverse_y,
                                   dMax_y_in = dMax_y,
                                   iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                   dMin_y_in =dMin_y,
                                   aLabel_legend_in=[sLabel_y],
                                   sLabel_y_in= sLabel_y,
                                   sVariable_in = sVariable )
    
   
    sVariable='qcharge'
    sLabel_y = r'Groundwater recharge (mm/s)'
    iReverse_y = 0
    iFlag_scientific_notation_in = 1
    dMin_y=None
    dMax_y=None
    elm_tsplot_variable_singlegrid(oCase,
                                   iReverse_y_in= iReverse_y,
                                   dMax_y_in = dMax_y,
                                   iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                   dMin_y_in =dMin_y,
                                   aLabel_legend_in=[sLabel_y],
                                   sLabel_y_in= sLabel_y,
                                   sVariable_in = sVariable )


    print('finished')
