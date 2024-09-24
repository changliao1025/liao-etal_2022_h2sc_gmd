import numpy as np
from pyearth.system.define_global_variables import *

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.unstructured.singlecell.plot.elm_tsplot_variable_singlecell import elm_tsplot_variable_singlecell
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file

sDate = '20230401'
iYear_start = 2000
iYear_end = 2009

sModel = 'e3sm'
sRegion='k34'

iFlag_zwt = 1
iFlag_wt_slp = 1
iFlag_recharge = 1
iFlag_drainage = 1
iFlag_over= 1
iFlag_runoff = 0
iFlag_downslope = 0
iFlag_seepage = 0
iFlag_macropore = 0
iFlag_gage_height = 1
iFlag_water_table_hillslope = 0

for i in range(1,2,1):
    iCase_index = i

    sVariable = 'zwt'
    sLabel_y = r'Water table depth (m)'
    sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/e3sm.xml'
    sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/case.xml'

    aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
    print(aParameter_e3sm)
    oE3SM = pye3sm(aParameter_e3sm)
    aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,
                                                           iCase_index_in =  iCase_index,
                                                           iYear_start_in = iYear_start,
                                                           iYear_end_in = iYear_end,
                                                           iYear_subset_start_in = iYear_start,
                                                           iYear_subset_end_in = iYear_end,
                                                           sDate_in= sDate,
                                                           sModel_in = sModel,
                                                           sRegion_in = sRegion,
                                                           sVariable_in = sVariable )

    oCase = pycase(aParameter_case)

    #water table depth
    if iFlag_zwt == 1:
        sVariable = 'zwt'
        sLabel_y = r'Water table depth (m)'
        iReverse_y=1
        iFlag_scientific_notation_in=0
        dMin_y=0
        dMax_y=10
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in = iReverse_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       dMax_y_in = dMax_y,
                                       dMin_y_in = dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )




    if iFlag_runoff == 1:
        sVariable='qrunoff'
        sLabel_y = r'Total runoff (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_over == 1:
        sVariable='qover'
        sLabel_y = r'Surface runoff (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=0,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    #subsurface runoff
    if iFlag_drainage == 1:
        sVariable='qdrai'
        sLabel_y=r'Total subsurface runoff (mm/s)'
        iReverse_y=0
        iFlag_scientific_notation_in=1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=0,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )



    #sVariable='qcharge'
    if iFlag_recharge == 1:
        sLabel_y = r'Groundwater recharge (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_downslope == 1:
        sVariable='qdrai_downslope'
        sLabel_y = r'Sursurface drainage (downslope) (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_seepage == 1:
        sVariable='qdrai_seepage'
        sLabel_y = r'Sursurface drainage (seepage) (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_macropore== 1:
        sVariable='qdrai_macropore'
        sLabel_y = r'Sursurface drainage (macropore) (mm/s)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 1
        dMin_y=0
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_wt_slp == 1:
        sVariable='wt_slp'
        sLabel_y = r'Groundwater table slope'
        iReverse_y = 0
        iFlag_scientific_notation_in = 0
        dMin_y=None
        dMax_y=None
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_gage_height == 1:
        sVariable='gage_height'
        sLabel_y = r'Gage height (m)'
        iReverse_y = 0
        iFlag_scientific_notation_in = 0
        dMin_y=0
        dMax_y=5
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=0,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )

    if iFlag_water_table_hillslope == 1:
        sVariable='wt_hillslope'
        sLabel_y = r'Wate table depth (m)'
        iReverse_y = 1
        iFlag_scientific_notation_in = 0
        dMin_y=-40
        dMax_y=40
        elm_tsplot_variable_singlecell(oCase,
                                       iReverse_y_in= iReverse_y,
                                       dMax_y_in = dMax_y,
                                       iFlag_scientific_notation_in=iFlag_scientific_notation_in,
                                       iFlag_monthly_in=1,
                                       iFlag_daily_in=1,
                                       dMin_y_in =dMin_y,
                                       aLabel_legend_in=[sLabel_y],
                                       sLabel_y_in= sLabel_y,
                                       sVariable_in = sVariable )



    print('finished')
