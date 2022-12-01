import numpy as np

from pyearth.system.define_global_variables import *   
from pyearth.gis.gdal.read.gdal_read_geotiff_file import gdal_read_geotiff_file

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.structured.twod.plot.elm_scatterplot_variable_x_2d import elm_scatterplot_variable_x_2d
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
#read observation
sFilename_in = '/qfs/people/liao313/data/e3sm/amazon/elm/' + 'wtd_extract' + sExtension_tiff
    
aData_obs, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatial_reference = gdal_read_geotiff_file(sFilename_in)
dLon_min = dOriginX
dResolution_x= dPixelWidth
dLon_max = dLon_min + ncolumn * dPixelWidth
dLat_min =  dOriginY - nrow * dPixelWidth
dLat_max= dOriginY
aImage_extent= [dLon_min ,dLon_max , dLat_min ,  dLat_max]
#read simulated and apply average

sDate = '20220701'
#read slope from a case

iIndex_start=52
iIndex_end=52
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end

iFlag_scientific_notation_colorbar_in = 0

aVariable = ['sur_slp']
aFlag_scientific_notation_colorbar=[0,1,1,1]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'


aTitle= [ 'Water table depth' ]
aUnit = [r'Unit: m']

aConversion = [100]

sExtend='both'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
nvariable = len(aVariable)
ncase= len(aCase_index)
aData_y_in = aData_obs
dMin_x = 0
dMax_x= 20
dMin_y = 0
dMax_y = 10
sLabel_x = r'Surface slope (percent)'
sLabel_y = r'Water table depth (m)'
sLabel_legend = 'Observed WTD'
sVariable_y_in=r'obs_wtd'
for i in range(0,ncase):

    iCase_index =  aCase_index[i]
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
       
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

        oCase_x = pycase(aParameter_case)


        aLegend = list()
        sCase = "{:0d}".format(iCase_index - 50)
        sText = 'Case ' + sCase + ' - ' + 'Observation' 
        aLegend.append(sText)
        
    
    

        elm_scatterplot_variable_x_2d(oE3SM,\
                                         oCase_x,\
                                         aData_y_in, \
                                            sVariable_y_in,\
                                         iFlag_scientific_notation_x_in=0,\
                                         iFlag_scientific_notation_y_in=0,\
                                          iFlag_log_x_in= 0,\
                                        iFlag_log_y_in = 0,\
                                         dMin_x_in = dMin_x, \
                                         dMax_x_in = dMax_x, \
                                         dMin_y_in = dMin_y, \
                                         dMax_y_in = dMax_y, \
                                         dSpace_x_in = 2, \
                                         dSpace_y_in = 2, \
                                         sLabel_x_in = sLabel_x, \
                                         sLabel_y_in = sLabel_y,\
                                         sLabel_legend_in = sLabel_legend )
        pass

print('finished')




#read slope from a case