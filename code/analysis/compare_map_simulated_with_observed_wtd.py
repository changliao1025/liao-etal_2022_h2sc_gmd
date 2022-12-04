import numpy as np

from pyearth.system.define_global_variables import *   
from pyearth.gis.gdal.read.gdal_read_geotiff_file import gdal_read_geotiff_file


from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.structured.twod.map.elm_map_variable_difference_w_observation_2d import elm_map_variable_difference_w_observation_2d
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
#read observation
sFilename_in = '/qfs/people/liao313/data/e3sm/amazon/elm/' + 'wtd_extract' + sExtension_tiff
    
aData_obs, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatial_reference = gdal_read_geotiff_file(sFilename_in)
wtd_max = 45.0
aData_obs[np.where(aData_obs>wtd_max )]=wtd_max

dLon_min = dOriginX
dResolution_x= dPixelWidth
dLon_max = dLon_min + ncolumn * dPixelWidth
dLat_min =  dOriginY - nrow * dPixelWidth
dLat_max= dOriginY
aImage_extent= [dLon_min ,dLon_max , dLat_min ,  dLat_max]
#read simulated and apply average

sDate = '20220701'


iIndex_start=51
iIndex_end=58
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end

iFlag_scientific_notation_colorbar_in = 0

aVariable = ['ZWT']
aFlag_scientific_notation_colorbar=[0,1,1,1]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'


aTitle= [ 'Water table depth differences' ]
aUnit = [r'Unit: m']
aData_min = [-15]
aData_max = [15]
aConversion = [1]

sExtend='both'
sFormat_contour='%.2f'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
nvariable = len(aVariable)
ncase= len(aCase_index)
for i in range(0,ncase):

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

        oCase_x = pycase(aParameter_case)


        aLegend = list()
        sCase = "{:0d}".format(iCase_index - 50)
        sText = 'Case ' + sCase + ' - ' + 'Observation' 
        aLegend.append(sText)
        elm_map_variable_difference_w_observation_2d(oE3SM, oCase_x , aData_obs, dData_min_in=dData_min, dData_max_in=dData_max,\
            iFlag_scientific_notation_colorbar_in = iFlag_scientific_notation_colorbar, \
                sExtend_in = sExtend,\
                    sFormat_contour_in=sFormat_contour,\
            sUnit_in = sUnit,\
            sTitle_in=  sTitle,
            aLegend_in = aLegend )
        
        print(iCase_index)
        pass

print('finished')

