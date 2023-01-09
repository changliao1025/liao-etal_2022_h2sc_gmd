import numpy as np
from osgeo import gdal, osr #the default operator
from pyearth.system.define_global_variables import *   
from pyearth.gis.gdal.read.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.structured.twod.map.elm_map_variable_difference_w_observation_dc_2d import elm_map_variable_difference_w_observation_dc_2d
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pye3sm.elm.general.structured.twod.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d
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
iIndex_end=59
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
dData_min = -30
dData_max = 5
aConversion = [1]

sExtend='both'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
nvariable = len(aVariable)
ncase= len(aCase_index)
aData_all = list()

for i in range(ncase):

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
        oCase_y = pycase(aParameter_case)

        aData_out_y = elm_retrieve_variable_2d(  oCase_y,\
                iFlag_monthly_in =0,\
                    iFlag_annual_mean_in=1,\
                iFlag_annual_total_in= 0     )

        dummy_y = aData_out_y.pop()

        nan_index = np.where(aData_obs ==  -9999)

        dummy=dummy_y-aData_obs
        

        dummy_index = np.where(dummy > dData_max)
        dummy[dummy_index] = dData_max
        dummy_index = np.where(dummy < dData_min)
        dummy[dummy_index] = dData_min  
        dummy[nan_index] = -9999

        aData_all.append(dummy)
        
        pass

#get 

aPercentiles_in = np.arange(10, 90, 13)
aInterval = cgpercentiles(aData_all, aPercentiles_in, missing_value_in = -9999)  

#aInterval = np.linspace(dData_min, dData_max, num=7)
aColor = create_diverge_rgb_color_hex(  len(aInterval) +1  )
aColor.reverse()

for iCase_index in ( aCase_index ):
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
        dData_min = dData_min
        dData_max = dData_max
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
        elm_map_variable_difference_w_observation_dc_2d(oE3SM, oCase_x , aData_obs, dData_min_in=dData_min, dData_max_in=dData_max,\
            iFlag_scientific_notation_colorbar_in = iFlag_scientific_notation_colorbar, \
                  iFlag_annual_mean_in=1,\
                iFlag_annual_total_in= 0,\
                    aInterval_in = aInterval,\
                sExtend_in = sExtend,\
            sUnit_in = sUnit,\
            sTitle_in=  sTitle,
              aColor_in=aColor,\
            aLegend_in = aLegend )
        
        print(iCase_index)
        pass

print('finished')

