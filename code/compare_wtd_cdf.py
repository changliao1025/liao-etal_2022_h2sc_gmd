import os
from pathlib import Path
import datetime
import numpy as np
from netCDF4 import Dataset #read netcdf
from osgeo import gdal, osr #the default operator
from pyearth.system.define_global_variables import *   
from pyearth.gis.gdal.read.gdal_read_envi_file import gdal_read_envi_file_multiple_band
from pyearth.gis.gdal.read.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.write.gdal_write_geotiff_file import gdal_write_geotiff_file
from pye3sm.elm.grid.elm_extract_grid_latlon_from_mosart import elm_extract_grid_latlon_from_mosart
from pyearth.visual.histogram.cdf_plot_multiple_data import cdf_plot_multiple_data
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.structured.twod.map.elm_map_variable_difference_w_observation_2d import elm_map_variable_difference_w_observation_2d
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


iIndex_start=51
iIndex_end=57
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


#read slope

sPath_parent = str(Path(__file__).parents[1]) # data is located two dir's up
sFilename_initial = '/compyfs/liao313/e3sm_scratch/e3sm20220701052/run/e3sm20220701052.elm.h0.1981-01.nc'
sFilename = sFilename_initial
sVariable='sur_slp'
aDatasets = Dataset(sFilename)
for sKey, aValue in aDatasets.variables.items():
    if sVariable == sKey.lower():
                       
        aData_slope = (aValue[:]).data    
        break
    else:
        print(sKey.lower())

aData_slope = np.reshape(aData_slope, (nrow, ncolumn))
index_high = np.where(aData_slope > 0.1)
index_low = np.where(aData_slope < 0.02)

sFilename_in = '/qfs/people/liao313/data/e3sm/amazon/elm/' + 'wtd_extract' + sExtension_tiff
    
aData_obs, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatial_reference = gdal_read_geotiff_file(sFilename_in)


aTitle= [ 'Water table depth' ]
aUnit = [r'Unit: m']
aData_min = [-5]
aData_max = [5]
aConversion = [1]

sExtend='both'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
nvariable = len(aVariable)
aData_high = list()
aData_low = list()
aLabel_legend = list()

aData_high.append(aData_obs[index_high])
aData_low.append(aData_obs[index_low])
aLabel_legend.append('Observation')

for iCase_index in ( aCase_index ):
    sCase_index =  "{:0d}".format(iCase_index-50)
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable].lower()
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
        iMonth = 1
        subset_index_start = (iYear_subset_start - iYear_start) * 12 + iMonth-1
        subset_index_end = (iYear_subset_end + 1 - iYear_start) * 12 + iMonth-1
        subset_index = np.arange( subset_index_start,subset_index_end, 1 )
        dates=np.array(dates)
        dates_subset = dates[subset_index]
        nstress_subset= len(dates_subset)
        sWorkspace_analysis_case = oCase_in.sWorkspace_analysis_case
        sWorkspace_variable_dat = sWorkspace_analysis_case + slash + sVariable +  slash + 'dat'
        sFilename = sWorkspace_variable_dat + slash + sVariable  + sExtension_envi

        aData_all = gdal_read_envi_file_multiple_band(sFilename)
        aVariable_total = aData_all[0]
        aVariable_total_subset = aVariable_total[subset_index,:,:]


        sWorkspace_analysis_case_variable = sWorkspace_analysis_case + slash + sVariable
        if not os.path.exists(sWorkspace_analysis_case_variable):
            os.makedirs(sWorkspace_analysis_case_variable)

        sWorkspace_analysis_case_region = sWorkspace_analysis_case_variable + slash + 'map'
        if not os.path.exists(sWorkspace_analysis_case_region):
            os.makedirs(sWorkspace_analysis_case_region)
            pass

        
        aImage_x = np.mean(aVariable_total_subset, axis=0 )   
       
            
        sLabel_legend = 'Case ' + sCase_index
            
        aData_high.append(aImage_x[index_high])

        aData_low.append(aImage_x[index_low])

        aLabel_legend.append(sLabel_legend)
            
        
        
        print(iCase_index)
        pass


sLabel_x = 'Water table depth (m)'
sLabel_y = 'Likelihood'

sFilename_out= sPath_parent + '/' + 'figures' + '/' + 'high_slope_cdf.png'
dMin_x = 0
dMax_x =20
sTitle= 'Slope > 0.1'
cdf_plot_multiple_data(aData_high, \
    sFilename_out, \
    dMin_x_in = dMin_x, \
    dMax_x_in = dMax_x, \
    sLabel_x_in = sLabel_x, \
    sLabel_y_in = sLabel_y, \
    sTitle_in = sTitle,\
        aLabel_legend_in = aLabel_legend)

dMin_x = 0
dMax_x =20
sTitle= 'Slope < 0.02'
sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'low_slope_cdf.png'
cdf_plot_multiple_data(aData_low, \
    sFilename_out, \
    dMin_x_in = dMin_x, \
    dMax_x_in = dMax_x, \
    sLabel_x_in = sLabel_x, \
    sLabel_y_in = sLabel_y, \
    sTitle_in = sTitle,\
        aLabel_legend_in = aLabel_legend)

print('finished')

