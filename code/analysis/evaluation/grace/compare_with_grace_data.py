#this script should be run using Python 2.7.8 instead of Python 3
#module load python/2.7.8
#maybe I was wrong? 20200305 Chang Liao (chang.liao@pnnl.gov)

import os
import datetime
import numpy as np
from os.path import realpath
import calendar
import scipy.ndimage as ndimage
from pyearth.system.define_global_variables import * 

from pyearth.gis.gdal.read.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.toolbox.reader.text_reader_string import text_reader_string

from pyearth.gis.gdal.read.gdal_read_envi_file import gdal_read_envi_file_multiple_band
from pyearth.gis.location.find_index_by_latlon import find_index_by_latlon
from pyearth.gis.location.convert_lat_lon_range import convert_360_to_180

from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data
from pyearth.visual.color.create_qualitative_rgb_color_hex import create_qualitative_rgb_color_hex
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase


from pye3sm.elm.mesh.elm_retrieve_case_dimension_info import elm_retrieve_case_dimension_info
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file


sPath_parent = str(Path(__file__).parents[4]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_figure = realpath( sPath_parent +  '/figures/' )


iFlag_same_grid = 1
sModel = 'e3sm'
sRegion = 'amazon'
sDate = '20220701'

iFlag_pixel = 1


iYear_start = 2000
iYear_end = 2009
dConversion=1.0
#iCase_index_start = iIndex_start
##iCase_index_end = iIndex_end
#aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)

aCase_index = np.array([59, 61])
#iCase_index = 240

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)

aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                   iCase_index_in =  61 ,\
                                                   iYear_start_in = iYear_start, \
                                                   iYear_end_in = iYear_end,\
                                                    iYear_subset_start_in = 2004, \
                                                     iYear_subset_end_in = 2008, \
                                                    dConversion_in= dConversion,\
                                                   sDate_in= sDate,\
                                                   sModel_in = sModel, \
                                                       sRegion_in=sRegion)
oCase_x = pycase(aParameter_case)
aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                   iCase_index_in =  59 ,\
                                                   iYear_start_in = iYear_start, \
                                                   iYear_end_in = iYear_end,\
                                                    iYear_subset_start_in = 2004, \
                                                     iYear_subset_end_in = 2008, \
                                                    dConversion_in= dConversion,\
                                                   sDate_in= sDate,\
                                                   sModel_in = sModel, \
                                                       sRegion_in=sRegion)
oCase_y = pycase(aParameter_case)

sModel = oCase_x.sModel
sRegion = oCase_x.sRegion
iFlag_same_grid = oCase_x.iFlag_same_grid
iYear_start = oCase_x.iYear_start
iYear_end = oCase_x.iYear_end
iYear_subset_start = oCase_x.iYear_subset_start
iYear_subset_end = oCase_x.iYear_subset_end
sVariable = oCase_x.sVariable
sCase = oCase_x.sCase
sWorkspace_analysis_case = oCase_x.sWorkspace_analysis_case
sWorkspace_analysis_case_domain = sWorkspace_analysis_case + slash + 'tsplot_tws'
if not os.path.exists(sWorkspace_analysis_case_domain):
    os.makedirs(sWorkspace_analysis_case_domain)
    pass    
#new approach
aLon, aLat , aMask_ll= elm_retrieve_case_dimension_info(oCase_x)
#dimension
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
iRow_index_top = int( ( 90-0.5*dResolution_y  - dLat_max ) /  dResolution_y )
iRow_index_bot = int(( 90-0.5*dResolution_y  - dLat_min ) /  dResolution_y )
iColumn_index_left = int( (dLon_min - ( -180+0.5*dResolution_x   ) ) / dResolution_x)
iColumn_index_right = int( (dLon_max - ( -180+0.5*dResolution_x   ) ) / dResolution_x)

dLon = 299.777709960938 
dLat = -2.62420988082886 

dLon1 = convert_360_to_180(dLon)
row_index, column_index  = find_index_by_latlon(dLon1, dLat, dLon_min,dLat_max,dResolution_x, dResolution_y, iFlag_center_in= 1)


aDate = list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dSimulation = datetime.datetime(iYear, iMonth, 28)
        aDate.append( dSimulation )
nstress = nyear * nmonth
iMonth = 1
index_start = (iYear_subset_start - iYear_start)* 12 + iMonth - 1
index_end = (iYear_subset_end + 1 - iYear_start)* 12 + iMonth - 1
subset_index = np.arange(index_start , index_end , 1 )
aDate=np.array(aDate)
aDate_subset = aDate[subset_index]
nstress_subset= len(aDate_subset)
#read greace date list
sFilename_grace_lut = '/qfs/people/liao313/data/h2sc/global/auxiliary/grace/data_date.txt'
dummy = text_reader_string(sFilename_grace_lut, cDelimiter_in=',')
grace_date0 = dummy[:,0]
grace_date1 = dummy[:,1]
grace_date2 = dummy[:,2]
grace_date3 = dummy[:,3]
grace_date4 = dummy[:,4]
aData_grace=np.full( (nstress_subset, nrow_extract, ncolumn_extract), -9999, dtype=float )
sFilename_pre = 'GRD-3_'
sFilename_suf = '_GRAC_JPLEM_BA01_0600_LND_v03.tif'
sWorkspace_grace  = '/compyfs/liao313/00raw/grace/RL06/v03/JPL'
iStress = 1
for iYear in np.arange( iYear_subset_start, iYear_subset_end +1):
    for iMonth in np.arange(1, 13):
        sYear = "{:04d}".format(iYear)
        sYear2 = sYear[2:4]
        sMonth =  "{:02d}".format(iMonth)
        sMonth2 = calendar.month_abbr[iMonth]
        sDate  = sMonth2 + sYear2
        dummy_index = np.where( grace_date0 == sDate )
        sDummy = grace_date4[dummy_index][0]
        sDummy1 = sDummy.strip()
        sDummy2 = sYear + sDummy1[5:8] + '-' + sYear + sDummy1[-3:]
        sFilename_grace = sWorkspace_grace + slash \
            + sFilename_pre + sDummy2 + sFilename_suf
        print(sFilename_grace)
        #read the file
        aData_dummy0  = gdal_read_geotiff_file(sFilename_grace)
        aData_dummy0 = aData_dummy0[0]
        #resample the data because resolution is different
        aData_dummy1 = np.flip(aData_dummy0, 0) #flip the data first
        aData_dummy1 = np.roll(aData_dummy1, 180, axis=1) # right
        #because the dimensions, we use a simple way to resample
        #the grace resolution is 1.0 degree
        
        zoom_level = (1.0 / dResolution_x )
        aData_dummy2 = ndimage.zoom(aData_dummy1, zoom_level, order=0)
        aData_dummy2[np.where(aData_dummy2==-99999)] = -9999
        aData_dummy2[np.where(aData_dummy2==-9999)] = np.nan
        #extract 
        aData_dummy3 = aData_dummy2[ iRow_index_top:iRow_index_bot+1, iColumn_index_left:iColumn_index_right+1 ]
        
        aData_grace[iStress - 1,:,:]=aData_dummy3
        iStress=iStress+1


#read the stack data for each variable
aVariable = ['rain','snow', 'qdrai','qvegt',  'qvege','qsoil', 'qover']
nvariable = len(aVariable)
aVariable_tws = np.full(( nvariable, nstress_subset, nrow, ncolumn ),-9999, dtype=float)
#aVariable_begin_end =['tws_month_begin','tws_month_end']
for i in np.arange(1, nvariable+1):
    sVariable = aVariable[i-1]
    sWorkspace_variable_dat = sWorkspace_analysis_case + slash + sVariable +  slash + 'dat'
    sFilename = sWorkspace_variable_dat + slash + sVariable  + sExtension_envi
    aData_all = gdal_read_envi_file_multiple_band(sFilename)
    aVariable_total = aData_all[0]
    aVariable_tws[i-1, :, :, :] = aVariable_total[subset_index,:,:]
aVariable_all_x = np.full((nvariable, nstress_subset), -9999, dtype=float)
for i in np.arange(1, nvariable+1):
    sVariable = aVariable[i-1]
    aVariable_total_subset = aVariable_tws[i-1, :, : , :]
    aVariable0 = aVariable_total_subset.reshape(nstress_subset,nrow , ncolumn)
    aVariable2 = np.full(nstress_subset, -9999, dtype=float)
    for iStress in  np.arange(1,nstress_subset+1):
        dummy = aVariable0[iStress-1, :,:]
        if iFlag_pixel ==1:

            aVariable2[iStress-1] = dummy[ row_index, column_index]
            pass
        else:
            dummy1 = dummy[aMask_ul_index]
            aVariable2[iStress-1] = np.nanmean(dummy1)

    aVariable_all_x[i-1, :] = aVariable2
    
iStress = 1
aVariable_grace = np.full(nstress_subset, -9999, dtype=float)
aVariable_grace2 = np.full(nstress_subset, -9999, dtype=float)

for iStress in np.arange(1,nstress_subset+1):
    dummy = aData_grace[iStress-1, :,:]

    if iFlag_pixel ==1:
        aVariable_grace[iStress-1] = dummy[row_index, column_index]
    else:
        dummy1 = dummy[aMask_ul_index]
        dummy1[dummy1==-9999] = np.nan

        aVariable_grace[iStress-1] = np.nanmean(dummy1)
    #use regional mean instead of grid

iStress=1
for iYear in range(iYear_subset_start, iYear_subset_end +1,1):
    dummy_index = np.arange(  (iYear-iYear_subset_start )*12,  (iYear-iYear_subset_start )*12 + 12 ,1  )
    dummy_reference  = np.mean(  aVariable_grace[dummy_index] )
    for iMonth in np.arange(1, 13):

        aVariable_grace2[iStress-1] = aVariable_grace[iStress-1] - dummy_reference
        iStress=iStress+1



#be careful
aTWS_ts_x = aVariable_all_x[0,:]+ aVariable_all_x[1,:] \
    - (aVariable_all_x[2,:]+ aVariable_all_x[3,:] + aVariable_all_x[4,:])\
    - (aVariable_all_x[5,:]+ aVariable_all_x[6,:])

sWorkspace_analysis_case = oCase_y.sWorkspace_analysis_case
for i in np.arange(1, nvariable+1):
    sVariable = aVariable[i-1]
    sWorkspace_variable_dat = sWorkspace_analysis_case + slash + sVariable +  slash + 'dat'
    sFilename = sWorkspace_variable_dat + slash + sVariable  + sExtension_envi
    aData_all = gdal_read_envi_file_multiple_band(sFilename)
    aVariable_total = aData_all[0]
    aVariable_tws[i-1, :, :, :] = aVariable_total[subset_index,:,:]
aVariable_all_y = np.full((nvariable, nstress_subset), -9999, dtype=float)
for i in np.arange(1, nvariable+1):
    sVariable = aVariable[i-1]
    aVariable_total_subset = aVariable_tws[i-1, :, : , :]
    aVariable0 = aVariable_total_subset.reshape(nstress_subset,nrow , ncolumn)
    aVariable2 = np.full(nstress_subset, -9999, dtype=float)
    for iStress in  np.arange(1,nstress_subset+1):
        dummy = aVariable0[iStress-1, :,:]
        dummy1 = dummy[aMask_ul_index]
        aVariable2[iStress-1] = np.nanmean(dummy1)

    aVariable_all_y[i-1, :] = aVariable2

aTWS_ts_y = aVariable_all_y[0,:]+ aVariable_all_y[1,:] \
    - (aVariable_all_y[2,:]+ aVariable_all_y[3,:] + aVariable_all_y[4,:])\
    - (aVariable_all_y[5,:]+ aVariable_all_y[6,:])

#combined
aDate_ts = [aDate_subset, aDate_subset, aDate_subset ]


aData_ts = [aTWS_ts_x *31* 24*3600/1000, aTWS_ts_y *31* 24*3600/1000, aVariable_grace2]
sLabel_legend = sRegion.title()
sFilename_out = sWorkspace_analysis_case_domain + slash \
    + 'tws_tsplot_' + sRegion +'.png'
sFilename_out=sWorkspace_figure + slash \
                +  'tws_tsplots_grace_pixel' +'.png' 
sLabel_Y = r'Total water storage variation (m)'

sTitle='Amazon'

aColor = ['blue','red','black']
plot_time_series_data(aDate_ts, aData_ts,\
                          sFilename_out,\
                          ncolumn_in=3,\
                          dMax_y_in = 0.5,\
                          dMin_y_in = -0.5,\
                          sTitle_in = sTitle, \
                          sLabel_y_in= sLabel_Y,\
                          aLabel_legend_in = ['ELM-Default','ELM-HLGF', 'GRACE'], \
                          aMarker_in=['+','*','o'],\
                            aColor_in=aColor,\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)
print("finished")
