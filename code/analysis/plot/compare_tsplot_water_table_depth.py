#python built in package
import numpy as np
import datetime
from pathlib import Path
from os.path import realpath

#support package
from pyearth.system.define_global_variables import * 
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.color.create_qualitative_rgb_color_hex import create_qualitative_rgb_color_hex
from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data
from pyearth.visual.timeseries.plot_time_series_data_w_variation import plot_time_series_data_w_variation

#pye3sm package
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pye3sm.elm.grid.elm_retrieve_case_dimension_info import elm_retrieve_case_dimension_info
from pye3sm.elm.general.structured.twod.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d

# getting the data
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_figure = realpath( sPath_parent +  '/figures/' )

sDate = '20220701'
iIndex_start=51
iIndex_end=59
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end

iFlag_scientific_notation_colorbar_in = 0

aVariable = ['ZWT']
aFlag_scientific_notation_colorbar=[1]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'

iFlag_monthly = 1
iFlag_annual_mean = 0
iFlag_log = 0

iFlag_median = 0

aTitle= [ r'Water table depeth' ]
aUnit = [ r'm']
aData_min = [0]
aData_max = [8]
aConversion = [1]
sColormap = 'Spectral'



sExtend='both'

sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
aCase_index = np.arange(iCase_index_start, iCase_index_end + 1, 1)
ncase= len(aCase_index)
nvariable = len(aVariable)

aDate_all =      list()
aData_all =      list()    
aData_upper_all= list()
aData_lower_all= list()
aLegend = list()
for i in range(ncase):
    iCase_index =  aCase_index[i]
    sCase = "{:0d}".format(iCase_index - 50)
    sText = 'Case ' + sCase
    aLegend.append(sText)
    for iVariable in np.arange(0, nvariable):
        sVariable = aVariable[iVariable]
        sUnit = aUnit[iVariable]
        sTitle = aTitle[iVariable]
        sTitle = '' #no need for title as y label is used
        dMin_y_in = aData_min[iVariable]
        dMax_y_in = aData_max[iVariable]
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
        aLon, aLat, aMask =elm_retrieve_case_dimension_info(oCase_in)
        #dimension
        nrow = np.array(aMask).shape[0]
        ncolumn = np.array(aMask).shape[1]
        aMask = np.where(aMask==0)
        iYear_start = oCase_in.iYear_start
        iYear_end = oCase_in.iYear_end    
        iYear_subset_start = oCase_in.iYear_subset_start
        iYear_subset_end = oCase_in.iYear_subset_end
        dates = list()
        dates_year=list()
        nyear = iYear_end - iYear_start + 1
        for iYear in range(iYear_start, iYear_end + 1):
            dSimulation0 = datetime.datetime(iYear, 6, 30)
            dates_year.append( dSimulation0 )
            for iMonth in range(1,13):
                dSimulation = datetime.datetime(iYear, iMonth, 15)
                dates.append( dSimulation )

        nstress = nyear * nmonth
        #take the subset
        iMonth = 1
        subset_index_start = (iYear_subset_start - iYear_start) * 12 + iMonth-1
        subset_index_end = (iYear_subset_end + 1 - iYear_start) * 12 + iMonth-1
        subset_index = np.arange( subset_index_start,subset_index_end, 1 )

        dates=np.array(dates)
        dates_subset = dates[subset_index]
        nstress_subset= len(dates_subset)
        if iFlag_monthly ==1:
            aData_ret = elm_retrieve_variable_2d( oCase_in, iFlag_monthly_in = 1)        
            aDataTs = np.full(nstress_subset, -9999, dtype=float)
            aUpper=np.full(nstress_subset, -9999, dtype=float)
            aLower=np.full(nstress_subset, -9999, dtype=float)
            for i in np.arange(0, nstress_subset, 1):
                aImage = aData_ret[i]     
                dummy1 = np.reshape(aImage, (nrow, ncolumn))
                good_index = np.where(dummy1 != -9999)
                dummy1=dummy1[good_index]
                
                if iFlag_median ==1:
                    dRange = 25
                    aPercentiles_in = np.arange(50-dRange*0.5, 50+dRange*0.5 + 1, 25)
                    aInterval = cgpercentiles(dummy1, aPercentiles_in, missing_value_in = -9999) 
                    aDataTs[i] = np.median(dummy1) 
                    aUpper[i ] =aInterval[1]
                    aLower[i ] =aInterval[0]
                else:
                    dummy0 = np.std(dummy1)    
                    aDataTs[i] = np.mean(dummy1) 
                    aUpper[i ] =aDataTs[i] + 0.5*dummy0
                    aLower[i ] =aDataTs[i] - 0.5*dummy0

            if iFlag_log  == 1:
                aDataTs = np.log10(aDataTs)
                aUpper = np.log10(aUpper)
                aLower = np.log10(aLower)

                #set inf to min
                bad_index = np.where( np.isinf(  aDataTs) == True  )
                aDataTs[bad_index] = dMin_y_in            

            aDate_all.append(dates_subset)
            aData_all.append(aDataTs) 
            aData_upper_all.append(aUpper)
            aData_lower_all.append(aLower)

            

sFilename_out = sWorkspace_figure + slash \
                + sVariable.lower() + '_tsplots_monthly' +'.png'        

aColor_in = create_qualitative_rgb_color_hex(ncase)
sFormat_y='{:.1f}' 
sLabel_y=r'Water table depth (m)'
plot_time_series_data(aDate_all,
                          aData_all,\
                          sFilename_out,\
                          iReverse_y_in = 1, \
                          iFlag_log_in = 0,\
                          iFlag_scientific_notation_in=0,\
                          iFlag_miniplot_in = 1,\
                          ncolumn_in =int(5),\
                          dMax_y_in = dMax_y_in,\
                          dMin_y_in = dMin_y_in,\
                          dSpace_y_in = 2, \
                          sTitle_in = sTitle, \
                          sLabel_y_in= sLabel_y,\
                          sFormat_y_in= sFormat_y ,\
                          aLabel_legend_in = aLegend, \
                          aColor_in =aColor_in,\
                          aMarker_in = None,\
                          sLocation_legend_in = 'lower left' ,\
                          aLocation_legend_in = (0.0, 1.0),\
                          aLinestyle_in = None,\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)

sFilename_out = sWorkspace_figure + slash \
                + sVariable.lower() + '_tsplots_monthly_w_variation' +'.png' 
                
plot_time_series_data_w_variation(aDate_all,
                          aData_all,\
                          aData_upper_all,\
                          aData_lower_all,\
                          sFilename_out,\
                          iReverse_y_in = 1, \
                          iFlag_log_in = 0,\
                            iFlag_miniplot_in = 1,\
                          iFlag_scientific_notation_in=0,\
                          ncolumn_in = int(5),\
                          dMax_y_in = dMax_y_in,\
                          dMin_y_in = dMin_y_in,\
                          dSpace_y_in = 2, \
                          sTitle_in = sTitle, \
                          sLabel_y_in= sLabel_y,\
                          sFormat_y_in= sFormat_y ,\
                          aLabel_legend_in = aLegend, \
                          aColor_in =aColor_in,\
                          aMarker_in = None,\
                          sLocation_legend_in = 'lower left' ,\
                          aLocation_legend_in = (0.0, 1.0),\
                          aLinestyle_in = None,\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)
print('finished')
