#python built in package
import numpy as np
import datetime
from pathlib import Path
from os.path import realpath

#support package
from pyearth.system.define_global_variables import * 
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.color.create_qualitative_rgb_color_hex import create_qualitative_rgb_color_hex
from pyearth.visual.histogram.histogram_plot import histogram_plot

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
iIndex_end=58
#start loop
iCase_index_start = iIndex_start
iCase_index_end = iIndex_end

iFlag_scientific_notation_colorbar_in = 0

aVariable = ['QDRAI']
aFlag_scientific_notation_colorbar=[1]
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='amazon'
iFlag_monthly = 0
iFlag_annual_mean=0
iFlag_annual_total =1
iFlag_log = 1 

iFlag_median = 0

aTitle= [ r'Subsurface runoff' ]
aUnit = [ r'Units: mm/s']
aData_min = [-8]
aData_max = [-3]
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

aData_all =      list() 
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
        dMin_x_in = aData_min[iVariable]
        dMax_x_in = aData_max[iVariable]
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
           
            for i in np.arange(0, nstress_subset, 1):
                aImage = aData_ret[i]     
                dummy1 = np.reshape(aImage, (nrow, ncolumn))
                good_index = np.where(dummy1 != -9999)  
                aData_monthly=dummy1[good_index]
              
            if iFlag_log  == 1:
                aData_monthly = np.log10(aData_monthly)
                #set inf to min
                bad_index = np.where( np.isinf(  aData_monthly) == True  )
                aData_monthly[bad_index] = dMin_x_in            

            aData_all.append(aData_monthly) 

        if iFlag_annual_total==1:
            aData_ret = elm_retrieve_variable_2d( oCase_in, iFlag_annual_total_in = 1)        
           
            for iYear in range(iYear_end, iYear_end + 1):
                aImage = aData_ret[iYear-iYear_start]
                sYear = "{:04d}".format(iYear)       
                dummy1 = np.reshape(aImage, (nrow, ncolumn))
                good_index = np.where(dummy1 != -9999)  
                aData_annual_total=dummy1[good_index]
              
            if iFlag_log  == 1:
                aData_annual_total = np.log10(aData_annual_total)
                #set inf to min
                bad_index = np.where( np.isinf(  aData_annual_total) == True  )
                aData_annual_total[bad_index] = dMin_x_in            

                bad_index = np.where(   aData_annual_total <  dMin_x_in  )
                aData_annual_total[bad_index] = dMin_x_in  

                bad_index = np.where(   aData_annual_total >  dMax_x_in  )
                aData_annual_total[bad_index] = dMax_x_in  

            aData_all.append(aData_annual_total) 


sFilename_out = sWorkspace_figure + slash \
                + sVariable.lower() + '_histogram_monthly' +'.png'      

aColor_in = create_qualitative_rgb_color_hex(ncase)
sFormat_x='{:.0f}' 
sLabel_x=r'Subsurface runoff (mm/s)'
sLabel_y='Density'

histogram_plot(  aData_all,\
                          sFilename_out,\
                          iFlag_log_in = 1,\
                          iFlag_scientific_notation_in=0,\
                          ncolumn_in =int(ncase/2),\
                          dMax_x_in = dMax_x_in,\
                          dMin_x_in = dMin_x_in,\
                          dSpace_x_in = 0.5, \
                          sTitle_in = sTitle, \
                          sLabel_x_in= sLabel_x,\
                            sLabel_y_in= sLabel_y,\
                          sFormat_x_in= sFormat_x ,\
                          aLabel_legend_in = aLegend, \
                          aColor_in =aColor_in,\
                          sLocation_legend_in = 'upper left' ,\
                          aLocation_legend_in = (0.0, 1.0),\
                          iSize_x_in = 12,\
                          iSize_y_in = 5)


                

print('finished')
