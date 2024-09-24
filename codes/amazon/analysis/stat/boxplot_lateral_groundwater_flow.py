from pathlib import Path
import numpy as np
from netCDF4 import Dataset #read netcdf
from pathlib import Path
from os.path import realpath
from pyearth.system.define_global_variables import * 
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pyearth.visual.color.create_qualitative_rgb_color_hex import create_qualitative_rgb_color_hex
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pye3sm.elm.general.structured.twod.retrieve.elm_retrieve_variable_2d import elm_retrieve_variable_2d

from pyearth.visual.boxplot.boxplot_data import boxplot_data

# getting the data
sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/' )
sWorkspace_figure = realpath( sPath_parent +  '/figures/' )

#read slope first

sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
sFilename_initial = '/compyfs/liao313/e3sm_scratch/e3sm20220701052/run/e3sm20220701052.elm.h0.1981-01.nc'
sFilename = sFilename_initial
sVariable='sur_slp'
aDatasets = Dataset(sFilename)
for sKey, aValue in aDatasets.variables.items():
    if sVariable == sKey.lower():
                       
        aData_ll = (aValue[:]).data    
        break
    else:
        print(sKey.lower())

nrow = 52
ncolumn = 58
aData_ll = np.reshape(aData_ll, (nrow, ncolumn))
aData_ul = np.flip(aData_ll, 0)  
aPercentiles_in = np.array([5,47.5,52.5, 95])
aData_ul[np.where(aData_ul == np.max(aData_ul))] = -9999

aInterval = cgpercentiles(aData_ul, aPercentiles_in, missing_value_in = -9999)  

aIndex_flat = np.where( (aData_ul <= aInterval[0]) & (aData_ul!=-9999) )

aIndex_moderate = np.where( (aData_ul < aInterval[2]) & ( aData_ul > aInterval[1] )& (aData_ul!=-9999) )

aIndex_mountain = np.where( (aData_ul >= aInterval[3]) & (aData_ul!=-9999) )



sDate = '20220701'
iFlag_log=0

iIndex_start=51
iIndex_end=59
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


aTitle= [ r'Subsurface runoff' ]
aUnit = [ r'Units: mm/s']
aData_min = [-1E-5]
aData_max = [1E-5]
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

aData_all = list()
aLabel_x_in=list()
for i in range(ncase):
    iCase_index =  aCase_index[i]
    sCase = "{:0d}".format(iCase_index - 50)
    sText = 'Case ' + sCase 
    aLabel_x_in.append(sText)
    
    
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
        aData_out_x = elm_retrieve_variable_2d(  oCase_x,\
                iFlag_monthly_in =0,\
                    iFlag_annual_mean_in=0,\
                iFlag_annual_total_in= 1     )
  
        dummy_x = aData_out_x.pop()

        nan_index = np.where(dummy_x ==  -9999)

        if iFlag_log  == 1:
            dummy_x = np.log10(dummy_x)            

        dummy_x[nan_index] = -9999
        dummy_flat = dummy_x[aIndex_flat]
        dummy_moderate = dummy_x[aIndex_moderate]
        dummy_mountain = dummy_x[aIndex_mountain]

        good_index = np.where( np.isinf(  dummy_flat) == False  )
        dummy_flat = dummy_flat[good_index] 

        good_index = np.where( np.isinf(  dummy_moderate) == False  )
        dummy_moderate = dummy_moderate[good_index] 

        good_index = np.where( np.isinf(  dummy_mountain) == False  )
        dummy_mountain = dummy_mountain[good_index]       
        
        print( np.mean(dummy_flat), np.std(dummy_flat) )
        print( np.mean(dummy_moderate), np.std(dummy_moderate) )
        print( np.mean(dummy_mountain), np.std(dummy_mountain) )

        aData_group=list()
        aData_group.append(dummy_flat)
        aData_group.append(dummy_moderate)
        aData_group.append(dummy_mountain)
        aData_all.append(aData_group)
        
        pass

#get 
aLabel_y_in=None

sFilename_out = sWorkspace_figure + slash \
                + sVariable.lower() + '_boxplot_annual_total' +'.png'      

aColor_in = create_qualitative_rgb_color_hex(3)

sFormat_x='{:.1f}' 
sLabel_y=r'Subsurface runoff (mm/s)'
aHatch=['/', '\\',  'x']

aLabel_legend_in =['Flat region' ,  'Moderate slope' ,  'Mountainous region']

boxplot_data(aData_all,\
                 aLabel_x_in,\
                 aLabel_y_in,\
                 sFilename_out,\
                 iDPI_in = None,\
                 ncolumn_in = 5,\
                 iSize_x_in = None, \
                 iSize_y_in = None, \
                 dMax_y_in = None, \
                 dMin_y_in = None, \
                 dSpace_y_in = None,\
                 aMarker_in = None,\
                 aColor_in = aColor_in,\
                 aHatch_in = aHatch,\
                 sLabel_y_in = sLabel_y, \
                 aLabel_legend_in = aLabel_legend_in,\
                 sLocation_legend_in = 'lower right' ,\
                 aLocation_legend_in = (1.0, 0.0),\
                 sFormat_y_in =None,\
                 sTitle_in = None)


print('finished')