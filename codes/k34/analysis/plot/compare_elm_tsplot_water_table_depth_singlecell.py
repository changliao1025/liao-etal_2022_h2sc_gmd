import os
import numpy as np
import netCDF4 as nc #read netcdf

from datetime import datetime
from pyearth.system.define_global_variables import *

from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.elm.general.unstructured.singlecell.plot.elm_tsplot_variable_singlecell import elm_tsplot_variable_singlecell
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file


sDate = '20240101'
iYear_start = 2000
iYear_end = 2009
sModel = 'e3sm'
sRegion='k34'
aTime = list()
aData = list()
aLabel_legend = list()
iFlag_first = 1

aDate = list()

nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dSimulation = datetime(iYear, iMonth, 15)
        aDate.append( dSimulation )
        pass
    
sVariable = 'zwt'   
sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/case.xml'
for i in range(2, 6,1):
    aLabel_legend.append('Case ' + str(i))
    aTime.append(aDate)
    iCase_index = i    

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
                                                           sRegion_in = sRegion,
                                                           sVariable_in = sVariable )
    #print(aParameter_case)
    oCase = pycase(aParameter_case)  
    sWorkspace_simulation_case_run = oCase.sWorkspace_simulation_case_run
    sCase = oCase.sCase
    iYear_subset_start = oCase.iYear_subset_start
    iYear_subset_end = oCase.iYear_subset_end
    iMonth = 1
    index_start = (iYear_subset_start - iYear_start)* 12 + iMonth - 1
    index_end = (iYear_subset_end + 1 - iYear_start)* 12 + iMonth - 1
    subset_index = np.arange(index_start , index_end , 1 )
    aDate=np.array(aDate)
    aDate_subset = aDate[subset_index]
    nstress_subset= len(aDate_subset)
    #retrieve zwt 

    aData_out = np.full(nstress_subset,missing_value, dtype=float)

    iStress=1
    for iYear in range(iYear_start, iYear_end + 1):
        sYear = "{:04d}".format(iYear) #str(iYear).zfill(4)

        for iMonth in range(iMonth_start, iMonth_end + 1):
            sMonth = str(iMonth).zfill(2)

            sDummy = '.elm.h0.' + sYear + '-' + sMonth + sExtension_netcdf
            sFilename = sWorkspace_simulation_case_run + slash + sCase + sDummy

            #read before modification

            if os.path.exists(sFilename):
                print("Yep, I can read that file: " + sFilename)
                pass
            else:
                print(sFilename)
                print("Nope, the path doesn't reach your file. Go research filepath in python")
                quit()
                pass

            aDatasets = nc.Dataset(sFilename)

            #read the actual data
            for sKey, aValue in aDatasets.variables.items():
                if sVariable == sKey.lower():
                    #for attrname in aValue.ncattrs():
                    #print("{} -- {}".format(attrname, getattr(aValue, attrname)))
                    aData_tmp = (aValue[:]).data
                    #print(aData)
                    missing_value1 = np.max(aData_tmp)
                    aData_out[iStress-1]= aData_tmp
                    iStress= iStress + 1
                    pass

                else:
                    pass

    aData.append(aData_out)


#plot multiple case     
    
#repeat a numpy array n times


aTime = np.array(aTime)
aData = np.array(aData)

current_file_path = os.path.abspath(__file__)

current_file_directory = os.path.dirname(current_file_path)

sFilename_out = current_file_directory + slash + 'zwt.png'
iFlag_scientific_notation=0 
iReverse_y=1 
sTitle = 'Water table depth (m)'
sLabel_y = r'Water table depth (m)'


dMin_y=0
dMax_y=10


plot_time_series_data( aTime, aData,                      
                           sFilename_out,
                           iFlag_scientific_notation_in=iFlag_scientific_notation,
                           iReverse_y_in = iReverse_y,
                           sTitle_in = sTitle,
                           sLabel_y_in= sLabel_y,
                           dMax_y_in = None,
                           dMin_y_in =dMin_y,                           
                           aLabel_legend_in=aLabel_legend,
                           aLocation_legend_in = [1.0,0.0],
                           sLocation_legend_in='lower right',
                           iSize_x_in = 12,
                           iSize_y_in = 5)    
    
    
   
    


print('finished')
