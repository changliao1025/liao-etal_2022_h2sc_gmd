import os, sys
import shutil
import datetime
import glob
import copy
from netCDF4 import Dataset #read netcdf
import numpy as np
from pyearth.system.define_global_variables import *
from pyearth.toolbox.date.day_in_month import day_in_month

from pyearth.system.create_symlink import create_symlink

from pyearth.toolbox.data.beta.replace_variable_in_netcdf import replace_variable_in_netcdf

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pye3sm.atm.general.atm_retrieve_forcing_data_info import atm_retrieve_forcing_data_info
from pye3sm.elm.mesh.elm_retrieve_case_dimension_info import elm_retrieve_case_dimension_info
#the script is used to generate the forcing data that are the same with the parflow forcing data 
#it will replace existing forcing at certain cells with the parflow forcing

lat = -2.62420988082886 
lon = 299.777709960938 
sWorkspace_parflow ='/qfs/people/lili400/datashare/k34_updt_0205'

sWorkspace_prec = '/compyfs/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/Precip'
sWorkspace_solar = '/compyfs/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/Solar'
sWorkspace_tpqw = '/compyfs/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/TPHWL'



#find the orginal forcing data path using a case
#use the following two files as examples:
#sFilename_user_datm_prec = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1.Precip'
#sFilename_user_dlnd = '/qfs/people/liao313/data/e3sm/sag/mosart/dlnd.streams.txt.lnd.gpcc'

#this first file only change the path to the new prep data

#the second one is not used in coupled cases

#read parflow data

#set up group becaues elm use several file for different variables



def replace_solar_forcing(oCase_in):
    sFolder_out  = '/compyfs/liao313/00raw/prec/parflow/solar'
    isExist = os.path.exists(sFolder_out)
    if not isExist:    
       # Create a new directory because it does not exist
       os.makedirs(sFolder_out)
       
    sVariable_forcing_in ='solar'
    iYear_start_data = 1950
    iYear_end_data = 2010

    iYear_start_avail = 2002
    iYear_end_avail =2005

    sFolder_origin, aField, aFilename = atm_retrieve_forcing_data_info (oCase_in, sVariable_forcing_in)
    sField = aField[0]
   
    aLon, aLat , aMask_ll= elm_retrieve_case_dimension_info(oCase_in)
    aLon = np.flip(aLon, 0) 
    aLat = np.flip(aLat, 0) 
    aMask = np.flip(aMask_ll, 0) 
    nrow_extract, ncolumn_extract = aLon.shape
    #resolution
    dLon_min = np.min(aLon)
    dLon_max = np.max(aLon)
    dLat_min = np.min(aLat)
    dLat_max = np.max(aLat)
    dResolution_x = (dLon_max - dLon_min) / (ncolumn_extract-1)
    dResolution_y = (dLat_max - dLat_min) / (nrow_extract-1)
    aImage_extent =  [dLon_min- dResolution_x ,dLon_max + dResolution_x, dLat_min -dResolution_x,  dLat_max+dResolution_x]
  
    dResoultion_elm = dResolution_x
   

    iYear_start = oCase_in.iYear_start
    iYear_end = oCase_in.iYear_end

    dates = list()
    dates_year =list()
    nyear = iYear_end - iYear_start + 1

    for iYear in range(iYear_start, iYear_end + 1):
        dSimulation0 = datetime.datetime(iYear, 6, 30)
        dates_year.append( dSimulation0 )
        for iMonth in range(1,13):
            dSimulation = datetime.datetime(iYear, iMonth, 15)
            dates.append( dSimulation )

    for iYear in range(iYear_start_data, iYear_end_data + 1, 1):
        sYear =  "{:04d}".format(iYear)
        iYear_origin = iYear
        sYear_origin =  "{:04d}".format(iYear_origin)
        if iYear < iYear_start_avail or  iYear > iYear_end_avail:
            #if earlier than available data, we just copy it
            for iMonth in range(1, 12 + 1, 1):
                sMonth =  "{:02d}".format(iMonth)
                sDate = sYear + '-' + sMonth
                sDate_origin = sYear_origin + '-' + sMonth
                dummy = '*'+sDate_origin+'*'
                sRegex_origin = os.path.join( sFolder_origin, dummy )
           

                dom = day_in_month(iYear, iMonth, iFlag_leap_year_in = 0)
                nts = dom * 8 #3 hour temporal resolution

                for sFilename_old in glob.glob(sRegex_origin):
                        
                    sBasename = os.path.basename(sFilename_old)
                    sBasename_new = sBasename.replace(sYear, sYear_origin, 1)
                    sFilename_new = os.path.join(sFolder_out,sBasename_new )
                    #shutil.copyfile(sFilename_old, sFilename_new)
                    create_symlink(sFilename_old, sFilename_new)
                    pass
            continue
        else:
            #we need to replace 
            pass
        
        #get the file by year
        for iMonth in range(1, 12 + 1, 1):
            sMonth =  "{:02d}".format(iMonth)
            sDate = sYear + '-' + sMonth
            sDate_origin = sYear_origin + '-' + sMonth
            dummy = '*'+sDate_origin+'*'
            sRegex_origin = os.path.join( sFolder_origin, dummy )
           

            dom = day_in_month(iYear, iMonth, iFlag_leap_year_in = 0)
            nts = dom * 8 #3 hour temporal resolution

            for sFilename in glob.glob(sRegex_origin):
                aDatasets = Dataset(sFilename)
                for sKey, aValue in aDatasets.variables.items():
                    if (sKey == 'LONGXY'):                   
                        aLongitude_origin = (aValue[:]).data
                        continue
                    if (sKey == 'LATIXY'):                    
                        aLatitude_origin = (aValue[:]).data
                        continue
                    if (sKey == sField):                    
                        aData0 = (aValue[:]).data
                        continue  

                sFilename_old = sFilename    

           
            dLon_min = np.min(aLongitude_origin)
            dLon_max = np.max(aLongitude_origin)
            dLat_min = np.min(aLatitude_origin)
            dLat_max = np.max(aLatitude_origin)
            nrow_full, ncolumn_full = aLongitude_origin.shape
            dResolution_x = (dLon_max - dLon_min) / (ncolumn_full-1)
            dResolution_y = (dLat_max - dLat_min) / (nrow_full-1)


            row_index = int((dLat_max+0.5*(dResolution_y)-lat)/dResolution_y)
            column_index = int( (lon-dLon_min-0.5*(dResolution_x))/dResolution_x )             


            aData0 = np.flip(aData0, 1)

            #read new data
            dummy = sDate_origin+'*'
            sRegex=  os.path.join( sWorkspace_parflow, dummy )

            for sFilename in glob.glob(sRegex):
                aDatasets = Dataset(sFilename)
                for sKey, aValue in aDatasets.variables.items():
                   
                    if (sKey == 'FSDS'):                    
                        aData1 = (aValue[:]).data
                        break
                                
                aData_out_extract =  copy.copy(aData0)
                aData1 = np.reshape(aData1, aData1.size)

                for iStep in range(nts):
                    iStep_start = iStep * 3
                    iStep_end = (iStep+1) * 3 

                    dummy1 = np.mean(aData1[iStep_start: iStep_end])
                    dummy0 = aData0[iStep, :,:]

                    dummy = dummy0.reshape(nrow_full, ncolumn_full)
                    
                    dummy[row_index, column_index] = dummy1
              
                    aData_out_extract[iStep,:,:] = dummy

                #save to a new location
                aData_out_extract = np.flip(aData_out_extract, 1)

                sBasename = os.path.basename(sFilename_old)
                sBasename_new = sBasename.replace(sYear, sYear_origin, 1)
                sFilename_new = os.path.join(sFolder_out,sBasename_new )
                
                replace_variable_in_netcdf(sFilename_old, sFilename_new, aData_out_extract, sField)
                
                   

            pass
        pass
    

if __name__ == '__main__':
    sFilename_case_configuration = '/qfs/people/liao313/workspace/python/pye3sm/pye3sm/case.xml'
    iCase_index_y = 52
    iYear_start = 2000
    iYear_end = 2009
    sDate = '20220701'
    sModel = 'e3sm'
    sRegion='amazon'

    aParameter_case  = pye3sm_read_case_configuration_file(sFilename_case_configuration,\
                                                       iCase_index_in =  iCase_index_y ,\
                                                       iYear_start_in = iYear_start, \
                                                       iYear_end_in = iYear_end,\
                                                        iYear_subset_start_in = iYear_start, \
                                                         iYear_subset_end_in = iYear_end, \
                                           
                                                       sDate_in= sDate,\
                                                       sModel_in = sModel, \
                                                           sRegion_in=sRegion )
    oCase_y_in = pycase(aParameter_case)


    replace_solar_forcing(oCase_y_in)