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

sWorkspace_out = ''

#find the orginal forcing data path using a case
#use the following two files as examples:
#sFilename_user_datm_prec = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1.Precip'
#sFilename_user_dlnd = '/qfs/people/liao313/data/e3sm/sag/mosart/dlnd.streams.txt.lnd.gpcc'

#this first file only change the path to the new prep data

#the second one is not used in coupled cases

#read parflow data

#set up group becaues elm use several file for different variables



def replace_temp_forcing(oCase_in):
    sFolder_out  = '/compyfs/liao313/00raw/prec/parflow/tpqw'
    isExist = os.path.exists(sFolder_out)
    if not isExist:    
       # Create a new directory because it does not exist
       os.makedirs(sFolder_out)
    sVariable_forcing_in = 'tpqw'
    iYear_start_data = 1950
    iYear_end_data = 2010

    iYear_start_avail = 2002
    iYear_end_avail = 2005

    sFolder_origin, aField, aFilename = atm_retrieve_forcing_data_info (oCase_in, sVariable_forcing_in)
   
    sField1=aField[0]
    sField2=aField[1]
    sField3=aField[2]
    sField4=aField[3]
    sField5=aField[4]
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
                    if (sKey == sField1):                    
                        aData1 = (aValue[:]).data
                        continue  
                    if (sKey == sField2):                    
                        aData2 = (aValue[:]).data
                        continue  
                    if (sKey == sField3):                    
                        aData3 = (aValue[:]).data
                        continue  
                    if (sKey == sField4):                    
                        aData4 = (aValue[:]).data
                        continue  
                    if (sKey == sField5):                    
                        aData5 = (aValue[:]).data
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
           


            aData1 = np.flip(aData1, 1)
            aData2 = np.flip(aData2, 1)
            aData3 = np.flip(aData3, 1)
            aData4 = np.flip(aData4, 1)
            aData5 = np.flip(aData5, 1)

            #read new data
            dummy = sDate_origin+'*'
            sRegex=  os.path.join( sWorkspace_parflow, dummy )

            for sFilename in glob.glob(sRegex):
                aDatasets = Dataset(sFilename)
                for sKey, aValue in aDatasets.variables.items():
                   
                    if (sKey == 'TBOT'):                    
                        aData11 = (aValue[:]).data
                        continue
                    if (sKey == 'WIND'):                    
                        aData21 = (aValue[:]).data
                        continue
                    if (sKey == 'QBOT'):                    
                        aData31 = (aValue[:]).data
                        continue
                    if (sKey == 'PSRF'):                    
                        aData41 = (aValue[:]).data
                        continue
                    if (sKey == 'FLDS'):                    
                        aData51 = (aValue[:]).data
                        continue
                    
                                
                aData_out_extract1 =  copy.copy(aData1)
                aData_out_extract2 =  copy.copy(aData2)
                aData_out_extract3 =  copy.copy(aData3)
                aData_out_extract4 =  copy.copy(aData4)
                aData_out_extract5 =  copy.copy(aData5)

                aData11 = np.reshape(aData11, aData11.size)
                aData21 = np.reshape(aData21, aData21.size)
                aData31 = np.reshape(aData31, aData31.size)
                aData41 = np.reshape(aData41, aData41.size)
                aData51 = np.reshape(aData51, aData51.size)

                for iStep in range(nts):

                    dummy01 = aData1[iStep, :,:]
                    dummy02 = aData2[iStep, :,:]
                    dummy03 = aData3[iStep, :,:]
                    dummy04 = aData4[iStep, :,:]
                    dummy05 = aData5[iStep, :,:]

                    iStep_start = iStep * 3
                    iStep_end = (iStep+1) * 3 

                    dummy1 = np.mean(aData11[iStep_start: iStep_end])
                    dummy2 = np.mean(aData21[iStep_start: iStep_end])
                    dummy3 = np.mean(aData31[iStep_start: iStep_end])
                    dummy4 = np.mean(aData41[iStep_start: iStep_end])
                    dummy5 = np.mean(aData51[iStep_start: iStep_end])
                    

                    dummy = dummy01.reshape(nrow_full, ncolumn_full)       
                    dummy[row_index, column_index] = dummy1              
                    aData_out_extract1[iStep,:,:] = dummy

                    dummy = dummy02.reshape(nrow_full, ncolumn_full)       
                    dummy[row_index, column_index] = dummy2              
                    aData_out_extract2[iStep,:,:] = dummy

                    dummy = dummy03.reshape(nrow_full, ncolumn_full)       
                    dummy[row_index, column_index] = dummy3              
                    aData_out_extract3[iStep,:,:] = dummy

                    dummy = dummy04.reshape(nrow_full, ncolumn_full)       
                    dummy[row_index, column_index] = dummy4              
                    aData_out_extract4[iStep,:,:] = dummy

                    dummy = dummy05.reshape(nrow_full, ncolumn_full)       
                    dummy[row_index, column_index] = dummy5              
                    aData_out_extract5[iStep,:,:] = dummy

                #save to a new location
                aData_out_extract1 = np.flip(aData_out_extract1, 1)
                aData_out_extract2 = np.flip(aData_out_extract2, 1)
                aData_out_extract3 = np.flip(aData_out_extract3, 1)
                aData_out_extract4 = np.flip(aData_out_extract4, 1)
                aData_out_extract5 = np.flip(aData_out_extract5, 1)

                sBasename = os.path.basename(sFilename_old)
                sBasename_new = sBasename.replace(sYear, sYear_origin, 1)
                sFilename_new = os.path.join(sFolder_out,sBasename_new )

                aData_out_extract= [aData_out_extract1, aData_out_extract2, aData_out_extract3, aData_out_extract4, aData_out_extract5 ]
                
                
                replace_variable_in_netcdf(sFilename_old, sFilename_new, aData_out_extract, aField)
                
                   

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


    replace_temp_forcing(oCase_y_in)