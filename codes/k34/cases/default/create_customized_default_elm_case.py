import os
from pathlib import Path
import numpy as np
import datetime
from shutil import copyfile
from pyearth.system.define_global_variables import *
from pyearth.toolbox.data.beta.add_variable_to_netcdf import add_multiple_variable_to_netcdf

from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pye3sm.case.e3sm_create_case import e3sm_create_case
from pye3sm.case.e3sm_choose_res_and_compset import e3sm_choose_res_and_compset
from pye3sm.elm.mesh.elm_create_customized_domain import elm_create_customized_domain



sModel = 'e3sm'
iFlag_region = 0


sRegion = 'k34'
#the lat/lon only used when in single grid case
#k34 site
dLatitude = -2.6091
dLongitude = -60.2093

dResolution = 0.5 #this is predefined by the mosart grid and elm grid

aMask=np.full((1,1), 1, dtype=np.int32)

#case index and date
iCase = 1

sDate = '20230401'
sDate_spinup = '20230401'

#e3sm components setting

#atm
iFlag_atm = 0
iFlag_datm = 1
iFlag_create_atm_grid = 0

iFlag_replace_datm_forcing=1
iFlag_replace_dlnd_forcing=0
iFlag_replace_drof_forcing=0

#elm
iFlag_lnd = 1
iFlag_dlnd =0
iFlag_create_lnd_grid = 1

#mosart
iFlag_rof = 0
iFlag_drof= 0
iFlag_create_rof_grid = 0

#order of grid

iFlag_rof_lnd_atm = 0 #rof is on
iFlag_lnd_atm_rof = 1 #rof is off

nrow = 1
ncolumn = 1

#set compset name
res, compset = e3sm_choose_res_and_compset(iFlag_atm_in=iFlag_atm, iFlag_datm_in=iFlag_datm,
                                           iFlag_lnd_in = iFlag_lnd, iFlag_dlnd_in = iFlag_dlnd,
                                           iFlag_rof_in= iFlag_rof, iFlag_drof_in= iFlag_drof)

iFlag_2d_to_1d = 0
iFlag_create_case = 1
iFlag_submit_case = 0

iFlag_default = 1
iFlag_debug = 0 #is this a debug run
iFlag_branch = 0
iFlag_initial = 1 #use restart file as initial
iFlag_lnd_spinup = 0 #is this a spinup run
iFlag_short = 0 #do you run it on short queue
iFlag_continue = 0 #is this a continue run
iFlag_resubmit = 0 #is this a resubmit
iFlag_optimal_parameter = 0


#for a single grid case, we can create this file on the fly
sPath = os.path.dirname(os.path.realpath(__file__))
pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)
sCase_spinup =  sModel + sDate_spinup + "{:03d}".format(0)

sWorkspace_scratch = '/compyfs/liao313'

#prepare a ELM namelist based on your input
sWorkspace_region = sWorkspace_scratch + slash + '04model' + slash + sModel + slash + sRegion + slash \
    + 'cases'

sWorkspace_region1 = sWorkspace_scratch + slash + '04model' + slash + sModel + slash + sRegion + slash \
    + 'cases_aux'
if not os.path.exists(sWorkspace_region):
    Path(sWorkspace_region).mkdir(parents=True, exist_ok=True)

if not os.path.exists(sWorkspace_region1):
    Path(sWorkspace_region1).mkdir(parents=True, exist_ok=True)

#some pre-defined files
sFilename_elm_surface_data_default='/compyfs/inputdata/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr2010_c191025.nc'
sFilename_elm_domain_default='/compyfs/inputdata/share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'

sFilename_dlnd_stream = '/qfs/people/liao313/data/e3sm/sag/mosart/dlnd.streams.txt.lnd.gpcc'

sFilename_initial = '/compyfs/liao313/e3sm_scratch/' + sRegion + slash  + sCase_spinup + '/run/'  + sCase_spinup +  '.elm.rh0.1980-01-01-00000.nc'
#check if the file exists
if not os.path.exists(sFilename_initial):
    print('Error: the initial file does not exist:', sFilename_initial)
    exit()



sFilename_user_datm_prec = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.Precip'
sFilename_user_datm_solar = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.Solar'
sFilename_user_datm_temp = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.TPQW'

sFilename_atm_domain=None
sFilename_lnd_domain=None
sFilename_datm_namelist=None
sFilename_rof_domain=None
sFilename_rof_namelist=None
sFilename_rof_parameter =None


sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/case.xml'

sCIME_directory ='/qfs/people/liao313/workspace/fortran/e3sm/E3SM/cime/scripts'

#why is this one needed?
sFilename_configuration = '/people/liao313/workspace/python/pye3sm/pye3sm/elm/mesh/elm_sparse_grid.cfg'

sCase_date = sDate + "{:03d}".format(iCase)
sCase = sModel + sDate + "{:03d}".format(iCase)

sWorkspace_region2 = sWorkspace_region1 + slash + sCase
if not os.path.exists(sWorkspace_region2):
    Path(sWorkspace_region2).mkdir(parents=True, exist_ok=True)

#prepare grid


if iFlag_lnd_atm_rof ==1:
    #lnd
    if iFlag_lnd ==1:
        if iFlag_create_lnd_grid ==1:
            #maybe single grid
            #aLon aLat should be used for a list of location
            aLon =np.array([dLongitude])
            aLat =np.array([dLatitude])
            sFilename_lon_lat_in = sWorkspace_region2 + slash + 'elm_' + sCase_date +'.txt'
            ofs = open(sFilename_lon_lat_in, 'w')
            ngrid = 1
            sGrid =  "{:0d}".format( ngrid)
            sLine = sGrid + '\n'
            ofs.write(sLine)
            for i in range(ngrid):
                dLongitude = aLon[i]
                dLatitude = aLat[i]
                sLine = "{:0f}".format( dLongitude ) + ' ' +  "{:0f}".format( dLatitude) + '\n'
                ofs.write(sLine)
                ofs.close()
            pass
        else:

            pass
        pass
    else:
        if iFlag_dlnd ==1:
            pass
        

if iFlag_create_case ==1:       

    if iFlag_lnd_atm_rof ==1:
        #elm component
        if iFlag_lnd ==1:
            sFilename_lnd_namelist = sWorkspace_region2 + slash + 'user_nl_elm_' + sCase_date
            sFilename_elm_surface_data_out = sWorkspace_region2 + slash + 'elm_surfdata_' + sCase_date + '.nc'
            sFilename_lnd_domain_out = sWorkspace_region2 + slash +  'elm_domain_' + sCase_date + '.nc'
            elm_create_customized_domain( aLon, aLat, aMask, dResolution, dResolution,
                                          sFilename_configuration,
                                          sFilename_elm_surface_data_default,
                                          sFilename_elm_domain_default,
                                          sFilename_elm_surface_data_out,
                                          sFilename_lnd_domain_out)
            sFilename_lnd_domain = sFilename_lnd_domain_out

            
            ofs = open(sFilename_lnd_namelist, 'w')
            sCommand_out = "fsurdat = " + "'" \
                + sFilename_elm_surface_data_out + "'" + '\n'
            ofs.write(sCommand_out)
            sCommand_out = "flndtopo = " + "'" \
                + sFilename_elm_surface_data_out  + "'" + '\n'
            ofs.write(sCommand_out)

            sLine = "hist_fincl1 = 'QOVER', 'QDRAI', 'QRUNOFF', 'ZWT', 'QCHARGE' "  + '\n'
            ofs.write(sLine)
            
            #this is a case that use existing restart file
            #be careful with the filename!!!
            sLine = "finidat = " + "'"+ sFilename_initial +"'" + '\n'
            ofs.write(sLine)
            ofs.close()

            if iFlag_create_atm_grid==1:
                pass
            else:
                sFilename_atm_domain = sFilename_lnd_domain
                pass
        else:   
            pass

        #atm component
        if iFlag_atm == 1:
            pass
        else:
            if iFlag_datm ==1:                

                sFilename_atm_domain = sFilename_lnd_domain
                pass
        
        #rof component        
        sFilename_drof_namelist=None
        pass  


    aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration ,
                                                          iFlag_debug_in = iFlag_debug,
                                                          iFlag_branch_in = iFlag_branch,
                                                          iFlag_continue_in = iFlag_continue,
                                                          iFlag_resubmit_in = iFlag_resubmit,
                                                          iFlag_short_in = iFlag_short ,
                                                          RES_in =res,
                                                          COMPSET_in = compset ,
                                                          sCIME_directory_in = sCIME_directory)
    oE3SM = pye3sm(aParameter_e3sm)

    
    iYear_start = 1980
    iYear_end = 2009

    aParameter_case = pye3sm_read_case_configuration_file(sFilename_case_configuration,
                                                              iFlag_replace_datm_forcing_in = iFlag_replace_datm_forcing,
                                                              iFlag_replace_dlnd_forcing_in = iFlag_replace_dlnd_forcing,
                                                              iFlag_replace_drof_forcing_in = iFlag_replace_drof_forcing,
                                                              iFlag_debug_case_in = 0,
                                                              iFlag_atm_in = iFlag_atm, iFlag_datm_in = iFlag_datm,
                                                              iFlag_lnd_in= iFlag_lnd, iFlag_dlnd_in= iFlag_dlnd,
                                                              iFlag_lnd_spinup_in = iFlag_lnd_spinup,
                                                              iFlag_rof_in= iFlag_rof, iFlag_drof_in= iFlag_drof,
                                                              iYear_start_in = iYear_start,
                                                              iYear_end_in = iYear_end,
                                                              iYear_data_datm_end_in = 2009,
                                                              iYear_data_datm_start_in = 1980  ,
                                                              iCase_index_in = iCase,
                                                              sDate_in = sDate,
                                                              sModel_in = sModel,
                                                              sRegion_in = sRegion,
                                                              sFilename_atm_domain_in = sFilename_atm_domain,
                                                              sFilename_datm_namelist_in = sFilename_datm_namelist ,
                                                              sFilename_user_datm_prec_in= sFilename_user_datm_prec,
                                                              sFilename_user_datm_temp_in= sFilename_user_datm_temp,
                                                              sFilename_user_datm_solar_in= sFilename_user_datm_solar,
                                                              sFilename_lnd_namelist_in = sFilename_lnd_namelist,
                                                              sFilename_lnd_domain_in = sFilename_lnd_domain,
                                                              sWorkspace_scratch_in =   sWorkspace_scratch )
   
    
    #print(aParameter_case)

    oCase = pycase(aParameter_case)

    e3sm_create_case(oE3SM, oCase )
