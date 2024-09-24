import os
from pathlib import Path
import shutil

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
from pye3sm.mosart.mesh.structured.mosart_extract_cell_by_coordinates import mosart_extract_cell_by_coordinates
from pye3sm.mosart.mesh.structured.mosart_extract_variables_for_elm import mosart_extract_variables_for_elm

iFlag_resubmit = 1
nSubmit = 2
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
iCase = 2
dMacropore_fraction = 0.1
sDate = '20240401'

aHillslope_width = [700*1445, 1035*671, 30*671, 3850*33, 1035*671,1035*671,1035*671,1035*671]
aHillslope_length = [635, 978, 1367,5000, 978,978,978,978]
aHillslope_slope = [0.08, 0.07,0.04,0.03, 0.07, 0.07, 0.07, 0.07]
dHillslope_width = aHillslope_width[iCase-2]
dHillslope_length = aHillslope_length[iCase-2]
dHillslope_slope = aHillslope_slope[iCase-2]


sDate_spinup = '20230401'

#e3sm components setting

#atm
iFlag_atm = 0
iFlag_datm = 1
iFlag_create_atm_grid = 0

iFlag_replace_datm_forcing=1
iFlag_replace_dlnd_forcing=0
iFlag_replace_drof_forcing=1

#elm
iFlag_lnd = 1
iFlag_dlnd =0
iFlag_create_lnd_grid = 1

#mosart
iFlag_rof = 0
iFlag_drof = 1
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

nTask = 1
iFlag_2d_to_1d = 0
iFlag_create_case = 1
iFlag_submit_case = 0

iFlag_default = 0
iFlag_debug = 0 #is this a debug run
iFlag_branch = 0
iFlag_initial = 1 #use restart file as initial
iFlag_lnd_spinup = 0 #is this a spinup run
iFlag_short = 0 #do you run it on short queue
iFlag_continue = 0 #is this a continue run

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

sFilename_initial = '/compyfs/liao313/e3sm_scratch/' + sRegion + slash  + sCase_spinup + '/run/'  + sCase_spinup +  '.elm.rh0.1980-01-01-00000.nc'
#check if the file exists
if not os.path.exists(sFilename_initial):
    print('Error: the initial file does not exist:', sFilename_initial)
    #exit()

sFilename_mosart_parameter_default = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
sFilename_mosart_parameter_extracted = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/data/mosart_parameter_extracted.nc'
mosart_extract_cell_by_coordinates(sFilename_mosart_parameter_default, dLongitude, dLatitude, sFilename_mosart_out=sFilename_mosart_parameter_extracted)


sFilename_user_datm_prec_origin = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.Precip'
sFilename_user_datm_solar_origin = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.Solar'
sFilename_user_datm_temp_origin = '/compyfs/liao313/04model/e3sm/amazon/user_datm.streams.txt.CLMGSWP3v1_parflow.TPQW'

sFilename_user_drof_gage_height_origin= '/compyfs/liao313/04model/e3sm/amazon/user_drof.streams.txt.MOSART.gage_height'



sFilename_atm_domain=None
sFilename_lnd_domain=None
sFilename_datm_namelist=None
sFilename_rof_domain=None
sFilename_rof_namelist=None
sFilename_rof_parameter =None


sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/case.xml'


sCIME_directory ='/qfs/people/liao313/workspace/fortran/e3sm/E3SM_DROF/cime/scripts'

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
        #

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

            sLine = "use_h2sc = .true." + '\n'
            ofs.write(sLine)

            sMacropore_fraction = "{:0.4f}".format( dMacropore_fraction)
            sLine = "macropore_fraction = " +  sMacropore_fraction  + '\n'
            ofs.write(sLine)

            sHillslope_width = "{:0.4f}".format( dHillslope_width)
            sHillslope_length = "{:0.4f}".format( dHillslope_length)
            sHillslope_slope = "{:0.4f}".format( dHillslope_slope)
            sLine = "hillslope_width = " +  sHillslope_width  + '\n'
            ofs.write(sLine)
            sLine = "hillslope_length =" +  sHillslope_length + '\n'
            ofs.write(sLine)
            sLine = "hillslope_slope = " +  sHillslope_slope + '\n'
            ofs.write(sLine)


            sLine = 'hist_empty_htapes = .true.' + '\n'
            ofs.write(sLine)
            sLine = "hist_fincl1 = 'QOVER', 'QDRAI', 'QRUNOFF', 'ZWT', 'QCHARGE','wt_slp','wt_hillslope','RAIN','SNOW','QSOIL', 'QVEGE','QVEGT', 'QDRAI_downslope','QDRAI_seepage','QDRAI_macropore' "  + '\n'
            ofs.write(sLine)

            sLine = "hist_fincl2 = 'QOVER', 'QDRAI', 'QRUNOFF', 'ZWT', 'wt_slp','wt_hillslope', 'QDRAI_downslope','QDRAI_seepage', 'QDRAI_macropore','gage_height' "  + '\n'
            ofs.write(sLine)

            sLine = 'hist_nhtfrq = 0,-24 '+ '\n'
            ofs.write(sLine)

            sLine = 'hist_mfilt = 1,365 '+ '\n'
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
            if iFlag_dlnd ==1:
                #should we use user_nl_dlnd?
                sFilename_lnd_domain=sFilename_rof_domain
                sFilename_lnd_namelist = sWorkspace_region2 + slash + 'user_nl_dlnd_' + sCase_date
                ofs = open(sFilename_lnd_namelist, 'w')
                sLine = 'dtlimit=2.0e0' + '\n'
                ofs.write(sLine)
                ofs.close()

            pass

        #atm component
        if iFlag_atm == 1:
            pass
        else:
            if iFlag_datm ==1:
                if iFlag_lnd_spinup ==1:
                    #this is a case for spin up
                    sFilename_datm_namelist = sWorkspace_region2 + slash + 'user_nl_datm_' + sCase_date
                    ofs = open(sFilename_datm_namelist, 'w')
                    sLine = 'taxmode = "cycle", "cycle", "cycle"' + '\n'
                    ofs.write(sLine)
                    ofs.close()
                    pass
                else:
                    pass

                sFilename_atm_domain = sFilename_lnd_domain
                pass
            #rof component

        #river component
        if iFlag_drof == 1:
            sFilename_drof_namelist = sWorkspace_region2 + slash + 'user_nl_drof_' + sCase_date
            ofs = open(sFilename_drof_namelist, 'w')
            sLine = 'mapalgo = "nn"' + '\n'
            ofs.write(sLine)
            #opt_elevprof = 1
            ofs.close()
            sFilename_rof_domain = sFilename_lnd_domain
            pass
        pass




    #add elevation profile into surface data
    if iFlag_lnd == 1:
        #make a copy first
        sFilename_orginal = sWorkspace_region2 + slash + 'elm_surfdata_' + sCase_date + '_original.nc'
        dest = copyfile(sFilename_elm_surface_data_out, sFilename_orginal)
        #remove the old file
        os.remove(sFilename_elm_surface_data_out)
        sFilename_old=sFilename_orginal
        sFilename_new  = sWorkspace_region2 + slash + 'elm_surfdata_' + sCase_date + '_elevation_profile.nc'
        aVariable_all=['ele0', 'ele1','ele2', 'ele3','ele4','ele5','ele6', 'ele7','ele8', 'ele9','ele10']
        aElevation_profile = mosart_extract_variables_for_elm(sFilename_mosart_parameter_extracted, aVariable_all)
        aUnit_all= ['m', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm','m']
        aDimension= [nrow, ncolumn]
        nElev=11
        aDimension_all= list()
        for i in range(nElev):
            aDimension_all.append( aDimension)
            pass
        add_multiple_variable_to_netcdf(sFilename_old, sFilename_new, aElevation_profile, aVariable_all, aUnit_all,  aDimension_all)
        sFilename_old=sFilename_new
        sFilename_new = sFilename_elm_surface_data_out
        aVariable_all = ['gxr','rdep','hslp', 'rlen']
        aUnit_all =['m', 'm-1','unitless','m']
        aDimension_all=[aDimension,aDimension,aDimension, aDimension ]
        aVariable_mosart = mosart_extract_variables_for_elm(sFilename_mosart_parameter_extracted, aVariable_all)
        add_multiple_variable_to_netcdf(sFilename_old, sFilename_new, aVariable_mosart, aVariable_all, aUnit_all,  aDimension_all)
        sFilename_elm_surface_data_out =  sFilename_new


    aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration ,
                                                          iFlag_debug_in = iFlag_debug,
                                                          iFlag_branch_in = iFlag_branch,
                                                          iFlag_continue_in = iFlag_continue,
                                                          iFlag_resubmit_in = iFlag_resubmit,
                                                          iFlag_short_in = iFlag_short ,
                                                          nTask_in= nTask,
                                                          nSubmit_in= nSubmit,
                                                          RES_in =res,
                                                          COMPSET_in = compset ,
                                                          sCIME_directory_in = sCIME_directory)
    oE3SM = pye3sm(aParameter_e3sm)
    iYear_start = 1980
    iYear_end = 2009
    aParameter_case = pye3sm_read_case_configuration_file(sFilename_case_configuration,   iCase_index_in = iCase,   sDate_in = sDate,  sModel_in = sModel,        sRegion_in = sRegion)


    #print(aParameter_case)

    oCase = pycase(aParameter_case)
    sWorkspace_output = oCase.sWorkspace_case_aux

    sFilename_user_datm_prec = sWorkspace_output + slash + 'user_datm.streams.txt.CLMGSWP3v1_parflow.Precip'
    #if not os.path.exists(sFilename_user_datm_prec):
    shutil.copyfile(sFilename_user_datm_prec_origin, sFilename_user_datm_prec)

    sFilename_user_datm_temp = sWorkspace_output + slash + 'user_datm.streams.txt.CLMGSWP3v1_parflow.TPQW'
    #if not os.path.exists(sFilename_user_datm_temp):
    shutil.copyfile(sFilename_user_datm_temp_origin, sFilename_user_datm_temp)

    sFilename_user_datm_solar = sWorkspace_output + slash + 'user_datm.streams.txt.CLMGSWP3v1_parflow.Solar'
    #if not os.path.exists(sFilename_user_datm_solar):
    shutil.copyfile(sFilename_user_datm_solar_origin, sFilename_user_datm_solar)

    sFilename_user_drof_gage_height = sWorkspace_output + '/dlnd.streams.txt.lnd.gpcc'
    #if not os.path.exists(sFilename_user_drof_gage_height):
    shutil.copyfile(sFilename_user_drof_gage_height_origin, sFilename_user_drof_gage_height)

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
                                                              iYear_end_in = iYear_start + 10,
                                                              iYear_data_datm_end_in = 2009,
                                                              iYear_data_datm_start_in = 1980  ,
                                                              iYear_data_dlnd_end_in = 2009,
                                                              iYear_data_dlnd_start_in = 1980  ,
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
                                                              sFilename_rof_namelist_in = sFilename_rof_namelist,
                                                              sFilename_rof_parameter_in = sFilename_rof_parameter,
                                                              sFilename_rof_domain_in = sFilename_rof_domain,
                                                              sFilename_drof_namelist_in = sFilename_drof_namelist,
                                                              sFilename_user_drof_gage_height_in= sFilename_user_drof_gage_height,
                                                              sWorkspace_scratch_in =   sWorkspace_scratch )

    oCase = pycase(aParameter_case)
    e3sm_create_case(oE3SM, oCase )
