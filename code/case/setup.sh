#!/bin/bash
iIndex=$1
echo $iIndex
width=3
sDate=e3sm20230101

sIndex=$(printf "%0*d" $width $iIndex)
sCasename="${sDate}${sIndex}"

sCaseFolder=/compyfs/liao313/04model/e3sm/site/cases/
sCasePath="${sCaseFolder}${sCasename}"

if [ -d $sCasePath ] 
then
    echo "${sCasePath} exist, and it will be removed."
    rm -r $sCasePath
else
    echo "${sCasePath} does not exist."
fi

sCaseBuild="/compyfs/liao313/e3sm_scratch/${sCasename}/bld"

if [ -d $sCaseBuild ] 
then
    echo "${sCaseBuild} exist, and it will be removed."
    rm -r $sCaseBuild
else
    echo "${sCaseBuild} does not exist."
fi

sCaseRun="/compyfs/liao313/e3sm_scratch/${sCasename}/run"
if [ -d $sCaseRun ] 
then
    echo "${sCaseRun} exist, and it will be removed."
    rm -r $sCaseRun
else
    echo "${sCaseRun} does not exist."
fi

echo $sCasename
echo $sCasePath

CIME_Case=/qfs/people/liao313/workspace/fortran/e3sm/E3SM/cime/scripts/

sCreate="./create_newcase --case $sCasePath --res ELM_USRDAT --compset IELMDROF --project e3sm"
echo $sCreate
cd $CIME_Case
./create_newcase --case $sCasePath --res ELM_USRDAT --compset IELMDROF --project e3sm

cd $sCasePath
./xmlchange JOB_WALLCLOCK_TIME=6:00:00
./xmlchange JOB_QUEUE=slurm --force
./xmlchange NTASKS=1
./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=1980-01-01,STOP_OPTION=nyears,STOP_N=30
./xmlchange DATM_CLMNCEP_YR_START=1980
./xmlchange DATM_CLMNCEP_YR_END=2009
./xmlchange DATM_CLMNCEP_YR_ALIGN=1980
./xmlchange DROF_MOSART_YR_START=1980
./xmlchange DROF_MOSART_YR_END=2009
./xmlchange DROF_MOSART_YR_ALIGN=1980
./xmlchange DATM_MODE=CLMGSWP3v1
./xmlchange DROF_MODE=MOSART
./xmlchange --file env_run.xml --id DOUT_S -val FALSE
./xmlchange --file env_run.xml --id INFO_DBUG -val 2
#./xmlchange ATM_DOMAIN_FILE=elm_domain_${sCasename}.nc
./xmlchange LND_DOMAIN_FILE=elm_domain_${sCasename}.nc
./xmlchange ROF_DOMAIN_FILE=elm_domain_${sCasename}.nc
#./xmlchange ATM_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm${sCasename}
./xmlchange LND_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm${sCasename}
./xmlchange ROF_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm${sCasename}
./xmlchange ELM_USRDAT_NAME=site
./case.setup

