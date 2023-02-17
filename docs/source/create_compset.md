#############
New compset
#############

In order to run a single-cell H2SC model, we need to set up a special compset.

This compset has DROF, which is data MOSART, and serves as the boundary condition for ELM.

An existing MOSART simulation output can be used to extract the stream file.


The first example:

https://github.com/E3SM-Project/E3SM/compare/master...donghuix:E3SM:donghuix/ocn-lnd/one-way-coupling#diff-114250fbefb5fc4258e488a61018e7a94e9e9693374e9221b6a928977a4c9a42


Issues in this example:

1. Naming issue::

```
{
    <entry id="strm_domdir"  skip_default_entry="true">
        <type>char</type>
        <category>streams</category>
        <group>streams_file</group>
        <desc>Stream domain file directory.</desc>
        <values>
          <value>$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_ann_rx1"      >$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_ann_ais00_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_ann_ais45_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_ann_ais55_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_iaf_rx1"      >$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_iaf_ais00_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_iaf_ais45_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.diatren_iaf_ais55_rx1">$DIN_LOC_ROOT/lnd/dlnd7/RX1</value>
          <value stream="rof.iaf_jra_1p5"          >$DIN_LOC_ROOT/lnd/dlnd7/JRA55</value>
          <value stream="rof.iaf_jra.*"            >$DIN_LOC_ROOT/lnd/dlnd7/JRA55</value>
          <value stream="rof.ryf*"                 >$DIN_LOC_ROOT/lnd/dlnd7/JRA55</value>
          <value stream="rof.mosart"               >$DIN_LOC_ROOT/rof/drof/MOSART</value>
          <value stream="rof.cplhist"              >null</value>
        </values>
      </entry>
}
```

    In this case, the DROF uses the land folder and dataset. By definition, the DROF mode should be variable produced by RTM instead of ELM.


After the modification, the new compset can be used as:

A test path is /compyfs/liao313/04model/e3sm/site/cases

res = MOS_USRDAT

compset= IELMDROF

New compset::

```
{
  <compset>
    <alias>IELMDROF</alias>
    <lname>2000_DATM%QIA_ELM%SP_SICE_SOCN_DROF%MOSART_SGLC_SWAV</lname>
  </compset>
}
```


./create_newcase --case /compyfs/liao313/04model/e3sm/site/cases/e3sm20230101001  --res ELM_USRDAT --compset IELMDROF --project e3sm

For comparison
./create_newcase --case /compyfs/liao313/04model/e3sm/amazon/cases/e3sm20230101001  --res ELMMOS_USRDAT --compset IELM --project e3sm


https://github.com/donghuix/E3SM/tree/donghuix/ocn-lnd/one-way-coupling


Prepare the amazon mosart result as a stream file

We will use the Ming Pan runoff as the template to generate the stream files.

We will use 2000-2009 as a 10 years forcing data.

We will use river gage height as the variable.

We will use a case from the earlier simulation to extract the mosart variable.
To be consistent, we will use the default ELM_MOSART, which is case 11/12. However, due to the output frequency, we will use case 13, which is 63.

Variables to be extracted: (use template from mingpan )

Example:
netcdf ming_daily_2019 {
dimensions:
	lat = 360 ;
	lon = 720 ;
	time = 365 ;
variables:
	double QDRAI(time, lat, lon) ;
		QDRAI:standard_name = "subsurface runoff" ;
		QDRAI:units = "mm/s" ;
	double QOVER(time, lat, lon) ;
		QOVER:standard_name = "surface runoff" ;
		QOVER:units = "mm/s" ;
	double QRUNOFF(time, lat, lon) ;
		QRUNOFF:standard_name = "total runoff" ;
		QRUNOFF:units = "mm/s" ;
	double lat(lat) ;
		lat:standard_name = "latitude" ;
		lat:long_name = "latitude" ;
		lat:units = "degrees_north" ;
		lat:axis = "Y" ;
	double lon(lon) ;
		lon:standard_name = "longitude" ;
		lon:long_name = "longitude" ;
		lon:units = "degrees_east" ;
		lon:axis = "X" ;
	double time(time) ;
		time:standard_name = "time" ;
		time:calendar = "noleap" ;
		time:units = "days since 2019-01-01 00:00:00" ;
		time:axis = "T" ;

// global attributes:
		:Created_by = "xudo627" ;
		:Created_on = "Wed Nov 25 14:32:19 2020 " ;
}

lon
lat
area

main channel water level

Stream file need domain file?

  <filePath>
     /compyfs/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516
  </filePath>
  <fileNames>
     domain.lnd.360x720_gswp3.0v1.c170606.nc
  </fileNames>


Add stream info into the compset. Code example from `pye3sm`:

        if iFlag_replace_dlnd_forcing==1:
            sLine = 'cp ' + sFilename_user_dlnd + ' ./user_dlnd.streams.txt.lnd.gpcc' + '\n'
            sLine = sLine.lstrip()
            ofs.write(sLine) 
        else:
            pass

An example job file:
        #!/bin/bash
        #SBATCH -A e3sm
        #SBATCH -p short
        #SBATCH -t 1:00:00
        #SBATCH -N 1
        #SBATCH -n 1
        #SBATCH -J create_case
        #SBATCH -o stdout.out
        #SBATCH -e stderr.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=chang.liao@pnnl.gov
        cd $SLURM_SUBMIT_DIR
        module purge
        sCasename=/compyfs/liao313/04model/e3sm/amazon/cases//        e3sm20220701062
        cd $sCasename
        ./xmlchange JOB_WALLCLOCK_TIME=6:00:00
        ./xmlchange JOB_QUEUE=slurm --force
        ./xmlchange NTASKS=-3
        ./xmlchange RUN_TYPE=startup
        ./xmlchange RUN_STARTDATE=1980-01-01,STOP_OPTION=nyears,        STOP_N=30
        ./xmlchange DATM_CLMNCEP_YR_START=1980
        ./xmlchange DATM_CLMNCEP_YR_END=2009
        ./xmlchange DATM_CLMNCEP_YR_ALIGN=1980
        ./xmlchange DATM_MODE=CLMGSWP3v1
        ./xmlchange -file env_run.xml -id DOUT_S -val FALSE
        ./xmlchange -file env_run.xml -id INFO_DBUG -val 2
        ./xmlchange ATM_DOMAIN_FILE=elm_domain_20220701062.nc
        ./xmlchange LND_DOMAIN_FILE=elm_domain_20220701062.nc
        ./xmlchange ATM_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm20220701062
        ./xmlchange LND_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm20220701062
        ./xmlchange ELM_USRDAT_NAME=amazon
        cp /compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm20220701062/user_nl_rtm_20220701062 ./user_nl_mosart
        cp /compyfs/liao313/04model/e3sm/amazon/cases_aux/e3sm20220701062/user_nl_elm_20220701062 ./user_nl_elm
        ./case.setup
        ./preview_namelists
        ./case.build
        ./case.submit

Since our new drof compset will also have stream data file, we will add corresponding xml configuration into E3SM.

        ./xmlchange DROF_MOSART_YR_START=1980
        ./xmlchange DROF_MOSART_YR_END=2009
        ./xmlchange DROF_MOSART_YR_ALIGN=1980
        ./xmlchange ROF_DOMAIN_FILE=elm_domain_20220701062.nc        
        ./xmlchange ROF_DOMAIN_PATH=/compyfs/liao313/04model/e3sm/        

  
A bash file will be used to test the configuration files.

In a single cell case, there is no need to use a different DATM forcing data, therefore, there is no need to use the domain file for DATM.
However, we still need to use a domain file for the DROF.
