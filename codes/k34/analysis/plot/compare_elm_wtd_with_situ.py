import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyearth.toolbox.date.dt2cal import dt2cal
from pyearth.system.define_global_variables import *
import matplotlib as mpl
from pyearth.toolbox.date.day_in_month import day_in_month
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.create_qualitative_rgb_color_hex import create_qualitative_rgb_color_hex
from datetime import datetime
import netCDF4 as nc #read netcdf
from pye3sm.shared.e3sm import pye3sm
from pye3sm.shared.case import pycase
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file

iFlag_obs=1
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
#deafault has no slope, one sinlge wtd

#hexwatershed case
aHillslope_width = [700, 1035, 30, 3850, 1035,1035,1035,1035]
aHillslope_length = [635, 978, 1367,5000, 978,978,978,978]
aHillslope_slope = [0.08, 0.07,0.04,0.03, 0.07, 0.07, 0.07, 0.07]

aElevation_obs = np.array([59, 59, 60, 61, 60, 87,87, 81, 101, 96, 101,101, 101 ])
#aElevation = np.array([59, 59, 60, 61, 60, 87,87, 81, 101, 96, 101,101, 101 ])

aElevation_mosart = np.array([38, 90, 95, 99, 103, 107,111, 116, 120,125, 418 ])

dElevation_diff = 59-38
#mosart based
dHillslope_width = 1035
dArea = 3087.741393 * 1.0E6
dHillslope_length = np.sqrt(dArea) / 2.0 #from cell size
dHillslope_slope = (aElevation_mosart[10]-aElevation_mosart[0])/dHillslope_length
print(dHillslope_slope)

aHillslope_length.append(dHillslope_length)
aHillslope_slope.append(dHillslope_slope)

iYear_start = 2001
iYear_end = 2009

aTime = list()
aData = list()
aLabel_legend = list()
iFlag_first = 1


sFilename_e3sm_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/e3sm.xml'
sFilename_case_configuration = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/codes/k34/case.xml'


aDate_host=list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dom = day_in_month(iYear, iMonth)
        for iDay in range(1, dom+1):
            dSimulation = datetime(iYear, iMonth, iDay)
            aDate_host.append( dSimulation )
            pass
aDate_host=np.array(aDate_host)
nobs_host = len(aDate_host)
aData_host = np.full( (13,nobs_host), np.nan, dtype=float)


if iFlag_obs ==1:
    sWorkspace_data='/qfs/people/liao313/data'
    sModel = 'h2sc'
    sRegion='global'
    sWorkspace_auxiliary = sWorkspace_data + slash + sModel + slash + sRegion + slash + 'auxiliary'
    sFilename = sWorkspace_auxiliary + slash + 'situ' + slash + 'INPA-LBA_WellData_2001_2016.xlsx'
    xl = pd.ExcelFile(sFilename)
    aSheet = xl.sheet_names  # see all sheet names 14, last one is summary


    # we skip some data because language is not in english
    #4 of them not used
    aFlag = np.full( 13, 0, dtype=int )
    for iSheet in np.arange(1,14, 1):
        sSheet = aSheet[iSheet-1]
        if sSheet == 'PZ_PT-07':
            continue
        if sSheet == 'PP03':
            continue
        if sSheet == 'PP2':
            continue
        if sSheet == 'PP01':
            continue


        aFlag[iSheet -1 ] =1
        df = pd.read_excel(sFilename,
                           sheet_name=sSheet,
                           header=None,
                           skiprows=range(5),
                           usecols='A,E')
        df.columns = ['Date','WTD']
        dummy1 = df['Date']
        dummy2 = np.array(dummy1)
        dummy3 = dt2cal(dummy2)
        #aDate_obs = pd.to_datetime(np.array(dummy1))
        nobs =len(dummy3)
        aDate_obs = list()
        for iObs in range(nobs):
            d1=dummy3[iObs,0]
            d2=dummy3[iObs,1]
            d3= dummy3[iObs,2]
            dummy4= datetime(d1,d2 , d3 )
            aDate_obs.append( dummy4 )
            pass

        aDate_obs= np.array(aDate_obs)
        dummy5 = df['WTD']
        aWTD_obs_dummy = np.array(dummy5)  # mg/l

        dummy_index = aDate_obs-aDate_host[0]
        dummy_index1 = [x.days for x in dummy_index ]
        dummy_index1 = np.array(dummy_index1)
        dummy_index2 = np.where( dummy_index1 < nobs_host )

        dummy_obs= aWTD_obs_dummy[dummy_index2]
        dummy_index3 = dummy_index1[dummy_index2]
        aData_host[iSheet-1, dummy_index3 ] = dummy_obs
        pass

    #remove unused sites
    dummy = np.where( aFlag == 1)
    aElevation1= aElevation_obs[dummy]

    aElevation2 , indices = np.unique(aElevation1, return_index=True)
    #aOrder  = np.argsort(aElevation_obs)
    #aElevation_sort = np.sort(aElevation2)
    dummy2 = dummy[0][indices]
    aData_host2 = aData_host[dummy2, :  ]
    #get the average
    aWTD_obs = np.nanmean(aData_host2, axis=1)
    aWT_obs = aElevation2 - aWTD_obs

    dslp = (np.max(aElevation2) - np.min(aElevation2)) / 850.0

    x = ( aElevation2- np.min(aElevation2)  ) / dslp  #np.array([5]) #aElevation2
    y0 = aElevation2
    y1 = aWT_obs

x2 = np.array([0, dHillslope_length])
y2 = aElevation_mosart[0] + x2 * dHillslope_slope

iYear_start = 2000
iYear_end = 2009
iMonth_start=1
iMonth_end=12
aDate = list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(iMonth_start,iMonth_end+1):
        dSimulation = datetime(iYear, iMonth, 15)
        aDate.append( dSimulation )
        pass

nCase = 9
nData = 10 + 3

#pick a year
iYear = 2009
sModel = 'e3sm'
sRegion='k34'

aX_all = list()
aY_all = list()

aLabel_legend.append('In situ land surface')
aLabel_legend.append('In situ water table observation')
aLabel_legend.append('MOSART elevation profile')

#add default
aLabel_legend.append('Case 1 (Default)')
x_default = [0, dHillslope_length]

aParameter_e3sm = pye3sm_read_e3sm_configuration_file(sFilename_e3sm_configuration)
print(aParameter_e3sm)
oE3SM = pye3sm(aParameter_e3sm)
#read default
sDate = '20230401'
iCase_index = 1
sVariable = 'ZWT'
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
month_index = np.arange(2, nstress_subset, 12) #wet or dry
iStress=1
for iYear in range(iYear_start, iYear_end + 1):
    sYear = "{:04d}".format(iYear) #str(iYear).zfill(4)
    for iMonth in range(iMonth_start, iMonth_end + 1):
        sMonth = str(iMonth).zfill(2)
        sDummy = '.elm.h0.' + sYear + '-' + sMonth + sExtension_netcdf
        sFilename = sWorkspace_simulation_case_run + slash + sCase + sDummy
        #read before modification
        if os.path.exists(sFilename):
            #print("Yep, I can read that file: " + sFilename)
            pass
        else:
            print(sFilename)
            print("Nope, the path doesn't reach your file. Go research filepath in python")
            quit()
            pass
        aDatasets = nc.Dataset(sFilename)
        #read the actual data
        for sKey, aValue in aDatasets.variables.items():
            if sVariable.lower() == sKey.lower():
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

y_mean = aElevation_mosart[0] + dElevation_diff - np.mean(aData_out[month_index])
y_default = [y_mean,  y_mean ]

aX_all.append(x)
aY_all.append(y0) #land surface

aX_all.append(x)
aY_all.append(y1) #obs

aX_all.append(x2)
aY_all.append(y2+dElevation_diff) #mosart

aX_all.append(x_default)
aY_all.append(y_default)
sDate = '20240401'

aData_out_wt = np.full(nstress_subset,missing_value, dtype=float)

for i in range(2, 11, 1):
    iCase_index = i
    x = [0, aHillslope_length[i-2]]

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

        for iMonth in range(1, iMonth_end + 1):
            sMonth = str(iMonth).zfill(2)

            sDummy = '.elm.h0.' + sYear + '-' + sMonth + sExtension_netcdf
            sFilename = sWorkspace_simulation_case_run + slash + sCase + sDummy

            #read before modification

            if os.path.exists(sFilename):
                #print("Yep, I can read that file: " + sFilename)
                pass
            else:
                print(sFilename)
                print("Nope, the path doesn't reach your file. Go research filepath in python")
                quit()
                pass

            aDatasets = nc.Dataset(sFilename)

            #read the actual data
            for sKey, aValue in aDatasets.variables.items():
                if 'wt_hillslope' == sKey.lower():
                    #for attrname in aValue.ncattrs():
                    #print("{} -- {}".format(attrname, getattr(aValue, attrname)))
                    aData_tmp = (aValue[:]).data
                    #print(aData)
                    missing_value1 = np.max(aData_tmp)
                    aData_out[iStress-1]= aData_tmp
                    pass

                if 'wt_slp' == sKey.lower():
                    #for attrname in aValue.ncattrs():
                    #print("{} -- {}".format(attrname, getattr(aValue, attrname)))
                    aData_tmp = (aValue[:]).data
                    #print(aData)
                    missing_value1 = np.max(aData_tmp)
                    aData_out_wt[iStress-1]= aData_tmp

            iStress= iStress + 1

    y_mean = np.mean(aData_out)
    dSlope_wt = np.mean(aData_out_wt)

    #only use the specific month from the year

    y_mean = np.mean(aData_out[month_index])
    dSlope_wt = np.mean(aData_out_wt[month_index])
    y0= y_mean + dElevation_diff
    print(iCase_index, y_mean, dElevation_diff, 59-y0, dSlope_wt)
    y1=y0+ dSlope_wt * aHillslope_length[i-2]

    y = [y0, y1]
    aX_all.append(x)
    aY_all.append(y)

    aLabel_legend.append('Case ' + str(i))
    pass

#add mosart

#create the plot
iFlag = 1
iDPI = 150
iSize_x = 10
iSize_y = 8
nData_half = 5
#aColor_in0 = create_diverge_rgb_color_hex( nData_half)
aColor_in0 = create_qualitative_rgb_color_hex( 5)

aLinestyle0 = ['solid', 'dashed']
aMarker0 = ['o', 'x', 's']
#combine color and line style to create full list
aColor = list()
aLinestyle = list()
aMarker = list()
aLinethickness = list()

aColor.append('blue')
aLinestyle.append(aLinestyle0[0])
aMarker.append(aMarker0[1])
aLinethickness.append(1)

aColor.append('blue')
aLinestyle.append(aLinestyle0[1])
aMarker.append(aMarker0[0])
aLinethickness.append(1)

aColor.append('blue')
aLinestyle.append(aLinestyle0[1])
aMarker.append(aMarker0[1])
aLinethickness.append(1)

for i in range(nData_half):
    aColor.append(aColor_in0[i])
    aColor.append(aColor_in0[i])

    aLinestyle.append(aLinestyle0[0])
    aLinestyle.append(aLinestyle0[1])

    aMarker.append(aMarker0[0])
    aMarker.append(aMarker0[1])
    aLinethickness.append(1.5)
    aLinethickness.append(1.5)
    pass

dMax_x = np.max(aHillslope_length)
sFilename_out = 'wtd_hillslope_dry.png'

fig = plt.figure(dpi=iDPI)
fig.set_figwidth(iSize_x)
fig.set_figheight(iSize_y)


left, width = 0.1, 0.8
bottom, height = 0.2, 0.4
rect_full = [left, bottom, width, height]
ax_full = plt.axes(rect_full)

dMin_mini_x = 0
dMax_mini_x = 80
dMin_mini_y = 52
dMax_mini_y = 65
dX_mini = 0.125
dY_mini = 0.30
width_mini = 0.25
heigh_mini = 0.20
rect_mini = [dX_mini, dY_mini, width_mini, heigh_mini]
ax_mini = plt.axes(rect_mini)

ax_all = [ax_full, ax_mini]

for iax in range(len(ax_all)):
    ax = ax_all[iax]
    if iax == 0:
        aLegend_artist = []
        aLabel = []

    for i in np.arange(1, nData+1):
        x1 = aX_all[i-1]
        y1 = aY_all[i-1]
        #print(x1)
        #print(y1)
        ax.plot(x1, y1,
                color=aColor[i-1], linestyle=aLinestyle[i-1],
                        marker=aMarker[i-1], linewidth=aLinethickness[i-1],
                label=aLabel_legend[i-1])

    if iax == 0:
        ax.axis('on')
        ax.grid(which='major', color='grey', linestyle='--', axis='y')

        # ax.set_aspect(dRatio)  #this one set the y / x ratio

        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=10)

        ax.set_xmargin(0.05)
        ax.set_ymargin(0.15)


        sLabel_x ='Length (m)'
        sLabel_y = 'Elevation (m)'
        sTitle = 'Water table along the hillslope in February'

        ax.set_xlabel(sLabel_x, fontsize=12)
        ax.set_ylabel(sLabel_y, fontsize=12)
        ax.set_title(sTitle, loc='center', fontsize=12)

        ax.set_xlim(0, 1000)
        ax.set_ylim(50, 110)
        #ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(sFormat_y))
        ax.yaxis.set_major_locator(mpl.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.legend(loc='upper left',
          fontsize=8,
          ncol=5)

    else:
        ax.set_xlim(dMin_mini_x, dMax_mini_x)
        ax.set_ylim(dMin_mini_y, dMax_mini_y)
        ax.yaxis.set_major_locator(mpl.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

plt.savefig(sFilename_out, bbox_inches='tight')
plt.close('all')
plt.clf()

