import numpy as np

from datetime import datetime
import netCDF4 as nc
from pyearth.toolbox.date.leap_year import leap_year 
from pyearth.toolbox.date.day_in_month import day_in_month
from pyearth.toolbox.data.convert_time_series_daily_to_monthly import convert_time_series_daily_to_monthly
from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data
#use the mingpan runoff data to compare with the simulated total runoff

sFolder_mingpan_runoff = '/compyfs/liao313/00raw/mingpan_runoff/amazon' #0.05 degree

#we will focus on the time period from 2000 to 2009
iYear_start = 2000
iYear_end = 2009

aDate_daily = list()
aDate_monthly=list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dSimulation = datetime(iYear, iMonth, 15)
        aDate_monthly.append( dSimulation )

        if iMonth == 2:
            dayinmonth = 28
        else:
            dayinmonth = day_in_month(iYear, iMonth)

        for iDay in range(1, dayinmonth+1):
            dSimulation = datetime(iYear, iMonth, iDay)
            aDate_daily.append( dSimulation )
        pass

#use a for loop read the mingpan runoff


#use the first year data to retrieve the domain information
sYear = "{:04d}".format(iYear_start)
sFilename = sFolder_mingpan_runoff+'/' + 'ming_daily_' + sYear + '.nc'
pData_domain = nc.Dataset(sFilename)
#read the longitude and latitude
aLon = pData_domain.variables['lon'][:]
aLat = pData_domain.variables['lat'][:]

dLongitude_min = np.min(aLon)
dLatitude_max = np.max(aLat)

#get the dimension of the domain
ncolumn = len(aLon)
nrow = len(aLat)
#get the resolution of the data in x and y direction
dResolution_x = aLon[1] - aLon[0]
dResolution_y = aLat[1] - aLat[0]
#close the file
pData_domain.close()

#calculate the row and index based on the longitude and latitude
#k34 site
dLatitude = -2.6091
dLongitude = -60.2093

#calculate the row and column index
iRow = int(( dLatitude_max - dLatitude )/dResolution_y)
iColumn = int((dLongitude - aLon[0])/dResolution_x)

aRunoff_ts_daily = np.zeros(( iYear_end-iYear_start+1, 365))
aRunoff_ts_monthly = np.zeros(( iYear_end-iYear_start+1, 12))
for iYear in range(iYear_start, iYear_end+1):
    if leap_year(iYear):
        iDay_end = 30
    else:
        iDay_end = 31
    sYear = "{:04d}".format(iYear)
    sFilename = sFolder_mingpan_runoff + '/' +'ming_daily_' + sYear + '.nc'

    #read the runoff data
    pData_runoff = nc.Dataset(sFilename)
    #read the runoff data
    aRunoff = pData_runoff.variables['QOVER'][:,:,:]
    #flip the data in the y direction
    aRunoff = np.flip(aRunoff, axis=1)
    #check the time series dimension
    nstep = aRunoff.shape[0]
    if nstep != 365:
        aRunoff = aRunoff[0:365,:,:]

    dummy_daily = aRunoff[:,iRow,iColumn]

    aRunoff_ts_daily[iYear-iYear_start,:] = dummy_daily


    #convert daily to monthly

    dummy_monthly = convert_time_series_daily_to_monthly(dummy_daily, 
                                iYear, 1, 1, 
                                iYear, 12, iDay_end , sType_in = 'mean'  )
    
    aRunoff_ts_monthly[iYear-iYear_start,:] = dummy_monthly

    pData_runoff.close()
    pass


#convert the daily 2d to 1d
aRunoff_ts_daily = np.reshape(aRunoff_ts_daily, (iYear_end-iYear_start+1)*365)
aRunoff_ts_monthly = np.reshape(aRunoff_ts_monthly, (iYear_end-iYear_start+1)*12)

#save the runoff data
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/preprocess/mingpan_runoff_daily.txt'

np.savetxt(sFilename_out, aRunoff_ts_daily, delimiter=',')
#make a plot


aTime = np.array(aDate_daily)
aData = np.array(aRunoff_ts_daily)

sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/preprocess/mingpan_runoff_daily.png'

plot_time_series_data( [aTime], [aData],                      
                           sFilename_out,
                           iFlag_scientific_notation_in=1,
                           iReverse_y_in = 0,
                           sTitle_in = '',
                           sLabel_y_in= 'Runoff',
                           
                           dMin_y_in =0,                  
                           aLocation_legend_in = [1.0,0.0],
                           sLocation_legend_in='lower right',
                           iSize_x_in = 12,
                           iSize_y_in = 5)  


#save the runoff data
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/preprocess/mingpan_runoff_monthly.txt'

np.savetxt(sFilename_out, aRunoff_ts_monthly, delimiter=',')
#make a plot


aTime = np.array(aDate_monthly)
aData = np.array(aRunoff_ts_monthly)

sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/preprocess/mingpan_runoff_monthly.png'

plot_time_series_data( [aTime], [aData],                      
                           sFilename_out,
                           iFlag_scientific_notation_in=1,
                           iReverse_y_in = 0,
                           sTitle_in = '',
                           sLabel_y_in= 'Runoff',
                           
                           dMin_y_in =0,                  
                           aLocation_legend_in = [1.0,0.0],
                           sLocation_legend_in='lower right',
                           iSize_x_in = 12,
                           iSize_y_in = 5)  
