import numpy as np
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from datetime import datetime
import netCDF4 as nc

from pyearth.toolbox.date.day_in_month import day_in_month

from pyearth.visual.timeseries.plot_time_series_data import plot_time_series_data
#only total runoff can be evaluated

#read mingpan's data
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/preprocess/mingpan_runoff_monthly.txt'


dummy = text_reader_string(sFilename_out)
aRunoff_ts_mingpan =dummy[:,0].astype(np.float64)

#read the simulated runoff
iSize_x = 12
iSize_y = 8
sDate = '20230401'
iYear_start = 2000
iYear_end = 2009

aDate_daily = list()

aDate_monthly = list()

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

aDate_monthly=np.array(aDate_monthly)
nstress_monthly = len(aDate_monthly)
aDate_daily=np.array(aDate_daily)
nstress_daily = len(aDate_daily)


aDate_all=list()
aTotal_runoff_all = list()


#add mingpan

aDate_all.append(aDate_monthly)
aTotal_runoff_all.append(aRunoff_ts_mingpan)

for i in range(2,8,1):
    sCase = 'e3sm'+ sDate + format(i, '03d')
    nData = 4
    sFilename_total_runoff = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qrunoff/qrunoff_monthly_tsplot.csv'
    dummy = text_reader_string(sFilename_total_runoff)
    aTotal_runoff =dummy[:,0].astype(np.float64)

    aDate_all.append(aDate_monthly)
    aTotal_runoff_all.append(aTotal_runoff)


sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/code/k34/analysis/evaluation/runoff/mingpan_runoff_monthly.png'

plot_time_series_data( aDate_all, aTotal_runoff_all,                      
                           sFilename_out,
                           iFlag_scientific_notation_in=1,
                           iReverse_y_in = 0,
                           sTitle_in = '',
                           sLabel_y_in= 'Total runoff',
                           dMin_y_in =0,                  
                           aLocation_legend_in = [1.0,0.0],
                           sLocation_legend_in='lower right',
                           iSize_x_in = 12,
                           iSize_y_in = 5)  





  
    
