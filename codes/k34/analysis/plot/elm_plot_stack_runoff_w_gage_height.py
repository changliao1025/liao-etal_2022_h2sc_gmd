

from datetime import datetime
import numpy as np


from pyearth.toolbox.date.day_in_month import day_in_month

import matplotlib.pyplot as plt


from pyearth.system.define_global_variables import *
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.formatter import  MathTextSciFormatter
iDPI = 300

iSize_x = 12
iSize_y = 8

iFlag_default = 0
sDate = '20240401'
#sDate = '20230401'
iYear_start = 2000
iYear_end = 2009
aLabel_y_in = [ 'Subsurface runoff (downslope)', 'Subsurface runoff (seepage)', 'Subsurface runoff (macropore)','Overland runoff']
#aLabel_y_in = [ 'Subsurface runoff','Overland runoff']
sLabel_y = 'Runoff (mm/s)'

sTitle = 'Runoff partition (macropore: 0.10)'

#read 3 runoff data
#read gage height data

iYear_start = 2000
iYear_end = 2009
iFlag_monthly = 1
iFlag_daily = 1
aDate_monthly = list()
aDate_daily = list()
nyear = iYear_end - iYear_start + 1
for iYear in range(iYear_start, iYear_end + 1):
    for iMonth in range(1,13):
        dSimulation = datetime(iYear, iMonth, 15)
        aDate_monthly.append( dSimulation )
        iDay_end = day_in_month(iYear, iMonth)
        if iMonth == 2:
            iDay_end = 28 #always use 28 days for February
        for iDay in range(1,iDay_end+1):
            aDate_daily.append(datetime(iYear,iMonth,iDay))
        pass

aDate_monthly=np.array(aDate_monthly)
nstress_monthly = len(aDate_monthly)
aDate_daily=np.array(aDate_daily)
nstress_daily = len(aDate_daily)

month_index = np.arange(2+6, nstress_monthly, 12)

for i in range(3,4,1):
    sCase = 'e3sm'+ sDate + format(i, '03d')
    if iFlag_default ==1:
        nData = 2
    else:
        nData = 4
    sFilename_overland = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qover/qover_monthly_tsplot.csv'
    sFilename_downslope = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qdrai_downslope/qdrai_downslope_monthly_tsplot.csv'
    sFilename_seepage = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qdrai_seepage/qdrai_seepage_monthly_tsplot.csv'
    sFilename_macropore = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qdrai_macropore/qdrai_macropore_monthly_tsplot.csv'
    sFilename_gage_height = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/gage_height/gage_height_daily_tsplot.csv'
    sFilename_subrunoff = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+'/qdrai/qdrai_monthly_tsplot.csv'
    sFilename_out = '/compyfs/liao313/04model/e3sm/k34/analysis/'+sCase+ '/runoff_partition'+'.png'

    if iFlag_default ==1:
        aQover = text_reader_string(sFilename_overland)[:,0].astype(np.float64)
        aSubrunoff = text_reader_string(sFilename_subrunoff)[:,0].astype(np.float64)
        sFilename_gage_height = '/compyfs/liao313/04model/e3sm/k34/analysis/'+'e3sm20240401003' +'/gage_height/gage_height_daily_tsplot.csv'
    else:
        aQover = text_reader_string(sFilename_overland)[:,0].astype(np.float64)
        aQdownslope = text_reader_string(sFilename_downslope)[:,0].astype(np.float64)
        aQseepage = text_reader_string(sFilename_seepage)[:,0].astype(np.float64)
        aQmacropore = text_reader_string(sFilename_macropore)[:,0].astype(np.float64)

    aGage_height = text_reader_string(sFilename_gage_height)[:,0].astype(np.float64)


    if iFlag_default ==1:
        aColor = create_diverge_rgb_color_hex( 3 ) # include the gage on the right side
        aData = [aSubrunoff, aQover]
        print( np.sum(aQover) / (np.sum(aSubrunoff) + np.sum(aQover) ) )
        print( np.sum(aQover[month_index]) / (np.sum(aSubrunoff[month_index]) + np.sum(aQover[month_index]) ) )
    else:
        aColor = create_diverge_rgb_color_hex( 5 ) # include the gage on the right side
        aData = [aQdownslope, aQseepage, aQmacropore, aQover]
        print( np.sum(aQover) / (np.sum(aQdownslope) + np.sum(aQseepage) + np.sum(aQmacropore) + np.sum(aQover)  ) )
        print( np.sum(aQmacropore) / (np.sum(aQdownslope) + np.sum(aQseepage) + np.sum(aQmacropore) + np.sum(aQover)  ) )

        print( np.sum(aQover[month_index]) / (np.sum(aQdownslope[month_index]) + np.sum(aQseepage[month_index]) + np.sum(aQmacropore[month_index]) + np.sum(aQover[month_index])  ) )
        print( np.sum(aQmacropore[month_index]) / (np.sum(aQdownslope[month_index]) + np.sum(aQseepage[month_index]) + np.sum(aQmacropore[month_index]) + np.sum(aQover[month_index])  ) )

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(sLabel_y)
    ax.set_title(sTitle)

    total_width = 0.7
    width = 17.0
    bottom = np.zeros(nstress_monthly)



    for i in np.arange(1, nData + 1, 1):
        data1 = aData[i-1]
        rects = ax.bar( aDate_monthly, data1, width, label= aLabel_y_in[i-1], color = aColor[i-1], bottom = bottom )
        bottom = bottom + data1
        pass

    #set ax y axis limie
    ax.set_ylim(0, 1.E-4)

    #plot gage height
    ax2 = ax.twinx()

    ax2.set_ylabel('Gage height (m)')
    ax2.plot(aDate_daily, aGage_height, color = 'black', linewidth = 1.0, linestyle = '--', #marker = "." , fillstyle = 'none',
             label = 'Gage height')


    ax.yaxis.set_major_formatter( MathTextSciFormatter("%1.2e"))

    aLocation_legend=(1.0,1.0)
    sLocation_legend = "upper right"
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    handles.insert(0, handles2[0])
    labels.insert(0, labels2[0])

    ax.legend(handles[::-1], labels[::-1],
              bbox_to_anchor=aLocation_legend,
              loc=sLocation_legend,
              fontsize=12,
              ncol= 2)


    #sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/runoff_partition'+'.png'
    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close('all')
    plt.clf()
    print('finished')