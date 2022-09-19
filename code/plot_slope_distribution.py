#we will use a surface data to analyze the surface slope distribution
import os, sys, stat
from pathlib import Path
from netCDF4 import Dataset #read netcdf
from pyearth.visual.histogram.histogram_w_cdf_plot import histogram_w_cdf_plot

sPath_parent = str(Path(__file__).parents[1]) # data is located two dir's up
sFilename_initial = '/compyfs/liao313/e3sm_scratch/e3sm20220701052/run/e3sm20220701052.elm.h0.1981-01.nc'
sFilename = sFilename_initial
sVariable='sur_slp'
aDatasets = Dataset(sFilename)
for sKey, aValue in aDatasets.variables.items():
    if sVariable == sKey.lower():
                       
        aData_ll = (aValue[:]).data    
        break
    else:
        print(sKey.lower())

#use pyearth for visualization
sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'slope_distribution.png'
sLabel_x = 'Surface slope (percent)'
sLabel_y = 'Likelihood'
sTitle = 'Slope distribution'
sLegend = 'CDF'
histogram_w_cdf_plot( aData_ll,\
            sFilename_out, \
            iSize_x_in = 12, \
            iSize_y_in = 5, \
            iDPI_in = 150, \
            dMin_x_in = 0, \
            dMax_x_in = 0.5, \
            dSpace_x_in = 0.01, \
            sLabel_x_in = sLabel_x, \
            sLabel_y_in = sLabel_y, \
            sTitle_in = sTitle,\
                sLabel_legend_in= sLegend)


