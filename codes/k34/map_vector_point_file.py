import os
import datetime
import textwrap
from collections import defaultdict
from urllib.error import URLError
import numpy as np
from osgeo import  osr, gdal, ogr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PathCollection
from matplotlib.collections import LineCollection
from matplotlib.path import Path
import cartopy as cpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.img_tiles import OSM
import cartopy.io.img_tiles as cimgt
import cartopy.io.ogc_clients as cvec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.toolbox.math.stat.remap import remap
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.visual.map.map_servers import StadiaStamen, EsriTerrain, EsriHydro, Stadia_terrain_images, Esri_terrain_images, Esri_hydro_images
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
import requests
from PIL import Image
from io import BytesIO
import shapely.geometry as sgeom

iYear_current = datetime.datetime.now().year
sYear = str(iYear_current)

# Function to fetch and process the tile image
# Function to fetch and process the tile image

# Function to process the composite image
def process_image(img):
    if img.shape[2] == 3:  # If the image is RGB
            # Add an alpha channel
        alpha_channel = np.ones((img.shape[0], img.shape[1], 1), dtype=img.dtype) * 255
        img = np.concatenate((img, alpha_channel), axis=2)

    # Make black pixels transparent
    black = (img[:, :, 0] == 0) & (img[:, :, 1] == 0) & (img[:, :, 2] == 0)
    img[black, 3] = 0  # Set alpha channel to 0 for black pixels

    return img

def map_vector_point_file(iFiletype_in,
                            sFilename_in,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFlag_color_in = None,
                            iFlag_colorbar_in=None,
                            iFlag_zebra_in=None,
                            iFlag_size_in=None,
                            iFont_size_in=None,
                            iFlag_discrete_in=None,
                            iFlag_filter_in = None,
                            iFlag_openstreetmap_in = None,
                            iFlag_openstreetmap_level_in = None,
                            iFlag_terrain_image_in = None,
                            sColormap_in=None,
                            sField_size_in = None,
                            sField_color_in=None,
                            sTitle_in=None,
                            iDPI_in=None,
                            iSize_x_in=None,
                            iSize_y_in=None,
                            dMissing_value_in=None,
                            dData_max_in=None,
                            dData_min_in=None,
                            sExtend_in=None,
                            sLocation_legend_in=None,
                            sFont_in=None,
                            sUnit_in=None,
                            aLegend_in=None,
                            aExtent_in=None,
                            pProjection_map_in=None,
                            pProjection_data_in = None):

    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSRS_geodetic = ccrs.Geodetic()

    iCount0 = 0
    iCount1 = 0

    if iFiletype_in == 1:  # geojson
        pDriver = ogr.GetDriverByName('GeoJSON')
    else:
        if iFiletype_in == 2:  # shapefile
            pDriver = ogr.GetDriverByName('Esri Shapefile')
        else:
            pDriver = ogr.GetDriverByName('Parquet')

    if os.path.exists(sFilename_in) is False:
        print('file does not exist')
        return

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    sFilename_out= sFilename_output_in


    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 100

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8

    if iFlag_size_in is not None:
        iFlag_size = iFlag_size_in
    else:
        iFlag_size = 0

    if iFlag_discrete_in is not None:
        iFlag_discrete = iFlag_discrete_in
    else:
        iFlag_discrete = 0

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if iFlag_filter_in is not None:
        iFlag_filter = iFlag_filter_in
    else:
        iFlag_filter = 0

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0


    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

    if dData_min_in is not None:
        iFlag_data_min = 1
        dData_min = dData_min_in
    else:
        iFlag_data_min = 0
        pass

    if dData_max_in is not None:
        iFlag_data_max = 1
        dData_max = dData_max_in
    else:
        iFlag_data_max = 0
        pass

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap = 'rainbow'

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title = 1
    else:
        iFlag_title = 0
        sTitle = ''

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend = 'max'

    if sLocation_legend_in is not None:
        sLocation_legend = sLocation_legend_in
    else:
        sLocation_legend = 'best'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit = ''

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    if iFlag_color == 1:
        #check color field
        if sField_color_in  is not None:
            sField_color = sField_color_in
        else:
            print('Color field is not provided')
            return

    if iFlag_size == 1:
        if sField_size_in is not None:
            sField_size = sField_size_in
        else:
            print('Size field is not provided')
            return

    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'


    fig = plt.figure( dpi=iDPI)
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    lID = 0
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aValue_field_color = list()
    nFeature = pLayer.GetFeatureCount()
    aValue_field_color = list()
    aValue_field_size = list()

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_color == 1:
            field_value = pFeature.GetField(sField_color)
            dValue = float(field_value)
            aValue_field_color.append(dValue)

        if iFlag_size == 1:
            field_value = pFeature.GetField(sField_size)
            dValue = float(field_value)
            aValue_field_size.append(dValue)

        if sGeometry_type =='POINT':
            aCoords_gcs =   get_geometry_coordinates(pGeometry_in)
            #aCoords_gcs = aCoords_gcs[:,0:2]
            dLon_max = float(np.max( [dLon_max, np.max(aCoords_gcs[:,0])] ))
            dLon_min = float(np.min( [dLon_min, np.min(aCoords_gcs[:,0])] ))
            dLat_max = float(np.max( [dLat_max, np.max(aCoords_gcs[:,1])] ))
            dLat_min = float(np.min( [dLat_min, np.min(aCoords_gcs[:,1])] ))

    if iFlag_color == 1:
        aValue_field_color = np.array(aValue_field_color)
        if iFlag_discrete ==1:
            #convert to integer
            aValue_field_color = np.array(aValue_field_color, dtype=int)
            #reorder it
            aValue_field_color = np.unique(aValue_field_color)
            aValue_field_color = np.sort(aValue_field_color)
            nValue_field = len(aValue_field_color) #get unique values from the aValue_field_color

        if iFlag_data_min == 1:  # min is provided
            dValue_min = dData_min  # np.min(aValue_field_color)
        else:
            aValue_field_color = aValue_field_color[aValue_field_color != dMissing_value]
            dValue_min = np.min(aValue_field_color)

        if iFlag_data_max == 1:  # max is provided
            dValue_max = dData_max  # np.max(aValue_field_color)
        else:
            aValue_field_color = aValue_field_color[aValue_field_color != dMissing_value]
            dValue_max = np.max(aValue_field_color)

        aValue_field_color = np.clip(aValue_field_color, a_min=dValue_min, a_max=dValue_max)
        #print(dValue_min, dValue_max)
        if dValue_max == dValue_min:
            return

    if iFlag_size == 1:
        aValue_field_size = np.array(aValue_field_size)
        dValue_size_max = np.max(aValue_field_size)
        dValue_size_min = np.min(aValue_field_size)

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = cpl.crs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                                central_latitude = 0.50*(dLat_max+dLat_min), globe=None)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84

    ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=pProjection_map ) #request.crs
    ax.set_global()

    if iFlag_discrete ==1:
        aIndex = np.linspace(0,1,nValue_field)
        prng = np.random.RandomState(1234567890)
        prng.shuffle(aIndex)
        colors = plt.colormaps[sColormap](aIndex)
        pCmap = ListedColormap(colors)
    else:
        pCmap = plt.colormaps[sColormap]

    iSize_max = 150.0
    iSize_min = 50.0

    if aExtent_in is None:
        marginx = (dLon_max - dLon_min) / 20
        marginy = (dLat_max - dLat_min) / 20
        if (dLat_max + marginy)> 90:
            dLat_max = 90
        else:
            dLat_max = dLat_max + marginy
        if (dLat_min - marginy) < -90:
            dLat_min = -90
        else:
            dLat_min = dLat_min - marginy
        if (dLon_max + marginx) > 180:
            dLon_max = 180
        else:
            dLon_max = dLon_max + marginx
        if (dLon_min - marginx) < -180:
            dLon_min = -180
        else:
            dLon_min = dLon_min - marginx
        aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
    else:
        aExtent = aExtent_in

    print(aExtent)
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    minx,  maxx, miny, maxy = aExtent
    ax.coastlines(color='black', linewidth=1,resolution='10m')


    sRegion = 'k34'

    if sRegion == 'k34':
        sFilename_k34 = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/k34.geojson'
        sFilename_dem = '/qfs/people/liao313/workspace/python/subgrid/data/k34/input/k34-DEM.tif'
        dummy = gdal_read_geotiff_file(sFilename_dem)
        img_dem = dummy['dataOut']
        pixelWidth = dummy['pixelWidth']
        pixelHeight = dummy['pixelHeight']
        nrow = dummy['nrow']
        ncolumn = dummy['ncolumn']
        aGeotransform = dummy['geotransform']
        pProjection = dummy['projection']
        dMissing_value = dummy['missingValue']
        dXmin = aGeotransform[0]
        dYmax = aGeotransform[3]
        dXmax = aGeotransform[0] + aGeotransform[1] * ncolumn
        dYmin = aGeotransform[3] + aGeotransform[5] * nrow
        aImage_extent = [dXmin, dXmax, dYmin, dYmax]
        pSpatial_refernce = osr.SpatialReference(pProjection)
        if pSpatial_refernce.IsProjected():
            sCode = pSpatial_refernce.GetAttrValue('AUTHORITY', 1)
            #pProjection_dem = ccrs.epsg(sCode)
            pProjection_data = ccrs.epsg(sCode)
        else:
            pProjection_data = pSRS_wgs84

        cmap = plt.colormaps['terrain']
        dData_max = np.max(img_dem)
        dData_min = np.min(img_dem[np.where(img_dem != dMissing_value)])
        dummy_index = np.where(img_dem > dData_max)
        img_dem[dummy_index] = dData_max
        dummy_index = np.where(img_dem < dData_min)
        img_dem[dummy_index] = dData_min
        rasterplot = ax.imshow(img_dem, origin='upper', extent=aImage_extent,cmap=cmap,
                   transform=pProjection_data, alpha=0.8)

        rasterplot.set_clim(vmin=dData_min, vmax=dData_max)

        #plot the single point
        #sFilename_k34 = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/k34.geojson'
        #or we can just lon lat
        dLongitude = -60.2093
        dLatitude = -2.6091
        #plot the point using the marker star

        ax.plot(dLongitude, dLatitude, 'ro', markersize=10, transform=pSRS_wgs84 )
        #add a label using anotate with some shift in x and y
        ax.annotate('K34', (dLongitude, dLatitude), transform=pSRS_wgs84,
                    textcoords="offset points", xytext=(25,0), ha='center')
    else:
        pass


    if iFlag_filter == 1:
        pLayer.SetSpatialFilterRect(minx, miny, maxx, maxy)

    try:
        dAlpha = 1.0
        if iFlag_openstreetmap_in is not None and iFlag_openstreetmap_in == 1:
            from cartopy.io.img_tiles import OSM
            if iFlag_openstreetmap_level_in is not None:
                iFlag_openstreetmap_level = iFlag_openstreetmap_level_in
            else:
                iFlag_openstreetmap_level = 9
                pass

            osm_tiles = OSM()
            #Add the OSM image to the map
            ax.add_image(osm_tiles, iFlag_openstreetmap_level, alpha=0.5)
            sLicense_info = "© OpenStreetMap contributors "+ sYear + "." + " Distributed under the Open Data Commons Open Database License (ODbL) v1.0."
            sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
            ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                    color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

            #we also need to set transparency for the image to be added
            dAlpha = 0.9
        else:
            if iFlag_terrain_image_in is not None and iFlag_terrain_image_in == 1:
                if iFlag_openstreetmap_level_in is not None:
                    iFlag_zoom_level = iFlag_openstreetmap_level_in
                else:
                    iFlag_zoom_level = 8
                    pass
                iFlag_zoom_level_hydro = 8
                #satellite_tiles = cimgt.Stamen('terrain-background')
                stamen_terrain = StadiaStamen("terrain-background")
                #ax.add_image(stamen_terrain, iFlag_zoom_level) #, alpha=0.5
                # Add elevation data using ESRI tiles
                esri_terrain = EsriTerrain()
                esri_hydro = EsriHydro()
                ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
                multi_poly = esri_hydro.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
                target_domain = multi_poly.geoms[0]
                _, aExtent_stamen, _ = stamen_terrain.image_for_domain(target_domain, iFlag_zoom_level)
                #_, aExtent_terrain, _ = esri_terrain.image_for_domain(target_domain, iFlag_zoom_level)
                _, aExtent_hydro, _ = esri_hydro.image_for_domain(target_domain, iFlag_zoom_level_hydro)
                img_stadia_terrain = Stadia_terrain_images(aExtent, iFlag_zoom_level)
                #img_esri_terrain  = Esri_terrain_images(aExtent, iFlag_zoom_level)
                img_eari_hydro  = Esri_hydro_images(aExtent, iFlag_zoom_level_hydro)
                #ax.imshow(img_stadia_terrain,  extent=aExtent_stamen, transform=stamen_terrain.crs, alpha=0.8)
                #ax.imshow(img_esri_terrain,  extent=aExtent_terrain, transform=esri_terrain.crs)
                ax.imshow(img_eari_hydro,  extent=aExtent_hydro, transform=esri_hydro.crs)
                #add the license information
                sLicense_info = "© Stamen Design, under a Creative Commons Attribution (CC BY 3.0) license."
                #sLicense_info = "© Esri Hydro Reference Overlay"
                sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
                ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                        color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
                dAlpha = 0.8
            else:
                pass
    except URLError as e:
        print('No internet connection')
        dAlpha = 1.0

    #ax.add_feature(cfeature.LAKES, alpha=0.5)
    #rivers = cfeature.NaturalEarthFeature(
    #    category='physical',
    #    name='rivers_lake_centerlines',
    #    scale='10m',
    #    edgecolor='blue',
    #    facecolor='none',
    #    linewidth=1.0)
    #ax.add_feature(rivers)


    aPoint_x = list()
    aPoint_y = list()
    aColor = list()
    aSize = list()
    aMarker = list()

    iFlag_special = 1
    aName = list()
    aName.append('PZ_PT-07')
    aName.append('PP03')
    aName.append('PP02')
    aName.append('PP01')

    if iFlag_special ==1:
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        aPolygon = list()
        custom_dash_pattern = (0, (5, 3))  # 5 points on, 10 points off

    aLabel = list()
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        sName = pFeature.GetField('well')

        if sName in aName:
            continue

        aLabel.append(sName)
        if sGeometry_type =='POINT':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            aPoint_x.append(aCoords_gcs[0,0])
            aPoint_y.append(aCoords_gcs[0,1])
            if iFlag_color ==1:
                pass

            else:
                color = 'blue'
            aColor.append(color)
            aSize.append(50)

            iThickness = 1.0
            lID = lID + 1

    print(iCount0, iCount1)
    scatter = ax.scatter(aPoint_x, aPoint_y, c=aColor, s=aSize, alpha=dAlpha, transform=pProjection_data)

    x_offset = 0.1  # Adjust as needed
    y_offset = 0.0  # Adjust as needed

    #for i in range(len(aX_selected)):
    #    ax.text(aPoint_x[i] + x_offset, aPoint_y[i] + y_offset, aName[i], fontsize=10, color='red', transform=pProjection_data)

    # Add labels for each point
    for i, label in enumerate(aLabel):
        ax.annotate(label, (aPoint_x[i], aPoint_y[i]), transform=pProjection_data,
                    textcoords="offset points", xytext=(15,0), ha='center')

    ax.set_extent(aExtent, crs = pSRS_wgs84)
    #gridline
    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}



    if iFlag_zebra ==1:
        ax.set_xticks(np.arange(minx, maxx+(maxx-minx)/11, (maxx-minx)/10))
        ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
        ax.set_axis_off()

    if  iFlag_title == 1:
        sTitle_wrapped = "\n".join(textwrap.wrap(sTitle, width=50))
        ax.set_title( sTitle_wrapped )

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size-2 )
    else:
        pass

    if iFlag_colorbar == 1:
        fig.canvas.draw()
        ax_pos = ax.get_position()
        ax_cb = fig.add_axes([ax_pos.x1+0.06, ax_pos.y0, 0.02, ax_pos.height])
        if iFlag_scientific_notation_colorbar == 1:
            formatter = OOMFormatter(fformat="%1.1e")
            cb = plt.colorbar(rasterplot, cax=ax_cb,
                              extend=sExtend, format=formatter)
        else:
            formatter = OOMFormatter(fformat="%1.1f")
            cb = plt.colorbar(rasterplot, cax=ax_cb,
                              extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 3
        cb.ax.set_ylabel(sUnit, rotation=90)
        cb.ax.get_yaxis().set_label_position('left')
        cb.ax.tick_params(labelsize=6)

    iFlag_colorbar = 0
    if iFlag_colorbar == 1:
        fig.canvas.draw()
        # Section 2
        ax_pos = ax.get_position() # get the original position
        #use this ax to set the colorbar ax position
        ax_cb = fig.add_axes([ax_pos.x1+0.06, ax_pos.y0, 0.02, ax_pos.height])
        if iFlag_scientific_notation_colorbar == 1:
            formatter = OOMFormatter(fformat="%1.1e")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
        else:
            if iFlag_discrete ==1:
                formatter = mpl.ticker.FuncFormatter(lambda x, pos: "{:.0f}".format(x))
                cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
            else:
                formatter = OOMFormatter(fformat="%1.2f")
                cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 3
        cb.ax.set_ylabel(sUnit, rotation=90, fontsize=iFont_size-2)
        cb.ax.get_yaxis().set_label_position('left')
        cb.ax.tick_params(labelsize=iFont_size-2)


    if iFlag_zebra ==1:
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent, crs = pSRS_wgs84)
    pDataset = pLayer = pFeature  = None

    if sFilename_output_in is None:
        plt.show()
    else:
        #remove it if exists
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == '.png':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            if sExtension == '.pdf':
                plt.savefig(sFilename_out, bbox_inches='tight')
            else:
                plt.savefig(sFilename_out, bbox_inches='tight', format='ps')
        plt.close('all')
        plt.clf()

    print('Finished plotting point map')
