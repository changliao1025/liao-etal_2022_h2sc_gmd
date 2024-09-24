import os, sys
from osgeo import ogr, gdal, osr
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
sPath_project = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd'
sys.path.append(sPath_project)
from codes.k34.map_vector_point_file import map_vector_point_file
sPath_project = '/qfs/people/liao313/workspace/python/subgrid'
sys.path.append(sPath_project)
from subgrid.algorithms.auxiliary.create_box_from_longitude_latitude import create_box_from_longitude_latitude
from pyearth.visual.map.raster.map_raster_file import map_raster_file

dLongitude = -60.2093
dLatitude = -2.6091
dResolution_x = 0.5
dResolution_y = 0.5
dMin_lon = -180

aExtent, aCoordinates_out = create_box_from_longitude_latitude(dLongitude, dLatitude, dResolution_x, dResolution_y)
print(aExtent)
#aExtent = [-60.21 ,-60.20,  -2.61 , -2.60]
iFiletype_in = 1
sFilename_in = ''

sFilename_dem = '/qfs/people/liao313/workspace/python/subgrid/data/k34/input/k34-DEM.tif'
sFilename_out = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/dem.png'
#map_raster_file(sFilename_dem, dData_min_in=10, iDPI_in=100, iFlag_zebra_in =1,
#                sFilename_output_in=sFilename_out,
#                aExtent_in=aExtent)


#create a point geojson file for the site
sFilename_k34 = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/k34.geojson'
pDriver = ogr.GetDriverByName('GeoJSON')
pDS = pDriver.CreateDataSource(sFilename_k34)
pLayer = pDS.CreateLayer('k34', None, ogr.wkbPoint)
pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal))
pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal))
pLayerDefn = pLayer.GetLayerDefn()
pFeature = ogr.Feature(pLayerDefn)
pFeature.SetField('id', 1)
pFeature.SetField('longitude', dLongitude)
pFeature.SetField('latitude', dLatitude)
pPoint = ogr.Geometry(ogr.wkbPoint)
pPoint.AddPoint(dLongitude, dLatitude)
pFeature.SetGeometry(pPoint)
pLayer.CreateFeature(pFeature)
pDS = None

sFilename_site = '/qfs/people/liao313/data/e3sm/k34/vector/well_site.geojson'

#now prepare the site measurement sites
image_size = [1000, 1000]
dpi = 150
scale_denominator = calculate_scale_denominator(aExtent, image_size)
pSrc = osr.SpatialReference()
pSrc.ImportFromEPSG(3857) # mercator
pProjection = pSrc.ExportToWkt()
iFlag_openstreetmap_level = calculate_zoom_level(scale_denominator, pProjection, dpi=dpi)
print(iFlag_openstreetmap_level)
sFilename_output_in = '/qfs/people/liao313/workspace/python/liao-etal_2022_h2sc_gmd/figures/study_area_full.png'
sUnit = 'Unit: m'
map_vector_point_file(iFiletype_in,
                            sFilename_site,
                            sFilename_output_in=sFilename_output_in,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFlag_color_in = None,
                            iFlag_colorbar_in=1,
                            iFlag_zebra_in=1,
                            iFlag_size_in=None,
                            iFont_size_in=None,
                            iFlag_discrete_in=None,
                            iFlag_filter_in = None,
                            iFlag_openstreetmap_in = 0,
                            iFlag_openstreetmap_level_in = 9,
                            iFlag_terrain_image_in = 0,
                            sColormap_in=None,
                            sField_size_in = None,
                            sField_color_in=None,
                            sTitle_in=None,
                            iDPI_in=dpi,
                            iSize_x_in=None,
                            iSize_y_in=None,
                            dMissing_value_in=None,
                            dData_max_in=None,
                            dData_min_in=None,
                            sExtend_in=None,
                            sLocation_legend_in=None,
                            sFont_in=None,
                            sUnit_in=sUnit,
                            aLegend_in=None,
                            aExtent_in=aExtent,
                            pProjection_map_in=None,
                            pProjection_data_in = None)