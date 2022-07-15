# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:34:26 2020

@author: Gabin
"""


import geopandas
import pandas
import shapely.geometry, shapely.wkt
from osgeo import osr
import os
import time
import numpy as np
import pyproj


def MERAfiltering(pathToClimateData, pathToAgristatsQuadratsData):
    
    tic = time.time()
    
    pandasQuadrats = geopandas.read_file('{}/agriclimate_quadrats_Ireland.shp'.format(pathToAgristatsQuadratsData))
    
    prj_file = open('{}/agriclimate_quadrats_Ireland.prj'.format(pathToAgristatsQuadratsData), 'r')
    prj_txt = prj_file.read()
    
    srs = osr.SpatialReference()
    srs.ImportFromESRI([prj_txt])
    
    proj4Quadrats = srs.ExportToProj4()
    
    # IRENET95 / Irish Transverse Mercator, epsg:2157
    pandasQuadrats = geopandas.GeoDataFrame(pandasQuadrats, crs=proj4Quadrats, geometry='geometry')
    pandasQuadrats = pandasQuadrats.to_crs('epsg:2157')
    
    
    pandasFilteredAgristats = pandas.read_csv('{}/AgristatsSimplified.txt'.format(pathToAgristatsQuadratsData), sep=';')
    pandasFilteredAgristats['geometry'] = pandasFilteredAgristats['geometry'].apply(shapely.wkt.loads)
    pandasFilteredAgristats = geopandas.GeoDataFrame(pandasFilteredAgristats, crs='epsg:2157', geometry='geometry')
    
    
    prj_txt = """GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]"""
    
    srs = osr.SpatialReference()
    srs.ImportFromESRI([prj_txt])
    
    proj4Mera = srs.ExportToProj4()
    
    
    interestPoly = list(pandasFilteredAgristats.geometry.values)
    
    for polyQ in pandasQuadrats.geometry:
        
        intersect = []
        
        for polyA in pandasFilteredAgristats.geometry:
            if polyA.intersects(polyQ):
                intersect.append(polyA)
        
        diff = polyQ
        for polyI in intersect:
            diff = diff.difference(polyI)
            
        interestPoly.append(diff)
        
    pandasInterest = geopandas.GeoDataFrame(crs='epsg:2157', geometry=interestPoly)
    
    
    minx = min(pandasInterest.geometry.bounds.minx)
    miny = min(pandasInterest.geometry.bounds.miny)
    maxx = max(pandasInterest.geometry.bounds.maxx)
    maxy = max(pandasInterest.geometry.bounds.maxy)
    
    polyInterest = shapely.geometry.Polygon([[minx,maxy], [maxx,maxy], [maxx,miny], [minx,miny]])
    
    transformer = pyproj.Transformer.from_proj(pyproj.Proj('epsg:2157'), pyproj.Proj(proj4Mera))
    polyProjected = shapely.ops.transform(transformer.transform, polyInterest)
    
    boundaries = list(polyProjected.bounds)
    
    # MÉRA just uses positive longitudes
    if boundaries[0] < 0:
        boundaries[0] += 360
    
    if boundaries[2] < 0:
        boundaries[2] += 360
    
    
    columnsName = ['Latitude', 'Longitude', 'Value', 'dataDate', 'dataTime', 'validityDate', 'validityTime']
    columnsType = {'Latitude': float, 'Longitude': float, 'Value': float, 'dataDate': int, 'dataTime': int,
                   'validityDate': int, 'validityTime': int}
    
    pathResults = '{}/SimplifiedClimateFiles'.format(pathToClimateData)
    
    try:
        os.mkdir(pathResults)
    except FileExistsError:
        pass
    
    dirs = os.listdir(pathToClimateData)
    
    tac = time.time()
    print ('Time for preprocess: {} seconds'.format(round(tac-tic, 3)))
    
    for fileName in dirs:
        
        if fileName[-4::] != '.txt':
            continue
        if (fileName[:11:] != 'TotalPrecip') & (fileName[:4:] != 'Temp'):
            continue
        
        tic = time.time()
        print ('\nTreatment of file {}'.format(fileName))
        
        fileName = '{}/{}'.format(pathToClimateData, fileName)
        
        pandasClimate = pandas.read_csv(fileName, delim_whitespace=True, names=columnsName,
                                       na_values=['Latitude,', 'Longitude,', 'Value,', 'dataDate,', 'dataTime,', 'validityDate,', 'validityTime'])
        pandasClimate = pandasClimate.dropna()
        pandasClimate = pandasClimate.astype(columnsType)
        
        pandasClimate = pandasClimate.reset_index(drop=True)
        
        pandasClimateOne = pandasClimate[pandasClimate.validityDate == pandasClimate.validityDate[1]]
        
        nbDays = len(set(pandasClimate.validityDate))
        stepDays = len(pandasClimate[pandasClimate.validityDate == pandasClimate.validityDate[1]])
        
        
        pandasClimateOne = pandasClimateOne[(pandasClimateOne.Longitude > boundaries[0] - 0.1) & (pandasClimateOne.Longitude < boundaries[2] + 0.1)
                                           & (pandasClimateOne.Latitude > boundaries[1] - 0.1) & (pandasClimateOne.Latitude < boundaries[3] + 0.1)]
        
        
        geometry = [shapely.geometry.Point(lonlat) for lonlat in zip(pandasClimateOne.Longitude, pandasClimateOne.Latitude)]
        geoClimateOne = geopandas.GeoDataFrame(pandasClimateOne, crs=proj4Mera, geometry=geometry)
        
        geoClimateOne = geoClimateOne.to_crs('epsg:2157')
        
        
        geoClimateOne = geopandas.sjoin(geoClimateOne, pandasInterest, op='within') 
        geoClimateOne = geoClimateOne[geoClimateOne.index_right != 0]
        geoClimateOne = geoClimateOne.drop(columns=['index_right'])
        
        
        indexes = np.array(geoClimateOne.index.values)
        
        for i in range (1, nbDays):
            indexes = np.concatenate((indexes, np.array(geoClimateOne.index.values) + (i * stepDays)))
        
        pandasFilteredClimate = pandasClimate.iloc[indexes]
        
        pandasFilteredClimate = pandasFilteredClimate.drop(['dataDate', 'dataTime', 'validityTime'], axis='columns')
        
        if fileName.split('/')[-1][:4:] == 'Temp':
            pandasFilteredClimate = pandasFilteredClimate.groupby(['validityDate', 'Latitude', 'Longitude']).Value.mean().reset_index()
        else:
            pandasFilteredClimate = pandasFilteredClimate.reindex(columns = ['validityDate', 'Latitude', 'Longitude', 'Value']).reset_index(drop=True)
      
        geometry = [shapely.geometry.Point(lonlat) for lonlat in zip(pandasFilteredClimate.Longitude, pandasFilteredClimate.Latitude)]
        geoFilteredClimate = geopandas.GeoDataFrame(pandasFilteredClimate, crs=proj4Mera, geometry=geometry)
        
        geoFilteredClimate = geoFilteredClimate.to_crs('epsg:2157')        
        
        
        geoFilteredClimate.to_csv('{}/{}_Simplified.txt'.format(pathResults, fileName.split('/')[-1][:-4:]), index=False)
        
        
        tac = time.time()
        print ('Time for processing this file: {} seconds'.format(round(tac-tic, 3)))


if __name__ == "__main__":
    
    MERAfiltering(u'D:\Jubault\Documents\Cours\ENSG\Ing2\Stage\Données\Climate',
                  u'D:\Jubault\Documents\Cours\ENSG\Ing2\Stage\Données\Map\Quadrats\Quadrats')