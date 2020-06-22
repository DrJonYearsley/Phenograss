# -*- coding: utf-8 -*-

"""
    ENSG / UCD - Internship - May 2020
    ---
    Python.py : Used script for the first developments.
    ---
    G. Jubault
    jubault.gabin@gmail.com
"""



import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import copy
import os
import time



def browseFiles(folderPath):
    """
        Recover *d_sq.txt file paths in a list with a certain depth :
        [cluster[square[year]]]. This file contains MODIS data.
        
        :param folderPath: path folder containing all data
        :type folderPath: str
        
        :return: paths' list to desired data
        :rtype: list(str, ..., str)
    """    
    filePaths = []
    
    # Browse d_sq.txt files
    for cluster in os.listdir(folderPath):
        if cluster[:-1:] == 'Cluster':
            c = []
            for square in os.listdir('{}/{}'.format(folderPath,cluster)):
                if square[:-1:] == 'Square':
                    s = []
                    for year in os.listdir('{}/{}/{}'.format(folderPath,cluster,square)):
                        if year[:-2:] == '20':
                            y = []
                            for file in os.listdir('{}/{}/{}/{}'.format(folderPath,cluster,square,year)):
                                if file[-8::] == "d_sq.txt":
                                    filePath = '{}/{}/{}/{}/{}'.format(folderPath,cluster,square,year,file)
                                    y.append(filePath)
                            s.append(y)
                    c.append(s)
            filePaths.append(c)
    
    return (filePaths)



def prepareModisDates(modisPath):
    """
        Recover *r.file.date.txt file dates in a chronological sorted list.
        This file contains MODIS passages over Ireland.
        
        :param modisPath: path to the file from browseFiles()
        :type modisPath: str
        
        :return: list of each MODIS passages
        :rtype: list(str, datetime.datetime, ...)
    """    
    dates = []
    
    datePath = '/'.join(modisPath.split('/')[:-1:])
    
    for file in os.listdir(datePath):
        if file[-15::] == 'r.file.date.txt':
            datePath += '/' + file
            dateFile = open(datePath, 'r')
            
            for line in dateFile.readlines():
                if len(dates) == 0:
                    dates.append('date')
                else:
                    line = line.replace('\n', '')
                    line = line.split(' ')
                    dates.append(datetime.strptime(line[1], '%Y-%m-%d'))
    
    dates[1::].sort()
    
    return dates



def prepareQuadratsData(quadratsCoordPath):
    """
        Read a .txt file to recover ITM coordinates of quadrats : QuadratsCoordinates.txt
        on my computer. These translations were made with https://irish.gridreferencefinder.com/.
        
        :param quadratsCoordPath: path to the file
        :type quadratsCoordPath: str
        
        :return: list with a header containing quadrat's id and x and y coordinates in ITM
        :rtype: list(list, ..., list)
    """
    quadratsFile = open(quadratsCoordPath, 'r')
    
    quadratsData = []
    
    for line in quadratsFile.readlines():
        
        line = line.replace('\n', '')
        line = line.split(' ')
        
        if len(quadratsData) != 0:
            line[0] = int(line[0])
            line[1] = float(line[1])
            line[2] = float(line[2])
        
        quadratsData.append(line)
    
    quadratsFile.close()
    
    return quadratsData



def prepareModisData(modisPath):
    """
        Read a .txt file to recover MODIS data pre-processed with phenograss_process_MODIS_v2.R
        by Jon Yearsley. The output is sorted by chronological order.
        
        :param modisPath: path to the file from browseFiles()
        :type modisPath: str
        
        :return: array containing MODIS data
        :rtype: numpy.ndarray
    """
    modisFile = open(modisPath, 'r')
    modisData = []
    
    for line in modisFile.readlines():
        
        line = line.replace('"',  '')
        line = line.replace('\n', '')
        line = line.split(' ')
        
        if len(modisData) == 0:
            line.insert(0, 'id')
        else:
            line[0] = int(line[0])
            line[1] = float(line[1])
            line[2] = float(line[2])
            line[3] = float(line[3])
            line[4] = float(line[4])
            line[6] = int(line[6])
            line[7] = int(line[7])
            line[8] = int(line[8])
            line[9] = float(line[9])
            line[10] = float(line[10])
            line[11] = int(line[11])
            line[12] = datetime.strptime(line[12], '%Y-%m-%d')
            
        modisData.append(line)
        
    modisFile.close()
    
    header = modisData[0]
    modisData = sorted (modisData[1::], key=lambda day: day[8])
    modisData = np.insert(modisData, 0, header, axis=0)
    
    return modisData



def coords(modisData, path, quadratsData, rewrite=False):
    """
        Use MODIS data and quadrats coordinates to plot coordinates of both system.
        It permits to verify the surface area of the work zone, the centering on the principal
        point and roughly the projection.
        
        :param modisData: MODIS data from prepareModisData()
        :param path: name of the path to the .txt MODIS data
        :param quadratsData: quadrats coordinates from prepareQuadratsData()
        :param rewrite: define if plots already saved will be overwritten or not
        :type modisData: numpy.ndarray
        :type path: str
        :type quadratsData: list(list, ..., list)
        :type rewrite: bool
    """
    # Quadrat number search
    for name in (path.split('/')[-1]).split('_'):
        if 'square' in name:
            numQuadrat = int(name.replace('square',''))
            break
    
    # Creation of folders to store results
    path = "/".join(path.split('/')[:-1:])
    
    path += '/Python'
    createFolder(path)
    
    path += '/coordinates'
    createFolder(path)
    
    # Treatment
    pathFig = '{}/{}{}'.format(path,'coord','.png')
    
    if not rewrite:
        if os.path.isfile(pathFig):
            return
    
    transpose = modisData[1::].T
    xITM, yITM, xMODIS, yMODIS = transpose[1], transpose[2], transpose[3], transpose[4]
    
    xStatsITM = [min(xITM), max(xITM)]
    yStatsITM = [min(yITM), max(yITM)]
    xStatsMODIS = [min(xMODIS), max(xMODIS)]
    yStatsMODIS = [min(yMODIS), max(yMODIS)]    
    
    # Graphics
    fig = plt.figure(figsize=(12,9), dpi=200)
    
    plt.suptitle('Measured coordinates', size='xx-large', weight='bold', y=0.79)
    
    plt.subplot(121)
    plt.plot(xITM, yITM, 'b.')
    plt.plot(quadratsData[numQuadrat][1], quadratsData[numQuadrat][2], 'mo', markersize=10)
    plt.xlabel('xITM', size='large', style='italic')
    plt.ylabel('yITM', size='large', style='italic')
    plt.grid()
    
    plt.gca().set_aspect(aspect='equal')
    axisStats = plt.axis()
    
    plt.plot([axisStats[0],axisStats[1]], [yStatsITM[0],yStatsITM[0]], 'r-')
    plt.plot([axisStats[0],axisStats[1]], [yStatsITM[1],yStatsITM[1]], 'r-')
    plt.plot([xStatsITM[0],xStatsITM[0]], [axisStats[2],axisStats[3]], 'r-')
    plt.plot([xStatsITM[1],xStatsITM[1]], [axisStats[2],axisStats[3]], 'r-')
    plt.text(axisStats[0], axisStats[3]-200, str(round(xStatsITM[0])) + ' ; '
                                         + str(round(yStatsITM[1])), weight='bold',
                                         bbox=dict(facecolor='red', alpha=0.5))
    plt.text(axisStats[1]-3000, axisStats[2], str(round(xStatsITM[1])) + ' ; '
                                         + str(round(yStatsITM[0])), weight='bold',
                                         bbox=dict(facecolor='red', alpha=0.5))
    
    plt.subplot(122)
    plt.plot(xMODIS, yMODIS, 'b.')
    plt.xlabel('xMODIS', size='large', style='italic')
    plt.ylabel('yMODIS', size='large', style='italic', rotation=-90, labelpad=15)
    plt.grid()
    
    plt.gca().yaxis.set_ticks_position('right')
    plt.gca().yaxis.set_label_position('right')
    plt.gca().set_aspect(aspect='equal')
    axisStats = plt.axis()
    
    plt.plot([axisStats[0],axisStats[1]], [yStatsMODIS[0],yStatsMODIS[0]], 'r-')
    plt.plot([axisStats[0],axisStats[1]], [yStatsMODIS[1],yStatsMODIS[1]], 'r-')
    plt.plot([xStatsMODIS[0],xStatsMODIS[0]], [axisStats[2],axisStats[3]], 'r-')
    plt.plot([xStatsMODIS[1],xStatsMODIS[1]], [axisStats[2],axisStats[3]], 'r-')
    plt.text(axisStats[0], axisStats[3]-200, str(round(xStatsMODIS[0])) + ' ; '
                                         + str(round(yStatsMODIS[1])), weight='bold',
                                         bbox=dict(facecolor='red', alpha=0.5))
    plt.text(axisStats[1]-3000, axisStats[2], str(round(xStatsMODIS[1])) + ' ; '
                                         + str(round(yStatsMODIS[0])), weight='bold',
                                         bbox=dict(facecolor='red', alpha=0.5))
    
    fig.savefig(pathFig)
    
    plt.close(fig)



def shadeCoords(modisData, path, rewrite='False'):
    """
        Use MODIS data to plot coordinates with 'rainbow' shade of NDVI and EVI measures.
        The aim is to try to identify parcels with similarities.
                
        :param modisData: MODIS data from prepareModisData()
        :param path: name of the path to the .txt MODIS data
        :param rewrite: define if plots already saved will be overwritten or not
        :type modisData: numpy.ndarray
        :type path: str
        :type rewrite: bool
    """
    # Creation of folders to store results
    path = "/".join(path.split('/')[:-1:])
    
    path += '/Python'
    createFolder(path)
    
    path += '/coordinates'
    createFolder(path)
    
    # Treatment
    transpose = modisData[1::].T
    doy, evi, ndvi = transpose[8], transpose[9], transpose[10]
    xMODIS, yMODIS = transpose[3], transpose[4]
    xMin, xMax, yMin, yMax = min(xMODIS), max(xMODIS), min(yMODIS), max(yMODIS)    
    
    doyShade = {}
    
    # doyShade = {day in the year: [appearance indices, ...], ...}
    for i in range (len(doy)):
        key = doy[i]
        
        if key in doyShade:
            indexes = doyShade[key]
            indexes.append(i)
            doyShade[key] = indexes
        else:
            doyShade[key] = [i]
    
    # Keep days with the most measures
    doyShade = sorted(doyShade.items(), key=lambda num: len(num[1]), reverse=True)[:20:]
    
    for i in range (len(doyShade)):
        
        pathFig = '{}/{}{}'.format(path,'Shade-{}'.format(str(doyShade[i][0])),'.png')
        
        if not rewrite:
            if os.path.isfile(pathFig):
                return
        
        x, y, n, e = [], [], [], []
        
        for index in doyShade[i][1]:
            x.append(xMODIS[index])
            y.append(yMODIS[index])
            n.append(ndvi[index])
            e.append(evi[index])
        
        # Graphics
        cm = plt.cm.get_cmap('jet')
        
        fig, axes = plt.subplots(1, 2, figsize=(12,9), dpi=200)
        plt.suptitle("Nuanced coordinates to identify neighborhood", size='xx-large', weight='bold', y=0.80)
        
        plt.subplot(121)
        plt.scatter(x, y, c=n, marker='.', s=150, edgecolor='none', cmap=cm, vmin=0, vmax=1)
        plt.title('NDVI', size='x-large', weight='bold')
        plt.xlabel('xMODIS', size='large', style='italic')
        plt.ylabel('yMODIS', size='large', style='italic')
        plt.xlim(xMin-200, xMax+200)
        plt.ylim(yMin-200, yMax+200)
        plt.grid()
        
        plt.gca().set_aspect(aspect='equal')
        
        plt.subplot(122)
        plt.scatter(x, y, c=e, marker='.', s=150, edgecolor='none', cmap=cm, vmin=0, vmax=1)
        plt.title('EVI', size='x-large', weight='bold')
        plt.xlabel('xMODIS', size='large', style='italic')
        plt.ylabel('yMODIS', size='large', style='italic', rotation=-90, labelpad=15)
        plt.xlim(xMin-200, xMax+200)
        plt.ylim(yMin-200, yMax+200)
        plt.grid()
        
        plt.gca().set_aspect(aspect='equal')
        plt.gca().yaxis.set_ticks_position('right')
        plt.gca().yaxis.set_label_position('right') 
        
        sm = plt.cm.ScalarMappable(cmap=cm)
        fig.colorbar(sm, ax=axes.ravel().tolist(), orientation='vertical', fraction=0.52, shrink=0.5)
        
        fig.savefig(pathFig)
        
        plt.close(fig)



def simpleGraph(modisData, path, method, rewrite=False):
    """
        Use MODIS data to plot NDVI and EVI measures in a year.
        There are graphics with every measures beside one with the average of measures
        on each days.
        
        :param modisData: MODIS data from prepareModisData()
        :param path: name of the path to the .txt MODIS data
        :param method:
        :param rewrite: define if plots already saved will be overwritten or not
        :type modisData: numpy.ndarray
        :type path: str
        :type method:
        :type rewrite: bool
    """
    # Creation of folders to store results
    if method == 'year':
        path = "/".join(path.split('/')[:-1:])
    
    elif method == 'square':
        path = "/".join(path[0].split('/')[:-2:])
    
    path += '/Python'
    createFolder(path)
    
    path += '/simpleGraph'
    createFolder(path)
    
    # Treatment
    if method == 'year':
        nbYears = 1
        pathFig = '{}/{}{}'.format(path,'onYear','.png')
    
    elif method == 'square':
        modisData, nbYears = squareData(modisData)
        startYear = datetime.strftime(modisData[1][12], '%Y')
        endYear = datetime.strftime(modisData[-1][12], '%Y')
        
        pathFig = '{}/{}-{}{}'.format(path,startYear,endYear,'.png')
    
    if not rewrite:
        if os.path.isfile(pathFig):
            return
    
    transpose = modisData[1::].T
    doy, evi, ndvi = transpose[8], transpose[9], transpose[10]
    doyAverage, eviAverage, ndviAverage = sorted(set(doy)), [], []
    eviMax, eviMin, ndviMax, ndviMin = [], [], [], []
    
    # count = {day in the year: number of occurrences, ...}
    count = {}.fromkeys(doyAverage,0)
    for value in doy:
        count[value] += 1
    
    index = 0
    for key,value in count.items():
        
        # 50 is an arbitrary value to filter days without enough values to avoid noise
        if value > 50:
            eviSum, ndviSum = 0, 0
            
            for i in range (value):
                i += index
                
                eviSum += evi[i]
                ndviSum += ndvi[i]
            
            eviAverage.append(eviSum/value)
            ndviAverage.append(ndviSum/value)
            
            eviMax.append(max(evi[index:value+index:]))
            eviMin.append(min(evi[index:value+index:]))
            
            ndviMax.append(max(ndvi[index:value+index:]))
            ndviMin.append(min(ndvi[index:value+index:]))
        else:
            doyAverage.remove(key)
            
        index += value
    
    # Graphics
    fig = plt.figure(figsize=(12,9), dpi=200)
    
    plt.suptitle("Quadrat's MODIS vegetation indices", size='xx-large', weight='bold', y=0.96)
    
    plt.subplot(221)
    plt.title('All measures', size='x-large', weight='bold')
    plt.plot(doy, ndvi, 'b.', markersize=3)
    plt.xlim(0,365*nbYears)
    plt.ylabel('NDVI', size='large', style='italic')
    plt.ylim(0,1)
    plt.grid()
    
    plt.subplot(223)
    plt.plot(doy, evi, 'b.', markersize=3)
    plt.xlabel("Year's day", size='large', style='italic')
    plt.xlim(0,365*nbYears)
    plt.ylabel('EVI', size='large', style='italic')
    plt.ylim(0,1)
    plt.grid()
    
    plt.subplot(222)
    plt.title('Stats on points with more than 50 values', size='x-large', weight='bold')
    plt.plot(doyAverage, ndviMax, 'r:o', markersize=3, linewidth=1, label='max')
    plt.plot(doyAverage, ndviAverage, 'g:o', markersize=3, linewidth=1, label='average')
    plt.plot(doyAverage, ndviMin, 'b:o', markersize=3, linewidth=1, label='min')
    plt.legend(loc='lower right', framealpha=0.8)
    plt.xlim(0,365*nbYears)
    plt.ylim(0,1)
    plt.grid()
    
    plt.subplot(224)
    plt.plot(doyAverage, eviMax, 'r:o', markersize=3, linewidth=1, label='max')
    plt.plot(doyAverage, eviAverage, 'g:o', markersize=3, linewidth=1, label='average')
    plt.plot(doyAverage, eviMin, 'b:o', markersize=3, linewidth=1, label='min')
    plt.legend(loc='upper right', framealpha=0.8)
    plt.xlabel("Year's day", size='large', style='italic')
    plt.xlim(0,365*nbYears)
    plt.ylim(0,1)
    plt.grid()
    
    fig.savefig(pathFig)
    
    plt.close(fig)



def perPixel(modisData, path, method, rewrite=False):
    """
        Use MODIS data to plot NDVI and EVI time series for each pixel.
        Only the 20 pixels the most measured are saved.
        
        :param modisData: MODIS data from prepareModisData()
        :param path: name of the path to the .txt MODIS data
        :param method:
        :param rewrite: define if plots already saved will be overwritten or not
        :type modisData: numpy.ndarray
        :type path: str
        :type method:
        :type rewrite: bool
    """
    # Creation of folders to store results
    if method == 'year':
        path = "/".join(path.split('/')[:-1:])
    
    elif method == 'square':
        path = "/".join(path[0].split('/')[:-2:])
    
    path += '/Python'
    createFolder(path)
    
    path += '/perPixel'
    createFolder(path)
    
    # Treatment
    if method == 'year':
        nbYears = 1
    
    elif method == 'square':
        modisData, nbYears = squareData(modisData)
    
    transpose = modisData[1::].T
    xMODIS, yMODIS = transpose[3], transpose[4]
    coordMODIS = np.array([xMODIS,yMODIS]).T
    
    location = {}
    
    # location = {(xMODIS,yMODIS): [appearance indices, ...], ...}
    for i in range (len(coordMODIS)):
        key = (coordMODIS[i,0], coordMODIS[i,1])
        
        if key in location:
            indexes = location[key]
            indexes.append(i)
            location[key] = indexes
        else:
            location[key] = [i]
    
    # Keep coordinates with the most measures
    location = sorted(location.items(), key=lambda num: len(num[1]), reverse=True)[:20:]
    
    for count in location:
        
        key, value = count[0], count[1]
        
        pathFig = '{}/{}_{}_{}{}'.format(path,str(len(value)),str(key[0]),str(key[1]),'.png')
        
        if not rewrite:
            if os.path.isfile(pathFig):
                continue
        
        doy, evi, ndvi = [], [], []
        
        for index in value:
            doy.append(transpose[8][index])
            evi.append(transpose[9][index])
            ndvi.append(transpose[10][index])
        
#        dNdvi, dEvi, dDoy = [], [], []
#        for i in range (len(ndvi)-1):
#            dDoy.append((doy[i] + doy[i+1]) / 2)
#            dNdvi.append((ndvi[i+1] - ndvi[i]) / (doy[i+1] - doy[i]))
#            dEvi.append((evi[i+1] - evi[i]) / (doy[i+1] - doy[i]))
        
        # Graphics
        fig = plt.figure(figsize=(12,9), dpi=200)
        
        plt.suptitle("Pixel's time series", size='xx-large', weight='bold', y=0.92)
        
        plt.subplot(211)
        plt.plot(doy, ndvi, 'b:o')
        plt.xlim(0,nbYears*365)
        plt.ylabel('NDVI', size='large', style='italic')
        plt.ylim(0,1)
        plt.grid()
        
        plt.subplot(212)
        plt.plot(doy, evi, 'b:o')
        plt.xlabel("Year's day", size='large', style='italic')
        plt.xlim(0,nbYears*365)
        plt.ylabel('EVI', size='large', style='italic')
        plt.ylim(0,1)
        plt.grid()
        
        fig.savefig(pathFig)
        
        plt.close(fig)



def createFolder(path):
    if not os.path.isdir(path):
        os.makedirs(path)



def squareData(modisData):
    copyModisData = copy.deepcopy(modisData)
    for data in copyModisData:
        try:
            if len(globalModisData) > 0:
                data = data[1::]
                data[:,8] = data[:,8] + (nbYears * 365)
                globalModisData = np.concatenate((globalModisData, data), axis=0)
                nbYears += 1
        except NameError:
            globalModisData = data
            nbYears = 1
    
    return globalModisData, nbYears



if __name__ == '__main__':
    
    tic = time.time()
    
        
    folderPath = "Données/Archive"
    
    quadratsCoordPath = folderPath + '/QuadratsCoordinates.txt'
    quadratsCoordData = prepareQuadratsData(quadratsCoordPath)
    
    filePaths = browseFiles(folderPath)
    
    
    for cluster in filePaths:
        clusterModisData, clusterModisDates, clusterFilePaths = [], [], []
        for square in cluster:
            squareModisData, squareModisDates, squareFilePaths = [], [], []
            for year in square:
                filePath = year[0]
                
                # per year
                modisData = prepareModisData(filePath)
                modisDates = prepareModisDates(filePath)
#                coords(modisData, filePath, quadratsCoordData, True)
#                shadeCoords(modisData, filePath, True)
                perPixel(modisData, filePath, 'year', True)
                simpleGraph(modisData, filePath, 'year', True)
#                break
#            break
#        break
    
                # per square
                squareFilePaths.append(filePath)
                squareModisData.append(modisData)
                squareModisDates.append(modisDates)
#            perPixel(squareModisData, squareFilePaths, 'square', True)
#            simpleGraph(squareModisData, squareFilePaths, 'square', True)
#            break
#        break
    
#            # per cluster
#            clusterFilePaths.append(squareFilePaths)
#            clusterModisData.append(modisData)
#            clusterModisDates.append(modisDates)
    
    
    tac = time.time()
    print ("\nExecution time : {} seconds.".format(round(tac-tic,3)))