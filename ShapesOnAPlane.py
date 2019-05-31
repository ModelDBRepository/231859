# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 02:39:46 2015

TODO:
* properly design plotting
* Normalize according to number of cells in bin
* take out edge cases?

@author: Oddman
"""

import json
import pickle as pkl
from SkinBrian import *
import numpy as np
import bisect
import math
import radarPlot
from pylab import *
import os

def spikeSort(spikeSet):
    outTimes = []
    outCells = []
    for cellSpikes in spikeSet:
        cellId = cellSpikes[0]
        for spike in cellSpikes[1]:
            i = bisect.bisect_right(outTimes,spike)
            outTimes.insert(i,spike)
            outCells.insert(i,cellId)
    return [outCells,outTimes]
    
class ShapesOnAPlane:
    def __init__(self,targetDir,coincidenceWindow = 2.,maxDistBin = 32,
                 test=False):
        self.coincidenceWindow = coincidenceWindow
        self.maxDistBin = maxDistBin
        self.test = test
        self.targetDir = targetDir
        pathToMetadatafile = self.targetDir+'metadata.json'
        metadataFile = open(pathToMetadatafile,'rb')
        metadata = json.load(metadataFile)
        ## group runIds which are the same, except for dummy number:
        self.dummylessRunParamRes = [] 
        self.dummylessRunParams = []
        self.runIdGroups = []  
        for runId, runMetadata in metadata.items():
            runMetadata['networkParameters'].__delitem__('networkDummySeed')
            runMetadataRe = []
            runMetadataRe.append(sorted(runMetadata['networkParameters'].items()))
            runMetadataRe.append(sorted(runMetadata['cellModelParameters'].items()))
            runMetadataRe.append(sorted(runMetadata['cacheDir']))
            if runMetadataRe in self.dummylessRunParamRes:
                thisRunIndex = self.dummylessRunParamRes.index(runMetadataRe)
                self.runIdGroups[thisRunIndex].append(runId)
            else:
                self.dummylessRunParams.append(runMetadata)
                self.dummylessRunParamRes.append(runMetadataRe)
                self.runIdGroups.append([runId])
        metadataFile.close()
        print self.runIdGroups
        self.networkTopologies = [None for i in self.runIdGroups]

            
    def run(self):
        '''High-level function to manage all the stuff to do.'''
        self.buildNetworkTopologies()
        ## Check whether there is already a processed version.
        self.processCoincidenceDetection()
        self.processCoincidenceMatrices()

    def buildNetworkTopologies(self):
        print 'Building network topologies.'
        for runGroupI,runIds in enumerate(self.runIdGroups):
            runGroupMetadata = self.dummylessRunParams[runGroupI]
            networkType = runGroupMetadata['networkParameters']['networkType']
            self.networkTopologies[runGroupI] = NetworkTopology(networkType)
            cacheDir = runGroupMetadata['cacheDir']
            ## Build the network topology (ideally, load it from cache.)
            self.networkTopologies[runGroupI] = NetworkTopology(networkType)
            networkTopologyFileName = runGroupMetadata['networkTopologyFileName']
            self.networkTopologies[runGroupI].loadFromFile(cacheDir,networkTopologyFileName)

        
    def processCoincidenceDetection(self):
        print 'Performing coincidence detection.'
        self.mats = [None for i in self.runIdGroups]
        for runGroupI,runIds in enumerate(self.runIdGroups):
            cacheFilePath = self.targetDir + 'coincidenceMatrix_{0}.npy'.format(runGroupI)
            if os.path.exists(cacheFilePath):
                print 'Loading coincidence detection from file.'
                with open(cacheFilePath,'rb') as cacheFile:
                    self.mats[runGroupI] = np.load(cacheFile)
            else:
                runGroupMetadata = self.dummylessRunParams[runGroupI]
                circumference = runGroupMetadata['networkParameters']['length_circumference'][1]
                length = runGroupMetadata['networkParameters']['length_circumference'][0]
                cellAmount = len(self.networkTopologies[runGroupI].cells)
                self.mats[runGroupI] = np.zeros((cellAmount,cellAmount),dtype = 'int64')
                for runId in runIds:
                    print 'CoincidenceDetecting run {0}'.format(runId)
                    spikeFileBasePath = self.targetDir+'spikefile_{0}'.format(runId)
                    pathToSpikeZipfile = spikeFileBasePath+'.zip'
                    pathToSpikeFile = spikeFileBasePath+'.json'
                    spikeZipFile = zipfile.ZipFile(pathToSpikeZipfile,'r')
                    spikeZipFile.extract(spikeZipFile.filelist[0],path=self.targetDir)
                    spikeFile = open(pathToSpikeFile,'rb')
                    spikes = json.load(spikeFile)
                    spikeFile.close()
                    os.remove(pathToSpikeFile)
                    [spikeCells,spikeTimes] = spikeSort(spikes)
                    maxDist = (length**2+circumference**2)**.5
                    ## sliding window
                    while len(spikeTimes) > 0:
                        spikeTime = spikeTimes.pop(0)
                        cellIdx = spikeCells.pop(0)
                        toCells = spikeCells[:bisect.bisect_left(spikeTimes,spikeTime + self.coincidenceWindow)]
                        for toCell in toCells:
                            i,j = sorted([cellIdx,toCell])
                            self.mats[runGroupI][i,j] += 1
                with open(cacheFilePath,'wb') as cacheFile:
                    print 'Creating coicidence detection cache file.'
                    np.save(cacheFile,self.mats[runGroupI])
            
    def processCoincidenceMatrices(self):
        self.binnedSectorList = [None for i in self.mats]
        for runGroupI,mat in enumerate(self.mats):
            runGroupMetadata = self.dummylessRunParams[runGroupI]
            circumference = runGroupMetadata['networkParameters']['length_circumference'][1]
            length = runGroupMetadata['networkParameters']['length_circumference'][0]
            testset = [y for y in np.ndenumerate(mat)]
            binnedSectors = {}
            networkTopology = self.networkTopologies[runGroupI]
            for dist in range(max([length,circumference])-1):
                for angle in range(min([(dist+1)*3,self.maxDistBin])):
                    binnedSectors[(dist,angle)] = 0
            [self.handlePair(a,b,n,networkTopology,binnedSectors) for ((a,b),n) in testset if n != 0]
#                    [distance,(angleBin,angle)] = self.handlePair(a,b,n,networkTopology,binnedSectors)
#                    print [distance,(angleBin,angle)]
            self.binnedSectorList[runGroupI] = binnedSectors
        
    def handlePair(self,i,j,n,networkTopology,binnedSectors):
        fromCell = networkTopology.cells[i]
        toCell = networkTopology.cells[j]
        fromCellX = fromCell.flatCoordinates[0]
        fromCellY = fromCell.flatCoordinates[1]
        toCellX = toCell.flatCoordinates[0]
        toCellY = toCell.flatCoordinates[1]
        xDistance = toCellX - fromCellX
        yDistance = toCellY - fromCellY
        circumference = networkTopology.networkParameters['length_circumference'][1]
#        length = runGroupMetadata.networkParameters['length_circumference'][0]
#        print 'distances: '+repr([xDistance,yDistance])
        distances = []
        angles = []
#        for winding in [0]:
        for winding in range(-1,2,1):
            distances.append(np.sqrt((yDistance+(winding*circumference))**2 + xDistance**2))
            angles.append(math.atan2(yDistance+(winding*circumference),xDistance))
        minIdx = distances.index(min(distances))
        distance = int(distances[minIdx]) 
        distanceAngleFactor = min([distance * 3,self.maxDistBin]) * .9999999
        angleBin = int(((angles[minIdx] + math.pi * .5) / math.pi) * distanceAngleFactor)
#        if test:
#            binnedSectors[(distance,angleBin)] = angleBin
#        else:
        binnedSectors[(distance,angleBin)] += n
        return [distance,(angleBin,angles[minIdx])]
        
def shapeOnAPlane(targetDir,runId,coincidenceWindow=2.,maxDistBin = 32,test=False):
    ## open spikefile
    pathToSpikefile = targetDir+'spikefile_{0}.json'.format(runId)
    spikeFile = open(pathToSpikefile,'rb')
    spikes = json.load(spikeFile)
    spikeFile.close()
    pathToMetadatafile = targetDir+'metadata.json'
    metadataFile = open(pathToMetadatafile,'rb')
    metadata = json.load(metadataFile)
    print metadata[runId]
    metadataFile.close()
    circumference = metadata[runId]['networkParameters']['circumference'] 
    length = metadata[runId]['networkParameters']['length']
    networkType = metadata[runId]['networkParameters']['networkType']
    networkTopologyFilePath = metadata[runId]['networkTopologyFilePath']
    networkTopology = NetworkTopology(networkType)
    assert os.path.exists(networkTopologyFilePath)
    networkTopology.loadFromFile(networkTopologyFilePath)
    cellAmount = len(networkTopology.cells)
    mat = np.zeros((cellAmount,cellAmount),dtype = 'int64')
    [spikeCells,spikeTimes] = spikeSort(spikes)
    maxDist = (length**2+circumference**2)**.5
    
    ## sliding window
    while len(spikeTimes) > 0:
        spikeTime = spikeTimes.pop(0)
        cellIdx = spikeCells.pop(0)
#        print spikeTimes[:10]
#        k.append(bisect.bisect_left(spikeTimes,spikeTime + interval))
        toCells = spikeCells[:bisect.bisect_left(spikeTimes,spikeTime + coincidenceWindow)]
        for toCell in toCells:
            i,j = sorted([cellIdx,toCell])
            mat[i,j] += 1
    print mat
    def handlePair(i,j,n):
        fromCell = networkTopology.cells[i]
        toCell = networkTopology.cells[j]
        fromCellX = fromCell.flatCoordinates[0]
        fromCellY = fromCell.flatCoordinates[1]
        toCellX = toCell.flatCoordinates[0]
        toCellY = toCell.flatCoordinates[1]
        xDistance = toCellX - fromCellX
        yDistance = toCellY - fromCellY
#        print 'distances: '+repr([xDistance,yDistance])
        distances = []
        angles = []
        for winding in [0]:
#        for winding in range(-1,2,1):
            distances.append(np.sqrt((yDistance+(winding*circumference))**2 + xDistance**2))
            angles.append(math.atan2(yDistance+(winding*circumference),xDistance))
        minIdx = distances.index(min(distances))
        distance = int(distances[minIdx]) 
        distanceAngleFactor = min([distance * 3,maxDistBin]) * .9999999
        angleBin = int(((angles[minIdx] + math.pi * .5) / math.pi) * distanceAngleFactor)
        if test:
            binnedSectors[(distance,angleBin)] = angleBin
        else:
            binnedSectors[(distance,angleBin)] += n
        return [distance,(angleBin,angles[minIdx])]
        
    testset = [y for y in np.ndenumerate(mat)]
    binnedSectors = {}
    for dist in range(max([length,circumference])-1):
        for angle in range(min([(dist+1)*3,maxDistBin])):
            binnedSectors[(dist,angle)] = 0
    for ((a,b),n) in testset:
        if n != 0:
            [distance,(angleBin,angle)] = handlePair(a,b,n)
#            print [distance,(angleBin,angle)]
    return binnedSectors
    

if __name__ == '__main__':
    sop = ShapesOnAPlane('./output/20160201_210217/')
    sop.run()

#    binnedSectors = shapeOnAPlane('./output/20151014_100811/','0',test=False)
#    binnedSectors = shapeOnAPlane('./output/20150601_232153/','2',test=True)
#    print binnedSectors.items()[:10]
#    [xes,ys] = zip(*binnedSectors.keys())
    [xes,ys] = [[],[]]
    for key,n in binnedSectors.items():
        if key[0] == 5:
            xes.append(key[1])
            ys.append(n)
    scatter(xes,ys)
    show()
    n_distances = max(dist for (dist,ang) in binnedSectors.keys())
    rp = radarPlot.RadarPlot(n_distances,binnedSectors,0.1)
    rp.makeplot('rdr_plt')