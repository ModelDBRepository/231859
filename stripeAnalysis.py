# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 13:34:12 2015

@author: oddman
"""

from SkinBrian import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
#import bisect
import os
import collections
import zipfile
import pickle as pkl
import diagonalShifter
import copy
reload(diagonalShifter)

ROWDIST = .75**.5
RINGINESSTHRESHOLDS = [.5,1.,1.5,2.] ## * sd
MARKERS = ['o','v','^','8','s','*','h','x','+']
LINESTYLES = ['-','--']

SHIFTED_MATS = []

class StripeAnalysis:
    def __init__(
                self,targetDir, cacheDir, angles = 128, coincidenceWindow = 2.,
                groupDummies = False, useOutfile = False, getFromCache = True,
                fixedRingAmount = None, filters = []
    ):
        self.coincidenceWindow = coincidenceWindow
        self.useOutfile = useOutfile 
        self.targetDir = targetDir
        self.angles = angles
        self.filters = filters
        self.getFromCache = getFromCache
        if self.useOutfile:
            self.outFile = open(self.targetDir+'ring_aggregated.csv','w')
            self.outFile.write('')
            self.outFile.close()
            self.outFile = open(self.targetDir+'ring_aggregated.csv','a')
            self.outFile.write(';'.join(['ringiness','runIds','frame','ring','val'])+'\n')
        pathToMetadatafile = self.targetDir+'metadata.json'
        with open(pathToMetadatafile,'rb') as metadataFile:
            rawMetadata = json.load(metadataFile)
            ## metadata: turn length/circumference lists into tuples for set-ing
            for key in rawMetadata.keys():
                print rawMetadata
                rawMetadata[key][u'networkParameters'][u'length_circumference'] = tuple(rawMetadata[key][u'networkParameters'][u'length_circumference'])

            for filterAddress,filterValue in self.filters:
#                print rawMetadata['0'][filterAddress[0]][filterAddress[1]]
                print 'fV'
                print filterValue
                print [
                    val[filterAddress[0]][filterAddress[1]] for key,val 
                    in rawMetadata.items()]
                rawMetadata = dict(
                    (key,val) for key,val 
                    in rawMetadata.items()
                    if val[filterAddress[0]][filterAddress[1]] in filterValue)
                
            if len(rawMetadata) == 0:
                raise ValueError,'Filters leave no valid values.'
            self.metadata = rawMetadata
            
#            print self.metadata
        self.maxTime = self.metadata.values()[0]['runParameters']['t_end']
        self.mats = dict(
            (runId,dict(
                    (val,None) for val in fixedRingAmount
                )
            ) for runId,val in self.metadata.items()
        )
        self.devianceMats = dict(
            (runId,dict(
                    (val,None) for val in fixedRingAmount
                )
            ) for runId,val in self.metadata.items()
        )
        self.steppedMats = dict(
            (runId,dict(
                    (val,None) for val in fixedRingAmount
                )
            ) for runId,val in self.metadata.items()
        )
        self.timeBins = int((self.maxTime / self.coincidenceWindow) + 1)
        self.fixedRingAmount = fixedRingAmount
#        assert self.maxTime % self.coincidenceWindow == 0.0
        self.figs = []
        ## metadata: turn length/circumference lists into tuples for set-ing
        for key in self.metadata.keys():
            self.metadata[key][u'networkParameters'][u'length_circumference'] = tuple(self.metadata[key][u'networkParameters'][u'length_circumference'])
        
        self.cacheFilePath = self.targetDir+'stripeAnalysisCache.pkl'
        if os.path.exists(self.cacheFilePath) and self.getFromCache:
            ## try loading previous version.
            print "Loading cached stripe analysis data from: " + self.cacheFilePath
            self.loadFromFile(self.cacheFilePath)
        else:
            self.detectCoincidence(self.fixedRingAmount)
        
    def detectCoincidence(self,fixedRingAmount):
        '''fixedRingAmount is a list of values, each number representing
        the number of rings in the system. If it is None, every ring will be
        represented.
        If it is an integer, the number of rings to be analyzed should be 
        divisible by that integer. It is intended as a normalizing method for 
        an exponential body size parameter scan.
        '''
        for i,thisRunMetadata in sorted(self.metadata.items()):
            for thisFixedRingAmount in fixedRingAmount:
                if thisFixedRingAmount != None:
                    ## length = the number of rings in this run
                    length,circumference = thisRunMetadata[u'networkParameters'][u'length_circumference']
                    ## Make sure the fixedRingAMount is valid for this length
                    assert (length / float(thisFixedRingAmount)) % 1 == 0, 'Invalid fixedRingAmount/length combo - this body length is not divisible by fixedRingAmount'
                    self.processSingle(i,thisRunMetadata,thisFixedRingAmount)
                else:
                    self.processSingle(i,thisRunMetadata,1)
        if self.getFromCache:
            self.serialize(self.cacheFilePath) 

    def loadFromFile(self,cacheFilePath):
        with open(cacheFilePath,'rb') as openFile:
            [
                    self.mats,
                    self.devianceMats,
                    self.steppedMats
            ] = pkl.load(openFile) 

    def serialize(self,pathToCacheFile): 
        with open(pathToCacheFile,'wb') as openFile:
            pkl.dump(
                [
                    self.mats,
                    self.devianceMats,
                    self.steppedMats
                ],
                openFile
            )

    def processSingle(self,runId,thisRunMetadata,fixedRingAmount):
        '''Requires already loaded metadata.'''
        print runId
        
#        runId = int(runId)
        ## Get some relevant metadata.
        length = thisRunMetadata['networkParameters']['length_circumference'][0]
        ringsPerRingBin = int(length / float(fixedRingAmount))       
        networkType = thisRunMetadata['networkParameters']['networkType']
        networkTopologyFileName = thisRunMetadata['networkTopologyFileName']
        cacheDir = thisRunMetadata['cacheDir']
        ## Build the network topology (ideally, load it from cache.)
        thisNetworkTopology = NetworkTopology(networkType)
        thisNetworkTopology.generate(
            thisRunMetadata['networkParameters'], 
            cache = True, 
            cacheDir = cacheDir, 
            fileName = networkTopologyFileName,
            createCache = True
        )
        ## Set up the empty result matrix.
        rowBins = length / ringsPerRingBin
        self.mats[runId][fixedRingAmount] = np.zeros((rowBins,self.timeBins),dtype = 'int64')
        ## Start getting the actual data. Also involves unzipping.
        spikeFileBasePath = self.targetDir+'spikefile_{0}'.format(runId)
        pathToSpikeZipfile = spikeFileBasePath+'.zip'
        pathToSpikeFile = spikeFileBasePath+'.json'
        spikeZipFile = zipfile.ZipFile(pathToSpikeZipfile,'r')
        spikeZipFile.extract(spikeZipFile.filelist[0],path=self.targetDir)
        spikeFile = open(pathToSpikeFile,'rb')
        spikes = json.load(spikeFile)
        spikeFile.close()
        os.remove(pathToSpikeFile)
        ## Process the actual spikes.
        for cellIndex,spikeTrain in spikes:
#            timeBinIFactor = runGroupRunI
            cell = thisNetworkTopology.cells[cellIndex]
            row = int(cell.z/ROWDIST + .5)
            rowBin = int(row / float(ringsPerRingBin))
            for spike in spikeTrain:
                ## Put the spikes in the correct matrix cell one by one.
                spikeBin = int(spike / self.coincidenceWindow)
                self.mats[runId][fixedRingAmount][rowBin,spikeBin] += 1
        #======================================================================
        # Do some simple statistics.
        #======================================================================
        mean = np.sum(self.mats[runId][fixedRingAmount]) / float(np.product(self.mats[runId][fixedRingAmount].shape))
        deviance = (self.mats[runId][fixedRingAmount] - mean) ** 2
        self.devianceMats[runId][fixedRingAmount] = deviance
        BesselCorrectedNumerator = np.product(self.mats[runId][fixedRingAmount].shape) - 1
        sd = np.sum(deviance) / float(BesselCorrectedNumerator)
        binaryMats = [self.mats[runId][fixedRingAmount] > mean + sd * rt for rt in RINGINESSTHRESHOLDS]
        binMatsToSum = [np.expand_dims(a,2) for a in binaryMats]            
        binMatsToSum = np.concatenate(binMatsToSum,axis=2)
        self.steppedMats[runId][fixedRingAmount] = np.sum(binMatsToSum,axis=2)
        ## Optionally produce ;-separated output.
        if self.useOutfile == True:
            for rtI,bm in enumerate(self.binaryMats):
                rt = RINGINESSTHRESHOLDS[rtI]
                for (ring,frame),val in np.ndenumerate(bm):
                    self.outFile.write(';'.join([str(rt),''.join(runIds),str(frame),str(ring),str(int(val))])+'\n')
    
    def shiftAnalysisSingle(self,mat,savgol):
        print 'Savitzky-Golay filter?'
        print savgol
        if savgol:
            mat = signal.savgol_filter(mat,window_length = 7, polyorder = 3, axis = 1)
        shiftedMats = diagonalShifter.shiftMat(mat,self.angles)
        SHIFTED_MATS.append(shiftedMats)
        angleSignificance = diagonalShifter.angleSignificance(shiftedMats,product=False)
        return angleSignificance

    def aggregatedSignificancePlot(self):
        imshow(self.aggregatedSignificances,interpolation='nearest')#
#        for meep in self.aggregatedSignificances:
#            plot([i for i,x in enumerate(meep)],meep)  
        
    def performPerRunStats(self,runId,mats,fixedRingAmount,normalize,deviancePerRow):
        mat = mats[runId][fixedRingAmount]
        metadata = self.metadata[runId]
        rawMean = np.sum(mat) / float(np.product(mat.shape))
        multiMat = np.ones(mat.shape,dtype=np.float64)
        multiMat = multiMat * rawMean
        rowsWithData = np.sum(mat,axis=0) > 0
        runTime = np.sum(rowsWithData > 0) * self.coincidenceWindow ## in ms
        events = np.sum(mat)
        eventsPerCell = float(events) / runTime
        eventsPerCellPerSecond = eventsPerCell / (runTime*1000.) ## Hz
#        unnormed = copy.copy(mat)
        normed = mat / multiMat
        if normalize:
            mat = normed
        BesselCorrectedNumerator = np.product(mat.shape) - 1
        mean = np.sum(mat) / float(np.product(mat.shape))
        if not deviancePerRow:
            deviance = (mat - mean) ** 2
        else:
            meanDenominatorMat = (
                (np.ones(mat.shape[1]) * mat.shape[0])
            )
            meansPerTimebin = np.sum(mat,axis=0) / meanDenominatorMat
            meansPerTimebin = np.tile(meansPerTimebin,(mat.shape[0],1))
            deviancesPerTimebin = (mat - meansPerTimebin)**2
            deviance = np.sum(deviancesPerTimebin)
        sd = np.sqrt(np.sum(deviance) / float(BesselCorrectedNumerator))
        length,circ = metadata['networkParameters']['length_circumference']
        N = length * circ
        Lambda = .022 ## lambda sans cap is a keyword. 
        ## refractory = 20ms, delay = 2ms
        tau_r = .1 ## events per time unit - time unit is s; .1Hz so .1 events per s.
        theoreticalVariance = N * (1./(1.+Lambda*tau_r))**2*Lambda*tau_r
        M1 = np.sum(deviance) / theoreticalVariance
        
        expectedVariance = N * (1./(1.+eventsPerCellPerSecond))**2*eventsPerCellPerSecond
        M2_ = np.sum(deviance) / expectedVariance
        M2 = np.sum(deviance) / (expectedVariance*(1-expectedVariance)/N)
        print 'MikadoDistance: {0}\nNoise: {1}\nTotal number of cell-firings in sample: {2}\nMean: {3}\nDeviance: {4}\nSD: {5}\nVar(N_r): {6}\nM1: {7}\nM2_: {8},\nM2: {9}\n\n'.format(
            metadata['networkParameters']['excMikadoDistance'],
            metadata['cellModelParameters']['noise'],
            np.sum(mat),
            mean,
            np.sum(deviance),
            sd,
            theoreticalVariance,
            M1,
            M2_,
            M2
        )
#        print {'sd':sd,'VarObs':np.sum(deviance),'M1':M1,'M2_':M2_,'M2':M2}
        return {'sd':sd,'VarObs':np.sum(deviance),'M1':M1,'M2_':M2_,'M2':M2}

    def tableauCSV(self,stat,normalizeOverallSignificance,deviancePerRow):
        ## setup csv outfile
        outCSV = open(self.targetDir+'mikadoPatternedness.csv','w')
        outCSV.write('')
        outCSV.close()
        outCSV = open(self.targetDir+'mikadoPatternedness.csv','a')
        dimensionKeys = []
        runID, md = self.metadata.items()[0]
        for superkey in md.keys():
            if type(md[superkey]) == dict:
                for subkey in md[superkey]:
                    dimensionKeys.append((superkey,subkey))
                    catchHack = md[superkey][subkey]
                    
            else:
                dimensionKeys.append((superkey))
        outCSVHeaderLine = []
        for dimKey in dimensionKeys:
            if type(dimKey) == tuple:
                outCSVHeaderLine.append('_'.join(dimKey))
            else:
                outCSVHeaderLine.append(dimKey)
        outCSVHeaderLine.extend(['runID','sd8','sd16'])    
        outCSV.write('\t'.join(outCSVHeaderLine)+'\n')
        for runID, md in self.metadata.items():
            print 'tableauing run {0}'.format(runID)
            print md.items()
            sd16 = self.performPerRunStats(
                        runID, self.mats, 16,
                        normalize=normalizeOverallSignificance,
                        deviancePerRow = deviancePerRow
                    )['sd']
            sd8 = self.performPerRunStats(
                        runID, self.mats, 8,
                        normalize=normalizeOverallSignificance,
                        deviancePerRow = deviancePerRow
                    )['sd']
            outList1 = []
            outList2 = []
            for dimKey in dimensionKeys:
                if type(dimKey) == tuple:
                    outList1.append(str(md[dimKey[0]][dimKey[1]]))
                    outList2.append(str(md[dimKey[0]][dimKey[1]]))
                else:
                    outList1.append(str(md[dimKey]))
                    outList2.append(str(md[dimKey]))
            outList1.append(str(runID))
#            outList2.append(str(runID))
            outList1.append(str(sd8))
#            outList2.append(str(sd8))
            outList1.append(str(sd16))
#            outList2.append(str(sd16))
            outCSV.write('\t'.join(outList1)+'\n')
#            outCSV.write('\t'.join(outList2)+'\n')
                

        outCSV.close()
                
                
        
    def stripePlot(self,picTypes,measurement,stat,shiftAnalysis,shiftSavgol,savgol,
                   normalizeOverallSignificance,deviancePerRow,
                   errorbarDimension = None,
                    limitPic=None,
                    stripePlotScale = 'overall',
#                   xAxis=['cellModelParameters','noise'],
                    yAxis=['networkParameters','length_circumference'],
                    xAxis=['networkParameters','excMikadoDistance']

        ):
        if measurement == 'abs':
            mats = self.mats
        elif measurement == 'deviance':
            mats = self.devianceMats
        elif measurement == 'stepped':
            mats = self.steppedMats
        else:
            raise ValueError('Unknown measurement: {0}'.format(measurement))
        ## Create a location for image files.
        figOutputDirPath = self.targetDir+'figures/'
        if not os.path.exists(figOutputDirPath):
            os.mkdir(figOutputDirPath)
#        print self.metadata
        ## start Matplotlib work
        stripeFig = plt.figure(figsize = (18,12))
        self.figs.append(stripeFig)
        if shiftAnalysis:
            shiftFig = plt.figure(figsize = (18,12))
#        xAxisOrderMap = sorted([(md[xAxis[0]][xAxis[1]],int(i)) for i,md in self.metadata.items()])
        yAxisValues = sorted(
            list(set([(md[yAxis[0]][yAxis[1]]) for i,md in self.metadata.items()]))
        )
        xAxisValues = sorted(
            list(set([(md[xAxis[0]][xAxis[1]]) for i,md in self.metadata.items()]))
        )
        if errorbarDimension == 'perBinNorm':
            errorbarDimensionValues = ['overall','perBin']
        elif errorbarDimension == 'ringAmount':
            errorbarDimensionValues = self.fixedRingAmount
        elif errorbarDimension is not None:
            errorbarDimensionValues = sorted(
                list(set([(md[errorbarDimension[0]][errorbarDimension[1]]) for i,md in self.metadata.items()]))
                )
#            print errorbarDimensionValues
#        if errorbarDimension != 'ringAmount':
        defaultRingAmount = self.fixedRingAmount[0]
        aggregatedStatMat2 = np.zeros((len(yAxisValues),len(xAxisValues)),dtype=np.float32)
        aggregatedStatListofLists = [[[] for x in xAxisValues] for y in yAxisValues]
        if errorbarDimension is not None:
            errorbarStatListofLists = [[[[] for x in xAxisValues] for y in yAxisValues] for dim in errorbarDimensionValues]  
        else:
            errorbarStatListofLists = [[[[] for x in xAxisValues] for y in yAxisValues]]
#        instanceCounter = np.zeros((len(yAxisValues),len(xAxisValues)),dtype=np.float32)
#        print yAxisValues,xAxisValues
#        for xOrderI,(xAxisVal,i) in xAxisOrderMap:
            
        ## get the maximum number of events.
        imshowVmaxes = []
        for yIndex,yVal in enumerate(yAxisValues):
            thisYBinMetadata = dict((runId,md) for runId,md in self.metadata.items() if md[yAxis[0]][yAxis[1]] == yVal)
#            print   len(thisYBinMetadata),yVal
            for xIndex,xVal in enumerate(xAxisValues):
                theseRunIds = [runId for runId,md in thisYBinMetadata.items() if md[xAxis[0]][xAxis[1]] == xVal]
#                assert len(thisRunId) == 1, 'Aggregation not supported for stripePlot. Make sure every cell has only 1 run.'
                run0Id = theseRunIds[0]                
                mat = mats[run0Id].values()[0]
                if savgol == 'viz':
                    mat = signal.savgol_filter(mat,window_length = 5, polyorder = 2, axis = 0)                
                if limitPic != None:
                    vizMat = mat[:,limitPic[0]:limitPic[0]+limitPic[1]]
                else:
                    vizMat = mat
                imshowVmaxes.append(np.max(mat))
        imshowVmax = max(imshowVmaxes)
        ## Outer loop: y axis
        for yIndex,yVal in enumerate(yAxisValues):
            thisYBinMetadata = dict((runId,md) for runId,md in self.metadata.items() if md[yAxis[0]][yAxis[1]] == yVal)
#            print   len(thisYBinMetadata),yVal
            for xIndex,xVal in enumerate(xAxisValues):
                theseRunIds = [runId for runId,md in thisYBinMetadata.items() if md[xAxis[0]][xAxis[1]] == xVal]
#                assert len(thisRunId) == 1, 'Aggregation not supported for stripePlot. Make sure every cell has only 1 run.'
                run0Id = theseRunIds[0]                
                mat = mats[run0Id].values()[0]
                if savgol == 'viz':
                    mat = signal.savgol_filter(mat,window_length = 5, polyorder = 2, axis = 0)                
                if limitPic != None:
                    vizMat = mat[:,limitPic[0]:limitPic[0]+limitPic[1]]
                else:
                    vizMat = mat
                ## Do the overall significance thing
                if savgol == 'overall':
                    mat = signal.savgol_filter(mat,window_length = 5, polyorder = 2, axis = 0)   
#                print theseRunIds
                stats = [
                    self.performPerRunStats(
                        thisRunId, mats, defaultRingAmount,
                        normalize=normalizeOverallSignificance,
                        deviancePerRow = deviancePerRow
                    )[stat] 
                    for thisRunId in theseRunIds
                ]
#                aggregatedStatListofLists[yIndex][xIndex].append(stats)
                aggregatedStatListofLists[yIndex][xIndex] = stats
                if errorbarDimension == 'perBinNorm':
                    overallStat = [
                            self.performPerRunStats(
                                thisRunId, mats, defaultRingAmount,
                                normalize=normalizeOverallSignificance,
                                deviancePerRow = False
                            )[stat] 
                            for thisRunId in theseRunIds
                        ]
                    perBinStat = [
                            self.performPerRunStats(
                                thisRunId, mats, defaultRingAmount,
                                normalize=normalizeOverallSignificance,
                                deviancePerRow = True
                            )[stat] 
                            for thisRunId in theseRunIds
                        ]
                    errorbarStatListofLists[0][yIndex][xIndex] = overallStat
                    errorbarStatListofLists[1][yIndex][xIndex] = perBinStat
                elif errorbarDimension == 'ringAmount':
                    for ebI,errorbarDimensionValue in enumerate(errorbarDimensionValues):
                        errorbarStats = [
                            self.performPerRunStats(
                                thisRunId, mats, errorbarDimensionValue,
                                normalize=normalizeOverallSignificance,
                                deviancePerRow = deviancePerRow
                            )[stat] 
                            for thisRunId in theseRunIds
                        ]
                        errorbarStatListofLists[ebI][yIndex][xIndex] = errorbarStats

                elif errorbarDimension is not None:
                    for ebI,errorbarDimensionValue in enumerate(errorbarDimensionValues):
                        errorbarRunIds = [
                            runId 
                            for runId,md
                            in thisYBinMetadata.items() 
                            if (md[errorbarDimension[0]][errorbarDimension[1]] == errorbarDimensionValue) and
                                (runId in theseRunIds)
                        ]
                        errorbarStats = [
                            self.performPerRunStats(
                                thisRunId, mats, defaultRingAmount,
                                normalize=normalizeOverallSignificance,
                                deviancePerRow = deviancePerRow
                            )[stat] 
                            for thisRunId in errorbarRunIds
                        ]
                        errorbarStatListofLists[ebI][yIndex][xIndex] = errorbarStats
                else:
                    errorbarStatListofLists[0][yIndex][xIndex] = stats
#                print 'len(stats)'
#                print len(stats)
#                print 'stats'
#                print stats
                
#                print aggregatedStatListofLists
                aggregatedStatMat2[yIndex,xIndex] = stats[0]
#                print xIndex+1,yIndex+1
#                print (len(yAxisValues),len(xAxisValues),yIndex*len(xAxisValues)+xIndex+1)
                ax1 = stripeFig.add_subplot(len(yAxisValues),len(xAxisValues),yIndex*len(xAxisValues)+xIndex+1)
                if stripePlotScale is not 'overall':
                    imshowVmax = np.max(vizMat)
                strplt = ax1.imshow(
                    vizMat,interpolation='nearest',cmap=plt.cm.gray,
                    vmin = 0.,
                    vmax = imshowVmax)
                ax1.set_ylabel(ylabel='{0}'.format(yVal))
                ax1.set_xlabel(xlabel='{0}'.format(xVal))
                ax1.set_yticks([])
                if stripePlotScale is not 'overall':
                    ax1.set_title('Black = 0 spikes; white = {0} spikes'.format(imshowVmax))
                if shiftAnalysis:
                    ax1 = shiftFig.add_subplot(len(yAxisValues),len(xAxisValues),yIndex*len(xAxisValues)+xIndex+1)
                    ax1.plot(self.shiftAnalysisSingle(mat,shiftSavgol))#,cmap=plt.cm.gray)
                    ax1.set_ylabel(ylabel='{0}'.format(yVal))
                    ax1.set_xlabel(xlabel='{0}'.format(xVal))
#                cbar = plt.colorbar(strplt,ax=ax1,extendfrac = [.01,.01])
#                cbar.set_ylim(bottom = 0.)
        if stripePlotScale == 'overall':
            stripeFig.suptitle('Black = 0 spikes; white = {0} spikes'.format(imshowVmax))
        overallSigPlot = plt.figure(figsize = (len(yAxisValues),len(xAxisValues)))
        boxPlot = plt.figure(figsize = (len(yAxisValues)*2,len(xAxisValues)*2))
        errorbarPlot = plt.figure(figsize = (len(yAxisValues)*2,len(xAxisValues)*2))      
        yAxisLength = len(yAxisValues) 
#        print aggregatedStatListofLists
        errorbarLegendList = []
        errorbarLegendListNames = []
        for yIndex,yVal in enumerate(yAxisValues):

            subplotString = '{0}1{1}'.format(yAxisLength,yIndex+1)
            ax = boxPlot.add_subplot(subplotString)
            ax.boxplot(aggregatedStatListofLists[yIndex])
            ax.set_ylim(ymin=0.)
            print '....'
            ax.set_xticklabels(xAxisValues)
            ebx = errorbarPlot.add_subplot(subplotString)
            print aggregatedStatListofLists[yIndex]
            print 'boxplotinput'
            print np.average(aggregatedStatListofLists[yIndex],axis=1)
            print xAxisValues
            if errorbarDimension is not None:
                for dimIndex,dimVal in enumerate(errorbarDimensionValues):
                    curMarker = MARKERS[dimIndex]
                    curLinestyle = LINESTYLES[dimIndex%2]
                    eb = ebx.errorbar(
                        x = np.array(xAxisValues), 
                        y = np.average(errorbarStatListofLists[dimIndex][yIndex],axis=1),
                        yerr=std(errorbarStatListofLists[dimIndex][yIndex],axis=1),  
                        linestyle=curLinestyle,
                        marker = curMarker
                    )
                    errorbarLegendList.append(eb)
                    errorbarLegendListNames.append(str(dimVal))
            else:
                eb = ebx.errorbar(
                    x = np.array(xAxisValues), 
                    y = np.average(errorbarStatListofLists[0][yIndex],axis=1),
                    yerr=std(errorbarStatListofLists[0][yIndex],axis=1),  
                    fmt='-o'
                )
                
            ebx.set_ylim(bottom = 0.)
            ebx.set_xticks(xAxisValues)
            ebx.set_xlim(left = xAxisValues[0]-.1, right = xAxisValues[-1]+.1)
        
        ax1 = overallSigPlot.add_subplot('211')
#        ax2 = overallSigPlot.add_subplot('212')
#        aggregatedStatMat = np.array(aggregatedStatListofLists)
#        aggregatedStatMat = np.average(aggregatedStatMat,axis=2)
#        print aggregatedStatMat.shape
        ax1.imshow(aggregatedStatMat2,interpolation='nearest',cmap=plt.cm.gray)
#        ax2.imshow(aggregatedStatMat2,interpolation='nearest',cmap=plt.cm.gray)
        overallSigPlot.show()   
        filterNamePart = ''.join([
            '_'+'fVal'+str(i)+'_'+str(val[1]) for i,val in enumerate(self.filters)
        ])
        errorbarPlot.legend(
            errorbarLegendList,errorbarLegendListNames,'upper right')
        errorBarDimensionNamePart = str(errorbarDimension)
        if 'boxplot' in picTypes:
            boxPlot.savefig(figOutputDirPath+'BOXPLOT'+filterNamePart+'.png')
            boxPlot.savefig(figOutputDirPath+'BOXPLOT'+filterNamePart+'.eps')
            boxPlot.savefig(figOutputDirPath+'BOXPLOT'+filterNamePart+'.pdf')
            boxPlot.savefig(figOutputDirPath+'BOXPLOT'+filterNamePart+'.jpg')
            boxPlot.show()
        if 'stripeFig' in picTypes:
            stripeFig.savefig(figOutputDirPath+'stripeFig'+filterNamePart+'.png')
            stripeFig.savefig(figOutputDirPath+'stripeFig'+filterNamePart+'.eps')
            stripeFig.savefig(figOutputDirPath+'stripeFig'+filterNamePart+'.pdf')
            stripeFig.savefig(figOutputDirPath+'stripeFig'+filterNamePart+'.jpg')
            stripeFig.show()
        if 'errorbarPlot' in picTypes:
            errorbarPlot.savefig(figOutputDirPath+'errorBar'+errorBarDimensionNamePart+filterNamePart+'.png')
            errorbarPlot.savefig(figOutputDirPath+'errorBar'+errorBarDimensionNamePart+filterNamePart+'.eps')
            errorbarPlot.savefig(figOutputDirPath+'errorBar'+errorBarDimensionNamePart+filterNamePart+'.pdf')
            errorbarPlot.savefig(figOutputDirPath+'errorBar'+errorBarDimensionNamePart+filterNamePart+'.jpg')
            errorbarPlot.show()
        
if __name__ == '__main__':
    ## Initialize an analysis object.
    tiger = StripeAnalysis(
#            './output/supp3/',
#            './output/3x3x4/',
            './output/TS_incNone/',
            './topologyCaches/',
#            useOutfile = True,
            coincidenceWindow = 3.,
            angles = 37,
            getFromCache = True,
            fixedRingAmount = [16,8]
#            ,filters = [
#                [['networkParameters','transmissionSpeed'],[.01]]
#            ]
        )
    ## create tableauable output
    tiger.tableauCSV(stat='sd', normalizeOverallSignificance = True,
        deviancePerRow = True)
    
    ## Perform the stripey visualizations & an analysis visualization
    tiger.stripePlot(
        picTypes = ['boxplot','stripeFig', 'errorbarPlot'],
        measurement='abs',
        stat='sd',
        normalizeOverallSignificance = True,
        deviancePerRow = True,
        shiftAnalysis = False,
        shiftSavgol = False,
        savgol = None,
        limitPic = (150,150),
        stripePlotScale = 'overall'
        , errorbarDimension = ['networkParameters','excMikadoDensity']
#        , errorbarDimension = 'perBinNorm'
#        , errorbarDimension = 'ringAmount'
    )