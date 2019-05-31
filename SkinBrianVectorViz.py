# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 04:21:02 2015

@author: oddman
"""
import os,json, cProfile, subprocess
#from SkinBrian import *
from ShapesOnAPlane import *
from matplotlib.font_manager import FontProperties
from matplotlib.figure import SubplotParams
from matplotlib.patches import Ellipse, Arrow
from matplotlib.lines import Line2D
import multiprocessing
import numpy as np

RAMONYCAJALCOUNT = 15

class SkinBrainComic:
    def __init__(self,fromTime,toTime,timeStep,networkTopology,
                 spikes,vizOutputDirPath):
        self.maxMs = max([
            max(spikeTrain+[0]) 
            for cell_i,spikeTrain 
            in spikes.items()
        ])
        self.fromTime = fromTime
        self.toTime = toTime
        self.timeStep = float(timeStep)
        self.spikes = spikes
        self.vizOutputDirPath = vizOutputDirPath
        self.networkTopology = networkTopology
#        print spikeSort(self.spikes.items())
        print self.maxMs
        self.maxTime = int(self.maxMs) + 1
        self.alignEvents()
    
    def alignEvents(self):
        ## the format of spikes produced by SingleRunDataReader is a list of
        ## lists, spike times of cells.
        ## the time of the last spike in ms/binsize, rounded up
        intervalCount = int(self.maxTime / self.timeStep) + 1
        self.cellsFiringPerMoment = dict((i,[]) for i in range(intervalCount))
        ## for every spike, put the cellindex in the appropriate time-bin
        [
            [
                self.cellsFiringPerMoment[int(spikeTime / self.timeStep)].append(cellIndex)
                for spikeTime
                in thisCellSpikes
            ]
            for cellIndex,thisCellSpikes
            in self.spikes.items()
        ]
#        print self.cellsFiringPerMoment
#        self.linksfiringpermoment = dict((i,[]) for i in range(self.maxtime))
#        [
#            [
#                self.linksfiringpermoment[int(moment/self.binsize)].append(link)
#                for moment
#                in momentlist
#            ]
#            for (link,momentlist)
#            in self.thisrun.firinglinks.items()
#        ]
    def run(self,orientation,outputFormat,CPUCount,autoplay):
        firstInterval = int(self.fromTime / self.timeStep)
        lastInterval = int(self.toTime / self.timeStep)
        plottedIntervals = range(firstInterval,lastInterval)
        print 'Plotting {0} frames.'.format(len(plottedIntervals))
        movieVizOutputDirPath = self.vizOutputDirPath + 'movie/'
        if os.path.exists(movieVizOutputDirPath):
            pass
        else:
            os.mkdir(movieVizOutputDirPath)
#        pool = multiprocessing.Pool(processes = CPUCount)
        frameArgumentList = [
            [   
                intervalKey,outputFormat,movieVizOutputDirPath,
                self.networkTopology,self.cellsFiringPerMoment
            ] 
            for intervalKey in plottedIntervals
        ]
        for frameArguments in frameArgumentList:
            sfp = multiprocessing.Process(
                target = singleFrame, 
                args = frameArguments
            )
            sfp.start()
            sfp.join()
#        pool.map(frameMultiprocessHelper,poolArgumentList)
        singleFrame(
            intervalKey,outputFormat,
            self.vizOutputDirPath,self.networkTopology,
            self.cellsFiringPerMoment,cellSize = .25, 
            links = 'ramon_y_cajal')
        os.chdir(movieVizOutputDirPath)
#        intervals = len(plottedIntervals)
        if outputFormat == 'png':
            if os.path.exists('./skinBrainMovie.avi'):
                os.remove('./skinBrainMovie.avi')
            ffmpegCommandList = [
                    'ffmpeg',
                    '-f',
                    'image2',
                    '-framerate',
                    '30',
                    '-pattern_type',
                    'sequence',
                    '-start_number',
                    '{0}'.format(firstInterval),
                    '-r',
                    '12',
                    '-i',
                    '%05dvector_wires.png',
                    'skinBrainMovie.avi'
            ]
            subprocess.call(ffmpegCommandList)
            mplayerCommandList = ['gnome-mplayer','skinBrainMovie.avi']
            if autoplay:
                subprocess.call(mplayerCommandList)
        
def singleFrame(
        intervalKey,outputFormat,outputDirPath,networkTopology,
        cellsFiringPerMoment,cellSize=.75,links = False
    ):    
    """A function producing a single wire-frame 2d-representation
    of a skin brain at a certain time in the simulation. 
    """
    print intervalKey
    timestr = "%05d"%(intervalKey)
    if links == 'all':
        linkstr = 'links_'
    elif links == 'ramon_y_cajal':
        linkstr = 'ramon_y_cajal_links_'
    else:
        linkstr = ''
    graphName = timestr+linkstr+'vector_wires.'+outputFormat
    if os.path.exists(outputDirPath+graphName):
        return graphName
    interCellDistance = 5.
    inch = 25.4 ## mm
    ## size should be relative to length (x) and circumference (y).
    figXSize = (networkTopology.networkParameters['length_circumference'][0]*interCellDistance*.5*sqrt(2)) / inch
    figYSize = (networkTopology.networkParameters['length_circumference'][1]*interCellDistance) / inch
    fig = figure(figsize=(figXSize,figYSize), subplotpars=SubplotParams(left=0,right=1,top=1,bottom=0))
#        nodes = [i for i in enumerate(self.coords2d)]
#        print str(moment)
    ax = fig.add_subplot(111,frame_on=False,label=intervalKey)
    nodes = [(cell.overallIndex,cell.flatCoordinates) for cell in networkTopology.cells]      
    if links == 'all' or links == 'ramon_y_cajal':
        if links == 'ramon_y_cajal':
            inCells = np.random.choice(networkTopology.cells,RAMONYCAJALCOUNT)
        for link in networkTopology.links:
            (origin,target) = (
                link.fromCell,
                link.toCell
            )
            if links == 'all' or (links == 'ramon_y_cajal' and origin in inCells):
                thisline = Line2D(
                    *zip(*[list(origin.flatCoordinates),list(target.flatCoordinates)]),
                    linewidth=1.5,
                    solid_capstyle = 'butt')
                ax.add_artist(thisline)
                thisline.set_clip_box(ax.bbox)
                thisline.set_color([0.8,0.8,0.8])  
    #                    print self.colourdict[link]
    #            linkprune.__delitem__(link)
                thisline.set_zorder(1)
#            otherlines = [
#                Line2D(
#                    *zip(*[list(origin),list(target)]),
#                    linewidth=1.5,
#                    solid_capstyle = 'butt')
#                for (link,(origin,target)) 
#                in linkprune.items()]
#            for e in otherlines:
#                ax.add_artist(e)
#                e.set_clip_box(ax.bbox)
#                e.set_color(linkcolour)
#                e.set_zorder(0)
  
    ellipses = [(i,Ellipse(xy,cellSize,cellSize)) for (i,xy) in nodes]
    for (i,e) in ellipses:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        if i not in cellsFiringPerMoment[intervalKey]:
            e.set_color([0.9,0.9,0.9])
        else:
            e.set_color([0.5,0.0,0.0])
        if links == 'all':
             e.set_color([0.4,0.4,0.4])
        e.set_zorder(2)
    ax.tick_params(top=False,bottom=False,left=False,right=False)
    ax.set_xlim(-1, max([x for (i,(x,y)) in nodes]) + 1)
    ax.set_ylim(-1, max([y for (i,(x,y)) in nodes]) + 1)
    ax.set_yticklabels('')
    ax.set_xticklabels('')
#        linkstr = "_links" + repr(links)
    graphName = timestr+linkstr+'vector_wires.'+outputFormat
    savefig(outputDirPath+graphName,format=outputFormat)
    close('all')
#    fig.show()
    return graphName

def generateComic(
        fromTime,toTime,timeStep,targetFolder,topologyCachePath,runIndex,
        CPUCount,autoplay=True,show=True,outputFormat='png'
    ):
#    print os.listdir(targetFolder)
    metadataPath = targetFolder.strip('/') + '/metadata.json'
    vizOutputDirPath = targetFolder.strip('/') + '/viz_{0}/'.format(runIndex)
    if os.path.exists(vizOutputDirPath):
        pass
    else:
        print os.getcwd()
        os.mkdir(vizOutputDirPath)
    with open(metadataPath,'r') as metadataFile:
        metadata = json.load(metadataFile)
#        print self.metadata['0']['cellModelParameters']
    key = str(runIndex)
    networkParameters = metadata[key]['networkParameters']
    networkType = networkParameters['networkType']
    networkTopology = NetworkTopology(
        networkType,
        cacheDir = topologyCachePath
    )
    networkTopology.generate(
        networkParameters,
        cache=True
    )
#    print networkTopology.networkParameters
    spikeFilePath = targetFolder.strip('/') + '/spikefile_'+key+'.json'
    with open(spikeFilePath,'r') as spikeFile:
        spikes = dict(json.load(spikeFile))
    skinBrainComic = SkinBrainComic(fromTime,toTime,timeStep,networkTopology,spikes,vizOutputDirPath)
    skinBrainComic.run('vertical',outputFormat,CPUCount,autoplay)
    
        
if __name__ == '__main__':
    CPUCount = 1
#    thisRunFolder = './output/'+max(os.listdir('./output/'))
    thisRunFolder = './output/fullRun'
    print thisRunFolder
    generateComic(0,1000,2,thisRunFolder,'./topologyCaches/',1930,CPUCount,autoplay=False,outputFormat='png')