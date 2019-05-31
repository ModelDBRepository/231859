# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 00:22:34 2014

TODO:
* parallelize

@author: Oddman
"""

## import builtins whole
import os, collections, itertools, json, time, shutil, logging, zipfile
## ... or with an alias
import cPickle as pkl
## ... or partly
from multiprocessing import Pool
## Import brian
from brian2 import *
## Import packaged dependencies
from scipy.spatial import distance
#import visual as vp
## import local script dependencies
#import povexport

#==============================================================================
# Globals
#==============================================================================
FILE_LOG_LEVEL = logging.DEBUG
CONSOLE_LOG_LEVEL = logging.DEBUG
PARALLEL = True
CPUS = 2
codegen.target = 'cython'
#==============================================================================
# Define namedtuples
#==============================================================================
Link = collections.namedtuple(
            'Link',
            [
                'fromCell',
                'toCell',
                'fromCellIndex',
                'toCellIndex',
                'linkType',
                'linkWeight',
                'rhosum' ## actual distance over elongations
            ]
)

#def characterizer(
#        runProfile,
#    ):
#    x = 
##    print x
#    print json.dumps(x)

#==============================================================================
# Define functions
#==============================================================================

def dictifyArray(indices,spikeTimes,totalCells):
    '''Turns a pair of arrays containing spiketimes and indices into a dict
    with indices for keys and their spikes as values.'''
    targetDict = dict((i,[]) for i in range(totalCells))
#    print [indices,spikeTimes,totalCells]
    for i in range(len(indices)):
        targetDict[indices[i]].append(spikeTimes[i])
    return targetDict

def productProfiler(multiprofileDict):
    return None
    

def spawnTriangularGridTubularCells(
        length,
        circumference,
        excWeight,
        pacemakerType = None,
        excMikadoDensity = 0.,
        excMikadoDistance = None,
        inhibMikadoDensity = 0.,
        inhibMikadoDistance = None,
        inhibWeight = 0.
    ):
    cellTypes = {}
    cellTypeCount = 0
    inCellTypeIndices = {}
    c = 0
    l = 0
    nCells = length * circumference
    assert excMikadoDensity + inhibMikadoDensity <= 1.
    for i in range(nCells):
        if pacemakerType == 'singleCellExc' and i == 0:
            cellType = 'activeExcitatory'
#            inhib = False
            weight = excWeight
        elif pacemakerType == 'ringExc' and i < circumference:
            cellType = 'activeExcitatory'
#            inhib = False
            weight = excWeight
        else:
            cellType = 'passiveExcitatory'
            weight = excWeight
        mikadoRoll = rand()
        if mikadoRoll < excMikadoDensity:
            connectivity = 'excMikado'
            mikadoOrientation = (rand()*pi*2) - pi#(random.randint(-36,35) / 36.) * pi
            mikadoDistance = excMikadoDistance            
            weight = excWeight
        elif mikadoRoll < excMikadoDensity + inhibMikadoDensity:
            connectivity = 'inhibMikado'
            mikadoOrientation = (rand()*pi*2) - pi#(random.randint(-36,35) / 36.) * math.pi
            mikadoDistance = inhibMikadoDistance            
            weight = inhibWeight
        else:
            connectivity = None
            mikadoOrientation = None
            weight = excWeight
            mikadoDistance = None
        circumferenceIndex = c
        lengthIndex = l
        c += 1
        if c == circumference:
            c = 0
            l += 1
        try:
            cellTypeIndex = cellTypes[cellType]
        except KeyError:
            cellTypes[cellType] = cellTypeCount
            cellTypeIndex = cellTypes[cellType]
            cellTypeCount += 1
        try:
            inCellTypeIndices[cellType] += 1
        except KeyError:
            inCellTypeIndices[cellType] = 0
        radius = circumference / (2 * pi)  
        angleProp = (circumferenceIndex / float(circumference))
        angle = angleProp * 2 * pi
        if lengthIndex % 2 == 1:
            angle += (2 * pi) / (circumference * 2)
        x = sin(angle) * radius
        y = cos(angle) * radius
        z = lengthIndex * sqrt(.75)
        coordinates = (x,y,z)
        flatX = lengthIndex * sqrt(.75)
        flatY = circumferenceIndex
        if lengthIndex % 2 == 1:
            flatY += sqrt(.75) * .5 
        flatCoordinates = (flatX,flatY)
        cell = Cell(
                i,
                weight,
                coordinates,
                flatCoordinates = flatCoordinates,
                cellType = cellType,
                connectivity = connectivity, 
                mikadoOrientation = mikadoOrientation,
                mikadoDistance = mikadoDistance,
                cellTypeIndex = cellTypeIndex, 
                inCellTypeIndex = inCellTypeIndices[cellType]
            )
        yield cell

def modelRunWorker(cellModelParameters,networkParameters,_runParameters,_runCntr,_topologyCachePath,_outputDir,_useCache,_metadataPickler):
    '''Multiprocessing worker function to execute a SkinBrain model.'''
    thisRunMetadata = collections.OrderedDict()          
#    cellModelParameters = _runProfile[0]._asdict()
#    print cellModelParameters.keys()
#    networkParameters = _runProfile[1]._asdict()
    cellModelType = cellModelParameters['cellModelType']            
    networkType = networkParameters['networkType']
    ## set up some metadata for ease of file identification
    thisRunMetadata['cellModelParameters'] = cellModelParameters
    thisRunMetadata['networkParameters'] = networkParameters     
    thisRunMetadata['runParameters'] = _runParameters         
    networkTopology = NetworkTopology(
        networkType,
        cacheDir = _topologyCachePath
    )
    networkTopology.generate(networkParameters,cache=_useCache)
    thisRunMetadata['networkTopologyFileName'] = networkTopology.fileName
    thisRunMetadata['cacheDir'] = networkTopology.cacheDir
    cellModelImplementation = BrianCellMechanismImplementation(
        networkTopology,
        cellModelType,
        cellModelParameters,
        logger = None
    )
    json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')
    spikeFilePath = _outputDir + 'spikefile_%s.json'%_runCntr
    spikeZipFilePath = _outputDir + 'spikefile_%s.zip'%_runCntr
    indices,spikeTimes = cellModelImplementation.run(_runParameters)#.items()
    runOut = dictifyArray(indices,spikeTimes,len(networkTopology.cells))
    runOut = dict((key,np.array(val)) for key,val in runOut.items())
    with open(spikeFilePath,'wb') as spikeFile:
        json.dump(
            [
                (i,list(spikeArray*1000))
                for i,spikeArray 
                in runOut.items()
            ],
            spikeFile
        )
    with zipfile.ZipFile(spikeZipFilePath,'w',compression=zipfile.ZIP_DEFLATED) as zipSpikeFile:
        zipSpikeFile.write(spikeFilePath,'spikefile_%s.json'%_runCntr)
    os.remove(spikeFilePath)
    print thisRunMetadata
#    successFile = open(_outputDir+'last_successful_run_index','w')
#    successFile.write(str(_runCntr))
#    successFile.close()
    _metadataPickler.dumpMetadata(_runCntr,thisRunMetadata)
    return _runCntr
    

#==============================================================================
# Define classes
#==============================================================================
#class MetadataContainer(collections.OrderedDict):
#    def setDumpFilename(self,dumpFilename):
#        self.dumpFilename = dumpFilename
#        
#    def updateMetadata(self,key,val):
#        self.__setitem__(key,val)
#        self.dumpJSON()
#        
#    def dumpJSON(self):
#        with open(self.dumpFilename,'wb') as metadatafile:
#            metadatafile.write('')
#        with open(self.dumpFilename,'wb') as metadatafile:
#            json.dump(self,metadatafile,indent=4)
            
class MetadataPickler:
    def __init__(self,dumpFilenameBase):
        self.dumpFilenameBase = dumpFilenameBase
       
    def dumpMetadata(self,runCntr,thisRunMetadata):
        runNum = ('0000'+str(runCntr))[-5:]
        with open(self.dumpFilenameBase+runNum+'.pkl','wb') as outfile:
            pkl.dump((runCntr,thisRunMetadata),outfile)
        
#    def dumpJSON(self):
#        with open(self.dumpFilename,'wb') as metadatafile:
#            metadatafile.write('')
#        with open(self.dumpFilename,'wb') as metadatafile:
#            json.dump(self,metadatafile,indent=4)
        

class Cell:
    def __init__(self, 
         overallIndex, weight, coordinates, flatCoordinates = None, cellType = None, 
         inhib = False, connectivity = None, mikadoOrientation = None,
         mikadoDistance = None, cellTypeIndex = None, inCellTypeIndex = None
    ):
        self.overallIndex = overallIndex
        self.coordinates = coordinates
        [self.x,self.y,self.z] = self.coordinates
        self.flatCoordinates = flatCoordinates
        self.cellType = cellType
        self.connectivity = connectivity
        self.mikadoOrientation = mikadoOrientation
        self.mikadoDistance = mikadoDistance
        self.cellTypeIndex = cellTypeIndex
        self.inCellTypeIndex = inCellTypeIndex
        self.weight = weight
#        for k,v in self.__dict__.items():
#            print [k,v]


class NetworkTopology:
    def __init__(self,networkType, cacheDir = './'):
        self.networkType = networkType
        self.cacheDir = cacheDir.rstrip('/') + '/'
        
    def loadFromFile(self,cacheDir,fileName):
        print "Loading cached network topology from: " + cacheDir + fileName+'.zip'
        self.pathToCacheFile = cacheDir+fileName+'.zip'
        with zipfile.ZipFile(cacheDir+fileName+'.zip','r') as openZip:
            openFile = openZip.open(fileName+'.pkl')
            [
                    loadedNetworkType,
                    self.networkParameters,
                    self.cells,
                    self.links
            ] = pkl.load(openFile) 
            openFile.close()
        assert loadedNetworkType == self.networkType
        
    def generate(self,networkParameters, cache = True, cacheDir = None, fileName = None,createCache = True):
        self.networkParameters = networkParameters
        if fileName == None:
            self.fileName = [
                (k[:4] + str(v)).replace('.','') 
                for k,v in sorted(networkParameters.items())
            ]
            self.fileName = (
                'NetImp_' + self.networkType + '_' + 
                '_'.join(self.fileName)
            )
            assert len(self.fileName) <= 255
        else: self.fileName = fileName
        currentCacheDir = [cacheDir,self.cacheDir][cacheDir == None]
            
        if cache and os.path.exists(currentCacheDir + self.fileName+'.zip'):
            ## try loading previous version.
            self.loadFromFile(currentCacheDir,self.fileName)
        else:
            ## no cache or no cache found? run it :)
             self.buildTopology()
             ## and then build a cache file
             if createCache:
                 self.serialize()
          
    def buildTopology(self):
        if self.networkType == 'triangularGridTubular':
            ## turn the network parameters into object parameters of this
            ## NetworkTopologyImplementation.
            self.__dict__.update(self.networkParameters)
            [self.cells, self.links] = self.generateTriangularGridTubularNetwork()

                
    def generateTriangularGridTubularNetwork(self):
        #======================================================================
        # First, generate the nodes / cells.
        #======================================================================
        cellBuildingBlocks = []
        overallIndex = 0
        seed(self.networkDummySeed)
        self.length = self.length_circumference[0]
        self.circumference = self.length_circumference[1]
        cells = [
            i for i in 
            spawnTriangularGridTubularCells(
                self.length, self.circumference, self.excSynapseWeight,
                pacemakerType = self.pacemakerStructure,                
                excMikadoDensity = self.excMikadoDensity, 
                excMikadoDistance = self.excMikadoDistance,
                inhibMikadoDensity = self.inhibMikadoDensity,
                inhibMikadoDistance = self.inhibMikadoDistance,
                inhibWeight = self.inhibSynapseWeight
            )
        ]
        #======================================================================
        # Nodes done? Do the links / synapses.
        #======================================================================
        links = []
        for fromCellIndex,fromCell in enumerate(cells):
            for toCellIndex,toCell in enumerate(cells):
                autapse = fromCellIndex == toCellIndex     ## no autapses
                doubleCheck = toCellIndex < fromCellIndex  ## no double checking
                if not autapse and not doubleCheck:
                    rowDist = absolute(fromCell.z - toCell.z)
                    if rowDist <= self.connectionDistance:
                        linkDistance = distance.euclidean(
                                fromCell.coordinates,toCell.coordinates)
                        if linkDistance <= self.connectionDistance:
                            links.append(
                                Link(
                                    fromCell,
                                    toCell,
                                    fromCellIndex,
                                    toCellIndex,
                                    'pancake',
                                    fromCell.weight,
                                    1 ## distance measure
                                )
                            )
                            links.append(
                                Link(
                                    toCell,
                                    fromCell,
                                    toCellIndex,
                                    fromCellIndex,
                                    'pancake',
                                    toCell.weight,
                                    1 ## distance measure
                                )
                            )
                    #==========================================================
                    # Add the mikado-links.
                    #==========================================================   
                    if (    fromCell.connectivity in ['excMikado','inhibMikado'] and 
                            toCell.connectivity in ['excMikado','inhibMikado']
                    ):
                        mikadoLinks = self.mikadoIze(
                            fromCell,
                            toCell)
                        if mikadoLinks != None:
                            links.extend(mikadoLinks)
        return [cells,links]
        
    def serialize(self): 
        if not os.path.exists(self.cacheDir):
            os.mkdir(self.cacheDir)
        self.pathToCacheFile = self.cacheDir + self.fileName+'.zip'
        with open(self.cacheDir + self.fileName+'.pkl','wb') as openFile:
            pkl.dump(
                [
                    self.networkType,
                    self.networkParameters,
                    self.cells,
                    self.links
                ],
                openFile
            )
        with zipfile.ZipFile(self.pathToCacheFile,'w',compression=zipfile.ZIP_DEFLATED) as zipCache:
            zipCache.write(self.cacheDir + self.fileName+'.pkl',self.fileName+'.pkl')
        os.remove(self.cacheDir + self.fileName+'.pkl')
        
        
            
    def mikadoIze(self,fromCell,toCell):
        weight = False
        maxWinding=int(max([fromCell.mikadoDistance,toCell.mikadoDistance])/self.circumference)+1
        if fromCell.overallIndex < toCell.overallIndex:
            pass
        else:
            fromCell,toCell = toCell,fromCell
        fromCellDistance = fromCell.mikadoDistance
        toCellDistance = toCell.mikadoDistance#[self.excMikadoDistance,self.inhibMikadoDistance][toCell.connectivity == 'inhibMikado']
        totalDistance = fromCellDistance + toCellDistance
        fromCellX = fromCell.flatCoordinates[0]
        fromCellY = fromCell.flatCoordinates[1]
        fromCellPhi = fromCell.mikadoOrientation
        toCellX = toCell.flatCoordinates[0]
        toCellY = toCell.flatCoordinates[1]
        toCellPhi = toCell.mikadoOrientation
        xDistance = toCellX - fromCellX
        yDistance = toCellY - fromCellY
#        print 'distances: '+repr([xDistance,yDistance])
        for winding in range(-maxWinding,maxWinding+1,1):
            distance = sqrt((yDistance+(winding*self.circumference))**2 + xDistance**2)
            phid = arctan2(yDistance+(winding*self.circumference),xDistance)
#            phid = math.atan2(xDistance,yDistance+(winding*self.circumference))
            ## Only ever do anything if distance < 2 * mikadoDistance
            if distance < totalDistance:
                 ##analyse degenerate cases separately
                if (    (toCellPhi == (fromCellPhi - pi)) or 
                        (toCellPhi == (fromCellPhi + pi)) or
                        (toCellPhi == fromCellPhi) 
                ):
#                    print 'degenerate case'
                    if (    (phid == fromCellPhi - pi) or 
                            (phid == fromCellPhi + pi)
                    ):
                        if (    (   (phid == toCellPhi - pi) or 
                                    (phid == toCellPhi + pi)
                                ) and 
                                distance < fromCellDistance 
                        ):
                            weight = True
                            cache_rhosum = distance	
                    elif (phid == fromCellPhi):
                        if (distance < self.excMikadoDistance):
                            weight = True
                            cache_rhosum = distance
                        elif (  (   (phid == toCellPhi - pi) or 
                                    (phid == toCellPhi + pi)
                                ) and distance < totalDistance
                        ):
                            weight = True
                            cache_rhosum = distance
                else:
#                    print 'non-degenerate case.'
                    rho1 = (
                        (distance * sin(phid - toCellPhi)) 
                        / 
                        (sin(toCellPhi - fromCellPhi))
                    )
                    rho2 = (
                        (distance * sin(phid - fromCellPhi)) 
                        / 
                        (sin(toCellPhi - fromCellPhi))
                    )
#                    print 'rhos: '+repr([rho1,rho2])
                    if (    (rho1 < fromCellDistance) 
                            and 
                            (rho2 < toCellDistance)
                        ):
                        ## Check whether crossing is on the cylinder
                        if(-0.5 <= fromCellX+rho1*cos(toCellPhi) <= (self.length-1)*sqrt(0.75)+0.5): 
                            weight = True
                            cache_rhosum = abs(rho1)+abs(rho2)	
        if weight != False:      
#            print 'some result.'  
            return [Link(
                    fromCell,
                    toCell,
                    fromCell.overallIndex,
                    toCell.overallIndex,
                    'mikado',
                    fromCell.weight,
                    cache_rhosum
            ),
                Link(
                    toCell,
                    fromCell,
                    toCell.overallIndex,
                    fromCell.overallIndex,
                    'mikado',
                    toCell.weight,
                    cache_rhosum
            ),
            ]
#        else:
##            print 'No dice.'
#            pass
        
        
class BrianCellMechanismImplementation:
    def __init__(
            self,
            networkTopology,
            cellMechanismType,
            cellMechanismParameters,
            logger = None
    ):
        self.cellMechanismType = cellMechanismType
        self.networkTopology = networkTopology
        self.cellMechanismParameters = cellMechanismParameters
        if logger == None:
            self.logger = logging.getLogger('SkinBrian_2')
        else:
            self.logger = logger
            
    def run(self,runParameters):
        if self.cellMechanismType == 'uniformHH':
            spikeMonitor = self.uniformHH(self.cellMechanismParameters, runParameters)
        if self.cellMechanismType == 'IF':
            spikeMonitor = self.IF(self.cellMechanismParameters, runParameters)
#        print spikeMonitor.t
#        print 
        return [np.array(spikeMonitor.i),np.array(spikeMonitor.t)]

    #==========================================================================
    # HH functions
    #==========================================================================
    def uniformHH(self,cellMechanismParameters,uniformHHRunParameters):
        ## clear up the environment.
        clear(erase=True)
        k = MagicNetwork(verbose=False)
        k.reinit()
        # test whether parameters get passed
#        print cellMechanismParameters['cellModelDummySeed']
        area=1257*umetre**2
        Cm=(1.*ufarad*cm**-2)*area ## == 12.57 pF
        gl=(3e-5*siemens*cm**-2)*area
        El=-54.3*mV
        EK=-77*mV
        ENa=50*mV
        g_na=(120*msiemens*cm**-2)*area
        g_kd=(36*msiemens*cm**-2)*area
        VT=-63*mV
        # Time constants
        taue=2*ms
        taui=10*ms
        # Reversal potentials
        Ee=0*mV
        Ei=-80*mV
        we=6*nS # excitatory synaptic weight (voltage)
        wi=67*nS # inhibitory synaptic weight
        # The model
        eqs=Equations('''
dv/dt = (gl*(El-v)+ge*(Ee-v)+gi*(Ei-v)-g_na*(m*m*m)*h*(v-ENa)-g_kd*(n*n*n*n)*(v-EK))/Cm : volt
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1
dge/dt = -ge*(1./taue) : siemens
dgi/dt = -gi*(1./taui) : siemens
alpham = 0.32*(mV**-1)*(13*mV-v+VT)/(exp((13*mV-v+VT)/(4*mV))-1.)/ms : Hz
betam = 0.28*(mV**-1)*(v-VT-40*mV)/(exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
alphah = 0.128*exp((17*mV-v+VT)/(18*mV))/ms : Hz
betah = 4./(1+exp((40*mV-v+VT)/(5*mV)))/ms : Hz
alphan = 0.032*(mV**-1)*(15*mV-v+VT)/(exp((15*mV-v+VT)/(5*mV))-1.)/ms : Hz
betan = .5*exp((10*mV-v+VT)/(40*mV))/ms : Hz
        ''')
        ## set up the cells in Brian
        ## (uniform, so only one group)
        self.cellGroup = NeuronGroup(
                len(self.networkTopology.cells), 
                model = eqs,
                threshold = EmpiricalThreshold(
                        threshold=-20*mV,
                        refractory=3*ms
                ),
                implicit=True,
                freeze=True
        )            
        ## set up the connections
        ## starts empty, so this is nothing yet
        self.cellConnections = Connection(
                self.cellGroup, 
                self.cellGroup, 
                'ge', 
#                weight=we, 
                delay = 3.0 * ms
            )  
        ## fill in the connections according to the topology
        for link in self.networkTopology.links:
            self.cellConnections[link.fromCellIndex, link.toCellIndex] = (
                we * link.linkWeight
            )
        ## set the starting values to something random
        self.cellGroup.v = El+(randn(len(self.cellGroup))*5-5)*mV
        self.cellGroup.ge = (randn(len(self.cellGroup))*1.5+4)*10.*nS
        self.cellGroup.gi = (randn(len(self.cellGroup))*12+20)*10.*nS 
        ## set up a poisson noise generator and connect it
        poissonNoiseGenerator = PoissonGroup(
            len(self.networkTopology.cells),
            rates = self.cellMechanismParameters['noise'] * hertz, 
        )
        poissonConnections = Connection(
                poissonNoiseGenerator, ## from
                self.cellGroup,        ## to 
                'ge'                   ## acts on
        )
        poissonConnections.connect_one_to_one(
                poissonNoiseGenerator,
                self.cellGroup,
                weight = we * self.networkTopology.networkParameters['noiseWeight']
        )
        ## set up the spike monitor
        self.spikeMonitor = SpikeMonitor(self.cellGroup)   
        
        ## set up the custom end condition, where we end when the spike count
        ## reaches a certain number.
        ## (The simulation still terminates then t_end is reached; it's 
        ## whichever condition is reached first.)
        @network_operation
        def spikeCountStop():
            maxCount = uniformHHRunParameters['maxSpikeCountPerCell'] * len(self.networkTopology.cells)
            if self.spikeMonitor.nspikes >= maxCount:
                self.network.stop()
        ## bung everything into a network
        network = Network(
                spikeCountStop,
                self.cellGroup,
                self.cellConnections,
                self.spikeMonitor,
                poissonNoiseGenerator,
                poissonConnections
        )
        ## run the network & return the spike train
        self.logger.debug('Start actual running.')
        network.run(uniformHHRunParameters['t_end']*ms)
        return self.spikeMonitor
#        raster_plot(self.spikeMonitor)
#        show() 

    def IF(self,cellMechanismParameters,IFRunParameters):
        '''
        Simple but functional.
        '''
        self.logger.debug('starting IF run')
        ## clear up the environment.
#        clear(erase=True)
#        k.remove()
#        k = MagicNetwork()
#        k.reinit()
        ## mechanism itself
        tau = 10 * msecond
        Vt = 0 * mvolt
        Vr = -10 * mvolt
        El = 1 * mvolt
#        Ac = 0
        baseW = 10 * mvolt
        cellsPerType = dict()
        for cell in self.networkTopology.cells:
#            print cell.cellType
            try:           
                cellsPerType[cell.cellType].append(cell)
            except KeyError:
                cellsPerType[cell.cellType] = [cell]
        self.logger.debug('{0} cell objects initialized'.format(len(self.networkTopology.cells)))
        cellGroupDefs = [(len(v), k) for k,v in cellsPerType.items()]
        model = '''
            dV/dt = Ac*-(V-El)/tau : volt (unless refractory)
            Ac : 1
        '''
        
#        '''
#            dv/dt  = (ge+gi-(v-El))/taum : volt (unless refractory)
#            dge/dt = -ge/taue : volt
#            dgi/dt = -gi/taui : volt
#        '''
        cellGroups = []
        typeGroupDict = {}
        allCells = NeuronGroup(
                N=sum([a for a,b in cellGroupDefs]), 
                model = model, 
                threshold = 'V > Vt', 
                reset = 'V = Vr',
                refractory = 20. * ms
            )
        
        for cell in self.networkTopology.cells:
            if cell.cellType == 'activeExcitatory': 
                cI = cell.overallIndex
                allCells[cI:cI+1].Ac = 1. ## can only slice, hence a slice of 1
        allConns = Synapses(
                allCells, 
                allCells, 
                model = 'w : volt',
                on_pre =  'V_post += w'
        )
        self.logger.debug('Synapses set up.')
#        allConns.w = psp
        
#        allConns.w = psp
#        allConns.connect_full(
#                allCells, 
#                allCells, 
#                weight = 0*psp, 
#                delay = 6.
#        )
        self.logger.debug('Setting up {0} links (synapses)'.format(len(self.networkTopology.links)))
        fromList = []
        toList = []
        for link in self.networkTopology.links:
#            print link
            fromList.append(link.fromCell.overallIndex)
            toList.append(link.toCell.overallIndex)
#            allConns.connect(link.fromCell.overallIndex, link.toCell.overallIndex)
#        print [fromList,toList]
        allConns.connect(i=fromList,j=toList)
        allConns.w = 1.01 * baseW
#        allConns.delay = 2.*msecond
        if self.networkTopology.transmissionSpeed == None:
            allConns.delay = 2.*msecond
        else:
            for link in self.networkTopology.links:
                i = link.fromCell.overallIndex
                j = link.toCell.overallIndex
                thisDelay = (
                    self.networkTopology.cellSize *
                    link.rhosum
                ) / self.networkTopology.transmissionSpeed
                allConns.delay[i,j] = 2.*msecond + thisDelay *second
#        for link in self.networkTopology.links:
#            lid = (link.fromCell.overallIndex,link.toCell.overallIndex)
#            allConns[lid].w = baseW * link.linkWeight
        self.logger.debug('Synapses connected.')
        poissonInput = PoissonInput(
            allCells, target_var = 'V', N=1,#len(self.networkTopology.cells), 
            rate=self.cellMechanismParameters['noise']*Hz,
            weight=self.networkTopology.networkParameters['noiseWeight'] * baseW
        )
        ## set up the spike monitor
        self.spikeMonitor = SpikeMonitor(allCells)
        ## set up a single cell monitor
        self.cellMonitor = StateMonitor(allCells,'V',record=[0, 1,2,3,4,5])
        
        ## set up the custom end condition, where we end when the spike count
        ## reaches a certain number.
        ## (The simulation still terminates then t_end is reached; it's 
        ## whichever condition is reached first.)
        @network_operation
        def spikeCountStop():
            if self.spikeMonitor.num_spikes >= IFRunParameters['maxSpikeCountPerCell'] * len(self.networkTopology.cells):
                self.network.stop()
        ## bung everything into a network
        self.network = Network(
                spikeCountStop,
                self.spikeMonitor,
                self.cellMonitor,
                poissonInput,
                allCells,
                allConns
        )
#        for i in cellGroups + connGroups:
#            self.network.add(i)
        ## run the network & return the spike train
        self.logger.debug('Start actual running.')
        self.network.run(IFRunParameters['t_end']*ms)
        subplot(211)
        plot(self.spikeMonitor.t/ms, self.spikeMonitor.i, '.k')
#        raster_plot(self.spikeMonitor)
        subplot(212)
        plot(self.cellMonitor.t / ms, self.cellMonitor.V[0] / mV)
        plot(self.cellMonitor.t / ms, self.cellMonitor.V[1] / mV)
        print 'should show'
        show()
        return self.spikeMonitor

#    def IF1(self,cellMechanismParameters,IFRunParameters):
#        '''
#        Simple but functional.
#        '''
#        ## clear up the environment.
#        clear(erase=True)
#        k = MagicNetwork(verbose=False)
#        k.reinit()
#        ## mechanism itself
#        tau = 10 * msecond
#        Vt = 0 * mvolt
#        Vr = -10 * mvolt
#        El = 1 * mvolt
#        Ac = 0
#        psp = 10 * mvolt
#        cellsPerType = dict()
#        for cell in self.networkTopology.cells:
##            print cell.cellType
#            try:           
#                cellsPerType[cell.cellType].append(cell)
#            except KeyError:
#                cellsPerType[cell.cellType] = [cell]
#        cellGroupDefs = [(len(v), k) for k,v in cellsPerType.items()]
#        model = '''
#            dV/dt = Ac*-(V-El)/tau : volt
#            Ac : 1
#        '''
#        cellGroups = []
#        typeGroupDict = {}
#        allCells = NeuronGroup(
#                N=sum([a for a,b in cellGroupDefs]), 
#                model = model, 
#                threshold = Vt, 
#                reset = Vr,
#                refractory = 20. * ms
#            )
#        
#        for cell in self.networkTopology.cells:
#            if cell.cellType == 'activeExcitatory': 
#                allCells[cell.overallIndex].Ac = 1.
#        allConns = DelayConnection(
#                allCells, 
#                allCells, 
#                'V',
#                delay = 3.0*msecond,
#                maxDelay = 6.*msecond
#        )
#        allConns.connect_full(
#                allCells, 
#                allCells, 
#                weight = 0*psp, 
#                delay = 6.
#        )
#        for link in self.networkTopology.links:
##            print link
#            allConns[link.fromCell.overallIndex, link.toCell.overallIndex] = (
#                psp * link.linkWeight
#            )
#        poissonInput = PoissonInput(
#            allCells, N=1,#len(self.networkTopology.cells), 
#            rate=self.cellMechanismParameters['noise'],
#            weight=self.networkTopology.networkParameters['noiseWeight'] * psp,
#            state='V'
#        )
#        ## set up the spike monitor
#        self.spikeMonitor = SpikeMonitor(allCells)
#        ## set up a single cell monitor
#        self.cellMonitor = StateMonitor(allCells,'V',record=[0, 1])
#        
#        ## set up the custom end condition, where we end when the spike count
#        ## reaches a certain number.
#        ## (The simulation still terminates then t_end is reached; it's 
#        ## whichever condition is reached first.)
#        @network_operation
#        def spikeCountStop():
#            if self.spikeMonitor.nspikes >= IFRunParameters['maxSpikeCountPerCell'] * len(self.networkTopology.cells):
#                self.network.stop()
#        ## bung everything into a network
#        self.network = Network(
#                spikeCountStop,
#                self.spikeMonitor,
#                self.cellMonitor,
#                poissonInput,
#                allCells,
#                allConns
#        )
##        for i in cellGroups + connGroups:
##            self.network.add(i)
#        ## run the network & return the spike train
#        self.network.run(IFRunParameters['t_end']*ms)
#        subplot(211)
#        raster_plot(self.spikeMonitor)
#        subplot(212)
#        plot(self.cellMonitor.times / ms, self.cellMonitor[0] / mV)
#        plot(self.cellMonitor.times / ms, self.cellMonitor[1] / mV)
#        show()
#        return self.spikeMonitor

class NetworkRun3dVisualizer:
    def __init__(self,networkTopology,spikes = None, cellRadius = .5,
                 cylinderRadius = .01, stepSize = 2, showLinks = False, 
                 flatCamera = False, opacity = .8, output = False, 
                 targetFolder = None):
        self.name = networkTopology.fileName
        nodes = networkTopology.cells
        links = networkTopology.links
#        print len(nodes)
        zippedCoordinates = zip(*[node.coordinates for node in nodes])
        [maxX,maxY,maxZ] = [max(coord) for coord in zippedCoordinates]
        [minX,minY,minZ] = [min(coord) for coord in zippedCoordinates]
        print [ [maxX,maxY,maxZ], [minX,minY,minZ]]
        zippedCoordinates = zip(*[node.coordinates for node in nodes])
        [maxX,maxY,maxZ] = [max(coord) for coord in zippedCoordinates]
        [minX,minY,minZ] = [min(coord) for coord in zippedCoordinates]
        print [ [maxX,maxY,maxZ], [minX,minY,minZ]]
        ## nobody cares about x; we want to know the deltaY / deltaZ proportion
        deltaY = maxY - minY
        deltaZ = maxZ - minZ
        [propY,propZ] = [i/float(max([deltaY,deltaZ])) for i in [deltaY,deltaZ]]
        print [propY,propZ]
        self.stepSize = stepSize
        self.cellRadius = cellRadius
        self.scene = vp.scene
        self.scene.width = 1000 #int(1000*propZ)
        self.scene.height = 1000 #int(1000*propY)
        self.scene.autocenter = False
        self.scene.background = vp.vector(1,1,1)##vp.vector(.8,.8,.8)        
        self.scene.forward = vp.vector(-1,0.,0.)
        self.scene.center = vp.vector(0,0,deltaZ*.5)
        self.scene.range = max([deltaY,deltaZ])
        self.activeNodeColour = vp.vector(1.,0.,0.)
        self.inactiveNodeColour = vp.vector(1.,1.,1.)
        self.opacity = opacity
        self.output = output
        self.targetFolder = targetFolder
        self.nodeDict = dict(
            (i,
            vp.sphere(
                pos = (vp.vector(node.coordinates)),
                opacity = opacity,
                radius = self.cellRadius))
            for i,node
            in enumerate(nodes)
        )
        links = [(link[0].overallIndex, link[1].overallIndex) for link in links]
        self.linkdict = dict()
        if showLinks == True:
            for link in links:
                self.linkdict[link] = vp.cylinder(
                    pos = self.nodeDict[list(link)[0]].pos,
                    axis = self.nodeDict[list(link)[1]].pos - self.nodeDict[list(link)[0]].pos,
                    opacity = opacity,
                    radius = cylinderRadius)
                    
    def alignEvents(self,spikes):
        ## the format of spikes produced by brian is a list of
        ## arrays, spike times of cells.
        ## the time of the last spike in ms/binsize, rounded up
        self.cellsfiringpermoment = dict((i,[]) for i in range(self.maxtime))
        ## for every spike, put the cellindex in the appropriate time-bin
#        for cellindex,spikeTrain in spikes.items():
#            print max(spikeTrain)
#            for spiketime in spikeTrain:
#                print max()
#                self.cellsfiringpermoment[int(spiketime * 1000.)].append(cellindex)       
        [
            [
                self.cellsfiringpermoment[int(spiketime)].append(cellindex)
                for spiketime
                in spikeTrain
            ]
            for cellindex,spikeTrain
            in spikes.items()
        ]
        print self.cellsfiringpermoment

    def step(self):
        print 'curtime (ms): ' + str(self.curtime)
        self.adjustColor(self.curtime)
        self.curtime += self.stepSize
        if self.output == True:
            vp.scene.visible = False
#            vp.scene.title = 
            vp.scene.visible = True
            ## take a picture
            povexport.export(display=vp.scene,filename=("viz_shot_" + str(self.curtime).rjust(3,'0')+'.pov'), include_list = ['colors.inc'],xy_ratio=1)
            povexport.export()

    def adjustColor(self,curTime):
#        print 'adjusting colors!'
#        print self.cellsfiringpermoment[curTime]
        ## reset stuff to inactive
        [
            node.__setattr__('color',self.inactiveNodeColour) 
            for node in self.nodeDict.values()
        ]
        [
            node.__setattr__('opacity',self.opacity) 
            for node in self.nodeDict.values()
        ]
#        print [self.nodeDict[cellidx] for cellidx in self.cellsfiringpermoment[curTime]]
        try:
            firings = [self.cellsfiringpermoment[t] for t in range(curTime,curTime+self.stepSize,1)]
            activeCells = []
            for f in firings:
                activeCells.extend(f)
        except KeyError:
           activeCells = []
#        print activeCells
        [
            self.nodeDict[cellidx].__setattr__('opacity',1.) 
            for cellidx in activeCells
        ]
        [
            self.nodeDict[cellidx].__setattr__('color',self.activeNodeColour) 
            for cellidx in activeCells
        ]
#        if self.links:
#            [link.__setattr__('color',self.linkcolour) for link in self.linkdict.values()]
#            [self.linkdict[link].__setattr__('color',self.colourdict[link]) for link in self.linksfiringpermoment[time]]

    def run(self,spikeFilePath,mintime=0,maxtime=None,stepSize = 10.):
        with open(spikeFilePath,'r') as spikeFile:
            spikes = dict(json.load(spikeFile))
#        try:
        self.maxMs = max([max(spikeTrain+[0]) for cell_i,spikeTrain in spikes.items()])
#        except ValueError:
#            self.maxMs = 0
#        self.maxtime = int(max([max(i) for i in spikes])/self.binsize + 1) 
        self.maxtime = int(self.maxMs) + 1
        print self.maxtime
        self.alignEvents(spikes)
#        if self.movie == False:
#            self.linegraph(mintime,maxtime,plotbinsize)
        self.curtime = mintime
        if maxtime != None and maxtime < self.maxtime:
             maxtime = maxtime
        else:
            maxtime = self.maxtime
        print str(self.curtime) + ' < ' + str(maxtime)
        while self.curtime < maxtime:
            time.sleep(.2)
#            vp.sleep(.2)
            self.step()
#        while True:
#            vp.sleep(.2)
#            self.scene.autoscale = False
    def MOOOVIE(self):
        pass
        ## C:\Program Files\ImageMagick-6.9.1-Q16

class ModelRunner:
    def __init__(self,dataDirectoryPath,topologyCachePath,pathToCaller):
        #======================================================================
        # Output and cache directories. Create them if they don't exist.
        #======================================================================
        if os.path.exists(dataDirectoryPath):
            pass
        else:
            os.mkdir(dataDirectoryPath)
        self.dataDirectoryPath = dataDirectoryPath
        if os.path.exists(topologyCachePath):
            pass
        else:
            os.mkdir(topologyCachePath)
        self.topologyCachePath = topologyCachePath
        assert os.path.exists(pathToCaller), "The model needs to know which script is calling it. Use os.path.abspath(__file__) to generate this argument." 
        self.pathToCaller = pathToCaller.replace('\\','/')
        self.pathToModule = os.path.abspath(__file__).replace('\\','/')
   
    def setup(self,
                  networkParameters,
                  cellModelParameters,
                  useCache,
                  continueFrom = None,
                  loglevel = logging.DEBUG
        ):
        '''Create the network topologies according to given parameters.
        Pass the distilled network parameters (which cell connects to which,
        and so on) to the cell mechanism implementation.
        '''
        self.useCache = useCache
        #======================================================================
        # Set up a folder for this run.
        #======================================================================
        self.runTitle = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        if continueFrom != None:
            self.runTitle = continueFrom
        self.outputDir = self.dataDirectoryPath.rstrip('/')+'/'+self.runTitle+'/'
        self.metadataPickler = MetadataPickler(self.outputDir + 'metadata')
        ## given that the moment, down to the second, is in the file name, it
        ## is reasonable to assume that the directory does not exist yet.
        if continueFrom == None:
            os.mkdir(self.outputDir)
        ## Copy both the script and the module to the target dir for reference.
        print (self.pathToCaller,self.pathToModule)##os.path.abspath(__file__)
        callerFilename = self.pathToCaller.split('/')[-1]
        shutil.copyfile(self.pathToCaller,self.outputDir + callerFilename)
        moduleFilename = self.pathToModule.split('/')[-1]
        if self.pathToCaller != self.pathToModule:
            shutil.copyfile(self.pathToModule,self.outputDir + moduleFilename)
#        print (self.pathToCaller,self.pathToModule)
        #==============================================================================
        # Setting up the logger
        #==============================================================================
        self.logger = logging.getLogger('SkinBrian_2')
        self.logger.setLevel(loglevel)
        ## and a formatter
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ## set file handler which records according to FILE_LOG_LEVEL
        ## (also create a new log file)
        logname = time.strftime(self.outputDir+"./%Y-%m-%d %H:%M:%S", time.gmtime()) + '.log'
        f = open(logname,'wb'); f.close()
        fh = logging.FileHandler(logname)
        fh.setLevel(FILE_LOG_LEVEL)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        ## set console handler which displays according to CONSOLE_LOG_LEVEL
        ch = logging.StreamHandler()
        ch.setLevel(CONSOLE_LOG_LEVEL)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.info('Logger initialized. Setting up ModelRunner...')
        #======================================================================
        # create a single network profile for every individual possible
        # instance
        #======================================================================
        networkParameterNames = networkParameters.keys()
        NetworkProfile = collections.namedtuple(
                'NetworkProfile',networkParameterNames)
        iterableNetworkParameters = [
            name 
            for name,value 
            in networkParameters.items() 
            if type(value) == list
        ]
        uniterableNetworkParameters = [
            (k,v) for k,v in networkParameters.items() 
            if k not in iterableNetworkParameters
        ]
        self.networkProfiles = [
            dict(
                zip(iterableNetworkParameters,profileProduct) + 
                uniterableNetworkParameters
            )
            for profileProduct
            in itertools.product(*[
                networkParameters[iterableName]
                for iterableName
                in iterableNetworkParameters
            ])
        ]
        self.logger.debug('List network profiles:')
        for prf in self.networkProfiles:
            self.logger.debug(str(prf))
#        print self.networkProfiles
        self.networkProfiles = [
            NetworkProfile(**npDict)
            for npDict in self.networkProfiles
        ]
        self.logger.debug('----------------------------------------- List end')
        #======================================================================
        # create a single cell model profile for every individual possible
        # instance
        #======================================================================
        cellModelParameterNames = cellModelParameters.keys()
        CellModelProfile = collections.namedtuple(
                'CellModelProfile',cellModelParameterNames)
        iterableCellModelParameters = [
            name 
            for name,value 
            in cellModelParameters.items() 
            if type(value) == list
        ]
        uniterableCellModelParameters = [
            (k,v) for k,v in cellModelParameters.items() 
            if k not in iterableCellModelParameters
        ]
        self.cellModelProfiles = [
            dict(
                zip(iterableCellModelParameters,profileProduct) + 
                uniterableCellModelParameters
            )
            for profileProduct
            in itertools.product(*[
                cellModelParameters[iterableName]
                for iterableName
                in iterableCellModelParameters
            ])
        ]
        self.cellModelProfiles = [
            CellModelProfile(**cmpDict)
            for cmpDict in self.cellModelProfiles
        ]
        self.runProfiles = [i for i in itertools.product(
            self.cellModelProfiles,self.networkProfiles)]
        
#        print self.networkProfiles[0]._asdict()
#        oneNetwork =     
#==============================================================================
#         testing code follows.
#==============================================================================
#        print self.networkProfiles[1]
#        networkTopology = NetworkTopology(
#                self.networkProfiles[1]._asdict()['networkType'],
#                cacheDir = self.topologyCachePath
#            )
#        networkTopology.generate(self.networkProfiles[1]._asdict(),cache=self.useCache)
#        viz = NetworkRun3dVisualizer(networkTopology)
        self.logger.info('ModelRunner setup (including topology generation) done.')
            
    def run(self,runParameters):
        '''
        Iterate over all the topologies and parameters and generate output.
        '''
        runCntr=0
#        metadata = collections.OrderedDict()
        processPool = Pool(processes=CPUS)
        doneRunPaths = self.getSuccessfulMetadataFilePaths()
        doneRunIndices = [int(path[-9:][:5]) for path in doneRunPaths]
        print doneRunIndices
        metadataResultObjects = []
        done = []
        for runProfile in self.runProfiles:
            if runCntr not in doneRunIndices:
                self.logger.info('Starting run {0}. RunProfile: {1}'.format(runCntr,str(runProfile)))
                cellModelParameters = runProfile[0]._asdict()
                networkParameters = runProfile[1]._asdict()            
                if not PARALLEL:
                    thisRunMetadata = modelRunWorker(cellModelParameters,networkParameters,runParameters,runCntr,self.topologyCachePath,self.outputDir,self.useCache,self.metadataPickler)
                    self.metadataPickler.dumpMetadata(runCntr,thisRunMetadata)
                else:
                    metadataResultObjects.append(
                        (
                            runCntr,
                            processPool.apply_async(
                                modelRunWorker,
                                (   
                                    cellModelParameters,networkParameters,runParameters,runCntr,
                                    self.topologyCachePath, self.outputDir,self.useCache,self.metadataPickler
                                 )
                             )
                         )
                    )
            runCntr+=1
        if PARALLEL:
            for rC,robj in metadataResultObjects:
                done.append(robj.get())
                print done
#        self.successfullyConcludedMetadata.dumpJSON()

        ## needed: rebuid all of the metadata.
        ## Read in all the metadata pkl files, put them in an OrderedDict and 
        ## JSON-sump the orderedDict.
        self.reconstitute()
        
    def getSuccessfulMetadataFilePaths(self):
        outputFiles = os.listdir(self.outputDir)
        metadatafiles = sorted([
            fn for fn in outputFiles 
            if (fn[-4:] == '.pkl' and fn[:8] == 'metadata')
        ])
        return metadatafiles
        
    def reconstitute(self):
        metadatafiles = self.getSuccessfulMetadataFilePaths()
        successfulRunMetadata = collections.OrderedDict()
        for metadatafilePath in metadatafiles:
            with open(self.outputDir+metadatafilePath,'rb') as metadatafile:
                runIndex,runMetadata = pkl.load(metadatafile)
            successfulRunMetadata[runIndex] = runMetadata
        
        with open(self.outputDir + 'metadata.json','wb') as metadatafile:
            json.dump(successfulRunMetadata,metadatafile,indent=4)
            

class VPVisualizer:
    def __init__(self,targetFolder,topologyCachePath):
#        os.listdir(targetFolder)
        self.targetFolder = targetFolder
        self.topologyCachePath = topologyCachePath
        metadataPath = self.targetFolder.strip('/') + '/metadata.json'
        with open(metadataPath,'r') as metadataFile:
            self.metadata = json.load(metadataFile)
#        print self.metadata['0']['cellModelParameters']
            
    def vizSingle(self,index,vizType='balls',output=False):
        key = str(index)
        networkParameters = self.metadata[key]['networkParameters']
        networkType = networkParameters['networkType']
        networkTopology = NetworkTopology(
                networkType,
                cacheDir = self.topologyCachePath
        )
        networkTopology.generate(
                networkParameters,
                cache=True
        )
        if vizType == 'balls':
            single3dViz = NetworkRun3dVisualizer(networkTopology,
                     cellRadius = .5, showLinks = False, opacity = .8,
                     output = output, targetFolder = self.targetFolder)
        elif vizType == 'links':
            single3dViz = NetworkRun3dVisualizer(networkTopology,
                     cellRadius = .1, showLinks = True, opacity = .4,
                     output = output, targetFolder = self.targetFolder)
        spikeFilePath = self.targetFolder.strip('/') + '/spikefile_'+str(index)+'.json'
        single3dViz.run(spikeFilePath)
    
    def vizAxes(self):
        pass
        
        
if __name__ == '__main__':
#    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
#    def __init__(self, **entries): 
#        self.__dict__.update(entries)    
    mr = ModelRunner('./output/','./topologyCaches/',os.path.abspath(__file__))
    randomTubularParameters = {
        'acceptableDistanceFactor': .45,
        'cellDensity': 10.,
        'length': [1,2,3,4], 
        'circumference': [.2, .3, .4]
    }
    triangularGridTubularParameters = {
        'networkType':'triangularGridTubular', ## mandatoryl
        'length_circumference': [(32,8),(64,16),(128,32),(256,64)],
#        'length': 256, 
#        'circumference': 64,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,1.5,2,2.5,3,3.5,4,6,10],
        'excMikadoDensity': [.02,.05,.1,.2,.5,1.],
#        'excMikadoDistance': 0.,
#        'excMikadoDensity': .0,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [4, 14, 25, 38, 44, 51, 60, 76, 83, 95],#, 106, 120, 127, 131, 150, 154, 165, 176, 190, 195],
#        'networkDummySeed': [4,],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }
    
    triangularGridTubularMaxtestParameters = {
        'networkType':'triangularGridTubular', ## mandatoryl
        'length_circumference': [(256,64)],
#        'length': 256, 
#        'circumference': 64,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [10],
        'excMikadoDensity': [1.],
#        'excMikadoDistance': 0.,
#        'excMikadoDensity': .0,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [83, 95],#, 106, 120, 127, 131, 150, 154, 165, 176, 190, 195],
#        'networkDummySeed': [4,],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }

    triangularGridTubularMikadolessParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32),(256,64)],
#        'length': 256, 
#        'circumference': 64,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
#        'excMikadoDistance': [1.5,2,2.5,3,3.5,4,6,10],
#        'excMikadoDensity': [.05,.1,.2,.3,.4,.5,1.],
        'excMikadoDistance': 0.,
        'excMikadoDensity': .0,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': 19,
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    } 
    
    triangularGridTubular5Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(128,32)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,1.5,2.,2.5,3.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [19],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    } 
   
    triangularGridTubular3Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(128,32)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [19],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }    
    
    triangularGridTubular3x4Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32),(256,64)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [25],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }
    
    triangularGridTubular3x3x4Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [90578, 30542, 44346, 2234],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }    
    
    triangularGridTubular3x3x2Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': [.1,.5],
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
        'transmissionSpeed': 1, ## m/s
        'cellSize': 50*10**-6, ## m
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': 2234,
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }    
    
    triangularGridTubularTransmissionSpeedtestParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(128,32)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': 10.,
        'excMikadoDensity': [.5],
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
        'transmissionSpeed': [.01,.1,1.,10.], ## m/s
        'cellSize': 50*10**-6, ## m
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [17,563,2234],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }    
    
    triangularGridTubularTS3x4Parameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(256,64)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'transmissionSpeed': [.01,.1,1.,10.], ## m/s
        'cellSize': 50*10**-6, ## m
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [25],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }
    
    triangularGridTubularTSParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32),(256,64)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,2.,4.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'transmissionSpeed': [.01,.1,1.,10.,None], ## m/s
        'cellSize': 50*10**-6, ## m
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [19],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }
    

    triangularGridTubularTSSmallParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(64,16)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.],
        'excMikadoDensity': .5,
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'transmissionSpeed': [.01,.1,1.,10.], ## m/s
        'cellSize': 50*10**-6, ## m
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [4, 14, 25, 38, 44, 51, 60, 76, 83, 95],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }
    
   
    triangularGridTubularMVSParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(32,8),(64,16),(128,32),(256,64)], 
#        'circumference': 32,
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': [0.,1.5,2,2.5,3,3.5,4,6,10],
        'excMikadoDensity': [.1,.5],
        'inhibMikadoDistance': None,
        'inhibMikadoDensity': 0.,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': 0., ## negative but small
        'networkDummySeed': [4, 14, 25, 38],
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    }    
    
    triangularGridTubularSingleParameters = {
        'networkType':'triangularGridTubular', ## mandatory
        'length_circumference': [(256,64)], 
        'connectionDistance': 1.1, ## connects to all cells within ~ circ
        'excMikadoDistance': 6,
        'excMikadoDensity': .5,
        'inhibMikadoDistance': 3.0,
        'inhibMikadoDensity': 0.0,
#        'inhibitoryDensity': 0.,
        'excSynapseWeight': 1.01, ## just above threshold
        'inhibSynapseWeight': -.2, ## negative but small
        'networkDummySeed': 19,
#        'cellTypeUniform': 'genericCell',
#        'pacemakerStructure':'singleCellExc',
#        'pacemakerStructure':'ringExc',
        'pacemakerStructure':None,
        'noiseWeight': 1.01
    } 
    uniformHHParameters = {
        'cellModelType': 'uniformHH', ## mandatory
        'noise': [0.1] ## in hertz (units defined in implementation method)
    }
    IFParameters = {
        'cellModelType': 'IF', ## mandatory
#        'noise': [.01,.1,1.] ## in hertz (units defined in implementation method)
        'noise': .1
    }
    genericRunParameters = {
        'maxSpikeCountPerCell': 100,
        't_end': 10000 ## in ms (units defined in implementation method)
    }
    localtime = time.time()
    mr.setup(
        triangularGridTubularTSParameters,
        IFParameters,
        useCache=False
        ,continueFrom = 'TS_incNone'
    )
    mr.run(genericRunParameters)
    lastRunFolder = './output/'+max(os.listdir('./output/'))
    print lastRunFolder
    print 'Run duration: {0} s'.format(time.time()-localtime)
#    viz = VPVisualizer('./output/20150203_210508','./topologyCaches/')
#    viz = VPVisualizer(lastRunFolder,'./topologyCaches/')
#    viz.vizSingle(0,vizType='balls',output=False)

#    q = [i for i in spawnCellBuildingBlock(4,20,'ringExc',.2,.2)]
    
    