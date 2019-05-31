# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:39:26 2016

@author: oddman
"""
from SkinBrian import *
import SkinBrianVectorViz
reload(SkinBrianVectorViz)
import stripeAnalysis
reload(stripeAnalysis)
import os

CPUCount = 3
thisRunFolder = './output/TS/'
#==============================================================================
# Comic part (manually point to the correct folder and unzip spikefiles)
#==============================================================================
os.chdir('/home/oddman/gitRepos/skinbrainspython/')
#SkinBrianVectorViz.generateComic(828,852,3,thisRunFolder,'./topologyCaches/',25,CPUCount,autoplay=False,outputFormat='pdf')
SkinBrianVectorViz.generateComic(500,900,3,thisRunFolder,'./topologyCaches/',25,CPUCount,autoplay=False,outputFormat='png')
os.chdir('/set correct absolute dir here/')

#==============================================================================
# Stripe plot part
#==============================================================================
tiger = stripeAnalysis.StripeAnalysis(
    './output/fullRun/',
    './topologyCaches/',
    coincidenceWindow = 3.,
    angles = 37,
    getFromCache = False,
    fixedRingAmount = [16]
    ,filters = [
        [['networkParameters','excMikadoDensity'],[.5]]
    ]

#    ,filters = [
#        [['networkParameters','excMikadoDensity'],[.5]],
#        [['networkParameters','length_circumference'],[(128,32)]],
#        ,
#        [['networkParameters','excMikadoDistance'],[0.,2.,4.]]
#    ]
)
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
#    limitPic = (150,80),
    stripePlotScale = 'indiv'
    , errorbarDimension =  ['networkParameters','excMikadoDensity'] # 'perBinNorm'#
)