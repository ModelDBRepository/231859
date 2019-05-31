# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 02:45:56 2016

@author: Oddman
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

def shiftMat(sourceMat,angles=None,savgol=False):
    '''This function takes a 2d array / matrix as its argument ("SourceMat")
    and produces a list of matrices where the rows (shape[0]) are shifted 
    diagonally. The items in the list are pairs of the square matrix turned
    into parallellograms with the triangles at the end cut off, like so:
    ______________               ______________      ______________
    |             |        (    /!          !  /     \  !          !\    )
    |             |     -> (   / !          ! /       \ !          ! \   )
    |_____________|        (  /__!__________!/   ,     \!__________!__\  )
    
    The index in the list of outputs is a function of the various possible 
    angles of shifting. The optional argument "angles" specifies the number of
    angle-subdivisions from vertical to 45°. When None, that number is equal 
    to the number of rows in the source matrix.
    '''
    ## Identify the number of rows in the source matrix.
    rowNum = sourceMat.shape[0]
    outList = []
    ## Assert that the source matrix is at least four times as long as it is 
    ## wide - owing to the total loss (2 * number of rows), one really does not
    ## want less than that.
    ## (When applying to stripe plots, remember that rings is on the number of
    ## rows and that time is on columns, so this should not be a problem.
    ## If it is, run longer simulations.)
    assert rowNum * 4 < sourceMat.shape[1], 'Source matrix too tall'
    if angles == None:
        ## angles not specified? Use the number of rows.
        angles = rowNum
    else:
        ## Divisor must be an integer. (We need iteration steps.)
        assert type(angles) == int, '"angles" must be an integer.'
    ## One more angle; the 0th angle is the base array, only truncated.
    angleIndices = range(angles+1)
    ## Define the target shape. 
    targetShape = list(sourceMat.shape)
    targetShape[1] = targetShape[1] - (2 * (rowNum - 1))
    print 'target shape: {0}'.format(targetShape)
    ## Define what the expanded array will look like.
    expandedTargetShape = targetShape[:]
    expandedTargetShape[1] *= angles
    ## To do the slanting, we chop the cells  of the source array into shares;
    ## every cell into the number of angles. This basically expands the array,
    ## making its rows angles times as long.
    ##
    ## Create the expanded array,  which will be reused for the various
    ## angle-options:
    sourceMatSplit = np.split(sourceMat,sourceMat.shape[1],axis=1)
    concat = []
    [concat.extend([mat/float(angles) for i in range(angles)]) for mat in sourceMatSplit]
    expandedMat = np.concatenate(concat,axis=1)
    if savgol == True:
        savgolWindowLength = angles*2
        savgolWindowLength += 1
        print 'herpderp'
        expandedMat = signal.savgol_filter(expandedMat,window_length = savgolWindowLength, polyorder = 3, axis = 1)
        plt.imshow(expandedMat,interpolation='nearest',cmap=plt.cm.gray)
    for angleIndex in angleIndices:
        ## For every possible angle, plus the case of no shift at all, do:
        ## create two matrices; left and right-slanting
        rightSlantExpanded = np.zeros(expandedTargetShape,dtype=expandedMat.dtype)
        leftSlantExpanded = np.zeros(expandedTargetShape,dtype=expandedMat.dtype)
        ## Then, for every row in the matrix, calculate the to-and-from
        ## indices for the expanded matrix to result in the correctly slanted
        ## output:
        for rowIndex in range(rowNum):
            shiftRightAmount = angleIndex*rowIndex
            shiftLeftAmount = (angles*(rowNum-1))-(angleIndex*rowIndex)
            ## defensive check to see whether the size is correct (numpy won't 
            ## give an error if the size is incorrect, it will 'broadcast')            
            assert np.shape(expandedMat[
                    rowIndex,
                    shiftRightAmount
                    :
                    expandedTargetShape[1]+shiftRightAmount
            ]) == np.shape(rightSlantExpanded[rowIndex,:])
            rightSlantExpanded[rowIndex,:] = expandedMat[
                    rowIndex,
                    shiftRightAmount
                    :
                    expandedTargetShape[1]+shiftRightAmount
            ]
            ## Same defensive business for the right-hand case
            assert np.shape(expandedMat[
                    rowIndex,
                    shiftLeftAmount
                    :
                    expandedTargetShape[1]+shiftLeftAmount
            ]) == np.shape(leftSlantExpanded[rowIndex,:])
            leftSlantExpanded[rowIndex,:] = expandedMat[
                    rowIndex,
                    shiftLeftAmount
                    :
                    expandedTargetShape[1]+shiftLeftAmount
            ]
        ## The cutting has been done. Now, to condense the matrices.
        leftSlant = np.zeros(shape=targetShape,dtype=expandedMat.dtype)
        rightSlant = np.zeros(shape=targetShape,dtype=expandedMat.dtype)
        for colNum in range(expandedTargetShape[1]):
            ## Add the expanded rows to the condensed matrix one by one
            targetColNum = colNum // angles
            leftSlant[:,targetColNum] += leftSlantExpanded[:,colNum]
            rightSlant[:,targetColNum] += rightSlantExpanded[:,colNum]
        ## For every angle, put in a right-slanted and left-slanted matrix.
        outList.append([rightSlant,leftSlant])
    return outList
    
def angleSignificance(listOfAnglePairs,product=False):
    angles = len(listOfAnglePairs)
    ## some defensive shape checking
    ## (should be no problem if shiftMat output is used)
    baseShape = listOfAnglePairs[0][0].shape
    base_dtype = listOfAnglePairs[0][0].dtype
    for pair in listOfAnglePairs:
        for mat in pair:
            assert mat.shape == baseShape 
    allSamples = []
    angleMatShape = np.copy(baseShape)
    angleMatShape[1] *= 2
    totalVec = np.zeros(shape = angleMatShape[1]*angles, dtype = base_dtype)
    for i,pair in enumerate(listOfAnglePairs):
        totalPair = np.concatenate(pair,axis=1)
        multiMat = np.ones(totalPair.shape,dtype=np.float64)
        multiMat = multiMat * np.mean(totalPair)
        normed = totalPair / multiMat
        normed -= np.mean(normed)
        normed *= normed
#        angleVec = np.zeros(shape = angleVecLen, dtype = base_dtype)
        ## For ease of computation,  stick the opposed angles together in a
        ## single vector.
#        angleVec[0:baseShape[1]] = np.sum(pair[0],axis = 0)
#        angleVec[baseShape[1]:] = np.sum(pair[1],axis = 0)
        if product:
            angleVec = normed + .1
            angleVec = np.product(angleVec)
        else:
            angleVec = np.sum(normed,axis=0)
        totalVec[i*angleMatShape[1]:i*angleMatShape[1]+angleMatShape[1]] = angleVec
        allSamples.append(angleVec)
    overallMean = np.mean(totalVec)
#    print totalVec
#    print overallMean
    BesselCorrectedNumerator = angleMatShape[1] - 1
    sdPerAngle = []
    for sample in allSamples:
        deviance = (sample - np.mean(sample)) ** 2
#        deviance = sample
        sd = np.sqrt(np.sum(deviance) / float(BesselCorrectedNumerator))
        sdPerAngle.append(sd)
#        sdPerAngle.append(np.sum(deviance))
#        thisAngleSummedPerFrame = [np.sum(mat,axis = 0) for mat in pair]
#        mean = 
    return sdPerAngle
    
#def angleSquaredErrorPerRow(listOfAnglePairs):
#    angles = len(listOfAnglePairs)
#    ## some defensive shape checking
#    ## (should be no problem if shiftMat output is used)
#    baseShape = listOfAnglePairs[0][0].shape
#    base_dtype = listOfAnglePairs[0][0].dtype
#    for pair in listOfAnglePairs:
#        for mat in pair:
#            assert mat.shape == baseShape 
#    ## Set up target data containers
#    allSamples = []
#    angleMatShape = np.copy(baseShape)
#    angleMatShape[1] *= 2
#    totalVec = np.zeros(shape = angleMatShape[1]*angles, dtype = base_dtype)
#    ## Fill the data containers
#    for i,pair in enumerate(listOfAnglePairs):
#        
#        totalPair = np.concatenate(pair,axis=1)
#        multiMat = np.ones(totalPair.shape,dtype=np.float64)
#        multiMat = multiMat * np.mean(totalPair)
#        normed = totalPair / multiMat
##        prodded = np.product(normed,axis=0)        
#        
#        
##        print i
##        angleVec = np.zeros(shape = angleVecLen, dtype = base_dtype)
#        ## For ease of computation,  stick the opposed angles together in a
#        ## single vector.
##        angleVec[0:baseShape[1]] = np.sum(pair[0],axis = 0)
##        angleVec[baseShape[1]:] = np.sum(pair[1],axis = 0)
#        totalVec[i*angleMatShape[1]:i*angleMatShape[1]+angleMatShape[1]] = normed
#        allSamples.append(normed)
#    overallMean = np.mean(totalVec)
##    print totalVec
##    print overallMean
#    BesselCorrectedNumerator = angleMatShape[1] - 1
#    sdPerAngle = []
#    for sample in allSamples:
#        deviance = (sample - np.mean(sample)) ** 2
#        sdPerAngle.append(np.sqrt(np.sum(deviance) / float(BesselCorrectedNumerator)))
#        
##        thisAngleSummedPerFrame = [np.sum(mat,axis = 0) for mat in pair]
##        mean = 
#    return sdPerAngle
#    
#def angleSignificanceProduct(listOfAnglePairs):
#    angles = len(listOfAnglePairs)
#    ## some defensive shape checking
#    ## (should be no problem if shiftMat output is used)
#    baseShape = listOfAnglePairs[0][0].shape
#    base_dtype = listOfAnglePairs[0][0].dtype
#    for pair in listOfAnglePairs:
#        for mat in pair:
#            assert mat.shape == baseShape 
#    allSamples = []
#    angleMatShape = np.copy(baseShape)
#    angleMatShape[1] *= 2
#    totalVec = np.zeros(shape = angleMatShape[1]*angles, dtype = base_dtype)
#    for i,pair in enumerate(listOfAnglePairs):
##        print pair
#        totalPair = np.concatenate(pair,axis=1)
##        print totalPair
##        print np.max(totalPair)
#        multiMat = np.ones(totalPair.shape,dtype=np.float64)
#        multiMat = multiMat * np.mean(totalPair)
##        print multiMat
#        normed = totalPair / multiMat
#        normed += .1
#        prodded = np.product(normed,axis=0)
#        
##        print prodded
#        totalVec[i*angleMatShape[1]:i*angleMatShape[1]+angleMatShape[1]] = prodded
#        allSamples.append(prodded)

        
#    for i,pair in enumerate(listOfAnglePairs):
##        print i
#        angleVec = np.zeros(shape = angleVecLen, dtype = base_dtype)
#        ## For ease of computation,  stick the opposed angles together in a
#        ## single vector.
#        angleVec[0:baseShape[1]] = np.sum(pair[0],axis = 0)
#        angleVec[baseShape[1]:] = np.sum(pair[1],axis = 0)
#        totalVec[i*angleVecLen:i*angleVecLen+angleVecLen] = angleVec
#        allSamples.append(angleVec)
    overallMean = np.mean(totalVec)
#    print totalVec
#    print overallMean
    BesselCorrectedNumerator = angleMatShape[1] - 1
    sdPerAngle = []
    for sample in allSamples:
        deviance = (sample - np.mean(sample)) ** 2
        sdPerAngle.append(np.sqrt(np.sum(deviance) / float(BesselCorrectedNumerator)))
        
#        thisAngleSummedPerFrame = [np.sum(mat,axis = 0) for mat in pair]
#        mean = 
    return sdPerAngle    
    
if __name__ == '__main__':
#==============================================================================
#     Only runs when called as script. This is essentially unit testing code.
#==============================================================================
    ## Create a sample array.
    testMat1 = np.array([[0,1,2,1],[1,2,1,0],[2,1,0,1],[1,0,1,2],[0,1,2,1],[1,2,1,0],[2,1,0,1],[1,0,1,2],[0,1,2,1],[1,2,1,0],[2,1,0,1],[1,0,1,2],[0,1,2,1],[1,2,1,0],[2,1,0,1],[1,0,1,2],[0,1,2,1],[1,2,1,0]],dtype=np.float64)
    testMat = np.array([[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],dtype=np.float64)
    testMat = np.transpose(testMat)
#==============================================================================
#     Case 1: no Savitzky-Golay filtering.
#==============================================================================
    ## Apply the shifting.
    shiftedMats = shiftMat(testMat,5)
    mats = [[testMat,testMat]] + shiftedMats
    matNo = len(mats)
    ## Display the result.
    stripeFig = plt.figure(figsize = (8,8))
    for i,mat in enumerate(mats):
        mat = mats[i]
#        ax1 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+1))
#        ax2 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+2))
        ax1 = stripeFig.add_subplot(matNo,2,i*2+1)
        ax2 = stripeFig.add_subplot(matNo,2,i*2+2)
        ax1.imshow(mat[0],interpolation='nearest',cmap=plt.cm.gray)
        ax2.imshow(mat[1],interpolation='nearest',cmap=plt.cm.gray)
    stripeFig.show()     
    ## Apply the variance analysis.
    print angleSignificance(shiftedMats,product=False)
#    shiftMat(testMat,savgol=True)
#==============================================================================
#     Case 2: Savitzky-Golay filtering on expanded matrix.
#==============================================================================
    ## Apply the shifting.
    shiftedMats = shiftMat(testMat,5,savgol=True)
    mats = [[testMat,testMat]] + shiftedMats
    matNo = len(mats)
    ## Display the result.
    stripeFig2 = plt.figure(figsize = (8,8))
    for i,mat in enumerate(mats):
        mat = mats[i]
#        ax1 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+1))
#        ax2 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+2))
        ax1 = stripeFig2.add_subplot(matNo,2,i*2+1)
        ax2 = stripeFig2.add_subplot(matNo,2,i*2+2)
        ax1.imshow(mat[0],interpolation='nearest',cmap=plt.cm.gray)
        ax2.imshow(mat[1],interpolation='nearest',cmap=plt.cm.gray)
    stripeFig2.show()     
#==============================================================================
#     Case 2: Savitzky-Golay filtering on input.
#==============================================================================
    ## Apply the shifting.
    testMat = signal.savgol_filter(testMat,window_length = 5, polyorder = 3, axis = 1)
    shiftedMats = shiftMat(testMat,5,savgol=False)
    mats = [[testMat,testMat]] + shiftedMats
    matNo = len(mats)
    ## Display the result.
    stripeFig3 = plt.figure(figsize = (8,8))
    for i,mat in enumerate(mats):
        mat = mats[i]
#        ax1 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+1))
#        ax2 = stripeFig.add_subplot('{0}2{1}'.format(matNo,i*2+2))
        ax1 = stripeFig3.add_subplot(matNo,2,i*2+1)
        ax2 = stripeFig3.add_subplot(matNo,2,i*2+2)
        ax1.imshow(mat[0],interpolation='nearest',cmap=plt.cm.gray)
        ax2.imshow(mat[1],interpolation='nearest',cmap=plt.cm.gray)
    stripeFig3.show()     
    
