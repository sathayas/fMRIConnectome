#
# ROINetUtil.py
#
# A collection of utility functions to construct and analyze ROI-based
# networks.
#
#

import os
import numpy as np
import nibabel as nib
import networkx as nx


def extract_ts(ffMRI, fAtlas, roiMin=-1, roiMax=-1):
    '''
    A function to extract the average ROI time series from 4D fMRI
    data.

    input parameters:
          ffMRI:    The file name of the 4D fMRI image. The fMRI data should
                    have already been normalized and preprocessed.
          fAtlas:   The file name of the atlas image defining different ROIs.
                    The atlas image is assumed to be in the same space as the 
                    fMRI data. In other words, it needs to be re-sliced to 
                    the fMRI data voxel size beforehand.
          roiMin:   The starting number for the ROI index. The default is -1, 
                    and it uses the smallest ROI value available in the atlas
                    image.
          roiMax:   The ending number for the ROI index. The default is -1, 
                    and it uses the largest ROI value available in the atlas
                    image.


    returns:
          roi_ts:   An array of the extracted time series. Rows correspond to
                    time points, the columns corresponds to ROIs. The ROIs are
                    in the same order as roi_ind.
          roi_ind:  A vector of ROI numbers, in the same order as the columns
                    of the roi_ts.

    '''

    # reading in the fMRI data
    fMRI = nib.load(ffMRI)
    datafMRI = fMRI.get_data()

    # reading in the brain atlas data
    atlas = nib.load(fAtlas)
    dataAtlas = np.flipud(np.squeeze(atlas.get_data()))

    # creating the roi indices
    if roiMin==-1:
        roiMin = np.min(np.unique(dataAtlas[dataAtlas>0]))
    if roiMax==-1:
        roiMax = np.max(np.unique(dataAtlas[dataAtlas>0]))
    roi_ind = np.arange(roiMin, roiMax+1)
    
    # preparing the output time series array
    nTime = datafMRI.shape[-1]
    roi_ts = np.zeros([nTime, len(roi_ind)])

    # for loop to calculate the time series
    for iTime in range(nTime):
        tmpfMRI = datafMRI[:,:,:,iTime]
        for i,iROI in enumerate(roi_ind):
            roi_ts[iTime, i] = np.nanmean(tmpfMRI[dataAtlas==iROI])

    # returning the results
    return roi_ts, roi_ind

