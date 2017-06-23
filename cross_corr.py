#
# cross_corr.py
#
# a collection of functions to calculate a cross correlation matrix
#

import os
import numpy as np
from scipy import stats, sparse
import nibabel as nib


def save_sparse_csr(filename,array):
    '''
    A function to save an array as a sparse matrix. It saves as 
    numpy array data.

    input parameters:
          filename:     The file name of the sparse matrix
          array:        The array to be saved as a sparse array

    returns:
          NONE

    '''

    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )



def calc_crosscorr(Y):
    '''
    A cross correlation calculation function.
    
    input parameters:
          Y:        An array of size T x V, where T is the number of
                    time points and V is the number of voxels.

    returns:
          Corr:     A correlation matrix of size V x V

    '''

    NScan = Y.shape[0]
    meanY = np.sum(Y,axis=0)/NScan
    varY = (NScan/(NScan-1)) * (np.sum(Y**2, axis=0)/NScan - meanY**2)
    tmpCov = (NScan/(NScan-1)) * (np.dot(Y.T, Y)/NScan - np.dot(np.matrix(meanY).T, np.matrix(meanY)))
    Corr = np.divide(tmpCov, np.dot(np.matrix(varY**0.5).T, np.matrix(varY**0.5)))
    return Corr



def pack_thresh(R, Th=0.3, PosOnly=0):
    '''
    A function to transform the correlation matrix to the uppder diagonal only,
    and threshold with a given parameter. The resulting file compresses very 
    compactly.

    input parameters:
          R:        A full correlation matrix of size V x V.
          Th:       A threshold to eliminate small correlation values.
                        For positive elements, R>Th
                        For negative elements, R<-Th (if PosOnly=0)
                    The default value of Th is 0.3.
          PosOnly:  A flag to indicate whether the thresholded correlation
                    matrix should only include positive values. 
                        PosOnly=0:    Both positive and negative correlations
                                      are retained.
                        PosOnly=1:    Only positive correlations are retained
                        PosOnly=-1:   Only negative correlations are retained
                                      (with R<-Th)
    returns:
          thrR_csr: The thresholded sparse correaltion matrix.
    '''
    
    # upper triangle only
    uR = np.triu(R,1)
    
    # thresholding
    if PosOnly==0:
        thuRpos = stats.threshold(uR, threshmin=Th, newval=0)
        thuRneg = stats.threshold(uR, threshmax=-Th, newval=0)
        thuR = thuRpos + thuRneg
    elif PosOnly==1:
        thuR = stats.threshold(uR, threshmin=Th, newval=0)
    elif PosOnly==-1:
        thuR = stats.threshold(uR, threshmax=-Th, newval=0)
    
    # sparse matrix
    thuR_csr = sparse.csr_matrix(thuR)
    return thuR_csr



def load_data(FeatDir):
    '''
    a function to load the data from the fMRI image data
    and to return a matrix Y of T x V, where T is tne number of scans
    and V is the number of within-mask voxels.

    input parameters:
          FeatDir:      The .feat directory containing the fMRI data that has
                        been normalized, resliced, band-pass filtered,
                        regressed, and motion-scrubbed.
       
    Returns:
          Y:            The T x V matrix of the fMRI data, where T is the
                        number of time points and V is the number of within-mask 
                        voxels.
          indMaskV:     The indices for the within-mask voxels.
    '''

    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fmask = os.path.join(RegDir, 'mask_fmri.nii.gz')
    ffmri = os.path.join(RegDir, 'func2standard_r_bp_reg_ms.nii.gz')
    # loading the image data
    img_data = nib.load(ffmri)
    X_data = img_data.get_data()
    img_mask = nib.load(fmask)
    X_mask = img_mask.get_data()
    NScan = X_data.shape[-1]
    # then creating a data matrix
    indMaskV = np.nonzero(X_mask)
    for iTime in range(NScan):
        tmpX_data = X_data[:,:,:,iTime]
        tmptrX_data = tmpX_data[indMaskV]
        if iTime==0:
            Y = np.array(tmptrX_data)
        else:
            Y = np.vstack((Y, tmptrX_data))
    # and returning the resulting matrix
    return Y, indMaskV


def mask_ts(x, ts):
    '''
    A function to mask the data matrix with a time series with zero-values.
    The data matrix is T x V, where T is the number of scans and V is the 
    number of withing-mask voxels. The time series ts, a 1D array of length T,
    may have some zeros. This function eliminates time points from the data
    matrix where ts=0. Moreover, the data matrix is multiplied non-zero values
    of ts.

    Input Parameters:
          x:        An array of size T x V, where T is the number of
                    time points and V is the number of voxels.
          ts:       A 1D array of length T. Some elements may be zero.

    Returns:
          y:        An array of size U x V, where U<T or U=T. This is
                    the data array masked by the time series ts. time 
                    points corresponding ts=0 are eliminated. The remaining
                    time points are weighted by the value of ts.
    '''
    # first, eliminate zero elements from x and ts
    tsNZ = ts[np.nonzero(ts)]
    xNZ = x[np.nonzero(ts)[0],:]
    # multiplying the data with ts by broadcasting
    y = xNZ * np.array([tsNZ]).T  
    # and returning the masked data
    return y


def run_crosscorr(FeatDir, Th=0.3, PosOnly=0, ts=[]):
    '''
    the wrapper function to calculate the cross correlaiton matrix
    and saves it. If a time series is provided as a 1D array, the 
    fMRI time series data is masked by that time series.

    input parameters:
          FeatDir:  The .feat directory containing the fMRI data that has
                    been normalized, resliced, band-pass filtered,
                    regressed, and motion-scrubbed.
                      **There is a hidden functionality associated with this
                        parameter. Instead of providing a .feat directory name,
                        the user can provide the directory where the correlation
                        information will be stored. The only catch is that this 
                        directory has to be somewhere under the .feat directory.
          Th:       A threshold to eliminate small correlation values.
                        For positive elements, R>Th
                        For negative elements, R<-Th (if PosOnly=0)
                    The default value of Th is 0.3.
          PosOnly:  A flag to indicate whether the thresholded correlation
                    matrix should only include positive values. 
                        PosOnly=0:    Both positive and negative correlations
                                      are retained.
                        PosOnly=1:    Only positive correlations are retained
                        PosOnly=-1:   Only negative correlations are retained
                                      (with R<-Th)
          ts:       A 1D array of length T. Some elements may be zero.


    returns:
          NONE

    output:
          This function creates a directory called CorrMat under the .feat directory.
          Within CorrMat directory, this function saves two files:
                    -corrmat_csr.npz     The correlation matrix, upper triangle only,
                                         thresholded, and sparse. In numpy array 
                                         format.
                    -corrmat_xyz.npz     The voxel indices for the correlation matrix.
                                         The order of the voxels is the same as the
                                         row / column of the correaltion matrix.
              **If a correlation matrix directory is provided as the FeatDir parameter
                instead of the actual .feat directory, then the correlation matrix
                information described above will be stored in that directory instead
                of the default CorrMat directory.

    '''
    # checking to see if FeatDir is actually a .feat directory
    tmpPath, tmpExt = os.path.splitext(FeatDir)
    if tmpExt != '.feat':
        CorrDir = FeatDir  # the input is the correlation info directory
        while True:  # finding the actual .feat directory, one directory up at a time
            FeatDir, tmpFName = os.path.split(FeatDir)
            tmpPath, tmpExt = os.path.splitext(FeatDir)
            if tmpExt == '.feat':
                break
    else:
        # the input is a .feat directory
        CorrDir = os.path.join(FeatDir, 'CorrMat') 
        
    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fCorrMat = os.path.join(CorrDir, 'corrmat_csr.npz')
    fCorrXYZ = os.path.join(CorrDir, 'corrmat_xyz.npz')
    # creating the output directory
    if os.path.isdir(CorrDir) == False:
        com_mkdir = 'mkdir ' + CorrDir
        res = os.system(com_mkdir)
    # loading the image data
    X, indMaskV = load_data(FeatDir)
    # calculating the cross correlation
    if len(ts)==0:
        R = calc_crosscorr(X)
    else:
        Y = mask_ts(X, ts)
        R = calc_crosscorr(Y)
    # compacting the correlation matrix
    thR = pack_thresh(R, Th, PosOnly)
    R = []
    # finally, saving the correlation matrix and voxel poitners
    save_sparse_csr(fCorrMat, thR)
    np.savez(fCorrXYZ, indMaskV=indMaskV)

    
    

