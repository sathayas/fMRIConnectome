#
# extract_mean.py
#
# a collection of functions to extract mean time courses from various
# mask images.
#

import os
import numpy as np
import nibabel as nib


def extract_TS(fData, fMask):
    '''
    a function to extract the mean time series with a mask.

    Input parameters:
          fData:   file name for the 4D fMRI data
          fMask:   file name for the mask image (brain, wm, csf, etc.)

    Returns:
          An array with the mean time course from the masked fMRI data.
    '''

    # first, loading the image data
    img_data = nib.load(fData)
    X_data = img_data.get_data()
    img_mask = nib.load(fMask)
    X_mask = img_mask.get_data()
    # then masking the data with the mask image
    mX_data = np.zeros_like(X_data)
    for iTime in range(X_data.shape[-1]):
        tmpX_data = X_data[:,:,:,iTime]
        mX_data[:,:,:,iTime] = tmpX_data * X_mask
    # finally calculating the mean
    meanTS = []
    for iTime in range(X_data.shape[-1]):
        tmpData = mX_data[:,:,:,iTime]
        tmpnzData = tmpData[np.nonzero(tmpData)]
        meanTS.append(tmpnzData.mean())
    return np.array(meanTS)


def extract_global(FeatDir):
    '''
    the wrapper that extracts time series for the brain, wm, and
    csf and save into a file.

    input parameter:
          FeatDir:      The .feat directory containing the high-res T1 image
                        normalized to the template space. It is assumed that
                        the reg directory under the .feat directory includes
                        fMRI data that has been band-pass filtered with the
                        _bp suffix.
    
    Returns:
          NONE

    Output:
          The extracted time series from the brain parenchyma, deep white
          matter, and CSF will be written as an T x 3 array (T time points)
          named PhysPar.npz in the reg directory under the .feat directory.
    '''

    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    ffmri = os.path.join(RegDir, 'func2standard_r_bp.nii.gz')
    fmask_brain = os.path.join(RegDir, 'highres2standard_seg_12_r.nii.gz')
    fmask_wm = os.path.join(RegDir, 'highres2standard_seg_2_ee_r.nii.gz')
    fmask_csf = os.path.join(RegDir, 'highres2standard_seg_0_r.nii.gz')
    fTS = os.path.join(RegDir, 'PhysPar.npz')
    # calling the time series extraction function
    TS_brain = extract_TS(ffmri, fmask_brain)
    TS_wm = extract_TS(ffmri, fmask_wm)
    TS_csf = extract_TS(ffmri, fmask_csf)
    # consolidating the time series and save
    PhysPar=np.vstack((TS_brain, TS_wm, TS_csf)).T
    np.savez(fTS,PhysPar=PhysPar)

