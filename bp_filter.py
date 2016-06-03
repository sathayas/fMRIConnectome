#
# bp_filter.py
#
# a program to apply a bandpass filter to a 4D fmri data.
#

import os
import numpy as np
import nibabel as nib
from scipy.signal import butter, filtfilt


def bandpass_fMRI(ffMRI, LowCutF=0.009, HiCutF=0.08):
    '''
    A band-pass filtering function. A Butterworth filter is used as a 
    band-pass filter.

    Input parameters:
          ffMRI:      The 4D fMRI time series to be filtered.
          LowCutF:    The lower cutoff frequency for the band-pass filter.
                      If not specified, the default is 0.009 Hz.
          HiCutF:     The upper cutoff frequency for the band-pass filter.
                      If not specified, the default is 0.08 Hz.

    Returns:
          NONE:


    Output:
          This function generates a band-pass filtered fMRI data, and saves
          it with the suffix _bp attached to the original fMRI data file name.

    '''

    # file name business
    WorkDir, fImg = os.path.split(os.path.abspath(ffMRI))
    tmpfname, tmpext = os.path.splitext(fImg)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    fout = os.path.join(WorkDir, tmpfname + '_bp.nii.gz')

    # loading the input image
    img = nib.load(ffMRI)
    img_hdr = img.header
    TR = img_hdr.get_zooms()[-1]   # getting the TR from the image data
    X = img.get_data()

    # setting up the Butterworth filter
    ButtOrder = 5
    NyqF = 0.5/TR
    Wn = [LowCutF/NyqF, HiCutF/NyqF]
    b, a = butter(ButtOrder, Wn, btype='band')

    # applying the band-pass filter
    filtX = np.float32(filtfilt(b, a, X, axis=3))
    # writing out the filtered image
    bpimg = nib.Nifti1Image(filtX, img.get_affine())
    nib.save(bpimg, fout)


def run_bp(FeatDir, LowCutF=0.009, HiCutF=0.08):
    '''
    A wrapper function to call band-pass filtering function of fMRI time series
    data. Details on band-pass filtering can be found in bandpass_fMRI function.

    Input parameters:
          FeatDir:    The feat directory containing the fMRI data
          LowCutF:    The lower cutoff frequency for the band-pass filter.
                      If not specified, the default is 0.009 Hz.
          HiCutF:     The upper cutoff frequency for the band-pass filter.
                      If not specified, the default is 0.08 Hz.

    Returns:
          NONE:


    Output:
          See bandpass_fMRI function.

    '''
    # file name business
    RegDir = os.path.join(FeatDir, 'reg')
    fIn = os.path.join(RegDir, 'func2standard_r.nii.gz') # input image

    # calling the bandpass filtering function
    bandpass_fMRI(fIn, LowCutF, HiCutF)

