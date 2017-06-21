#
# motion_scrub.py
#
# implements motion scrubbing
#

import os
import numpy as np
import nibabel as nib


def calc_fd(FeatDir):
    '''
    a function to calculates frame displacement (fd). 

    input parameters:
          FeatDir:      The .feat directory containing the high-res T1 image
                        normalized to the template space. The motion parameters
                        from mcflirt are assumed to reside under the mc
                        directory.

    returns:
          FD:           A time series for FD.

    outputs:
          This program also writes out the FD time series as
          -FD.npz (numpy array)
          -FD.par (text file)
          Under the mc directory.
    '''

    # directory and file names
    MCDir = os.path.join(FeatDir, 'mc')
    fMoPar = os.path.join(MCDir, 'prefiltered_func_data_mcf.par')
    fFD = os.path.join(MCDir, 'FD.npz')
    fFDtxt = os.path.join(MCDir, 'FD.par')
    # reading motion parameters
    MoPar = np.genfromtxt(fMoPar,
                          delimiter="  ",
                          missing_values=["NA"])
    NScan = MoPar.shape[0]
    # calculating FD (framewise displacement)
    dMoPar = np.diff(MoPar, axis=0)
    FD = np.sum(abs(dMoPar), axis=1)
    # writing out FD in an output file, just in case we need it
    np.savez(fFD, FD=FD, NScan=NScan)
    FDfile = open(fFDtxt, "w")
    for iFrame in range(NScan-1):
        FDfile.write("%.5f\n" % FD[iFrame])
    FDfile.close()
    return FD


def scrub_motion(FeatDir, FDTh=0.5):
    '''
    the function to remove volumes based on FD. The FD is calculated by
    calling calc_fd function.

    If the GLM model also exists, then the GLM model time series is also motion
    scrubbed. 

    input parameters:
          FeatDir:      The .feat directory containing the high-res T1 image
                        normalized to the template space. The motion parameters
                        from mcflirt are assumed to reside under the mc
                        directory.
          FDTh:         The cut off value for the FD. If the FD is greater than
                        FDTh at a certain time point, one time point prior
                        and two time points immedeately following are deleted.
                        The default FDTh is 0.5.
    
    returns:
          NONE
    
    outputs:
          This function writes out the following:
                -Motion scrubbed fMRI data (with _ms suffix)
                -Motion scrubbed GLM model, if exists. The motion scrubbed
                 GLM files are written as:
                      -GLM_model_ms.npz (numpy array)
                      -GLM_model_ms.txt (text file)
                -FD time series (see calc_fd function for more details)
                -Time series masks, with deleted time points indicated by
                 zeros.
                      -TimeMask.npz (numpy array)
                      -TimeMask.par (text file)
                 Both files are under the mc directory.

    '''
   

    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    MCDir = os.path.join(FeatDir, 'mc')
    ffmri = os.path.join(RegDir, 'func2standard_r_bp_reg.nii.gz')
    fout = os.path.join(RegDir, 'func2standard_r_bp_reg_ms.nii.gz')
    fGLMnpz = os.path.join(FeatDir, 'GLM_model.npz')
    fGLMnpz_out = os.path.join(FeatDir, 'GLM_model_ms.npz')
    fGLMtxt_out = os.path.join(FeatDir, 'GLM_model_ms.txt')
    ftmpBase = os.path.join(MCDir, 'vol')
    fTimeMask = os.path.join(MCDir, 'TimeMask.npz')
    fTimeMasktxt = os.path.join(MCDir, 'TimeMask.par')

    # if the GLM model exists, it is read
    if os.path.isfile(fGLMnpz):
        infile = np.load(fGLMnpz)
        X = infile['X']

    # calculating FD
    FD = calc_fd(FeatDir)
    NScan = len(FD)+1

    # vector of images to be deleted 
    TimeMask = np.ones(NScan)
    for iFrame in range(1, NScan):
        if FD[iFrame-1]>FDTh:
            TimeMask[iFrame-1] = 0 
            TimeMask[iFrame] = 0 
            if iFrame<NScan-1:
                TimeMask[iFrame+1] = 0
            if iFrame<NScan-2:
                TimeMask[iFrame+2] = 0

    # saving the time mask, just in case
    np.savez(fTimeMask, TimeMask=TimeMask)
    tmpFD = np.append(0, FD)
    TMfile = open(fTimeMasktxt, "w")
    for iFrame in range(NScan):
        TMfile.write("%.5f     " % tmpFD[iFrame])
        TMfile.write("%d\n" % TimeMask[iFrame])
    TMfile.close()

    # motion scrubbing the GLM
    if os.path.isfile(fGLMnpz):
        X = X[np.nonzero(TimeMask==1)[0], :]
        np.savez(fGLMnpz_out, X=X)
        f = open(fGLMtxt_out, 'w')
        for iRow in X:
            strRow = [str(i) for i in iRow]
            f.write('\t'.join(strRow))
            f.write('\n')
        f.close()

    # splitting the input image
    com_split = 'fslsplit ' + ffmri + ' ' + ftmpBase + ' -t'
    res = os.system(com_split)
    # removing volumes
    for iFrame in range(NScan):
        if TimeMask[iFrame]==0:
            com_rm = 'rm -f ' + ftmpBase + '%04d' % iFrame + '.nii.gz'
            res = os.system(com_rm)
    # re-constituting the fmri time series    
    com_merge =  'fslmerge -t ' + fout + ' ' + ftmpBase + '*'
    res = os.system(com_merge)
    # deleting temporary volumes
    com_rm = 'rm -f ' + ftmpBase + '*'
    res = os.system(com_rm)


    
