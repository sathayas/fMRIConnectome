#
# fmri_warp.py
#
# A collection of functions to call FSL to warp fMRI data to the
# standard space, and reslice it to an appropriate voxel size.
#
#

import os

def apply_warp(featDir):
    '''applying the warp parameter to the entire fMRI time series
    
    This function assumes that the motion corrected, filtered 4D fMRI 
    time series data resides in the .feat directory by default (i.e., 
    filtered_func_data.nii.gz). Then the warped fMRI data will be written
    to the reg directory under the .feat directory, with the name
    func2standard.nii.gz.

    This function can handle cases where
    (1). The structural image has been normalized to the template space with
         nonlinear and linear warp.
    (2). The structural image has been centered and re-oriented but still in 
         the native space, with 6 dof linear transformation only. 
    
    Input Parameters:
          featDir:    The .feat directory containing the warp parameters

    Returns:
          NONE
    
    '''
    # the directory where warped images reside
    RegDir = os.path.join(featDir, 'reg')
    # the warped T1 image in the standard space
    fwT1 = os.path.join(RegDir, 'highres2standard.nii.gz')
    # motion corrected 4D fMRI data
    mcfMRI = os.path.join(featDir, 'filtered_func_data.nii.gz')
    # output file name
    fout = os.path.join(RegDir, 'func2standard.nii.gz')
    # warp parameters for the fMRI data to the standard space
    fwarp = os.path.join(RegDir, 'example_func2standard_warp.nii.gz')
    # check if nonlinear warping has been run. 
    if os.path.isfile(fwarp):
        # applying the linear and non-linear warp
        #         NOTE: example_func2standard_warp file contains both
        #               linear and nonlinear warp, created by the 
        #               convertwarp command during FEAT.
        com_warp = 'applywarp --ref=' + fwT1 
        com_warp += ' --in=' + mcfMRI 
        com_warp += ' --out=' + fout
        com_warp += ' --warp=' + fwarp
        res = os.system(com_warp)
    else:
        # applying the 6 dof linear trasformation only
        fmat = os.path.join(RegDir, 'example_func2standard.mat')
        com_warp = 'flirt -in ' + mcfMRI
        com_warp += ' -ref ' + fwT1
        com_warp += ' -out ' + fout
        com_warp += ' -init ' + fmat
        com_warp += ' -applyxfm'
        res = os.system(com_warp)


def reslice_fmri(ffMRI, img_dim, vox_sz):
    '''a function to reslice the warped fMRI image to desired size.
    
    This function is necessary since warped fMRI has the same voxel size
    as the structural MRI. The resliced image is written with the same file
    name as the input file name with _r suffix attached at the end.
    
    Input parameters:
          ffMRI:        The file name of the 4D fMRI data to be resliced.
          img_dim:      A 3-element vector describing the number of voxles
                        in x, y, and z directions.
          vox_sz:       A 3-element vector describing the size of each voxel
                        in mm (x, y, and z sizes).
    
    Returns:
          NONE:

    '''

    # file name business first
    WorkDir, fImg = os.path.split(os.path.abspath(ffMRI))
    tmpfname, tmpext = os.path.splitext(fImg)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    # the fake header name
    ffakehdr = os.path.join(WorkDir, tmpfname + '_r_tmp')
    # the output 4D fMRI data
    fout = os.path.join(WorkDir, tmpfname + '_r')
    # the identity matrix (a la fsl)
    DirFSL = os.environ['FSLDIR']
    feye = os.path.join(DirFSL, 'etc/flirtsch/ident.mat')

    # then putting together the command to create the fake header
    com_hdr = 'fslcreatehd '
    for iDim in img_dim:
        com_hdr += ' ' + str(iDim)
    com_hdr += ' 1'
    for iSize in vox_sz:
        com_hdr += ' ' + str(iSize)
    com_hdr += ' 1 0 0 0 16 ' + ffakehdr
    # creating the faek header
    res = os.system(com_hdr)

    # then putting together the command for flirt for reslicing
    com_flirt = 'flirt -in ' + ffMRI
    com_flirt += ' -applyxfm -init ' + feye
    com_flirt += ' -out ' + fout
    com_flirt += ' -paddingsize 0.0'
    com_flirt += ' -interp nearestneighbour'
    com_flirt += ' -ref ' + ffakehdr
    # then calling flirt
    res = os.system(com_flirt)

    # finally removing the fake header
    com_rm = 'rm ' + ffakehdr + '.nii.gz'
    res = os.system(com_rm)



def run_warp(featDir, img_dim, vox_sz):
    '''a function to call the warp function and the reslice function

    This function calls apply_warp to warp fMRI data and 
    reslice_fmri to reslice the fMRI data into the desired dimension
    and voxel sizes.

    Input parameters:
          featDir:      The feat directory containing the results of
                        spatial normalization
          img_dim:      A 3-element vector describing the number of voxles
                        in x, y, and z directions.
          vox_sz:       A 3-element vector describing the size of each voxel
                        in mm (x, y, and z sizes).

    Returns:
          NONE

    '''

    # file name business -- warped fMRI data under the reg directory
    ffMRI = os.path.join(featDir, 'reg/func2standard.nii.gz')
    # calling the warp function
    apply_warp(featDir)
    # calling the reslice function
    reslice_fmri(ffMRI, img_dim, vox_sz)

