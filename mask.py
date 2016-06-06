#
# mask.py
#
# A collection of functions to generate mask images from segmented
# T1 images.
#

import os
import numpy as np
import nibabel as nib


def reslice_to_fMRI(fBrain, ffMRI):
    '''
    A function to reslice an image (same voxel size as the structural MRI)
    to the voxel size of the fMRI data.

    Input parameters:
          fBrain:    The file name of the high resolution image to be 
                     resliced.
          ffMRI:     The fMRI data. Its voxel size is used to reslice
                     the structural image.

    Returns:
          fOut:      The file name of the resliced image.

    Output:
          This function produces a resliced image with the suffix _r
          attached to the original file name.
    '''

    # output file name
    WorkDir, fImg = os.path.split(os.path.abspath(fBrain))
    tmpfname, tmpext = os.path.splitext(fImg)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    fOut = os.path.join(WorkDir, tmpfname+"_r.nii.gz")
    
    # the identity matrix (a la fsl)
    DirFSL = os.environ['FSLDIR']
    feye = os.path.join(DirFSL, 'etc/flirtsch/ident.mat')

    # actual reslicing using flirt functionality
    com_flirt = 'flirt -in ' + fBrain
    com_flirt += ' -applyxfm -init ' + feye
    com_flirt += ' -out ' + fOut
    com_flirt += ' -paddingsize 0.0'
    com_flirt += ' -interp nearestneighbour'
    com_flirt += ' -ref ' + ffMRI
    res = os.system(com_flirt)

    # returning the file name of the resliced
    return fOut



def mask_parenchyma(FeatDir):
    '''
    This function creates a parenchyma mask (GM + WM).

    Input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
    
    Returns:
          NONE:
    
    Output:
          This function generates a brain parenchyma mask image with the suffix 
          _seg12_r.
    '''
    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fseg1 = os.path.join(RegDir, 'highres2standard_seg_1.nii.gz')
    fseg2 = os.path.join(RegDir, 'highres2standard_seg_2.nii.gz')
    fout = os.path.join(RegDir, 'highres2standard_seg_12.nii.gz')
    frout = os.path.join(RegDir, 'highres2standard_seg_12_r.nii.gz')
    ffmri = os.path.join(RegDir, 'func2standard_r.nii.gz')

    # then first, calling fslmaths to add seg1 and seg2 images
    com_fslmaths = 'fslmaths ' + fseg1
    com_fslmaths += ' -add ' + fseg2
    com_fslmaths += ' ' + fout
    res = os.system(com_fslmaths)

    # then reslice the resulting mask image to the fMRI space
    res = reslice_to_fMRI(fout, ffmri)



def mask_wm(FeatDir):
    '''
    The function to create a deep white matter mask.
    
    This function erodes the white matter segmentation twice, then reslice
    the resulting image to the same voxel size as the fMRI data.


    Input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
    
    Returns:
          NONE:
    
    Output:
          This function generates a white matter mask image with the suffix 
          _seg2_ee_r.
    '''
    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fseg2 = os.path.join(RegDir, 'highres2standard_seg_2.nii.gz')
    fseg2_e = os.path.join(RegDir, 'highres2standard_seg_2_e.nii.gz')
    fseg2_ee = os.path.join(RegDir, 'highres2standard_seg_2_ee.nii.gz')
    fseg2_ee_r = os.path.join(RegDir, 'highres2standard_seg_2_ee_r.nii.gz')
    ffmri = os.path.join(RegDir, 'func2standard_r.nii.gz')

    # then first, calling fslmaths to erode the white matter image
    com_ero1 = 'fslmaths ' + fseg2
    com_ero1 += ' -ero ' + fseg2_e
    res = os.system(com_ero1)
    # calling fslmaths again to erode the white matter image
    com_ero2 = 'fslmaths ' + fseg2_e
    com_ero2 += ' -ero ' + fseg2_ee
    res = os.system(com_ero2)

    # then reslice the resulting mask image to the fMRI space
    res = reslice_to_fMRI(fseg2_ee, ffmri)



def mask_csf(FeatDir):
    '''
    The function to create a CSF mask based on segmented image. The
    CSF segmentation image is resliced to the fMRI voxel size.

    Input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
    
    Returns:
          NONE:
    
    Output:
          This function generates a CSF mask image with the suffix 
          _seg0_r.
    '''
    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fseg0 = os.path.join(RegDir, 'highres2standard_seg_0.nii.gz')
    fseg0_r = os.path.join(RegDir, 'highres2standard_seg_0_r.nii.gz')
    ffmri = os.path.join(RegDir, 'func2standard_r.nii.gz')

    # reslice the csf mask image to the fMRI space
    res = reslice_to_fMRI(fseg0, ffmri)



def mask_brain(FeatDir, fMask=""):
    '''
    This function creates a brain mask based on the subject's normalized
    T1_brain image. It can also incorporate an external mask.

    Input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
          fMask:      The file name for the binary mask image. If not provided,
                      then it will be omitted.
    
    Returns:
          NONE:
    
    Output:
          It produces a mask image callsed mask_brain.nii.gz under the
          reg directory of the FeatDir.
    '''
    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    fT1 = os.path.join(RegDir, 'highres2standard.nii.gz')
    fT1_bin = os.path.join(RegDir, 'highres2standard_bin.nii.gz')
    ffmri = os.path.join(RegDir, 'func2standard_r.nii.gz')
    fIntersect = os.path.join(RegDir, 'mask_brain.nii.gz')

    # threshold the T1 image
    com_fslmaths = 'fslmaths ' + fT1
    com_fslmaths += ' -bin ' + fT1_bin
    res = os.system(com_fslmaths)

    # reslicing the binarized T1 image to the fMRI space
    fT1_bin_r = reslice_to_fMRI(fT1_bin, ffmri)

    # if an external mask is provided
    if fMask!="":
        # first, copy the mask to the reg directory
        MaskDir, fMaskImg = os.path.split(os.path.abspath(fMask))
        fMaskCopy = os.path.join(RegDir, fMaskImg)
        com_cp = 'cp ' + os.path.join(MaskDir, fMaskImg) + ' ' + RegDir
        res = os.system(com_cp)

        # then reslice the mask copy to the voxel size of fMRI data
        fMaskCopy_r = reslice_to_fMRI(fMaskCopy, ffmri)

        # multiplying the resliced external mask and the brain mask
        # first, load mask images
        img_mask = nib.load(fMaskCopy_r)
        X_mask = img_mask.get_data()
        img_brain = nib.load(fT1_bin_r)
        X_brain = img_brain.get_data()
        # then multiply masks
        X_prod = X_mask * X_brain
        # then saving the new mask image
        maskimg = nib.Nifti1Image(X_prod, img_mask.get_affine())
        nib.save(maskimg, fIntersect)

    # if an external mask is not provided
    else:
        com_cp = 'cp ' + fT1_bin_r + ' ' + fIntersect
        res = os.system(com_cp)



def mask_fmri(FeatDir):
    '''
    This function creates a mask for fMRI data. The mask includes
    voxels with non-zero time course. The mask only includes the
    brain voxels identified by the mask_brain.nii.gz.

    Input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
    
    Returns:
          NONE:
    
    Output:
          It produces a mask image callsed mask_fmri.nii.gz under the
          reg directory of the FeatDir.
    '''

    # directory and file names
    RegDir = os.path.join(FeatDir, 'reg')
    ffmri = os.path.join(RegDir, 'func2standard_r.nii.gz')
    fbrain = os.path.join(RegDir, 'mask_brain.nii.gz')
    fmask = os.path.join(RegDir, 'mask_fmri.nii.gz')
    # loading the image data
    img_fmri = nib.load(ffmri)
    X_fmri = img_fmri.get_data()
    NScan = X_fmri.shape[3]
    img_mask = nib.load(fbrain)
    X_mask = img_mask.get_data()
    # then creating a data matrix of elements within-mask elements
    # the dimension of the matrix is T x V, where T is tne number of schans
    # and V is the number of within-mask voxels
    indMaskV = np.nonzero(X_mask)
    for iTime in range(NScan):
        tmpX_fmri = X_fmri[:,:,:,iTime]
        tmptrX_fmri = tmpX_fmri[indMaskV]
        if iTime==0:
            Y = np.array(tmptrX_fmri)
        else:
            Y = np.vstack((Y, tmptrX_fmri))
    # identifying voxels whose time course is zero 
    indZero = np.nonzero(np.sum(Y, axis=0) == 0)[0]
    indZeroV = [indMaskV[i][indZero] for i in range(3)]
    # deleting zero voxels from the mask data matrix
    trX_mask = np.array(X_mask)
    trX_mask[indZeroV] = 0
    # writing out the new mask image
    maskimg = nib.Nifti1Image(trX_mask, img_fmri.get_affine())
    nib.save(maskimg, fmask)


def generate_mask(FeatDir, fMask=""):
    '''
    A wrapper to call all the functions to generate masks.

    The following functions are called to create different masks.
    For details on individual mask creation, refer to the specific
    function.

        mask_parenchyma:      Brain parenchyma mask (i.e., GM + WM mask)
        mask_wm:              Deep white matter mask
        mask_csf:             CSF mask
        mask_brain:           The brain mask based on T1-weighted brain image.
                              An external mask can also be incoorporated in 
                              addition to the T1 image.
        mask_fmri:            The mask for fMRI data. Only those voxels with
                              non-zero time course are included. Also, only the
                              voxels that are in the brain mask (generated by
                              mask_brain) are included.

    input parameters:
          FeatDir:    The .feat directory where normalization and segmentation
                      results reside.
          fMask:      The file name for the binary mask image. If not provided,
                      then it will be omitted.
    
    Returns:
          NONE:
    
    Output:
          It produces a number of mask images. The most important one is 
          mask_fmri.nii.gz, which will be useful in fMRI connectivity analysis.
    '''
    
    # calling functions
    mask_parenchyma(FeatDir)
    mask_wm(FeatDir)
    mask_csf(FeatDir)
    mask_brain(FeatDir, fMask)
    mask_fmri(FeatDir)

