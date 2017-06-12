#
# fsl_feat_model_wrapper.py
#
# a collection of functions to call FSL's feat_model to construct
# regressors for fMRI first-level model
#

import os
import subprocess
import numpy as np
import nibabel as nib


def find_TR(ffMRI):
    '''
    a function to extract TR information from a 4D fMRI data
    
    Input Parameter:
          fFMRI:    The file name for 4D fMRI data
    
    Returns:
          TR:       TR (in second)
    '''
    img_data = nib.load(ffMRI)
    img_hdr = img_data.header
    TR = img_hdr.get_zooms()[-1]
    return TR

    
def find_nVoxels(ffMRI):
    '''
    a function to calculate the total number of voxels from a 4D fMRI data

    Input Parameter:
          fFMRI:    The file name for 4D fMRI data
    
    Returns:
          nVox:     Total number of voxels
    '''
    fslInfo = subprocess.check_output('fslinfo ' + ffMRI, shell=True)
    fslDimInfo = fslInfo.decode('utf-8').split()
    dimX = int(fslDimInfo[fslDimInfo.index('dim1') + 1])
    dimY = int(fslDimInfo[fslDimInfo.index('dim2') + 1])
    dimZ = int(fslDimInfo[fslDimInfo.index('dim3') + 1])
    dimT = int(fslDimInfo[fslDimInfo.index('dim4') + 1])
    return dimX*dimY*dimZ*dimT
     

def make_feat_model_design(fT1_brain, ffMRI, EVFileList, nVolDel=0):
    '''
    this function creates a design file (design.fsf) for feat
    
    Input Parameters:
          fT1_brain:    The brain extracted T1-weighted image
          ffMRI:        The 4D fMRI image
          EVFileList:   A list of EV timing files (in 3 column format).
                        Each element is a file name, with an absolute path.
          nVolDel:      The number of first volumes to be deleted.
                        The default is 0.
                        structural will not be normalized.
    
    Returns:
          fDes:         The design file name
    '''

    # some parameters
    fHiPass = 0.009 # the cut off freq for the high pass filtering
    sHiPass = 1/fHiPass

    # file name business to figure out the location of the temporary design file
    WorkDir, fImg = os.path.split(os.path.abspath(ffMRI))
    OutDir = os.path.join(WorkDir, 'design_stats') # temporary output directory
    if not os.path.exists(OutDir):  # if the output directory doesnt exist, make one
        os.makedirs(OutDir)
    tmpfname, tmpext = os.path.splitext(fImg)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    tmpfname = tmpfname + '_stats'
    fDesFile = os.path.join(OutDir, 'design_' + tmpfname + '.fsf')

    # MNI template location 
    DirFSL = os.environ['FSLDIR'] # getting the FSL directory from the environment variable
    fMNI = os.path.join(DirFSL, 'data/standard/MNI152_T1_2mm_brain')

    # First, opening the output file
    DesFile = open(fDesFile, "w")
    # then start writing some variables
    DesFile.write("set fmri(version) 6.00\n")
    DesFile.write("set fmri(inmelodic) 0\n")
    DesFile.write("set fmri(level) 1\n")
    DesFile.write("set fmri(analysis) 2\n")
    DesFile.write("set fmri(relative_yn) 0\n")
    DesFile.write("set fmri(help_yn) 1\n")
    DesFile.write("set fmri(featwatcher_yn) 0\n")
    DesFile.write("set fmri(sscleanup_yn) 0\n")
    # the output directory
    OutDir = os.path.join(os.path.abspath(WorkDir), tmpfname)
    DesFile.write("set fmri(outputdir) \"%s\"" % OutDir)
    DesFile.write("\n")
    # the TR
    TR = find_TR(ffMRI)
    DesFile.write("set fmri(tr) %.4f" % TR)
    DesFile.write("\n")
    # the number of volumes
    nVol = subprocess.check_output('fslnvols ' + ffMRI, shell=True)
    DesFile.write("set fmri(npts) %d" % int(nVol))
    DesFile.write("\n")
    # the number of volumes to be deleted
    DesFile.write("set fmri(ndelete) %d" % nVolDel)
    DesFile.write("\n")
    # more parameters
    DesFile.write("set fmri(tagfirst) 1\n")
    DesFile.write("set fmri(multiple) 1\n")
    DesFile.write("set fmri(inputtype) 2\n")
    DesFile.write("set fmri(filtering_yn) 0\n")
    DesFile.write("set fmri(brain_thresh) 10\n")
    DesFile.write("set fmri(critical_z) 5.3\n")
    DesFile.write("set fmri(noise) 0.66\n")
    DesFile.write("set fmri(noisear) 0.34\n")
    DesFile.write("set fmri(mc) 1\n")
    DesFile.write("set fmri(sh_yn) 0\n")
    DesFile.write("set fmri(regunwarp_yn) 0\n")
    DesFile.write("set fmri(dwell) 0.7\n")
    DesFile.write("set fmri(te) 35\n")  # not sure without DICOM
    DesFile.write("set fmri(signallossthresh) 10\n")
    DesFile.write("set fmri(unwarp_dir) y-\n")
    DesFile.write("set fmri(st) 0\n")
    DesFile.write("set fmri(st_file) \"\"\n")
    DesFile.write("set fmri(bet_yn) 1\n")
    DesFile.write("set fmri(smooth) 0.0\n")
    DesFile.write("set fmri(norm_yn) 0\n")
    DesFile.write("set fmri(perfsub_yn) 0\n")
    DesFile.write("set fmri(temphp_yn) 1\n")
    DesFile.write("set fmri(templp_yn) 0\n")
    DesFile.write("set fmri(melodic_yn) 0\n")
    # carry out the main stats
    DesFile.write("set fmri(stats_yn) 1\n")
    # no prewhitening
    DesFile.write("set fmri(prewhiten_yn) 0\n")
    # no extra variables
    DesFile.write("set fmri(motionevs) 0\n")
    DesFile.write("set fmri(motionevsbeta) \"\"\n")
    DesFile.write("set fmri(scriptevsbeta) \"\"\n")
    # some FLAME parameters
    DesFile.write("set fmri(robust_yn) 0\n")
    DesFile.write("set fmri(mixed_yn) 2\n")
    # number of EVs
    nEV = len(EVFileList)
    DesFile.write("set fmri(evs_orig) %d\n" % nEV)  # main effect
    DesFile.write("set fmri(evs_real) %d\n" % (2*nEV)) # with derivatives
    DesFile.write("set fmri(evs_vox) 0\n")
    # dummy contrast
    DesFile.write("set fmri(ncon_orig) 1\n")
    DesFile.write("set fmri(ncon_real) 1\n")
    DesFile.write("set fmri(nftests_orig) 0\n")
    DesFile.write("set fmri(nftests_real) 0\n")
    DesFile.write("set fmri(constcol) 0\n")
    # some dummy values for post-stats
    DesFile.write("set fmri(poststats_yn) 0\n")
    DesFile.write("set fmri(threshmask) \"\"\n")
    DesFile.write("set fmri(thresh) 3\n")
    DesFile.write("set fmri(prob_thresh) 0.05\n")
    DesFile.write("set fmri(z_thresh) 2.3\n")
    DesFile.write("set fmri(zdisplay) 0\n")
    DesFile.write("set fmri(zmin) 2\n")
    DesFile.write("set fmri(zmax) 8\n")
    DesFile.write("set fmri(rendertype) 1\n")
    DesFile.write("set fmri(bgimage) 1\n")
    DesFile.write("set fmri(tsplot_yn) 1\n")
    # registration and normalization are turned off
    DesFile.write("set fmri(reginitial_highres_yn) 0\n")
    DesFile.write("set fmri(reginitial_highres_search) 90\n")
    DesFile.write("set fmri(reginitial_highres_dof) 3\n")
    DesFile.write("set fmri(reghighres_yn) 1\n")
    DesFile.write("set fmri(reghighres_search) 90\n")
    DesFile.write("set fmri(reghighres_dof) 6\n")
    DesFile.write("set fmri(regstandard_yn) 0\n")
    DesFile.write("set fmri(alternateReference_yn) 0\n")
    # fsl's MNI template image
    DesFile.write("set fmri(regstandard) \"%s\"\n" % fMNI)
    # and more parameters ...
    DesFile.write("set fmri(regstandard_search) 90\n")
    DesFile.write("set fmri(regstandard_dof) 12\n")
    DesFile.write("set fmri(regstandard_nonlinear_yn) 1\n")
    DesFile.write("set fmri(regstandard_nonlinear_warpres) 10\n")
    # highpass filter cutoff
    DesFile.write("set fmri(paradigm_hp) %.2f\n" % sHiPass)
    # number of voxels
    nVox = find_nVoxels(ffMRI)
    DesFile.write("set fmri(totalVoxels) %d\n" % nVox)
    # and more parameters ...
    DesFile.write("set fmri(ncopeinputs) 0\n")
    # fmri time series
    DesFile.write("set feat_files(1) \"%s\"\n" % os.path.abspath(ffMRI))
    # and more parameters ...
    DesFile.write("set fmri(confoundevs) 0\n")
    # T1 brain extracted image
    DesFile.write("set highres_files(1) \"%s\"\n" % os.path.abspath(fT1_brain))
    # Writing out parameters for EVs
    for iEV in range(nEV):
        indEV = str(iEV + 1)
        nameEV = 'EV' + indEV
        DesFile.write("set fmri(evtitle" + indEV + ") \"" + nameEV + "\"\n")
        DesFile.write("set fmri(shape" + indEV + ") 3\n")
        DesFile.write("set fmri(convolve" + indEV + ") 3\n")
        DesFile.write("set fmri(convolve_phase" + indEV + ") 0\n")
        DesFile.write("set fmri(tempfilt_yn" + indEV + ") 1\n")
        DesFile.write("set fmri(deriv_yn" + indEV + ") 1\n")
        DesFile.write("set fmri(custom" + indEV + ") \"" + EVFileList[iEV]
                      + "\"\n")
        for jEV in range(nEV+1):
            DesFile.write("set fmri(ortho" + indEV + "." + str(jEV) + ") 0\n")
    # Parameters for a dummy contrast
    DesFile.write("set fmri(con_mode_old) orig\n")
    DesFile.write("set fmri(con_mode) orig\n")
    DesFile.write("set fmri(conpic_real.1) 0\n")
    DesFile.write("set fmri(conname_real.1) \"\"\n")
    # dummy contrast vector (real)
    for iCol in range(2*nEV):
        DesFile.write("set fmri(con_real1." + str(iCol+1) + ") 1\n")
    # more parameters for the dummy contrast
    DesFile.write("set fmri(conpic_orig.1) 0\n")
    DesFile.write("set fmri(conname_orig.1) \"\"\n")
    # dummy contrast vector (orig)
    for iCol in range(nEV):
        DesFile.write("set fmri(con_orig1." + str(iCol+1) + ") 1\n")
    # Other contrast related parameters
    DesFile.write("set fmri(conmask_zerothresh_yn) 0\n")
    DesFile.write("set fmri(conmask1_1) 0\n")
    # options not appearing on the GUI
    DesFile.write("set fmri(alternative_mask) \"\"\n")
    DesFile.write("set fmri(init_initial_highres) \"\"\n")
    DesFile.write("set fmri(init_highres) \"\"\n")
    DesFile.write("set fmri(init_standard) \"\"\n")
    DesFile.write("set fmri(overwrite_yn) 0\n")
    # finally closing the file
    DesFile.close()
    return fDesFile

##### START FROM HERE #####
def run_feat_model(fT1_brain, ffMRI, nVolDel=0, bNorm=True):
    '''
    The wrapper function to run feat to normalized T1 to MNI
 
   Input Parameters:
          fT1_brain:    The brain extracted T1-weighted image
          ffMRI:        The 4D fMRI image
          nVolDel:      The number of first volumes to be deleted.
                        The default is 0.
          bNorm:        The flag for normalization to the standard space. It can take
                        True or False. If True (default), then the structural will be
                        normalized to the standard template. If False, then the 
                        structural will not be normalized.
    
    Returns:
          DirFeat:      The output .feat directory name

    '''

    fDesFile = make_feat_design(fT1_brain, ffMRI, nVolDel, bNorm)
    com_feat = 'feat ' + fDesFile
    res = os.system(com_feat)
    
    # output directory name
    WorkDir, fImg = os.path.split(os.path.abspath(ffMRI))
    tmpfname, tmpext = os.path.splitext(fImg)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    DirFeat = os.path.join(os.path.abspath(WorkDir), tmpfname+'.feat')
    return DirFeat

