#!/phs/software/Python27/bin/python2.7

import sys
sys.path.append('/home/shayasak/Projects/Connectome/Development/')
from optparse import OptionParser
from MyPreproc import file_mover, feat_wrapper, fmri_warp
from MyPreproc import bp_filter, fast_wrapper, mask, extract_mean, regress
from MyPreproc import motion_scrub, cross_corr

### Parse parameters, store into variables: ###
optparser = OptionParser()
(options, args) = optparser.parse_args()

imgT1 = args[0]
imgT1_brain = args[1]
imgfMRI = args[2]
BaseDir = args[3]

# verifying the input
print 'imgT1 = ' + imgT1
print 'imgT1_brain = ' + imgT1_brain
print 'imgfMRI = ' + imgfMRI
print 'BaseDir = ' + BaseDir


# copying files to the working directory
file_mover.copy_img(imgT1, imgT1_brain, imgfMRI, BaseDir)

# TR=2.0 and the first 3 volumes are deleted
feat_wrapper.make_feat_design(BaseDir, 2.6, 3) 
feat_wrapper.run_feat(BaseDir)

# fMRI time series is warped to the standard space and re-sliced
fmri_warp.run_warp(BaseDir)

# band-pass filtering
bp_filter.run_bp(BaseDir, 2.6)

# segmentation
fast_wrapper.run_fast(BaseDir)

# generating masks
mask.generate_mask(BaseDir)

# extracting the mean time series
extract_mean.run_mean(BaseDir)

# running regression
regress.run_regression(BaseDir)

# motion scrubbing
motion_scrub.remove_vols(BaseDir)

#cross correlation matrix calculation
cross_corr.run_crosscorr(BaseDir)

