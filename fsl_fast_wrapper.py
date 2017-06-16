#
# fast_wrapper.py
#
# A wrapper function to call fast for segmentation in FSL.
#

import os

def run_fast(FeatDir):
    '''
    A wrapper function for FSL's FAST for tissue segmentation.

    Input parameters:
          FeatDir:      The .feat directory containing the high-res T1 image
                        normalized to the template space.
    
    Returns:
          NONE
    '''

    # file name business
    RegDir = os.path.join(FeatDir, 'reg')
    fIn = os.path.join(RegDir, 'highres2standard.nii.gz') # input image, normalized
    if os.path.isfile(fIn):
        # in case the normalized structural exists
        fOutBase = os.path.join(RegDir, 'highres2standard') # output base, normalized
    else:
        fIn = os.path.join(RegDir, 'highres.nii.gz')  # structural in native space
        fOutBase = os.path.join(RegDir, 'highres') # output base, normalized  
    # assemblying the command
    com_fast = 'fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g --nopve'
    com_fast += ' -o ' + fOutBase
    com_fast += ' ' + fIn
    # calling fast
    res = os.system(com_fast)
