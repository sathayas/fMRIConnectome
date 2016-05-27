#
# fsl_bet_wrapper.py
#
# A wrapper function to call BET from fsl for segmentation.
#

import os

def run_bet(fBrain, fThresh=0.3):
    #
    # A wapper function for FSL's BET.
    # Input Parameters:
    #     fBrain:     The file name for the high-res T1 image
    #     fThresh:    The fractional intensity threshold. Although FSL's
    #                 default is 0.5, this function uses 0.3 by default.
    #
    # 
    # file name business
    WorkDir, fHighRes = os.path.split(fBrain)
    tmpfname, tmpext = os.path.splitext(fHighRes)
    if tmpext == '.gz':
        # the extension is .nii.gz
        tmpfname, tmpext = os.path.splitext(tmpfname)
    fOutBrain = os.path.join(WorkDir, tmpfname+'_brain') # output file name
    # assemblying the command
    com_bet = 'bet ' + fBrain + ' ' + fOutBrain
    com_bet += ' -R -f ' + str(fThresh) + ' -g 0'
    # calling bet
    res = os.system(com_bet)
