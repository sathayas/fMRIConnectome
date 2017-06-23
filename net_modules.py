#
# Calculates network modules based on the Louvain method
#
# The modular parcellation of the largest connected component is calculated
# and written as an image, and module stats are saved in a separate .npz file.
#

import community
import networkx as nx
import numpy as np
import nibabel as nib

def run_Louvain(fNet, fMask, fOutImg, fOutInfo):
    '''
    A wrapper function for network community detection by the Louvain method.
    Only the largest connected component is parcellated into modules.
    
    input parameters:
          fNet:     the adjacency list filename for the network
          fMask:    the filename for the mask image. Its header
                    is used to create a modular parcellation image
          fOutImg:  the filename for the output image with modular
                    parcellation
          fOutInfo: the filename with information on modules and 
                    modularity.
    returns:
          NONE
    
    output:
          This function generates files recording modular parcellation.
               fOutImg:    Modular parcellation image
               fOutInfo:   Modular parcellation information as a numpy .npz file.
                           It includes:
                              Q:      The modularity Q
                              NMods:  The number of modules
                              ModID:  Module ID
                              NNodes: The number of nodes in a module. In the same
                                      order as ModID
    '''
    
    # loading the network data
    G = nx.read_adjlist(fNet, nodetype=int)
    # just the largest subgraph
    GC = max(nx.connected_component_subgraphs(G), key=len)
    # computing the best partition
    partition = community.best_partition(GC)
    # calculating the modularity
    Q = community.modularity(partition, GC)
    # converting the partition into arrays
    VoxInd = [int(i) for i in partition.keys()]
    ModInd = np.array(list(partition.values()))+1  # the module number starts with 1
    # calculating sizes of the modules
    NMods = np.max(ModInd)
    ModID = range(1,NMods+1)
    NNodes = []
    for iMod in ModID:
        tmpNNodes = len(np.nonzero(ModInd == iMod)[0])
        NNodes.append(tmpNNodes)
    # reading in the mask image header & data
    img_mask = nib.load(fMask)
    X_mask = img_mask.get_data()    
    # organizing the output
    Xout = np.zeros_like(X_mask)
    VoxXYZ = np.unravel_index(VoxInd, X_mask.shape)
    Xout[VoxXYZ] = ModInd
    # writing out the image
    modimg = nib.Nifti1Image(Xout, img_mask.get_affine())
    nib.save(modimg, fOutImg)
    # writing out module stats
    np.savez(fOutInfo, Q=Q, NMods=NMods, ModID=ModID, NNodes=NNodes)


def extract_mod_ts(fModImg, ffMRI):
    '''
    A function to extract module time course.

    Input Parameters:
          fModImg:    The modular parcellation image. Each module is denoted by
                      a distinct voxel value.
          ffMRI:      The fMRI time course data. It has to be in the same space as
                      the module image.

    Returns:
          tsMat:      A 2D array of size T x M, where T corresponds to the number of
                      time points, and M corresponds to the number of modules.
          modInd:     A 1D array of the module number, in the same order of columns
                      as tsMat.
    
    '''
    # first, reading in the image data
    img_mod = nib.load(fModImg)
    X_mod = img_mod.get_data().astype(int)
    img_fMRI = nib.load(ffMRI)
    X_fMRI = img_fMRI.get_data()
    # get dimension info from data
    nMod = X_mod.max()
    nT = X_fMRI.shape[-1]
    # initializing the output matrix
    tsMat = np.zeros((nT, nMod))
    modInd = np.zeros(nMod)
    # loop over modules
    for iMod in range(nMod):
        # loop over time points
        modInd[iMod] = iMod+1
        for iT in range(nT):
            tmpfMRI = X_fMRI[:,:,:,iT]
            modMean = tmpfMRI[X_mod==iMod+1].mean()
            tsMat[iT, iMod] = modMean
    # returning the time series and module index
    return tsMat, modInd

    
            
    
