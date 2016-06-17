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

