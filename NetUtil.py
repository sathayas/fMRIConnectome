#
# NetUtil.py
#
# A collection of utility functions, so that it doesn't have to be
# re-written over and over in separate programs
#
#

import os
import numpy as np
import nibabel as nib
import networkx as nx
import sys
from scipy import sparse



# functions to read and write sparse matrix
def save_sparse_csr(filename,array):
    '''
    A function to save an array as a sparse matrix. It saves as 
    numpy array data.

    input parameters:
          filename:     The file name of the sparse matrix
          array:        A sparse array to be saved

    returns:
          NONE

    '''

    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )

def load_sparse_csr(filename):
    '''
    A function to load a sparse matrix as a sparse array.

    input parameters:
          filename:     The file name of the sparse matrix in numpy array format

    returns:
          returns a sparse array read from the data.

    '''
    loader = np.load(filename)
    return sparse.csr_matrix((  loader['data'], 
                                loader['indices'], 
                                loader['indptr']),
                         shape = loader['shape'])



def load_corrmat(fCorr):
    '''
    A function to load the correlation matrix, with
    the main diagonal being zero-ed.
    
    This is a legend function. Use load_corrmat_sparse for correlation matrix saved
    as a sparse matrix.
    '''

    # then reading in the data
    infile = np.load(fCorr)
    if 'CorrMat' in infile:
        R = infile['CorrMat']
    else:
        R = infile['NMatUND']
    NNodes = R.shape[0]
    tR = R * (1 - np.identity(NNodes)) # main diagonal = 0
    # find if there is any row with all zeros
    maxRow = np.max(tR, axis=1)
    zeroRow = np.nonzero(maxRow==0)[0]
    # delete the zero rows / columns
    rtR = np.delete(tR, zeroRow, axis=0)
    crtR = np.delete(rtR, zeroRow, axis=1)
    # then returning
    return crtR



def load_corrmat_sparse(CorrDir, ffmri):
    '''
    A function to load the sparse correlation matrix, with
    the main diagonal being zero-ed.

    input parameters:
          CorrDir:    The directory name where the correaltion matrix files
                      reside. The directory should includ the following files.
                         corrmat_csr.npz     The correlation matrix, upper triangle 
                                             only, thresholded, and sparse. In numpy 
                                             array format.
                         corrmat_xyz.npz     The voxel indices for the correlation 
                                             matrix. The order of the voxels is the 
                                             same as the row / column of the 
                                             correaltion matrix.
          ffmri:      The file name of the fMRI data the correlation matrix is 
                      based on

    returns:
          crtR:       The dense correlation matrix with the main diagonal being zeros
                      and zero rows / colums deleted.
          voxInd:     The voxel indices corresponding to the rows / columns for the
                      correlation matix.

    '''

    # first, file name business
    fCorr = os.path.join(os.path.abspath(CorrDir), 'corrmat_csr.npz')
    fCorrXYZ = os.path.join(os.path.abspath(CorrDir), 'corrmat_xyz.npz')

    # loading the correlation matrix
    thuCorr = load_sparse_csr(fCorr)
    NNodes = thuCorr.shape[0]   # number of nodes

    # making an upper diagonal matrix to a full matrix
    tR = thuCorr + thuCorr.transpose()
    tR = tR.toarray()  # converting to a dense matrix
    thuCorr = [] # clearing the original correlation matrix to make more room.

    # find if there is any row with all zeros
    maxRow = np.max(tR, axis=1)
    zeroRow = np.nonzero(maxRow==0)[0]

    # delete the zero rows / columns
    rtR = np.delete(tR, zeroRow, axis=0)
    crtR = np.delete(rtR, zeroRow, axis=1)

    # next, loading the fMRI data
    img_data = nib.load(ffmri)
    X_data = img_data.get_data()
    ImDim = X_data.shape[:3]

    # loading the voxel coordinate data
    infile = np.load(fCorrXYZ)
    indMaskV = infile['indMaskV']
    rindMaskV = np.delete(indMaskV, zeroRow, axis=0)  # deleting the voxels with all R=0
    voxInd = np.ravel_multi_index(rindMaskV, ImDim)
    # then returning
    return crtR, voxInd



def net_builder_RankTh(R, NodeInd, d):
    '''
    a function to construct the network by the rank-based thresholding
     
    input parameters:
          R:         A dense correlation matrix array.
          NodeInd:   A list of nodes in the network.
          d:         The rank threshold for the rank-based thresholding.

    returns:
          G:         The resulting graph (networkX format)
          WorkR:     The correlation matrix with zeros for the elements
                     corresponding to the edges in G. This updated correlaton
                     matrix experdite the construction of networks if
                     such networks are constructed sequentially in the
                     increasing order of d.
    '''

    # first, initialize the graph
    G = nx.Graph()
    G.add_nodes_from(NodeInd)
    NNodes = R.shape[0]
    # then add edges
    WorkR = np.array(R)  # the working copy of R
    for iRank in range(d):
        I = np.arange(NNodes)
        J = np.argmax(WorkR, axis=1)
        # R has to be non-zero
        trI = [i for i in range(NNodes) if abs(WorkR[i, J[i]])>0]
        trJ = [J[i] for i in range(NNodes) if abs(WorkR[i, J[i]])>0]
        # adding connections (for R>0)
        Elist = np.vstack((NodeInd[trI], NodeInd[trJ])).T 
        G.add_edges_from(Elist)
        WorkR[trI, trJ] = 0  # clearing the correlation matrix
    # finally returning the resultant graph
    return G, WorkR



def net_builder_HardTh(R, NodeInd, K):
    '''
    a function to construct the network by the hard-thresholding.

    input parameters:
          R:         A dense correlation matrix array.
          NodeInd:   A list of nodes in the network.
          K:         The target K, the average connections at each node
    
    returns:
          G:         The resulting graph (networkX format)
          RTh:       The R-value of the threshold achieveing the target K.

    '''
    
    # first, initialize the graph
    G = nx.Graph()
    G.add_nodes_from(NodeInd)
    NNodes = R.shape[0]
    # upper triangle of the correlation matrix only
    I,J = np.triu_indices(NNodes,1)
    VecR = abs(R[I,J])   # good for both pos and neg correlations
    # the number of elements is too big, so we truncate it
    # first, find the appropriate threshold for R
    NthR = 0
    tmpRth = 0.95
    StepTh = 0.05
    while NthR<K*NNodes/2.0:
        tmpRth -= StepTh
        #print('Threshold = %.2f' % tmpRth)
        NthR = len(np.nonzero(VecR>tmpRth)[0])
    # second, truncate the variables
    IndVecR = np.nonzero(VecR>tmpRth)
    thVecR = VecR[IndVecR]
    thI = I[IndVecR]
    thJ = J[IndVecR]
    # sort the correlation values
    zipR = zip(thVecR, thI, thJ)
    zipsR = sorted(zipR, key = lambda t: t[0], reverse=True)
    sVecR, sI, sJ = zip(*zipsR)
    # then adding edges
    m = int(np.ceil(K*NNodes/2.0))  # the number of edges
    trI = np.array(sI[:m])
    trJ = np.array(sJ[:m])
    Elist = np.vstack((NodeInd[trI], NodeInd[trJ])).T
    G.add_edges_from(Elist)
    RTh = sVecR[m-1]  # the threshold
    # finally returning the resultant graph and the threshold
    return G, RTh


def net_builder_HardThE(R, NodeInd, E):
    '''
    a function to construct the network by the hard-thresholding
    with a number of edges specified.

    input parameters:
          R:         A dense correlation matrix array.
          NodeInd:   A list of nodes in the network.
          E:         The target E, the total number of edges
    
    returns:
          G:         The resulting graph (networkX format)
          RTh:       The R-value of the threshold achieveing the target E

    '''

    # first, initialize the graph
    G = nx.Graph()
    G.add_nodes_from(NodeInd)
    NNodes = R.shape[0]
    # upper triangle of the correlation matrix only
    I,J = np.triu_indices(NNodes,1)
    VecR = abs(R[I,J])  # for both pos and neg correlations
    # the number of elements is too big, so we truncate it
    # first, find the appropriate threshold for R
    NthR = 0
    tmpRth = 0.95
    StepTh = 0.05
    while NthR<E:
        tmpRth -= StepTh
        print('Threshold = %.2f' % tmpRth)
        NthR = len(np.nonzero(VecR>tmpRth)[0])
    # second, truncate the variables
    IndVecR = np.nonzero(VecR>tmpRth)
    thVecR = VecR[IndVecR]
    thI = I[IndVecR]
    thJ = J[IndVecR]
    # sort the correlation values
    zipR = zip(thVecR, thI, thJ)
    zipsR = sorted(zipR, key = lambda t: t[0], reverse=True)
    sVecR, sI, sJ = zip(*zipsR)
    # then adding edges
    m = E  # the target number of edges
    trI = np.array(sI[:m])
    trJ = np.array(sJ[:m])
    Elist = np.vstack((NodeInd[trI], NodeInd[trJ])).T
    G.add_edges_from(Elist)
    RTh = sVecR[m-1]  # the threshold
    # finally returning the resultant graph and the threshold
    return G, RTh
    

