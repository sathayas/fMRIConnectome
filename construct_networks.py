#
# A function to threshold a correlation matrix and write out the
# resulting adjacency list
#
import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import sys
sys.path.append('/study/connectome/Development/MyPreproc')
import NetUtil


def construct_RankTh(fCorr, ffMRI):
    #
    # a function to generate rank-based thresholding networks
    #
    #
    # some parameters
    Target_d = [3, 4, 5, 6, 8, 10, 15, 20, 30]
    # Output directory is relative to fCorr directory
    CorrDir, fCorrMat = os.path.split(fCorr)
    BaseDir, CorrDirName = os.path.split(CorrDir)
    OutBase = os.path.join(BaseDir, 'Adjlist')
    if not os.path.exists(OutBase):
        os.makedirs(OutBase)
    OutDir = os.path.join(OutBase, 'Network_RankTh')
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    # loading the correlation matrix
    R, NodeInd = NetUtil.load_corrmat_sparse(fCorr, ffMRI)
    # loop for generating rank-th networks
    for iTh in range(len(Target_d)):
        print "Generating a network with threshold d=" + str(Target_d[iTh])
        # generating the network
        if iTh==0:
            G, trR = NetUtil.net_builder_RankTh(R, NodeInd, Target_d[iTh])
            R = [] # releasing the memory
        else:
            # just generate the difference between the previous threshold.
            # then combine the resulting graphs
            deltaG, trR = NetUtil.net_builder_RankTh(trR, NodeInd, 
                                                Target_d[iTh]-Target_d[iTh-1])
            G = nx.compose(G, deltaG)
        # saving the network
        fNetFile = "Network_d" + str(Target_d[iTh]) + ".adjlist"
        fNet = os.path.join(OutDir,fNetFile)
        nx.write_adjlist(G, fNet)


def construct_HardTh(fCorr, ffMRI):
    #
    # a function to generate hard thresholding networks
    #
    #
    # some parameters
    Target_K = [10, 20, 30, 40, 50]
    # Output directory is relative to fCorr directory
    CorrDir, fCorrMat = os.path.split(fCorr)
    BaseDir, CorrDirName = os.path.split(CorrDir)
    OutBase = os.path.join(BaseDir, 'Adjlist')
    if not os.path.exists(OutBase):
        os.makedirs(OutBase)
    OutDir = os.path.join(OutBase, 'Network_HardTh')
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    # loading the correlation matrix
    R, NodeInd = NetUtil.load_corrmat_sparse(fCorr, ffMRI)
    # loop for generating hard-th networks
    for K in Target_K:
        print "Generating a network with threshold <k>=" + str(K)
        # generating the network
        G, RTh = NetUtil.net_builder_HardTh(R, NodeInd, K)
        # saving the network
        fNetFile = "Network_K" + str(K) + ".adjlist"
        fNet = os.path.join(OutDir,fNetFile)
        nx.write_adjlist(G, fNet)


def construct_HardThE(fCorr, ffMRI):
    #
    # a function to generate hard thresholding networks with the same number
    # of edges as rank-thresholded networks.
    #
    #
    # some parameters
    Target_d = [3, 4, 5, 6, 8, 10, 15, 20, 30]
    # Output directory is relative to fCorr directory
    CorrDir, fCorrMat = os.path.split(fCorr)
    BaseDir, CorrDirName = os.path.split(CorrDir)
    OutBase = os.path.join(BaseDir, 'Adjlist')
    if not os.path.exists(OutBase):
        os.makedirs(OutBase)
    OutDir = os.path.join(OutBase, 'Network_HardThE')
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    # directory where rank-th networks are
    RankDir = os.path.join(OutBase, 'Network_RankTh')
    # loading the correlation matrix
    R, NodeInd = NetUtil.load_corrmat_sparse(fCorr, ffMRI)
   # loop for generating hard-th networks
    for d in Target_d:
        print "Generating an equivalent hard thresholded network with d=" + str(d)
        # loading the rank thresholded network to determine the number of edges
        fdNetFile = "Network_d" + str(d) + ".adjlist"
        fdNet = os.path.join(RankDir,fdNetFile)
        tmpG = nx.read_adjlist(fdNet)
        E = len(tmpG.edges())
        # generating the network
        G, RTh = NetUtil.net_builder_HardThE(R, NodeInd, E)
        # saving the network
        fNetFile = "Network_EQd" + str(d) + ".adjlist"
        fNet = os.path.join(OutDir,fNetFile)
        nx.write_adjlist(G, fNet)


