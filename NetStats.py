#
# NetStats.py
#
# a collection of functions to calculate network statistics
#

import os
import numpy as np
import networkx as nx
import nibabel as nib


def eglob_node(G, xNode):
    '''
    A function to calculate the nodal global efficiency
    from a node.

    input parameters:
          G:      A graph in networkX format.
          xNode:  The node where the nodal global efficiency is calculated.
    
    returns:
          Eglob:  The nodal blobal efficiency at xNode.
    '''

    NNodes = len(G.nodes())
    Dx = list(nx.single_source_shortest_path_length(G, xNode).values())
    indZ = np.nonzero(np.array(Dx)==0)[0]
    nzDx = np.delete(Dx, indZ)
    if len(nzDx)>0:
        Eglob = (1.0/(NNodes-1.0)) * np.sum(1.0/nzDx)
    else:
        Eglob = 0
    # returning the nodal global efficiency
    return Eglob



def eglob_net(G):
    '''
    A function to calculate the network global efficiency
    
    input parameters:
          G:      A graph in networkX format.
    
    returns:
          Eglob:  The network wide average global efficiency
          Eglobi: Nodal global efficiency
          Nodes:  Nodes where nodal global efficency was calculated on. In the
                  same order as Eglobi.
    '''

    Nodes = G.nodes()
    if len(Nodes)>1:
        Eglobi = []
        nodecount = 1
        for iNode in Nodes:
            if (nodecount % 250)==0:
                print('Eglob:  Working on node: ' +str(nodecount))
            nodecount += 1
            tmpEglob = eglob_node(G, iNode)
            Eglobi.append(tmpEglob)
        Eglob = np.mean(Eglobi)
    else:
        Eglob = 0
        Eglobi = []
    return Eglob, Eglobi, Nodes



def calc_L(G):
    '''
    A function to calculate the average path lengths

    input parameters:
          G:      A graph in networkX format.
    
    returns:
          L:      The average path length for the largest connected component.
          GC:     The giant component size (in terms of number of nodes).
    '''
    
    subL = []
    subN = []
    subGs = list(nx.connected_component_subgraphs(G))
    indsubG = 1
    for H in subGs:
        subN.append(len(H.nodes()))
        if len(H.nodes())>1:
            print('L: Subgraph '+str(indsubG)+' with '+str(len(H.nodes()))+' nodes')
            subL.append(nx.average_shortest_path_length(H))
        else:
            subL.append(0)
        indsubG += 1
    # returning only the path length of the largest connected component
    iGC = np.argmax(subN)
    L = subL[iGC]
    GC = subN[iGC]
    return L, GC



def calc_C(G):
    '''
    A function to calculate the clustering coefficient

    input parameters:
          G:      A graph in networkX format.
    
    returns:
          C:      The average clustering coefficient for the entire network.
    '''

    Ci = list(nx.clustering(G).values())
    C = np.mean(Ci)
    return C



def calc_D(G):
    '''
    A function to calculate the diameter
    
    input parameters:
          G:      A graph in networkX format.
    
    returns:
          D:      The diameter, which is the largest diameter among all the sub
                  components.
    '''

    subD = []
    subGs = list(nx.connected_component_subgraphs(G))
    indsubG = 1
    for H in subGs:
        if len(H.nodes())>1:
            print('D: Subgraph '+str(indsubG)+' with '+str(len(H.nodes()))+' nodes')
            subD.append(nx.diameter(H))
        else:
            subD.append(0)
        indsubG += 1
    # returning the maximum diameter among all the sub components
    D = np.max(subD)
    return D


def calc_LDEglob_node(G, xNode):
    '''
    A function to calculate the nodal contribution of path length
    L, diameter D, and global efficiency Eglob from a node xNode.
    
    input parameters:
          G:      A graph in networkX format.
          xNode:  The node where the nodal global efficiency is calculated.
    
    returns:
          L:         The nodal path length at xNode.
          D:         The nodal diameter at xNode.
          EglobSum:  The nodal global efficiency.
    '''
         
    NNodes = len(G.nodes())
    Dx = list(nx.single_source_shortest_path_length(G, xNode).values())
    indZ = np.nonzero(np.array(Dx)==0)[0]
    nzDx = np.delete(Dx, indZ)
    if len(nzDx)>0:
        EglobSum = np.sum(1.0/nzDx)
        L = np.mean(nzDx)
        D = np.max(nzDx)
    else:
        EglobSum = 0
        L = 0
        D = 0
    # returning the nodal global efficiency
    return L, D, EglobSum


def calc_LDEglob_subnet(G):
    '''
    A function to calculate the path length L, diameter D, 
    and global efficiency Eglob from a node xNode within a
    connected component of a graph.

    This is a legacy code and its purpose is unclear
    '''

    Nodes = G.nodes()
    if len(Nodes)>1:
        EglobSum = []
        Li = []
        Di = []
        nodecount = 1
        for iNode in Nodes:
            if (nodecount % 250)==0:
                print('Calculating distance - Working on node: ' +str(nodecount))
            nodecount += 1
            tmpL, tmpD, tmpEglobSum = calc_LDEglob_node(G, iNode)
            Li.append(tmpL)
            Di.append(tmpD)
            EglobSum.append(tmpEglobSum)
        L = np.mean(Li)
        D = np.max(Di)
    else:
        L = 0
        D = 0
        EglobSum = [0] 
    return L, D, EglobSum, Nodes


def calc_LDEglob_net(G):
    '''
    This is a legacy function and its purpose is unclear
    (perhaps to speed up calculation by calculating L, D, and Eglob together?)
    '''
    # initializing the recorders
    subL = []
    subN = []
    subD = []
    subEglobSum = []
    subNodes = []
    # loop over connected subgraphs
    subGs = list(nx.connected_component_subgraphs(G))
    indsubG = 1
    for H in subGs:
        subN.append(len(H.nodes()))
        if len(H.nodes())>1:
            print('Subgraph '+str(indsubG)+' with '+str(len(H.nodes()))+' nodes')
        tmpL, tmpD, tmpEglobSum, tmpNodes = calc_LDEglob_subnet(H)
        subL.append(tmpL)
        subD.append(tmpD)
        subEglobSum += tmpEglobSum
        subNodes += tmpNodes
        indsubG += 1
    # organizing the output
    indGC = np.argmax(subN)
    L = subL[indGC]
    D = subD[indGC]
    GC = subN[indGC]
    NNodes = len(subNodes)
    tmpEglobi = list((1.0/(NNodes-1.0)) * np.array(subEglobSum))
    Eglob = np.mean(tmpEglobi)
    Nodes, Eglobi = sort_nodestat(subNodes, tmpEglobi)
    # returning the results
    return L, D, GC, Eglob, Eglobi, Nodes



def subgraph(G, xNode):
    ''''
    A function to extract a subgraph of a node xNode
    
    input parameters:
          G:      A graph in networkX format.
          xNode:  The node where the nodal global efficiency is calculated.

    returns:
          subG:   A subgraph of G, containing neighbors of xNode but not xNode
                  itself.
    '''
    subNodes = list(nx.all_neighbors(G, xNode))
    Edges = G.edges()
    subEdges = []       #create list of subgraph edges
    for iEdge in Edges:
        if (iEdge[0] in subNodes and iEdge[1] in subNodes):
            subEdges.append(iEdge)
    subG = nx.Graph()             # create subgraph
    subG.add_nodes_from(subNodes)    #populate subgraph with nodes
    subG.add_edges_from(subEdges)    # populate subgraph with edges
    return subG


def eloc_node(G, xNode):
    '''
    A function to calculate the nodal local efficiency
    from a node xNode.
    
    input parameters:
          G:      A graph in networkX format.
          xNode:  The node where the nodal global efficiency is calculated.

    returns:
          Eloc:   The nodal local efficiency at node xNode.
    '''

    subG = subgraph(G, xNode)
    #Eloc, tmpEloci, tmpNodes = eglob_net(subG)
    NNodes = len(subG.nodes())
    if NNodes>1:
        #Dij = nx.all_pairs_shortest_path_length(subG)
        Dij = nx.floyd_warshall(subG)
        D = [Dij[i].values() for i in subG.nodes()]
        cD = []
        for irow in D:
            cD += irow            
        NZD = np.array(cD)[np.nonzero(cD)]
        if len(NZD)>0:
            Eloc = (1.0/(NNodes*(NNodes-1.0))) * np.sum(1.0/NZD)
        else:
            Eloc = 0
    else:
        Eloc = 0
    return Eloc


def eloc_net(G):
    '''
    A function to calculate the network local efficiency
    
    input parameters:
          G:       A graph in networkX format.

    returns:
          Eloc:    The network-wide average local efficiency.
          sEloci:  The nodal local efficiency
          sNodes:  The nodes used in calculation of local efficiency. The same
                   order as the sEloci. Node numbers are sorted in the ascending
                   order.
    '''

    Nodes = G.nodes()
    Eloci = []
    nodecount = 1
    for iNode in Nodes:
        if (nodecount % 250)==0:
            print('Eloc:  Working on node: ' +str(nodecount))
        nodecount += 1
        tmpEloc = eloc_node(G, iNode)
        Eloci.append(tmpEloc)
    Eloc = np.mean(Eloci)
    # sorting the nodal local efficiency
    sNodes, sEloci = sort_nodestat(Nodes, Eloci)
    return Eloc, sEloci, sNodes



def GCSize(G):
    '''
    A function to caluclate the giant connected component size
    
    input parameters:
          G:       A graph in networkX format.

    returns:
          GC:      The giant component size (in terms of the number of nodes)
    '''

    cc = sorted(nx.connected_components(G), key = len, reverse=True)
    GC = len(cc[0])
    return GC



def degree_node(G):
    '''
    A function to calculate node degree

    input parameters:
          G:       A graph in networkX format.

    returns:
          K:       The average node degree.
          sKi:     Node degrees for individual nodes.
          sNodes:  The list of nodes, in the same order as sKi.
    '''

    Kinfo = G.degree()
    Nodes = list(Kinfo.keys())
    Ki = list(Kinfo.values())
    K = np.mean(Ki)
    # sorting the node degrees
    sNodes, sKi = sort_nodestat(Nodes, Ki)
    return K, sKi, sNodes



def sort_nodestat(NodeList, Stats):
    '''
    A function to sort node stats by ROI number
    
    This is a function used internally.
    '''
    iNode = [int(i) for i in NodeList]
    zipstat = zip(iNode, Stats)
    zipsstat = sorted(zipstat, key = lambda t: t[0])
    sNodeList, sStats = zip(*zipsstat)
    return sNodeList, sStats


def write_nodestat_nii(NodeList, Stats, fHDR, fOut):
    '''
    A function to write out node stats as an image
    
    input parameters:
          NodeList:    A list of nodes
          Stats:       A list of network stats, in the same order as NodeList
          fHDR:        An image whose header information will be used to write out
                       an image.
          fOut:        The file name for the output image

    returns:
          NONE
    
    output:
          This function generates an image with the file name specified by fOut.
    '''

    # loading the image data
    img_data = nib.load(fHDR)
    X_data = img_data.get_data()
    X_dim = X_data.shape
    # initializind the output image matrix
    X_out = np.zeros_like(X_data)
    X_out.dtype = 'float32'
    # converting node index to voxel coordinates
    NodeXYZ = np.unravel_index(NodeList, X_dim)
    X_out[NodeXYZ] = Stats
    # writing out the results
    OutImg = nib.Nifti1Image(X_out, img_data.get_affine())
    nib.save(OutImg, fOut)


