import h5py
import numpy as np
import crystallography as xtal
from Graph import GraphSet

def fuZqu(sym,roll=True):
    def fzQu(quat):
        quat = np.roll(quat,-1)
        quat = xtal.do2qu(quat)
        fZquat = sym.fzQu(quat)
        return xtal.qu2do(fZquat)
    return fzQu

def min_dif(x):
    def mindif(y):
        minimum_difference = np.min(np.sqrt((x-y)**2))
        return minimum_difference
    return mindif

def sum_min_dif(x):
    def summindif(y):
        return np.sum(list(map(min_dif(x),y)))
    return summindif

def load_dataset(filepath,symmetry=xtal.Symmetry('Cubic')):
    f = h5py.File(filename,'r')['DataContainers']
    try:
        data = f['ImageDataContainer']['CellFeatureData']
    except KeyError:
        data = f['SyntheticVolumeContainer']['CellFeatureData']
    neighbor_list = data['NeighborList']
    num_neighbors = data['NumNeighbors']
    misori_list = data['MisorientationList']
    quats = data['AvgQuats']
    diameters = data['EquivalentDiameters']
    surface_features = data['SurfaceFeatureList']
    centroids = data['Centroids']
    symmetry = symmetry

    return GraphSet(neighbor_list,num_neighbors,misori_list,quats,diameters,surface_features,centroids,symmetry)

def simultaneous_random_walk(Synthetic,Experimental,n=10):
    synthetic_nodes,synthetic_edges = Synthetic.generate_random_walk(n)
    experimental_nodes,experimental_edges = Experimental.generate_random_walk(n)

    if len(node_dif.shape == 1):
        node_dif = np.sum(list(map(mindif(synthetic_nodes),experimental_nodes)))
    else:
        node_dif = np.array(list(map(sum_min_dif(synthetic_edges),experimental_edges))

    if len(edge_dif.shape == 1):
        edge_dif = np.sum(list(map(mindif(synthetic_edges),experimental_edges)))
    else:
        edge_dif = np.array(list(map(sum_min_dif(synthetic_edges),experimental_edges)))
    
    try:
        return np.concatenate((node_dif,edge_dif),axis=0)
    except ValueError:
        return np.array([node_dif,edge_dif])

