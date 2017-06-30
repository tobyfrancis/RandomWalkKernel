import h5py
import numpy as np


class GraphSet:
	def __init__(self,neighbor_list,num_neighbors,misori_list,quats,diameters,surface_features,centroids,symmetry):
		self.neighbor_list = neighbor_list
		self.num_neighbors = num_neighbors
		self.misori_list = misori_list
		self.quats = quats
		self.diameters = diameters
		self.surface_features = surface_features
		self.centroids = centroids
		self.cluster_ids = clusterids
		self.symmetry = symmetry

		print('Computing Sigma3 Boundaries...')
		self.sigma_3_list = list(map(self.sigma_3_stage(),range(len(self.neighbor_list))))
		print('Computing Neighbor List Statistics...')
		self.node_indices = np.cumsum(num_neighbors+1)
		self.node_ids_sorted = np.sort(self.neighbor_list[self.node_indices])
		self.node_indices_sorted = np.argsort(self.neighbor_list[self.node_indices])
	
	def within_tolerance(self,a,b,tolerance):
		try:
			if not (a.shape == b.shape):
				raise Exception('Tolerance function only takes arrays of same shape')
			elif len(a.shape)>2:
				raise Exception('Tolerance function implemented only for up to 1D arrays')
		except AttributeError:
				pass
		try:
			return all(np.abs(a-b)<tolerance)
		except AttributeError:
			return np.abs(a-b) < tolerance
	
	def sigma_3_stage(self,angle_tolerance=2.0,ijk_tolerance=0.1,symmetry=self.symmetry):
		angle_tolerance = angle_tolerance * (np.pi/180)
		sigma3 = 60 * (np.pi/180)
		def is_sigma_3(edge_index):  
			quat_2 = self.quats[edge_index]
			search_index = np.searchsorted(self.node_indices,edge_index)
			if search_index == edge_index:
				return False
			start_index = self.node_indices[search_index-1]
			quat_1 = self.quats[start_index]
			disori = xtal.qu2ax(symmetry.disoQu(quat_1,quat_2))
			ax = 3*disori[0,:3]/np.sum(disori[0,:3])
			angle = disori[0,3]
			return self.within_tolerance(ax,np.ones(3),ijk_tolerance) and self.within_tolerance(angle,sigma3,angle_tolerance)
	return is_sigma_3
	
	def generate_random_twin_walk(self,n,in_cluster=False,exclude_surface=True):
		nodes_list = []
		node_features = []
		edge_features = []
		initialize = True

		edge_index = np.random.choice(self.sigma_3_list)
		edge = self.neighbor_list[edge_index]
		start_node_index = self.node_indices[np.searchsorted(self.node_indices, start_edge, side="left")-1]
		start_node = self.neighbor_list[start_node_index]
		if start_node not in self.surface_features and end_node not in self.surface_features:
			initialize = False
			
		node_features.append(self.diameters[start_node])
		edge_features.append(self.misori_list[edge_index])

		if in_cluster:
			cluster_id = self.cluster_ids[start_node]
			
		for i in range(n):
			start_node_index = self.node_indices_sorted[np.searchsorted(self.node_id_sorted,edge)]
			start_node_index = self.node_indices[start_node_index]
			start_node = self.neighbor_list[start_node_index]
			num_neighbors = self.num_neighbors[start_node]
			
			edge_choice_indices = np.arange(start_node_index + 1,start_node_index + num_neighbors + 1)
			edge_choices = self.neighbor_list[start_node_index + 1:start_node_index + num_neighbors + 1]
			
			probs = np.ones(len(edge_choices))/len(edge_choices)
			if in_cluster:
				probs = np.array(self.cluster_ids[edge_choices] == cluster_id,dtype=int)
			if exclude_surface:
				probs = probs and not np.array(self.surface_features[edge_choices])

			probs = probs/np.sum(probs)

			edge_index = np.choice(edge_choice_indices,p=probs)
			edge = self.neighbor_list[edge_index]
			node_features.append(self.diameters[start_node])
			if not i == n-1:
				edge_features.append(self.misori_list[edge_index])

		return np.array(node_features),np.array(edge_features)
	
	def generate_full_twin_network(self,n,in_cluster=True)
		#ToDo
		pass

	def vtk_random_twin_network(self,cutoff):
		#ToDo
		pass
	
		
		


# helper class for cluster connectivity calculations
class UndirectedGraph:
	# create a graph object from a Neighborlist
	def __init__(self,NeighborList,GrainFeatures):
		self.adjacency_matrix = self.construct_adjacency_matrix(NeighborList,GrainFeatures)
		self.edges = self.construct_edge_list(self.adjacency_matrix)
		self.num_nodes = self.adjacency_matrix.shape[0]
		self.edge_vectors = self.

	def construct_adjacency_matrix(NeighborList,GrainFeatures):

	# computes the minimum number of edges to connect each node to the root node
	def computeDistance(self, root):
		if self.numNodes is 0:
			return []
		distance = [-1 for i in range(self.numNodes)]
		distance[root] = 0
		changed = True
		while changed:
			changed = False
			for e in self.edges:
				d0 = distance[e[0]]
				d1 = distance[e[1]]
				if d0 < 0 and d1 >= 0: # e0 unassigned, e1 assigned
					distance[e[0]] = d1+1
					changed = True
				elif d1 < 0 and d0 >= 0: # e0 assigned, e1 unassigned
					distance[e[1]] = d0+1
					changed = True
				elif d0 >= 0 and d1 >= 0: # e0 assigned, e1 assigned
					if d0+1 < d1:
						distance[e[1]] = d0+1
						changed = True
					elif d1+1 < d0:
						distance[e[0]] = d1+1
						changed = True
				else: # e0 unassigned, e1 unassigned
					pass
		if min(distance) < 0:
			print('warning, isolated nodes exist')
		return distance

	def distanceBetween(self, n1, n2):
		return self.computeDistance(n1)[n2]

	# computes the root(s) to minimize the maximum distance from the root to any other node
	def findMinRoot(self):
		if self.numNodes is 0:
			return None
		root = [0]
		maxDist = max(self.computeDistance(0))
		for i in range(1, self.numNodes):
			dist = max(self.computeDistance(i))
			if dist < maxDist:
				root = [i]
				maxDist = dist
			elif dist == maxDist:
				root.append(i)
		return root

	# computes the root(s) to maximize the maximum distance from the root to any other node
	def findMaxRoot(self):
		if self.numNodes is 0:
			return None
		root = [0]
		maxDist = max(self.computeDistance(0))
		for i in range(1, self.numNodes):
			dist = max(self.computeDistance(i))
			if dist > maxDist:
				root = [i]
				maxDist = dist
			elif dist == maxDist:
				root.append(i)
		return root

	def maxDist(self):
		if self.numNodes is 0:
			return None
		return max(self.computeDistance(self.findMaxRoot()[0]))
