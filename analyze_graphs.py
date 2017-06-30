def compute_distance_matrix(adj_matrices,kernel):
    computed = []
    distance_matrix = np.zeros((len(adj_matrices),len(adj_matrices)))
    for i,adj_i in enumerate(adj_matrices):
        for j in range(i+1,len(adj_matrices)):
            if i ==j or (i,j) in computed or (j,i) in computed:
                pass
            else:
                adj_i = adj_matrices[j]
                distance = kernel._compute(adj_i,adj_j)
                distance_matrix[i,j] = distance
                distance_matrix[j,i] = distance
                computed.append((i,j))
    return distance_matrix


    