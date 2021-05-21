import numpy as np
from scipy.sparse import csc_matrix
from collections import defaultdict


def create_stochastic_matrix(graph):

    row = []
    col = []
    weights = []

    for protein_idx, neighbors in graph.items():
        # Add the neighbors to the matrix including the self-loop
        neighbors_with_self = neighbors.union({protein_idx})
        norm_weight = 1 / (len(neighbors_with_self))

        for neighbor in neighbors_with_self:
            col.append(protein_idx)
            row.append(neighbor)
            weights.append(norm_weight)

    return csc_matrix((np.array(weights), (np.array(row), np.array(col))),
                      shape=(len(graph), len(graph)))


def create_test_matrix():

    graph = {0:{1,2,3}, 1:{0,3}, 2:{0}, 3:{0,1}}
    return create_stochastic_matrix(graph)


def normalize(matrix):
    node_weight_sum_array = np.squeeze(np.asarray(matrix.sum(0)))
    rows, cols = matrix.nonzero()
    node_weight_sum_array_mat = csc_matrix(((1.0/node_weight_sum_array)[cols], (rows,cols)), shape=(matrix.shape))
    return matrix.multiply(node_weight_sum_array_mat)


def expand(matrix):
    matrix = matrix.dot(matrix)
    return normalize(matrix)


def inflate(matrix, inflation_param):
    matrix = matrix.power(inflation_param)
    return normalize(matrix)


def prune(matrix):
    node_max_weights_indices = np.squeeze(np.asarray(matrix.argmax(0)))
    node_max_weights = matrix.max(0).data
    values_to_keep_mask = (matrix > 1e-6)

    matrix = matrix.multiply(values_to_keep_mask)

    # The columns should not be empty - bring back one value that was the maximum before the pruning
    for col_idx in range(len(node_max_weights_indices)):
        max_val_row_idx = node_max_weights_indices[col_idx]
        if matrix[max_val_row_idx, col_idx] == 0:
            matrix[max_val_row_idx, col_idx] = node_max_weights[col_idx]

    return matrix


def is_converged(matrix):

    node_max_weight = matrix.max(0).data
    rows, cols = matrix.nonzero()
    node_max_weight_for_every_row = \
        csc_matrix((node_max_weight[cols], (rows,cols)), shape=(matrix.shape))

    # the matrix converges if for each columns there is only one value,
    # or all the values are equal (of almost equal)
    if ((node_max_weight_for_every_row - matrix) > 1e-6).sum() == 0:
        return True

    return False


def get_clusters(matrix, proteins):
    all_clusters_by_id = defaultdict(set)
    final_clusters_by_name = defaultdict(set)

    rows, cols = matrix.nonzero()
    for cluster_num, protein_id in zip(rows,cols):
        all_clusters_by_id[cluster_num].add(protein_id)

    cluster_id = 1
    for _, cluster in all_clusters_by_id.items():
        if len(cluster) >= 5:
            names_cluster = set(proteins[protein_id] for protein_id in cluster)
            final_clusters_by_name[cluster_id] = names_cluster
            cluster_id += 1

    return final_clusters_by_name


def create_clusters(network, inflation_param=2, steps=20):
    matrix = create_stochastic_matrix(network.graph)
    #matrix = create_test_matrix()

    for step in range(0, steps):
        print("start step ", step)

        matrix = expand(matrix)
        matrix = inflate(matrix, inflation_param)
        matrix = prune(matrix)

        print("end step ", step)

        if is_converged(matrix):
            break

    return get_clusters(matrix, network.proteins)