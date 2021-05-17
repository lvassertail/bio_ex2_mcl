import numpy as np
from scipy.sparse import csc_matrix


def create_stochastic_matrix(network):

    graph = network.graph
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


def expansion(matrix):
    pass


def inflation(matrix):
    pass


def pruning(matrix):
    pass


def is_converged(matrix):
    pass


def get_clusters(matrix):
    pass


def create_clusters(network):
    matrix = create_stochastic_matrix(network)

    for step in range(0,20):
        expansion(matrix)
        inflation(matrix)
        pruning(matrix)

        if is_converged(matrix):
            break

    return get_clusters(matrix)