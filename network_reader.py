
import csv
from collections import defaultdict
from network import Network


def remove_node_from_graph(graph, node):
    for neighbor in graph[node]:
        graph[neighbor].remove(node)
    graph.pop(node)


def keep_k_core_network(graph, k):

    nodes_to_check = set(graph.keys())
    while len(nodes_to_check) > 0:
        next_nodes_to_check = set()
        for node in nodes_to_check:
            if len(graph[node]) < k:
                next_nodes_to_check = next_nodes_to_check.union(graph[node])
                remove_node_from_graph(graph, node)
        nodes_to_check = next_nodes_to_check.intersection(graph.keys())


def create_network_from_graph(graph_by_names):
    proteins_id_to_name = defaultdict()
    proteins_name_to_id = defaultdict()
    graph_by_ids = defaultdict(set)

    next_protein_id = 0
    for protein_name in graph_by_names.keys():
        proteins_id_to_name[next_protein_id] = protein_name
        proteins_name_to_id[protein_name] = next_protein_id
        next_protein_id += 1

    for protein_name, neighbors in graph_by_names.items():
        protein_id = proteins_name_to_id[protein_name]
        for neighbor_name in neighbors:
            neighbor_id = proteins_name_to_id[neighbor_name]
            graph_by_ids[protein_id].add(neighbor_id)

    return Network(proteins_id_to_name, graph_by_ids)

def read_network(network_file_name):
    graph = defaultdict(set)

    with open(network_file_name, 'rt') as f:
        lines = list(csv.reader(f, delimiter='\t'))

        for node1, node2 in lines:
            # skipping any self-connected nodes
            if node1 == node2:
                continue
            graph[node1].add(node2)
            graph[node2].add(node1)

        keep_k_core_network(graph, 3)

        network = create_network_from_graph(graph)
        return network
