
from network_reader import read_network
from mcl import create_clusters

def main():
    network_graph = read_network("data/huri_symbol.tsv")
    clusters = create_clusters(network_graph.graph)

if __name__ == "__main__":
    main()