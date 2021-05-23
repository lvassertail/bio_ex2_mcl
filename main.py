import os

import numpy as np
from clusters_ps import compute_clusters_ps
from gaf_reader import read_gaf
from network_reader import read_network
from mcl import create_clusters


def save_clustered_proteins(clusters, filename):
    f = open(filename, "wt")
    f.write("protein name\tcluster\n")
    for cluster_id, proteins in clusters.items():
        for p in proteins:
            f.write(f"{p}\t{cluster_id}\n")


def save_clusters_p_values(clusters_ps, filename):
    f = open(filename, "wt")
    f.write("cluster\tp-value\tGO-term\n")
    for cluster_id, (go_term, p_value) in clusters_ps.items():
        p_value = -np.log(p_value)
        f.write(f"{cluster_id}\t{p_value}\t{go_term}\n")


def main():
    network_graph = read_network("data/huri_symbol.tsv")
    goa = read_gaf("data/goa_human.gaf")
    clusters = create_clusters(network_graph)
    clusters_ps = compute_clusters_ps(clusters, goa)

    save_clustered_proteins(clusters, "output/clustered_proteins.txt")
    save_clusters_p_values(clusters_ps, "output/clusters_P.txt")

    for inflation in [1.5,1.75,2.0,2.25,2.5,2.75,3,4,5]:
        clusters = create_clusters(network_graph, inflation_param=inflation)
        clusters_ps = compute_clusters_ps(clusters, goa)

        save_clustered_proteins(clusters, f"output/per_inflation/clustered_proteins_{inflation:.2f}.txt")
        save_clusters_p_values(clusters_ps, f"output/per_inflation/clusters_P_{inflation:.2f}.txt")

        print("Average cluster size", np.mean([len(c) for c in clusters.values()]))


if __name__ == "__main__":
    main()