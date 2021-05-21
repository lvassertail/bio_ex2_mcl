import scipy.stats as ss
from tqdm import tqdm


def compute_clusters_ps(predicted_clusters, goa_clusters):
    predicted_clusters = {a: p for a, p in predicted_clusters.items() if len(p) >= 3}
    goa_clusters = {a: p for a, p in goa_clusters.items() if len(p) >= 3}
    n_total_proteins = sum(len(p) for p in goa_clusters.values())

    top_p_values = {}

    for predict_cluster, predict_proteins in tqdm(predicted_clusters.items()):
        p_value = float('inf')
        top_goa_cluster = None
        for goa_cluster, goa_proteins in goa_clusters.items():
            n_goa_proteins = len(goa_proteins)
            n_predicted_proteins = len(predict_proteins)
            n_proteins_from_goa = len(goa_proteins.intersection(predict_proteins))
            goa_c_p_value = ss.hypergeom(n_total_proteins, n_goa_proteins, n_predicted_proteins).sf(n_proteins_from_goa)
            if goa_c_p_value < p_value:
                p_value = goa_c_p_value
                top_goa_cluster = goa_cluster
        top_p_values[predict_cluster] = (top_goa_cluster, p_value)

    return top_p_values
