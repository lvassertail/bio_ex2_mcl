from collections import defaultdict


def read_gaf(filepath):
    f = open(filepath, "rt")
    clusters = defaultdict(set)

    for l in f:
        if l.startswith('!'):
            continue
        values = l.split('\t')
        if values[8] != 'P':
            continue
        protein, annotation = values[2], values[4]
        clusters[annotation].add(protein)
    return clusters
