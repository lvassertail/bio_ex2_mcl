import matplotlib.pyplot as plt


def main():
    inf_values = ["1.50", "1.75", "2.00", "2.25", "2.50", "2.75", "3.00", "4.00", "5.00"]
    for inf in inf_values:
        f_p = open(f'output/per_inflation/clusters_P_{inf}.txt')
        f_c = open(f'output/per_inflation/clustered_proteins_{inf}.txt')
        clusters = [l.split('\t')[1] for i, l in enumerate(f_c) if i != 0]
        clusters_ps = [float(l.split('\t')[1]) for i, l in enumerate(f_p) if i != 0]
        clusters_ps.sort(reverse=True)
        plt.plot(clusters_ps, label=inf)
        print(f"Average size for {inf}: {len(clusters) / len(clusters_ps)}")

    plt.legend()
    # plt.show()
    plt.savefig('output/compare_inflate.pdf')


if __name__ == "__main__":
    main()