import matplotlib.pyplot as plt


def main():
    inf_values = ["1.50", "2.20", "2.90", "3.60", "4.30", "5.00"]
    for inf in inf_values:
        f = open(f'output/per_inflation/clusters_P_{inf}.txt')
        clusters_ps = [float(l.split('\t')[1]) for i, l in enumerate(f) if i != 0]
        clusters_ps.sort(reverse=True)
        plt.plot(clusters_ps, label=inf)

    plt.legend()
    # plt.show()
    plt.savefig('output/compare_inflate.png')


if __name__ == "__main__":
    main()