import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tool.core.sequence_analysis import SequenceAnalyzer

class PlotGenerator:
    @staticmethod
    def plot_difference_heatmap(wt_seq, variants, filename):
        plt.figure(figsize=(12,6), dpi=300)
        sns.heatmap(
            np.array([[1 if v[i]!=wt_seq[i] else 0 for i in range(len(wt_seq))] for v in variants]),
            cmap="Reds", cbar=True, linewidths=0.5, linecolor='gray',
            xticklabels=list(range(1,len(wt_seq)+1)),
            yticklabels=[f"Var {i+1}" for i in range(len(variants))],
            square=True,
            cbar_kws={"shrink": 0.8, "aspect": 20},
            vmin=0, vmax=1
        )
        plt.xticks(rotation=90, fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.title("Nucleotide Differences vs Wild-Type")
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)
        plt.close()

    @staticmethod
    def plot_hamming_heatmap(sequences, filename):
        n = len(sequences)
        matrix = np.zeros((n, n), dtype=int)
        
        for i in range(n):
            for j in range(n):
                matrix[i][j] = sum(c1 != c2 for c1, c2 in zip(sequences[i], sequences[j]))

        plt.figure(figsize=(10, 8), dpi=300)
        sns.heatmap(matrix, annot=True, fmt="d", cmap="coolwarm", 
                    xticklabels=[f"V{i+1}" for i in range(n)],
                    yticklabels=[f"V{i+1}" for i in range(n)],
                    square=True,
                    cbar_kws={"shrink": 0.8, "aspect": 20})
        plt.xticks(rotation=90, fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.title("Hamming Distance Between Variants")
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)
        plt.close()

    @staticmethod
    def plot_lmax_heatmap(sequences, filename):
        n = len(sequences)
        matrix = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                matrix[i][j] = SequenceAnalyzer.longest_common_substring(sequences[i], sequences[j])

        plt.figure(figsize=(10, 8), dpi=300)
        sns.heatmap(matrix, annot=True, fmt="d", cmap="YlGnBu", 
                    xticklabels=[f"V{i+1}" for i in range(n)],
                    yticklabels=[f"V{i+1}" for i in range(n)],
                    square=True,
                    cbar_kws={"shrink": 0.8, "aspect": 20})
        plt.xticks(rotation=90, fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.title("Lmax (Longest Shared Subsequence) Matrix")
        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)
        plt.close()