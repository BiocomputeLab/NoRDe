import numpy as np
from tool.core.sequence_analysis import SequenceAnalyzer

class DiversityAnalyzer:
    @staticmethod
    def compute_group_diversity(groups):
        all_seqs = [seq for group in groups for seq in group]
        seq_idx = {seq: idx for idx, seq in enumerate(all_seqs)}
        n = len(all_seqs)
        
        dist_matrix = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                dist_matrix[i, j] = sum(c1 != c2 for c1, c2 in zip(all_seqs[i], all_seqs[j]))
        
        intra_dists = []
        for group in groups:
            indices = [seq_idx[seq] for seq in group]
            for i in range(len(indices)):
                for j in range(i+1, len(indices)):
                    intra_dists.append(dist_matrix[indices[i], indices[j]])
        mean_intra = np.mean(intra_dists) if intra_dists else 0
        
        inter_dists = []
        for i, group1 in enumerate(groups):
            for j, group2 in enumerate(groups):
                if i >= j:
                    continue
                for s1 in group1:
                    for s2 in group2:
                        inter_dists.append(dist_matrix[seq_idx[s1], seq_idx[s2]])
        mean_inter = np.mean(inter_dists) if inter_dists else 0
        
        return mean_intra, mean_inter

    @staticmethod
    def compute_lmax_analysis(groups):
        all_seqs = [seq for group in groups for seq in group]
        seq_idx = {seq: i for i, seq in enumerate(all_seqs)}
        n = len(all_seqs)

        lmax_matrix = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j:
                    lmax_matrix[i][j] = SequenceAnalyzer.longest_common_substring(all_seqs[i], all_seqs[j])

        intra_lmax_vals = []
        for group in groups:
            indices = [seq_idx[seq] for seq in group]
            for i in range(len(indices)):
                for j in range(i+1, len(indices)):
                    intra_lmax_vals.append(lmax_matrix[indices[i]][indices[j]])
        mean_intra_lmax = np.mean(intra_lmax_vals) if intra_lmax_vals else 0

        inter_lmax_vals = []
        for i, group1 in enumerate(groups):
            for j, group2 in enumerate(groups):
                if i >= j:
                    continue
                for s1 in group1:
                    for s2 in group2:
                        inter_lmax_vals.append(lmax_matrix[seq_idx[s1], seq_idx[s2]])
        mean_inter_lmax = np.mean(inter_lmax_vals) if inter_lmax_vals else 0

        per_variant_lmax = np.max(lmax_matrix, axis=1)

        return {
            'mean_intra_lmax': mean_intra_lmax,
            'mean_inter_lmax': mean_inter_lmax,
            'per_variant_lmax': per_variant_lmax
        }