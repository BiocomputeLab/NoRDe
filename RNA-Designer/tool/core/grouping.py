import numpy as np
import random
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans
from sklearn.feature_extraction.text import CountVectorizer
from tool.core.sequence_analysis import SequenceAnalyzer
from tool.config import Config

class GroupGenerator:
    @staticmethod
    def generate_diverse_groups(variants, group_size, n_groups, wt_seq):
        seq_arrays = [[ord(c) for c in seq] for seq in variants]
        dist_matrix = pairwise_distances(seq_arrays, metric='hamming') * len(variants[0])

        selected_indices = set()
        groups = []

        for _ in range(n_groups):
            group = []
            available = list(set(range(len(variants))) - selected_indices)
            if not available:
                break

            if not group:
                wt_array = [ord(c) for c in wt_seq.replace("T", "U")]
                dists_to_wt = [np.sum(np.array(seq) != np.array(wt_array)) for seq in seq_arrays]
                first = max(available, key=lambda i: dists_to_wt[i])
            else:
                first = random.choice(available)
            
            group.append(first)
            selected_indices.add(first)

            while len(group) < group_size and len(selected_indices) < len(variants):
                remaining = list(set(range(len(variants))) - selected_indices)
                scores = []

                for idx in remaining:
                    min_dist = min(dist_matrix[idx][g] for g in group)
                    scores.append((min_dist, idx))
                
                _, best_idx = max(scores)
                group.append(best_idx)
                selected_indices.add(best_idx)
            
            groups.append([variants[i] for i in group])

        return groups

    @staticmethod
    def select_diverse_subset(variants, n):
        if n == 0:
            return []
        if len(variants) <= n:
            return variants
        if len(variants) < 100:
            return GroupGenerator.distance_based_selection(variants, n)
        else:
            return GroupGenerator.cluster_based_selection(variants, n)

    @staticmethod
    def distance_based_selection(variants, n):
        seq_arrays = [[ord(c) for c in seq] for seq in variants]
        dist_matrix = pairwise_distances(seq_arrays, metric='hamming') * len(variants[0])
        
        selected = []
        remaining_indices = set(range(len(variants)))
        
        first = max(remaining_indices, key=lambda i: np.sum(dist_matrix[i]))
        selected.append(first)
        remaining_indices.remove(first)
        
        while len(selected) < n and remaining_indices:
            best_score = -1
            best_idx = -1
            
            for idx in remaining_indices:
                min_dist = min(dist_matrix[idx][s] for s in selected)
                if min_dist > best_score:
                    best_score = min_dist
                    best_idx = idx
                    
            if best_idx == -1:
                break
                
            selected.append(best_idx)
            remaining_indices.remove(best_idx)
        
        return [variants[i] for i in selected]

    @staticmethod
    def cluster_based_selection(variants, n):
        vectorizer = CountVectorizer(analyzer='char', ngram_range=(1,3))
        X = vectorizer.fit_transform(variants)
        
        kmeans = KMeans(n_clusters=n, random_state=Config.RANDOM_SEED)
        clusters = kmeans.fit_predict(X)
        
        selected = []
        for cluster_id in range(n):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if len(cluster_indices) > 0:
                center = kmeans.cluster_centers_[cluster_id]
                distances = np.linalg.norm(X[cluster_indices].toarray() - center, axis=1)
                selected.append(variants[cluster_indices[np.argmin(distances)]])
        
        return selected