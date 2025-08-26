import RNA
import random
from functools import lru_cache

class InverseFolding:
    @staticmethod
    @lru_cache(maxsize=100)  # Cache frequent structures
    def generate(target_structure, num_variants, max_attempts=10000):
        variants = set()
        while len(variants) < num_variants and len(variants) < max_attempts:
            try:
                # Generate with random seed
                seed = "".join(random.choice("ACGU") for _ in range(len(target_structure)))
                seq, _ = RNA.inverse_fold(seed, target_structure)
                variants.add(seq)
            except RuntimeError:
                continue
        return list(variants)
    
    @staticmethod
    def get_seed_sequence(target_structure):
        """Get single high-quality sequence for mutation seeding"""
        return RNA.inverse_fold(
            "".join(random.choice("ACGU") for _ in range(len(target_structure))),
            target_structure
        )[0]