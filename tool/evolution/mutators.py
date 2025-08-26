import logging
from typing import List, Optional, Tuple, Any
import numpy as np
import random
import RNA
from tool.config import Config
from tool.utils.parallel import ParallelProcessor

logger = logging.getLogger(__name__)

class ConservationMutator:
    """Conservation-guided RNA sequence mutator."""
    
    def __init__(self, conservation_scores: Optional[List[float]] = None):
        """Initialize mutator with conservation scores.
        
        Args:
            conservation_scores: List of conservation scores per position
        """
        self.conservation = conservation_scores

    def generate(self, seed_seq: str, num_variants: int, target_structure: str, 
                max_attempts: int = 1000000) -> List[str]:
        """Generate variants using conservation-guided mutation.
        
        Args:
            seed_seq: Starting sequence
            num_variants: Number of variants to generate
            target_structure: Target secondary structure
            max_attempts: Maximum attempts for generation
            
        Returns:
            List of generated variant sequences
        """
        # Optimized pool size based on empirical success rates
        pool_size = min(max(num_variants * 50, 5000), max_attempts)
        args = [(seed_seq, self.conservation, target_structure) for _ in range(pool_size)]
        
        try:
            candidates = ParallelProcessor.parallel_map(self._generate_candidate, args)
        except Exception as e:
            logger.error(f"Error in parallel variant generation: {e}")
            return []
            
        # Filter valid and unique variants
        variants = []
        seen = set()
        for cand in candidates:
            if cand and cand not in seen:
                variants.append(cand)
                seen.add(cand)
            if len(variants) >= num_variants:
                break
                
        success_rate = len(variants) / pool_size * 100
        if len(variants) < num_variants:
            logger.warning(f"Only generated {len(variants)} variants after {pool_size} attempts (success rate: {success_rate:.2f}%)")
        else:
            logger.info(f"Successfully generated {len(variants)} variants (success rate: {success_rate:.2f}%)")
            
        return variants

    def _generate_candidate(self, args: Tuple[str, Optional[List[float]], str]) -> Optional[str]:
        """Generate a single candidate variant.
        
        Args:
            args: Tuple of (seed_seq, conservation, target_structure)
            
        Returns:
            Generated sequence if valid, None otherwise
        """
        seed_seq, conservation, target_structure = args
        
        try:
            if Config.MUTATION_STRATEGY == "fixed":
                n_mutations = Config.N_MUTATIONS
            else:
                weights = np.array(Config.MUTATION_WEIGHTS) / sum(Config.MUTATION_WEIGHTS)
                n_mutations = np.random.choice([1, 2, 3], p=weights)
                
            mutated = self._mutate_with_bias(seed_seq, n_mutations=n_mutations)
            struct, _ = RNA.fold(mutated)
            
            if struct == target_structure:
                return mutated
            return None
        except Exception as e:
            logger.debug(f"Error generating candidate: {e}")
            return None

    def _mutate_with_bias(self, seq: str, n_mutations: int = 20) -> str:
        """Mutate sequence with conservation bias.
        
        Args:
            seq: Input sequence
            n_mutations: Number of mutations to perform
            
        Returns:
            Mutated sequence
        """
        if not seq:
            return seq
            
        seq = list(seq)
        positions = range(len(seq))
        
        # Optimize mutation strategy based on sequence length
        if len(seq) <= 30:
            # For short sequences, use fewer mutations
            n_mutations = min(n_mutations, len(seq) // 3)
        else:
            # For longer sequences, use conservation bias
            n_mutations = min(n_mutations, len(seq) // 4)
        
        weights = np.ones(len(seq)) / len(seq)
        
        if self.conservation is not None:
            tolerance = 1 - np.array(self.conservation)
            tolerance = np.clip(tolerance, 0.1, 1.0)
            weights = (tolerance * Config.CONSERVATION_BIAS) + (np.ones(len(seq)) * (1 - Config.CONSERVATION_BIAS))
            weights /= weights.sum()
            
        try:
            mut_positions = np.random.choice(positions, size=n_mutations, replace=False, p=weights)
            for pos in mut_positions:
                # Optimize base choice for better structure compatibility
                current_base = seq[pos]
                alternatives = [b for b in 'ACGU' if b != current_base]
                
                # Prefer GC pairs for better structure stability
                if current_base in 'AU':
                    alternatives = ['G', 'C'] + [b for b in alternatives if b not in ['G', 'C']]
                
                seq[pos] = random.choice(alternatives)
        except ValueError as e:
            logger.error(f"Error in biased mutation: {e}")
            # Fallback to random mutation
            for _ in range(n_mutations):
                pos = random.randint(0, len(seq) - 1)
                seq[pos] = random.choice([b for b in 'ACGU' if b != seq[pos]])
                
        return ''.join(seq)