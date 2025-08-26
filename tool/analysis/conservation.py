import numpy as np
import RNA
from pathlib import Path
from tool.config import Config
from tool.core.sequence_analysis import SequenceAnalyzer
from tool.utils.parallel import ParallelProcessor

class ConservationAnalyzer:
    """Analyzes nucleotide conservation through mutation tolerance"""
    
    BASES = ['A', 'T', 'C', 'G']
    
    def __init__(self, sequence):
        self.seq = sequence
        self.struct, self.mfe = SequenceAnalyzer.run_rnafold(sequence)
        self.cache_file = Path(Config.OUTPUT_DIR) / "conservation_cache.npy"
        
    def analyze(self, attempts_per_pos=Config.CONSERVATION_ATTEMPTS):
        """
        Run conservation analysis with caching
        Returns: Normalized conservation scores (1=conserved, 0=tolerant)
        """
        if self.cache_file.exists():
            return np.load(self.cache_file)
            
        results = self._run_analysis(attempts_per_pos)
        np.save(self.cache_file, results['conservation'])
        return results['conservation']
    
    def _run_analysis(self, attempts_per_pos):
        """Core analysis workflow"""
        L = len(self.seq)
        mut_attempts = np.zeros((L, 3), dtype=float)  # 3 possible mutations per position
        mut_accepts = np.zeros((L, 3), dtype=float)
        
        # Parallelize position analysis
        results = ParallelProcessor.parallel_map(self._analyze_position, range(L))
        
        # Aggregate results
        for pos, (pos_attempts, pos_accepts) in enumerate(results):
            mut_attempts[pos] = pos_attempts
            mut_accepts[pos] = pos_accepts
            
        return self._calculate_scores(mut_attempts, mut_accepts)
    
    def _analyze_position(self, pos):
        """Analyze single position (runs in parallel)"""
        ref_base = self.seq[pos]
        alternatives = [b for b in self.BASES if b != ref_base]
        pos_attempts = np.zeros(3)
        pos_accepts = np.zeros(3)
        
        for j, alt_base in enumerate(alternatives):
            for _ in range(Config.CONSERVATION_ATTEMPTS):
                variant = self._create_variant(pos, alt_base)
                pos_attempts[j] += 1
                if self._is_valid_variant(variant):
                    pos_accepts[j] += 1
                    
        return pos_attempts, pos_accepts
    
    def _create_variant(self, pos, alt_base):
        """Generate single-point mutant"""
        return self.seq[:pos] + alt_base + self.seq[pos+1:]
    
    def _is_valid_variant(self, variant):
        """Check if variant meets all criteria"""
        struct, mfe = SequenceAnalyzer.run_rnafold(variant)
        return (struct == self.struct and
                passes_gc_filter(variant) and
                abs(mfe - self.mfe) <= Config.MFE_TOLERANCE and
                not SequenceAnalyzer.has_homopolymer(variant))
    
    def _calculate_scores(self, attempts, accepts):
        """Convert raw counts to conservation scores"""
        tolerance_matrix = np.divide(
            accepts, attempts,
            out=np.zeros_like(attempts),
            where=attempts != 0
        )
        return {
            'conservation': 1 - np.mean(tolerance_matrix, axis=1),
            'tolerance_matrix': tolerance_matrix,
            'attempts': attempts,
            'accepts': accepts
        }

# Helper functions (from original script)
def mutation_tolerance_analysis(seq, attempts=100):
    """Legacy interface from your original script"""
    analyzer = ConservationAnalyzer(seq)
    results = analyzer._run_analysis(attempts)
    return (results['conservation'], 
            results['tolerance_matrix'],
            results['attempts'], 
            results['accepts'])

def passes_gc_filter(seq, min_gc=0.3, max_gc=0.7):
    """Checks if the GC content of the sequence is within the allowed range."""
    gc_count = seq.count('G') + seq.count('C')
    gc_content = gc_count / len(seq)
    return min_gc <= gc_content <= max_gc