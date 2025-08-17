import logging
from typing import Tuple, Optional
import RNA
import numpy as np
from functools import lru_cache
from tool.utils.caching import CacheManager

logger = logging.getLogger(__name__)
cache_manager = CacheManager()

class SequenceAnalyzer:
    """Static class for RNA sequence analysis operations."""
    
    @staticmethod
    @lru_cache(maxsize=100000)
    def run_rnafold(sequence: str) -> Tuple[str, float]:
        """Run RNAfold to get secondary structure and MFE.
        
        Args:
            sequence: RNA sequence string
            
        Returns:
            Tuple of (structure, mfe) where structure is dot-bracket notation
        """
        cached = cache_manager.rnafold_cache.get(sequence)
        if cached is not None:
            return cached
        try:
            fc = RNA.fold_compound(sequence)
            result = fc.mfe()
            cache_manager.rnafold_cache.set(sequence, result)
            return result
        except (RuntimeError, ValueError) as e:
            logger.error(f"Error running RNAfold for sequence '{sequence}': {e}")
            return ("", 0.0)

    @staticmethod
    @lru_cache(maxsize=100000)
    def gc_content(seq: str) -> float:
        """Calculate GC content percentage of a sequence.
        
        Args:
            seq: DNA/RNA sequence string
            
        Returns:
            GC content as percentage (0-100)
        """
        cached = cache_manager.gc_content_cache.get(seq)
        if cached is not None:
            return cached
        
        if not seq:
            return 0.0
            
        result = (sum(1 for base in seq if base in "GC") / len(seq)) * 100
        cache_manager.gc_content_cache.set(seq, result)
        return result

    @staticmethod
    def has_homopolymer(seq: str, max_run: int = 3) -> bool:
        """Check if sequence contains homopolymer runs longer than max_run.
        
        Args:
            seq: DNA/RNA sequence string
            max_run: Maximum allowed consecutive identical bases
            
        Returns:
            True if homopolymer run exceeds max_run
        """
        if not seq:
            return False
            
        count = 1
        for i in range(1, len(seq)):
            count = count + 1 if seq[i] == seq[i-1] else 1
            if count > max_run: 
                return True
        return False

    @staticmethod
    def longest_common_substring(s1: str, s2: str) -> int:
        """Find length of longest common substring between two sequences.
        
        Args:
            s1: First sequence string
            s2: Second sequence string
            
        Returns:
            Length of longest common substring
        """
        if not s1 or not s2:
            return 0
            
        m, n = len(s1), len(s2)
        dp = [[0]*(n+1) for _ in range(m+1)]
        max_len = 0
        for i in range(m):
            for j in range(n):
                if s1[i] == s2[j]:
                    dp[i+1][j+1] = dp[i][j] + 1
                    max_len = max(max_len, dp[i+1][j+1])
        return max_len

    @staticmethod
    def structure_confidence(seq: str, target_structure: str) -> Tuple[float, float]:
        """Calculate structure probability and diversity for a sequence.
        
        Args:
            seq: RNA sequence string
            target_structure: Target structure in dot-bracket notation
            
        Returns:
            Tuple of (probability, diversity)
        """
        try:
            fc = RNA.fold_compound(seq)
            fc.pf()
            try:
                prob = fc.pr_structure(target_structure)
            except RuntimeError:
                prob = 0.0
            diversity = fc.mean_bp_distance()
            return prob, diversity
        except (RuntimeError, ValueError) as e:
            logger.error(f"Error calculating structure confidence for sequence '{seq}': {e}")
            return 0.0, 0.0