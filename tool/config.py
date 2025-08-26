import os
import logging
from pathlib import Path
from multiprocessing import cpu_count
import random
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# === RNA VARIANT DESIGN TOOLKIT CONFIGURATION ===
class Config:
    # --- Scaffold & Output ---
    WILD_TYPE_SCAFFOLD = "GTGAACTGCCGAGTAGGTAGCTGATAAC"  # Reference RNA sequence
    OUTPUT_DIR = "output"                                # Output directory for all results

    # --- Variant Generation ---
    TARGET_VARIANT_COUNT = 50                            # Number of variants to generate per run
    VARIANTS_PER_RUN_MULTIPLIER = 10                     # Multiplier for candidate pool size
    N_RUNS = 10                                          # Number of independent runs

    # --- Grouping ---
    GROUP_SIZE = 12                                      # Number of variants per group
    N_GROUPS = 1                                         # Number of groups to generate

    # --- Filtering & Constraints ---
    LMAX_THRESHOLD = 10                                  # Max allowed shared subsequence length (Lmax)
    GC_CONTENT_RANGE = (40, 60)                          # Allowed GC content percentage (min, max)
    MIN_STRUCTURE_PROB = 0.7                             # Minimum structure probability (ViennaRNA)
    MAX_STRUCTURE_DIVERSITY = 20                         # Maximum allowed structure diversity
    MFE_TOLERANCE = 5                                    # Allowed MFE deviation from wild-type (kcal/mol)

    # --- Parallelization ---
    NUM_CORES = max(1, cpu_count() - 1)                  # Number of CPU cores to use

    # --- Generation Method Control ---
    GENERATION_METHOD = "auto"                           # "auto" | "inverse" | "conservation"
    AUTO_SWITCH_LENGTH = 30                              # Length threshold for auto mode (≤ uses conservation-guided)
    MUTATION_STRATEGY = "fixed"                          # "fixed" | "random"

    # --- Conservation-Guided Mutation Parameters ---
    CONSERVATION_ATTEMPTS = 100                          # Attempts per position for conservation analysis
    CONSERVATION_BIAS = 0.7                              # [0-1] Strength of conservation bias in mutation selection
    N_MUTATIONS = 20                                     # Number of mutations per candidate (conservation-guided mode)
    MUTATION_WEIGHTS = [70, 25, 5]                       # Probabilities for 1, 2, or 3 mutations (if randomizing)

    # --- Inverse Folding Parameters ---
    INVERSE_MAX_ATTEMPTS = 10000                         # Max attempts for inverse folding per run

    # --- Reproducibility ---
    RANDOM_SEED = 1                                   # Set to an integer for reproducible runs

    @classmethod
    def setup(cls) -> None:
        """Create output directory if needed."""
        try:
            os.makedirs(cls.OUTPUT_DIR, exist_ok=True)
            logger.info(f"Output directory '{cls.OUTPUT_DIR}' ready")
        except OSError as e:
            logger.error(f"Error creating output directory '{cls.OUTPUT_DIR}': {e}")
            raise

    @classmethod
    def set_seed(cls, seed: int = None) -> None:
        """Set random seed for reproducible results.
        
        ⚠️  IMPORTANT REPRODUCIBILITY LIMITATION ⚠️
        
        ViennaRNA operations (structure prediction, inverse folding) are NOT 
        reproducible due to C library limitations. The seed only controls:
        
        ✅ Reproducible:
        - Mutation positions and base choices in variant generation
        - Clustering algorithms (K-means for grouping)
        - Random sampling operations
        - NumPy and Python random operations
        
        ❌ NOT Reproducible:
        - RNA structure prediction results (RNA.fold)
        - Inverse folding results (RNA.inverse_fold)
        - Final variant sequences (due to ViennaRNA dependency)
        
        For full reproducibility, consider:
        1. Caching ViennaRNA results
        2. Using alternative RNA tools with seeding support
        3. Documenting ViennaRNA version and system environment
        
        Args:
            seed: Integer seed value. If None, uses system random seed.
        """
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
            # Set scikit-learn random state for reproducible clustering
            try:
                from sklearn.utils import check_random_state
                check_random_state(seed)
            except ImportError:
                logger.warning("scikit-learn not available, clustering may not be reproducible")
            cls.RANDOM_SEED = seed
            logger.info(f"Random seed set to: {seed} (reproducible mode)")
            logger.warning("ViennaRNA operations remain non-reproducible due to C library limitations")
        else:
            logger.info("Using system random seed (non-reproducible mode)")



# --- Parameter Documentation ---
"""
Parameter Overview:

- Scaffold & Output:
    WILD_TYPE_SCAFFOLD: Reference RNA sequence for variant design.
    OUTPUT_DIR: Directory for all output files.

- Variant Generation:
    TARGET_VARIANT_COUNT: Number of variants to generate per run.
    VARIANTS_PER_RUN_MULTIPLIER: Multiplier for candidate pool size.
    N_RUNS: Number of independent runs for diversity.

- Grouping:
    GROUP_SIZE: Number of variants per group.
    N_GROUPS: Number of groups to generate.

- Filtering & Constraints:
    LMAX_THRESHOLD: Maximum allowed shared subsequence length (Lmax).
    GC_CONTENT_RANGE: Allowed GC content percentage (min, max).
    MIN_STRUCTURE_PROB: Minimum structure probability for acceptance.
    MAX_STRUCTURE_DIVERSITY: Maximum allowed structure diversity.
    MFE_TOLERANCE: Allowed deviation in minimum free energy (MFE).

- Parallelization:
    NUM_CORES: Number of CPU cores to use for parallel tasks.

- Generation Method Control:
    GENERATION_METHOD: "auto" (hybrid), "inverse", or "conservation".
    AUTO_SWITCH_LENGTH: Length threshold for auto mode. (Recommended ≤ 30)

- Conservation-Guided Mutation:
    CONSERVATION_ATTEMPTS: Attempts per position for conservation analysis.
    CONSERVATION_BIAS: Strength of conservation bias in mutation selection.
    N_MUTATIONS: Number of mutations per candidate.
    MUTATION_WEIGHTS: Probabilities for 1, 2, or 3 mutations (if randomizing).

- Inverse Folding:
    INVERSE_MAX_ATTEMPTS: Max attempts for inverse folding per run.

- Reproducibility:
    RANDOM_SEED: Set for reproducible results.

Mode-Specific Notes:
- Conservation-guided mode uses: N_MUTATIONS, CONSERVATION_ATTEMPTS, CONSERVATION_BIAS, MUTATION_WEIGHTS.
- Inverse folding mode uses: INVERSE_MAX_ATTEMPTS.
- Auto mode uses: AUTO_SWITCH_LENGTH to select between the above.

"""