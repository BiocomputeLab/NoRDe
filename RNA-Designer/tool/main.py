import logging
from typing import List, Tuple, Dict, Any
from tool.config import Config
from tool.core.sequence_analysis import SequenceAnalyzer
from tool.core.filtering import VariantFilter
from tool.core.grouping import GroupGenerator
from tool.visualization.plots import PlotGenerator
from tool.visualization.conservation import ConservationVisualizer
from tool.visualization.diversity_metrics import DiversityAnalyzer
from tool.utils.file_io import FileHandler
from tool.evolution.selector import VariantGenerator
from tool.analysis.conservation import ConservationAnalyzer

logger = logging.getLogger(__name__)

def aggregate_runs(n_runs: int = Config.N_RUNS, output_file: str = "output/aggregated_variants.fasta") -> List[str]:
    """Aggregate variants from multiple runs.
    
    Args:
        n_runs: Number of independent runs to perform
        output_file: Output FASTA file path
        
    Returns:
        List of final variant sequences
    """
    all_variants = set()
    
    for run in range(n_runs):
        logger.info(f"Starting Run {run+1}/{n_runs}")
        
        try:
            wt_struct, _ = SequenceAnalyzer.run_rnafold(Config.WILD_TYPE_SCAFFOLD)
            generator = VariantGenerator()
            raw_variants = generator.generate(wt_struct, Config.TARGET_VARIANT_COUNT*10)
            filtered = VariantFilter.filter_variants(raw_variants, Config.WILD_TYPE_SCAFFOLD, wt_struct)
            final = VariantFilter.enforce_lmax_filter(filtered, Config.WILD_TYPE_SCAFFOLD, Config.LMAX_THRESHOLD)
            
            all_variants.update(final[:Config.TARGET_VARIANT_COUNT])
            
            logger.info(f"Run {run+1} added {len(final[:Config.TARGET_VARIANT_COUNT])} variants")
            logger.info(f"Total unique variants so far: {len(all_variants)}")
        except Exception as e:
            logger.error(f"Error in run {run+1}: {e}")
            continue
    
    unique_variants = list(all_variants)
    actual_target = min(Config.TARGET_VARIANT_COUNT, len(unique_variants))
    
    if len(unique_variants) < Config.TARGET_VARIANT_COUNT:
        logger.warning(f"Only generated {len(unique_variants)} unique variants")
    
    final_variants = GroupGenerator.select_diverse_subset(unique_variants, actual_target)
    FileHandler.save_to_fasta(final_variants, output_file)
    return final_variants

def main() -> None:
    """Main workflow for RNA variant design and analysis."""
    Config.setup()
    
    # Set the seed if configured
    Config.set_seed(Config.RANDOM_SEED)
    
    # Warn about ViennaRNA reproducibility limitation
    logger.warning("⚠️  ViennaRNA operations (structure prediction, inverse folding) are NOT reproducible due to C library limitations")
    logger.info("Only mutation operations and clustering are seeded. Final variant sequences may vary between runs.")

    final_variants = aggregate_runs(n_runs=Config.N_RUNS)
    
    # --- Wild-type conservation analysis and visualization ---
    try:
        cons_analyzer = ConservationAnalyzer(Config.WILD_TYPE_SCAFFOLD)
        conservation_scores = cons_analyzer.analyze()

        # Heatmap of conservation
        ConservationVisualizer.plot_heatmap(
            Config.WILD_TYPE_SCAFFOLD,
            conservation_scores,
            "wildtype_conservation.png"
        )

        # Fold wild-type sequence to get secondary structure
        wt_struct, _ = SequenceAnalyzer.run_rnafold(Config.WILD_TYPE_SCAFFOLD)

        # Secondary structure plot colored by conservation scores
        ConservationVisualizer.plot_secondary_structure(
            sequence=Config.WILD_TYPE_SCAFFOLD,
            structure=wt_struct,
            conservation_scores=conservation_scores,
            filename="wildtype_secondary_structure_conservation"
        )
    except Exception as e:
        logger.error(f"Error in conservation analysis: {e}")

    # --- Group generation and further analysis ---
    try:
        groups = GroupGenerator.generate_diverse_groups(
            final_variants, 
            Config.GROUP_SIZE, 
            Config.N_GROUPS,
            Config.WILD_TYPE_SCAFFOLD
        )
        
        FileHandler.save_to_fasta(final_variants, "output/variants.fasta")
        FileHandler.save_groups_to_fasta(groups, "output/groups.fasta")
        
        PlotGenerator.plot_difference_heatmap(
            Config.WILD_TYPE_SCAFFOLD, 
            final_variants, 
            "output/variants_heatmap.png"
        )
        PlotGenerator.plot_hamming_heatmap(final_variants, "output/hamming_heatmap.png")
        PlotGenerator.plot_lmax_heatmap(final_variants, "output/lmax_heatmap.png")
        
        mean_intra, mean_inter = DiversityAnalyzer.compute_group_diversity(groups)
        logger.info(f"Mean Intra-group Hamming Distance: {mean_intra:.2f}")
        logger.info(f"Mean Inter-group Hamming Distance: {mean_inter:.2f}")
        
        lmax_results = DiversityAnalyzer.compute_lmax_analysis(groups)
        logger.info(f"Mean Intra-group Lmax: {lmax_results['mean_intra_lmax']:.2f}")
        logger.info(f"Mean Inter-group Lmax: {lmax_results['mean_inter_lmax']:.2f}")
        logger.info("Per-variant Lmax (max shared subsequence with any other variant):")
        for i, val in enumerate(lmax_results['per_variant_lmax']):
            logger.info(f"  Variant {i+1}: {val}")
    except Exception as e:
        logger.error(f"Error in group analysis: {e}")

if __name__ == "__main__":
    main()