from tool.core.sequence_analysis import SequenceAnalyzer
from tool.config import Config
from tool.utils.parallel import ParallelProcessor

class VariantFilter:
    @staticmethod
    def filter_variants(variants, wt_seq, wt_struct):
        results = ParallelProcessor.parallel_map(SequenceAnalyzer.run_rnafold, variants)

        variant_info = []
        for seq, (struct, mfe) in zip(variants, results):
            if struct != wt_struct or not (Config.GC_CONTENT_RANGE[0] <= SequenceAnalyzer.gc_content(seq) <= Config.GC_CONTENT_RANGE[1]) or SequenceAnalyzer.has_homopolymer(seq):
                continue
            
            prob, diversity = SequenceAnalyzer.structure_confidence(seq, wt_struct)
            if prob < Config.MIN_STRUCTURE_PROB or diversity > Config.MAX_STRUCTURE_DIVERSITY: 
                continue

            variant_info.append((seq, mfe, prob, diversity))

        variant_info.sort(key=lambda x: (-x[2], x[1]))
        top_variants = [seq for seq, _, _, _ in variant_info[:Config.TARGET_VARIANT_COUNT]]
        return top_variants

    @staticmethod
    def enforce_lmax_filter(variants, wt_seq, threshold=Config.LMAX_THRESHOLD):
        passed = []
        for candidate in variants:
            max_lmax = max(SequenceAnalyzer.longest_common_substring(candidate, seq) for seq in passed + [wt_seq])
            if max_lmax <= threshold:
                passed.append(candidate)
        return passed