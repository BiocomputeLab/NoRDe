from tool.config import Config
from tool.evolution.inverse_folding import InverseFolding
from tool.evolution.mutators import ConservationMutator

class VariantGenerator:
    def __init__(self, conservation_scores=None):
        self.inverse_folder = InverseFolding()
        self.mutator = ConservationMutator(conservation_scores)
        
    def generate(self, target_structure, num_variants):
        if self._use_inverse_folding(target_structure):
            return self.inverse_folder.generate(target_structure, num_variants)
        else:
            return self._generate_with_conservation(target_structure, num_variants)
    
    def _use_inverse_folding(self, target):
        if Config.GENERATION_METHOD == "inverse":
            return True
        if Config.GENERATION_METHOD == "conservation":
            return False
        # Auto mode:
        return len(target) > Config.AUTO_SWITCH_LENGTH
    
    def _generate_with_conservation(self, target, num_variants):
        seed = self.inverse_folder.get_seed_sequence(target)
        return self.mutator.generate(seed, num_variants, target)