"""
RNA Variant Design Toolkit

Core Submodules:
- core: Basic sequence analysis operations
- evolution: Variant generation strategies  
- analysis: Conservation and tolerance profiling
- visualization: Plotting and data presentation
"""

# Core exports
from tool.config import Config
from tool.main import main

# Submodule exports
from tool.evolution import VariantGenerator
from tool.analysis import ConservationAnalyzer

__all__ = [
    'Config',
    'main',
    'VariantGenerator',
    'ConservationAnalyzer'
]

# Package metadata
__version__ = "1.1.0"
__author__ = "Konstantinos Kariotis"
__email__ = "kariotis.konst@gmail.com"
__description__ = "A modular toolkit for RNA variant design and analysis"
__url__ = "https://github.com/kkariotis/rna-variant-tool"