import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import gridspec
import matplotlib.colors as mcolors
import RNA 
from pathlib import Path
from tool.config import Config

class ConservationVisualizer:
    """
    Handles all visualization for conservation analysis
    """
    
    @staticmethod
    def plot_heatmap(sequence, conservation_scores, filename=None):
        """
        Args:
            sequence: str - RNA sequence (e.g., "GTGAACTG...")
            conservation_scores: np.array - Your calculated scores
            filename: Optional output path
        Returns:
            matplotlib Figure object
        """
        fig = plt.figure(figsize=(len(sequence)*0.35, 2.5), dpi=300)
        gs = gridspec.GridSpec(1, 2, width_ratios=[25,1], wspace=0.05)
        
        # Main heatmap
        ax = plt.subplot(gs[0])
        sns.heatmap(
            [conservation_scores],
            cmap="plasma_r",
            cbar=False,
            xticklabels=list(range(1, len(sequence)+1)),
            yticklabels=[""],
            linewidths=0.3,
            linecolor='gray',
            ax=ax,
            vmin=0,
            vmax=1,
            square=True,
            cbar_kws={"shrink": 0.8, "aspect": 20}
        )
        ax.tick_params(axis='x', rotation=0, labelsize=8)
        ax.tick_params(axis='y', rotation=0, labelsize=8)
        ax.set_xlabel("Nucleotide Position", fontsize=10)
        ax.set_title("Nucleotide Conservation (Dark = Conserved)", fontsize=11, pad=10)
        
        # Colorbar
        ax_cb = plt.subplot(gs[1])
        norm = plt.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap="plasma_r", norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, cax=ax_cb)
        cbar.set_label("Conservation Score", fontsize=9)
        cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])
        cbar.set_ticklabels(['0', '0.25', '0.5', '0.75', '1'])
        cbar.ax.tick_params(labelsize=8)
        
        # Save or return
        if filename:
            output_path = Path(Config.OUTPUT_DIR) / filename
            plt.savefig(output_path.with_suffix('.png'), dpi=400, bbox_inches='tight', pad_inches=0.02)
            plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', pad_inches=0.02)
            plt.close()
            return str(output_path)
        return fig

    @staticmethod
    def plot_tolerance_matrix(sequence, tolerance_matrix, filename=None):
        """
        Visualizes position-specific mutation tolerance
        Args:
            tolerance_matrix: np.array from conservation analysis
        """
        fig, ax = plt.subplots(figsize=(len(sequence)*0.3, 3), dpi=300)
        
        # Convert to dataframe for better labeling
        positions = range(1, len(sequence)+1)
        bases = ['A','T','C','G']
        
        # Create annotated heatmap
        sns.heatmap(
            tolerance_matrix.T,
            annot=True,
            fmt=".2f",
            cmap="YlGnBu",
            xticklabels=positions,
            yticklabels=[b for b in bases if b != sequence[0]],
            ax=ax,
            square=True,
            cbar_kws={"shrink": 0.8, "aspect": 20}
        )
        ax.tick_params(axis='x', rotation=90, labelsize=8)
        ax.tick_params(axis='y', rotation=0, labelsize=8)
        ax.set_title("Position-Specific Mutation Tolerance")
        ax.set_xlabel("Nucleotide Position")
        ax.set_ylabel("Alternative Base")
        
        if filename:
            output_path = Path(Config.OUTPUT_DIR) / filename
            plt.savefig(output_path, dpi=400, bbox_inches='tight', pad_inches=0.02)
            plt.close()
            return str(output_path)
        return fig
    


    @staticmethod
    def plot_secondary_structure(sequence, structure, conservation_scores, filename=None):
        """
        Plot secondary structure with bases colored by conservation score.

        Args:
            sequence (str): RNA sequence.
            structure (str): Dot-bracket notation of RNA structure (ViennaRNA output).
            conservation_scores (np.ndarray): Scores (0-1) per base.
            filename (str, optional): If given, saves PNG and PDF.

        Integration in your package:
            - Generate structure via ViennaRNA:
                fc = RNA.fold_compound(sequence)
                structure, mfe = fc.mfe()
            - Call this function with sequence, structure, and scores.

        """
        # Generate coordinates using ViennaRNA's naview layout
        coords = RNA.naview_xy_coordinates(structure)

        fig, ax = plt.subplots(figsize=(8, 8), dpi=400)
        norm = mcolors.Normalize(vmin=0, vmax=1)
        # Use a slightly desaturated plasma_r colormap by blending with white
        base_cmap = plt.cm.plasma_r
        colors = base_cmap(np.linspace(0, 1, 256))
        white = np.array([1,1,1,1])
        desaturated_colors = colors * 0.85 + white * 0.15
        desaturated_cmap = mcolors.ListedColormap(desaturated_colors)

        # Draw bonds between paired bases
        pt = RNA.ptable(structure)
        pairs = [(i, pt[i]) for i in range(1, len(pt)) if pt[i] > i]
        for i, j in pairs:
            ax.plot(
                [coords[i-1].X, coords[j-1].X],
                [coords[i-1].Y, coords[j-1].Y],
                color='#A0A0A0', lw=0.7, zorder=1
            )

        # Draw backbone as a sequential line connecting each nucleotide (drawn underneath bonds and nucleotides)
        ax.plot(
            [coords[i].X for i in range(len(sequence))],
            [coords[i].Y for i in range(len(sequence))],
            color='black', lw=0.5, zorder=0
        )

        # Draw nucleotides
        for i, base in enumerate(sequence):
            ax.scatter(
                coords[i].X, coords[i].Y,
                c=[desaturated_cmap(norm(conservation_scores[i]))],
                s=380, edgecolor='black', linewidth=0.6, zorder=2
            )
            ax.text(coords[i].X, coords[i].Y, base,
                    ha='center', va='center', fontsize=11, fontweight='bold', zorder=3)

        ax.set_aspect('equal')
        ax.axis('off')

        sm = plt.cm.ScalarMappable(cmap=desaturated_cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Conservation Score", fontsize=10)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['0', '1'])
        cbar.outline.set_edgecolor('black')
        cbar.outline.set_linewidth(0.8)

        if filename:
            output_path = Path(Config.OUTPUT_DIR) / filename
            plt.savefig(output_path.with_suffix('.png'), dpi=400, bbox_inches='tight', pad_inches=0.01, facecolor='white')
            plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', pad_inches=0.01, facecolor='white')
            plt.close()
            return str(output_path)

        return fig