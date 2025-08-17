import os
from pathlib import Path

class FileHandler:
    @staticmethod
    def save_to_fasta(sequences, filename):
        with open(filename, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">scaffold_{i}\n{seq}\n")

    @staticmethod
    def save_groups_to_fasta(groups, filename):
        with open(filename, "w") as f:
            for g_num, group in enumerate(groups):
                for v_num, var in enumerate(group):
                    f.write(f">Group{g_num+1}_Var{v_num+1}\n{var}\n")