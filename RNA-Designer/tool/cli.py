import argparse
from tool.config import Config
from tool.main import main

def parse_args():
    parser = argparse.ArgumentParser(description="RNA Variant Design Toolkit CLI")
    parser.add_argument("--mode", choices=["auto", "inverse", "conservation"], default=Config.GENERATION_METHOD,
                        help="Variant generation mode")
    parser.add_argument("--n_runs", type=int, default=Config.N_RUNS, help="Number of runs")
    parser.add_argument("--output", type=str, default="output/variants.fasta", help="Output FASTA file")
    parser.add_argument("--group_size", type=int, default=Config.GROUP_SIZE, help="Group size")
    parser.add_argument("--n_groups", type=int, default=Config.N_GROUPS, help="Number of groups")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    return parser.parse_args()

def run_cli():
    args = parse_args()
    Config.GENERATION_METHOD = args.mode
    Config.N_RUNS = args.n_runs
    Config.GROUP_SIZE = args.group_size
    Config.N_GROUPS = args.n_groups
    Config.set_seed(args.seed)
    main()

if __name__ == "__main__":
    run_cli()