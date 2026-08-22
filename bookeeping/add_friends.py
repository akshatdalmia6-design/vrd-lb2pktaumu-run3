import sys
sys.path.append("../hadronic/kin_constraints")
from apply_LbFD_tauM_constraints import add_massreco_friend
from run3_samples_v1r5791 import samples
from itertools import product

from pathlib import Path
import glob

# Define friend configurations in a registry
FRIEND_CONFIGS = {
    "recomass": {
        "func": add_massreco_friend,
        "kwargs": {
            "isRSim": False,
            "selection": None,
            "pKmu_from_Lc": False,
            "is_MC": False,
            "keepBranches": False,
            "units_in_mum": False,
            "intermediatepart": "L1520",
            "is_run3": True
        }
    },
    # We can easily add more friends here later!
    # "bdt_score": {"func": add_bdt_friend, "kwargs": {...}},
}


def get_all_input_files(samples):
    """Generator to yield sample context and file path."""
    for sample_name, blocks in samples.items():
        if any(x in sample_name for x in ["b2oc", "psi", "rdlow"]):
            continue
        for block_path in blocks.values():
            for file_path in glob.glob(block_path):
                yield sample_name, file_path


# Main processing loop
saving_dir = Path("  ")
tree_names = [
    "Hlt2RD_LbToPKTauMu_TauTo3Pi_OS_DTF3",
    "Hlt2RD_LbToPKTauMu_TauTo3Pi_SS_DTF3"
]

for sample_name, in_file in get_all_input_files(samples):
    in_file_stem = Path(in_file).name.split('.')[0]

    for tree_name in tree_names:
        # Apply configurations
        for friend_type, config in FRIEND_CONFIGS.items():
            out_file = saving_dir / f"{friend_type}_{tree_name}_{in_file_stem}.root"

            print(f"Processing {in_file} -> {out_file}")

            config["func"](in_file, f"{tree_name}/DecayTree", str(out_file),
                           friend_type, **config["kwargs"])
