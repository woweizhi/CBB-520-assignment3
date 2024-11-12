from srcGroup10.util import *
import os

if __name__ == "__main__":
    input_dir = "./data"
    output_dir = "./results"
    output_prefix = "group10"

    if output_prefix:
        out_file_path = os.path.join(output_dir, output_prefix)
        if not os.path.isdir(out_file_path):
            os.mkdir(out_file_path)

    # Update the pdb_dir to your actual data directory
    pdb_dir = os.path.join(input_dir, "./test")
    output_file = os.path.join(out_file_path, "amino_acid_pair_analysis_v2.txt")

    # Optional: Validate input directory
    if not os.path.isdir(pdb_dir):
        logging.error(f"The specified PDB directory '{pdb_dir}' does not exist.")
    else:
        analyze_multiple_pdbs(
            pdb_dir=pdb_dir,
            output_file=output_file,
            randomizations=100,        # You can reduce this number for quicker testing
            sequence_gap=100,
            distance_threshold=10
        )
