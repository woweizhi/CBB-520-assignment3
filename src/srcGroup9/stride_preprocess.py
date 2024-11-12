import os
import subprocess

input_dir = "../../data/UP000002311_559292_YEAST_v4"
output_dir = "../../data/stride_output/"
stride_path = "/Users/mac/Downloads/stride/stride"  # specify the stride path in your own device
os.makedirs(output_dir, exist_ok=True)

for pdb_file in os.listdir(input_dir):
    if pdb_file.endswith(".pdb"):
        pdb_path = os.path.join(input_dir, pdb_file)
        output_file = os.path.join(output_dir, pdb_file[:-4] + '.txt')

        with open(output_file, 'w') as outfile:
            subprocess.run([stride_path, pdb_path], stdout=outfile)
        print(f"Processed: {pdb_file}")

print("All files processed!")