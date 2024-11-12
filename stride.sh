DIRECTORY="./data"
INPUT_DIR="${DIRECTORY}/UP000002311_559292_YEAST_v4"
OUTPUT_DIR="${DIRECTORY}/stride_output"

# make sure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Iterate over the files, calling stride each time
for FILE in "$INPUT_DIR"/*.pdb; do
    # We only need the file name
    BASENAME=$(basename "$FILE")
    # Call stride and create corresponding .out file
    # change the stride binary file directory to your own binary file address
    /Users/mac/Downloads/stride/stride -h "$FILE" > "${OUTPUT_DIR}/${BASENAME}.out"
done