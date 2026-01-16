## General helper functions: define input samples
from pathlib import Path

# Read the input directory and automatically discover
# input files with '.fastq.gz' extension.
INPUT_DIR = Path(config["input_directory"])
INPUT_FILES = list(INPUT_DIR.glob("*.fastq.gz"))
SAMPLES = [file.stem.replace(".fastq", "") for file in INPUT_FILES]
