import argparse
from timeit import default_timer as timer
from datetime import timedelta
from pathlib import Path, PurePath
import os
from lib.transform import make_reduce_file

def main():
    parser = argparse.ArgumentParser(description="Convert reduced file to plink")
    parser.add_argument("--matrix_folder", required=True, type=Path, help="Folder that has matrix-type SNP files.")
    parser.add_argument("--SNP_folder", required=True, type=Path, help="Folder that output reduced SNP files.")
    parser.add_argument("--skip_rows", type=int, default=9, help="How many rows in the matrix file before the IDs (default: 9)")
    parser.add_argument("--chunk_size", type=int, default=100000, help="How many snp will be transformed per chunk (default: 100000)")
    
    args = parser.parse_args()
    matrix_folder = PurePath(args.matrix_folder)
    reduce_folder = PurePath(args.SNP_folder)
    
    # Part 1
    def make_reduce_folder():
        # Generate output folder
        if not os.path.exists(reduce_folder):
            print("Part 1: Generate output folder")
            Path(reduce_folder).mkdir(parents=True, exist_ok=True)

    # Part 2
    def matrix_to_reduce():
        if os.path.exists(matrix_folder):
            make_reduce_file(matrix_folder, args.skip_rows, args.chunk_size, reduce_folder)

    # Main workflow
    make_reduce_folder()
    matrix_to_reduce()

if __name__ == "__main__":
    try:
        start = timer()
        main()
        end = timer()
        print("All tasks were done! This round spent " + str(timedelta(seconds=end-start)).split(".")[0] + ".")
    except Exception as e:
        print(e)