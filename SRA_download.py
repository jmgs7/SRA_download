"""
SRA_download.py
This script is an CLI for downloading fastq files from NCBI SRA using the SRA_download_lib.py library. 
Downloads fastq files for the list of samples in parallel using multiple download methods.
Requires joblib, kingfisher and Aspera Client  to be installed.
"""

from SRA_download_lib import download_fastq_parallel
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input_file")  # path to .txt file with sample ids
parser.add_argument("-O", "--output_dir")  # path to output directory
parser.add_argument(
    "-D", "--download_methods", default="ena-ascp ena-ftp aws-http prefetch"
)  # download methods (see Kingfisher documentation)
parser.add_argument(
    "-P", "--processes", type=int, default=10
)  # number of parallel downloads (recommended 10 as is the maximum simultaneuous request by NCBI)
parser.add_argument(
    "-M", "--use_max_processes", type=bool, default=False
)  # use maximum parallel downloads (occupy all cores)
args = vars(parser.parse_args())

file, out_dir, download_methods, processes, use_max_processes = (
    args["input_file"],
    args["output_dir"],
    args["download_methods"],
    args["processes"],
    args["use_max_processes"],
)

download_fastq_parallel(
    out_dir=out_dir,
    download_methods=download_methods,
    processes=processes,
    file=file,
    use_max_processes=use_max_processes,
)
