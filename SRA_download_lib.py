"""
SRA_download_lib.py
Contains functions for downloading data from the Sequence Read Archive (SRA).
"""

import os
import re


def get_GEO_info(GEO_id: str) -> dict:
    """
    Retrieve GEO information for the given GEO ID and return it as a dictionary.
    Requires the GEOparse package to be installed.
    """
    import GEOparse

    # Geoparse searchs for the accession number and generates a GSE object containing all the information.
    gse = GEOparse.get_GEO(GEO_id)
    return gse


def show_GSM_info(gse):
    """
    Showas GSM (sample) info.

    Args:
        gse: The GSE object containing GSM information.

    Returns:
        None
    """
    print()
    print("GSM (sample) info:")
    for gsm_name, gsm in gse.gsms.items():
        print("Name: ", gsm_name)
        print(
            "Metadata:",
        )
        for key, value in gsm.metadata.items():
            print(" - %s : %s" % (key, ", ".join(value)))
        print(
            "Table data:",
        )
        print(gsm.table.head())


def show_GPL_info(gse):
    """
    Function to display GPL (project) information.

    Args:
        gse: the GSE object containing GPL information

    Returns:
        None
    """
    print()
    print("GPL (project) info:")
    for gpl_name, gpl in gse.gpls.items():
        print("Name: ", gpl_name)
        print(
            "Metadata:",
        )
        for key, value in gpl.metadata.items():
            print(" - %s : %s" % (key, ", ".join(value)))
        print(
            "Table data:",
        )
        print(gpl.table.head())
        break


def list_to_dict(list: list, split: str = ":") -> dict:
    """
    Convert a list of strings into a dictionary by splitting each string at a specified character.

    Args:
        list: list of strings to convert into key-value pairs
        split: character to split each string in the list (default is ':')

    Returns:
        dict: dictionary containing key-value pairs from the input list
    """
    dict = {key: value for key, value in [pair.split(split) for pair in list]}
    return dict


def gsm_to_srr(gsm: str) -> str:
    """
    Convert a Gene Expression Omnibus (GEO) Sample accession (GSM) to Sequence Read Archive (SRA) run accession (SRR).
    Requires the pysradb package to be installed.

    Args:
        gsm (str): The Gene Expression Omnibus Sample accession (GSM) to be converted.

    Returns:
        str: The Sequence Read Archive run accession (SRR) corresponding to the input GSM.
    """

    from pysradb import SRAweb

    sra_db = SRAweb()
    sample_metadata = sra_db.sra_metadata(gsm)
    srr = sample_metadata.at[0, "run_accession"]
    return srr


def gsm_to_srp(gsm: str) -> str:
    """
    Convert a Gene Expression Omnibus (GEO) Sample accession (GSM) to Sequence Read Archive (SRA) project accession (SRP).
    Requires the pysradb package to be installed.

    Args:
        gsm (str): The Gene Expression Omnibus Sample accession (GSM) to be converted.

    Returns:
        srp: The Sequence Read Archive run accession (SRP) corresponding to the input GSM.
    """

    from pysradb import SRAweb

    sra_db = SRAweb()
    sample_metadata = sra_db.sra_metadata(gsm)
    srp = sample_metadata.at[0, "study_accession"]
    return srp


def download_fastq(
    out_dir: str,
    sample_ids: str | list = None,
    download_methods: str = "ena-ascp ena-ftp aws-http prefetch",
    file=None,
) -> None:
    """
    Downloads fastq files for the given sample ids using different download methods.
    Requires kingfisher and Aspera Client to be installed.

    Args:
        sample_ids (str | list): The ids of the samples to download. Can be a single id or a list of ids.
        out_dir (str, optional): The directory to save the downloaded files. Defaults to ".".
        download_methods (str, optional): The methods to use for downloading the files. Defaults to "ena-ascp ena-ftp aws-http prefetch".
        file = path to .txt file with sample ids. Each sample id should be on a separate line. If file, sample_list will be ignored.

    Returns:
        None
    """
    import subprocess

    if file:
        with open(file, "r") as file:
            sample_ids = [line.strip() for line in file.readlines()]

    if sample_ids.__class__ == list:
        sample_ids = " ".join(sample_ids)

    # In order to use current directory, --output-directory is not required.
    out_dir_parser = f"--output-directory {out_dir}" if out_dir else ""

    # Calls kingfisher to download fastq files.
    subprocess.call(
        f"kingfisher get --run-identifiers {sample_ids} {out_dir_parser} --download-methods {download_methods}",
        shell=True,
    )


def download_fastq_parallel(
    out_dir: str,
    sample_list: list = None,
    download_methods: str = "ena-ascp ena-ftp aws-http prefetch",
    processes: int = 10,
    file=None,
    use_max_processes=False,
) -> None:
    """
    Download fastq files for a list of samples in parallel using multiple download methods.
    Requires joblib, kingfisher and Aspera Client  to be installed.

    Parameters:
        sample_list (list): List of sample IDs to download fastq files for.
        out_dir (str): Output directory to save the downloaded files. Default is the current directory.
        download_methods (str): Methods to download the files (e.g., ena-ascp, ena-ftp, aws-http, prefetch).
        processes (int): Number of processes to use for parallel downloading. Default is 10.
        file = path to .txt file with sample ids. Each sample id should be on a separate line. If file, sample_list will be ignored.
        use_max_processes (bool): Flag to indicate whether to use the maximum available processes.

    Returns:
        None
    """

    if file:
        with open(file, "r") as file:
            sample_list = [line.strip() for line in file.readlines()]

    from joblib import Parallel, delayed

    if use_max_processes:
        processes = len(sample_list)

    Parallel(n_jobs=processes)(
        delayed(download_fastq)(out_dir, sample_id, download_methods)
        for sample_id in sample_list
    )


def download_GEO_dataset(
    GEO_id: str,
    out_dir: str,
    download_methods: str = "ena-ascp ena-ftp aws-http prefetch",
    processes: int = 10,
    use_max_processes=False,
) -> None:
    """
    Downloads a GEO dataset using the specified GEO ID.
    Requires: GEOparse, pysradb, joblib, Aspera Client and kingfisher.

    Args:
        GEO_id (str): The GEO ID of the dataset to be downloaded.
        out_dir (str, optional): The output directory for the downloaded files. Defaults to ".".
        download_methods (str, optional): The methods to be used for downloading. Defaults to "ena-ascp ena-ftp aws-http prefetch".
        processes (int, optional): The number of processes to be used for downloading. Defaults to 8.
        use_max_processes (bool): Flag to indicate whether to use the maximum available processes.

    Returns:
        None
    """

    gse = get_GEO_info(GEO_id)
    sample_ids = [gsm_to_srr(sample_id) for sample_id in gse.gsms.keys()]

    if use_max_processes:
        processes = len(sample_ids)

    download_fastq_parallel(
        sample_list=sample_ids,
        out_dir=out_dir,
        download_methods=download_methods,
        processes=processes,
    )
