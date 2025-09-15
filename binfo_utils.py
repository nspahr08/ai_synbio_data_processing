import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import os
import pandas as pd
import glob
import json


def count_reads(reads_file_path: str) -> int:
    """
    Counts the number of reads in a FASTQ or FASTA file (plain or gzipped).
    Automatically detects format by file extension.

    Args:
        reads_file_path (str): Path to the FASTQ or FASTA file.

    Returns:
        int: Total number of reads in the file, or -1 if an error occurs.
    """
    # Determine file type and compression
    file_lower = reads_file_path.lower()
    if file_lower.endswith('.gz'):
        open_func = lambda f: gzip.open(f, "rt")
        file_lower = file_lower[:-3]  # Remove .gz for type check
    else:
        open_func = lambda f: open(f, "r")

    # Detect format
    if file_lower.endswith('.fastq'):
        fmt = "fastq"
    elif file_lower.endswith('.fasta') or file_lower.endswith('.fa') or file_lower.endswith('.fna'):
        fmt = "fasta"
    else:
        print(f"Error: File extension not recognized for {reads_file_path}")
        return -1

    # Count records using SeqIO.parse
    try:
        with open_func(reads_file_path) as handle:
            return sum(1 for _ in SeqIO.parse(handle, fmt))
    except FileNotFoundError:
        print(f"Error: File not found at {reads_file_path}")
        return -1
    except Exception as e:
        print(f"An error occurred: {e}")
        return -1


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str(Seq(seq).reverse_complement())
    return complement


def parse_illumina_fastq_filename(filename: str) -> tuple:
    if '/' in filename:
        filename = os.path.basename(filename)
    parts = filename.split('_')
    if len(parts) < 4:
        raise ValueError("Filename does not contain enough parts to parse.")
    sample_name = parts[0]
    lane = parts[2]
    read = parts[3]
    return sample_name, lane, read


def create_manifest(folder: str, platform: str = 'illumina') -> pd.DataFrame:

    if platform == 'illumina':
        files = glob.glob(os.path.join(folder, '*.fastq*'))
        data = []
        for file in files:
            sample_name, _, read = parse_illumina_fastq_filename(file)
            if read == 'R1':
                data.append({
                    'sample_name': sample_name,
                    'fwd_fastq': file,
                    'rvs_fastq': None
                })
            elif read == 'R2':
                data.append({
                    'sample_name': sample_name,
                    'fwd_fastq': None,
                    'rvs_fastq': file
                })
        df = pd.DataFrame(data)
        manifest = df.groupby('sample_name', as_index=False).agg({
            'fwd_fastq': 'first',
            'rvs_fastq': 'first'
        })
        return manifest.sort_values('sample_name')

    elif platform == 'nanopore':
        data = []
        files = glob.glob(os.path.join(folder, '*.fastq*'))
        for file in files:
            sample_name = os.path.basename(file).split('.')[0].split('_')[0]
            data.append({'sample_name': sample_name, 'fastq': file}) 
        manifest = pd.DataFrame(data)

        return manifest.sort_values('sample_name')

    else:
        raise ValueError("Unsupported platform. Use 'illumina' or 'nanopore'.")


def write_df_to_fasta(df, seq_col, fasta_path, header_col=None):
    """
    Write DNA sequences from a DataFrame column to a fasta file.
    If header_col is provided, use its values as headers; otherwise use row index.
    """
    base_dir = os.path.dirname(fasta_path)
    if base_dir and not os.path.exists(base_dir):
        raise FileNotFoundError(
            f"Error: Directory '{base_dir}' does not exist for fasta output."
            )
    with open(fasta_path, 'w') as fasta_out:
        for idx, row in df.iterrows():
            header = str(row[header_col]) if header_col else str(idx)
            sequence = row[seq_col]
            fasta_out.write(f'>{header}\n{sequence}\n')

    return fasta_path


def read_fasta_into_df(fasta_path:str) -> pd.DataFrame:
    """
    Reads a fasta file and returns a DataFrame with columns 'seq_id' and 'sequence'.
    """
    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        records.append({'seq_id': record.id, 'sequence': str(record.seq)})
    return pd.DataFrame(records)


def json_to_dict(json_path: str) -> dict:
    """
    Reads a JSON file and returns its contents as a dictionary.
    """
    try:
        with open(json_path, 'r') as file:
            data = json.load(file)
        return data
    except FileNotFoundError:
        print(f"Error: {json_path} not found.")
        return {}

    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_path}.")
        return {}