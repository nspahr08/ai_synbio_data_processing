import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re
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


def parse_illumina_fastq_filename(filepath: str) -> tuple[str, str]:
    """
    Parse an Illumina fastq filename to extract sample name and read direction.
    Handles both gzipped (.fastq.gz) and non-gzipped (.fastq) files.
    
    Args:
        filepath (str): Path to the fastq file
            Examples: 
            - 'Data/Intensities/BaseCalls/SampleName_SampleNameContinued_S1_L001_R1_001.fastq.gz'
            - 'Data/Intensities/BaseCalls/SampleName_SampleNameContinued_S1_L001_R1_001.fastq'
    
    Returns:
        tuple[str, str]: (sample_name, read)
            Example: ('SampleName_SampleNameContinued', 'R1')
            
    Raises:
        ValueError: If the filename doesn't match expected Illumina format
    """
    
    # Get just the filename from the path
    filename = os.path.basename(filepath)
    
    # Pattern matches:
    # - Any characters up to _S\d+ (sample name)
    # - _S\d+ (sample number)
    # - _L\d{3} (lane number)
    # - _(R[12]) (read direction)
    # - _\d{3} (set number)
    # - \.fastq(\.gz)? (file extension, optional gz)
    pattern = r'^(.+)_S\d+_L\d{3}_(R[12])_\d{3}\w*.fastq(?:\.gz)?$'
    
    match = re.match(pattern, filename)
    if not match:
        raise ValueError(f"Filename {filename} doesn't match expected Illumina format")
    
    sample_name, read = match.groups()
    return sample_name, read


def create_manifest(folder: str, platform: str = 'illumina') -> pd.DataFrame:

    if platform == 'illumina':
        files = glob.glob(os.path.join(folder, '*.fastq*'))
        data = []
        for file in files:
            sample_name, read = parse_illumina_fastq_filename(file)
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

    elif platform == 'plasmidsaurus_hybrid':
        data = []
        files = glob.glob(os.path.join(folder, '*.fastq.gz'))
        sample_files = {}
        for file in files:
            basename = os.path.basename(file)
            sample_name = basename.replace('_illumina_R1.fastq.gz', '')\
                .replace('_illumina_R2.fastq.gz', '')\
                .replace('_nanopore.fastq.gz', '')\
                .replace('_illumina_R1_trimmed.fastq.gz', '')\
                .replace('_illumina_R2_trimmed.fastq.gz', '')\
                .replace('_nanopore_filtered.fastq.gz', '')
            if sample_name not in sample_files:
                sample_files[sample_name] = {'sample_name': sample_name}
            if file.endswith('_illumina_R1.fastq.gz') or file.endswith('_illumina_R1_trimmed.fastq.gz'):
                sample_files[sample_name]['fwd_fastq'] = file
            elif file.endswith('_illumina_R2.fastq.gz') or file.endswith('_illumina_R2_trimmed.fastq.gz'):
                sample_files[sample_name]['rvs_fastq'] = file
            elif file.endswith('_nanopore.fastq.gz') or file.endswith('_nanopore_filtered.fastq.gz'):
                sample_files[sample_name]['nanopore_fastq'] = file
        data = list(sample_files.values())
        manifest = pd.DataFrame(data)
        return manifest.sort_values('sample_name')
    
    else:
        raise ValueError("Unsupported platform. Use 'illumina', 'nanopore', or 'plasmidsaurus_hybrid'.")


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


def count_string_in_genbank(genbank_file: str, search_string: str) -> int:
    """
    Count occurrences of a string in a GenBank file, avoiding double-counting 
    between CDS and gene features that refer to the same gene.

    Args:
        genbank_file (str): Path to the GenBank file
        search_string (str): String to search for

    Returns:
        int: Number of unique occurrences of the string
    """
    count = 0
    counted_genes = set()  # Keep track of genes we've already counted

    for record in SeqIO.parse(genbank_file, 'genbank'):
        # First count occurrences in the sequence itself (case-insensitive)
        sequence = str(record.seq).upper()
        search_upper = search_string.upper()
        count += sequence.count(search_upper)

        # Process features
        for feature in record.features:
            # Skip if this is a gene feature and we've already counted its CDS
            if feature.type == 'gene':
                gene_id = feature.qualifiers.get('gene', [''])[0]
                if gene_id in counted_genes:
                    continue

            # Count occurrences in qualifiers
            for qualifiers in feature.qualifiers.values():
                for value in qualifiers:
                    # Handle string qualifiers only
                    if isinstance(value, str):
                        count += value.upper().count(search_upper)

            # If this is a CDS, add its gene to the counted set
            if feature.type == 'CDS':
                gene_id = feature.qualifiers.get('gene', [''])[0]
                if gene_id:
                    counted_genes.add(gene_id)

    return count


def convert_gbk_to_gff3(genbank_file, output_file, feature_type='gene'):
    with open(output_file, 'w') as out:
        # Write GFF3 header
        out.write('##gff-version 3\n')
        
        for record in SeqIO.parse(genbank_file, 'genbank'):
            for feature in record.features:
                if feature.type == feature_type:
                    strand = '+' if feature.strand == 1 else '-'
                    # GFF is 1-based
                    start = feature.location.start.position + 1
                    end = feature.location.end.position
                    
                    # Get feature ID, trying different qualifiers
                    feature_id = (
                        feature.qualifiers.get('locus_tag', [''])[0] or
                        feature.qualifiers.get('gene', [''])[0] or
                        feature.qualifiers.get('label', [''])[0] or
                        f"{record.id}_{start}_{end}"
                    )
                    
                    # Format GFF3 line
                    fields = [
                        record.id, '.', feature.type,
                        str(start), str(end), '.',
                        strand, '.', f"ID={feature_id}"
                    ]
                    gff_line = '\t'.join(fields) + '\n'
                    out.write(gff_line)
    
    return output_file
