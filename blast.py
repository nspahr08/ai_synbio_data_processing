import subprocess
import gzip
from binfo_utils import read_fasta_into_df, write_df_to_fasta, reverse_complement
import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional


def make_blast_db(fasta: str, blast_db: str, db_type: str = 'nucl') -> str:
    """Create a BLAST database from a FASTA file."""
    blast_db_dir = os.path.dirname(blast_db)
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Error: FASTA file '{fasta}' does not exist.")
    if not os.path.exists(blast_db_dir):
        raise FileNotFoundError(
            f"Error: Directory '{blast_db_dir}' does not exist for BLAST db."
            )

    cmd = [
        'makeblastdb',
        '-in', fasta,
        '-out', blast_db,
        # '-metadata_output_prefix', db_location,
        '-dbtype', db_type
    ]
    subprocess.run(cmd, check=True, cwd=blast_db_dir)
    return blast_db


def blast_reads(
        fastq, db, output, temp_dir='/storage/nspahr/tmp/',
         outfmt=None, num_threads=64, evalue=100
        ) -> tuple:
    if outfmt is None:
        outfmt = 'stitle sseqid qseqid pident nident qlen length mismatch gapopen qstart qend sstart send evalue'
    # From fastq files, create temp fasta file
    fasta = os.path.join(
        temp_dir, os.path.basename(fastq).replace('.fastq.gz', '.fasta')
        )
    with gzip.open(fastq, "rt") as fastq_in, open(fasta, "w") as fasta_out:
        SeqIO.convert(fastq_in, "fastq", fasta_out, "fasta")

    try:
        # Create blast results from fasta
        try:
            blast_cmd = [
                'blastn',
                '-db', db,
                '-query', fasta,
                '-out', output,
                '-num_threads', str(num_threads),
                '-evalue', str(evalue),
                '-outfmt', f'6 {outfmt}'
            ]
            print(' '.join(blast_cmd))
            subprocess.run(blast_cmd, check=True, cwd=os.path.dirname(output))
        except:
            os.remove(output)
            
        results = pd.read_csv(output, sep='\t', names=outfmt.split(' '))
    
        # Add read seq and overwrite tsv
        get_seqs = read_fasta_into_df(fasta)
        results = pd.merge(
            get_seqs, results, left_on='seq_id', right_on='qseqid', how='inner'
            )
        results.to_csv(output, sep='\t', index=False)
        print(f"BLAST results saved to {output}")
    
        # Delete temp fasta file
        os.remove(fasta)
    
        return (output, results)
    except:
        os.remove(fasta)
        raise


def get_insert_boundaries(query_positions: list) -> tuple:
    """Given a list of query positions, return the second and third positions in sorted order."""
    if len(query_positions) < 4:
        raise ValueError("List must contain exactly four positions.")
    ordered = sorted(query_positions)
    return ordered[1], ordered[2]


def calc_distance_between_flanks(query_positions: list) -> int:
    """Calculate the distance between the two flanking regions."""
    if len(query_positions) < 4:
        raise ValueError("List must contain exactly four positions.")
    ordered = sorted(query_positions)
    return ordered[2]-ordered[1]


def categorize_insertion_multiples(values, dec_to_nearest_int=0.1):
    categories = []
    for v in values:
        rounded = int(round(v))
        if abs(v - rounded) <= dec_to_nearest_int:
            categories.append(rounded)
        elif 0 < v < 1:
            categories.append('0<x<1')
        else:
            categories.append('other\nnon-int')
    return pd.Series(categories)


def process_flank_blast_results(blast_result_file: str, min_matching=None, hits_to_fasta=True) -> tuple:

    print(f"Processing {blast_result_file}")

    output_dir = os.path.dirname(blast_result_file)
    results_both_file = os.path.join(
        output_dir, str(os.path.basename(blast_result_file)).replace('.tsv', '_both_flanks.tsv'))
    results_one_file = os.path.join(
        output_dir, str(os.path.basename(blast_result_file)).replace('.tsv', '_one_flank.tsv'))

    # Check if file exists
    if not os.path.exists(blast_result_file) or os.path.getsize(blast_result_file) == 0:
        print(f"Warning: {blast_result_file} missing or empty.")
        return ((None, None), (None, None))

    blast_result_df = pd.read_csv(blast_result_file, sep='\t')

    # Calculate fraction of matching bases over length of the target/subject sequence.
    blast_result_df['pident/all'] = blast_result_df['nident']/300

    if min_matching is not None:
        # Only keep hits that have pident/all >= min_matching.
        blast_result_df = blast_result_df.loc[blast_result_df['pident/all'] >= min_matching]
        if len(blast_result_df) == 0:
            print(f"No valid hits in {blast_result_file} after filtering for min_matching={min_matching}.")
            return ((None, None), (None, None))

    # For each set of hits per query-subject, calculate the hit with the highest pident/all. Filter out all hits that don't have max_pident/all
    blast_result_df['max_pident/all'] = blast_result_df.groupby(['qseqid', 'sseqid'])['pident/all'].transform("max")
    blast_result_df = blast_result_df.loc[blast_result_df['pident/all'] == blast_result_df['max_pident/all']].drop_duplicates(['qseqid', 'sseqid'])

    if len(blast_result_df) == 0:
        print(f"No valid hits in {blast_result_file} after filtering.")
        return ((None, None), (None, None))

    # Split blast_result_df in two dataframes: those that have hits to both flanks, and those that have hits to only one flank
    flank_counts = blast_result_df.groupby('qseqid')['sseqid'].apply(lambda x: set(x))
    valid_qseqids = flank_counts[flank_counts.apply(lambda x: {'left_flank', 'right_flank'}.issubset(x))].index
    results_both = blast_result_df[blast_result_df['qseqid'].isin(valid_qseqids)]
    results_one = blast_result_df[~blast_result_df['qseqid'].isin(valid_qseqids)]

    if len(results_both) == 0:
        print(f"No reads with both flanks in {blast_result_file}.")
        results_both = None
    if len(results_one) == 0:
        print(f"No reads with one flank in {blast_result_file}.")
        results_one = None

    # Write dfs to file, write seqs to fasta (for manual inspection)
    if results_both is not None:
        results_both.to_csv(results_both_file, sep='\t', index=False)
        if hits_to_fasta:
            hits_fasta = results_both_file.replace('.tsv', '.fasta')
            write_df_to_fasta(results_both, 'sequence', hits_fasta, header_col='qseqid')
    else:
        results_both_file = None
    if results_one is not None:
        results_one.to_csv(results_one_file, sep='\t', index=False)
        if hits_to_fasta:
            hits_fasta = results_one_file.replace('.tsv', '.fasta')
            write_df_to_fasta(results_one, 'sequence', hits_fasta, header_col='qseqid')
    else:
        results_one_file = None

    return ((results_both_file, results_both), (results_one_file, results_one))


def count_features_in_blast_results(blast_result_file: str, features: dict):

    if not blast_result_file or not os.path.exists(blast_result_file):
        print(f"File does not exist.")
        return None

    results = pd.read_csv(blast_result_file, sep='\t')

    for name, seq in features.items():
        seq_rc = reverse_complement(seq)
        results[name] = results.apply(
            lambda x: x['sequence'].count(seq) + x['sequence'].count(seq_rc),
            axis=1
        )

    return results


def plot_feature_counts(counts_df: pd.DataFrame, features: dict,
                        out_pdf: Optional[str] = None,
                        greaterorequal: bool = True) -> None:

    if len(counts_df) == 0:
        print("No data to plot.")
        return None

    if not all(x in counts_df.columns for x in list(features.keys())):
        print("Error: Feature names not found in DataFrame columns")
        return None

    if greaterorequal:
        symbol = '>='
    else:
        symbol = ''

    # Create PDF with multiple pages, one for each feature
    if out_pdf:
        pdf = PdfPages(out_pdf)
    
    for tag in features.keys():
        feature_counts = counts_df[tag].value_counts().sort_index()
        feature_counts.index = feature_counts.index.map(str)
        feature_counts.index = [symbol + i for i in feature_counts.index]

        fig, ax = plt.subplots()
        plot = ax.bar(feature_counts.index, feature_counts.values)
        ax.set_ylabel(tag + ' count')
        ax.bar_label(plot, label_type='center')
        plt.show()

        if out_pdf:
            pdf.savefig(fig)
            plt.close(fig)

    # Close the PDF file if we created one
    if out_pdf:
        pdf.close()

    return None


def excise_insertion_seq(
        blast_result_file: str,
        expect_amplicon_len: int,
        out_pdf: Optional[str] = None):
    """
    Extract insertion sequences from blast results and analyze amplicon counts.
    
    Args:
        blast_result_file: Path to blast results file
        expect_amplicon_len: Expected length of one amplicon unit
        out_pdf: Path to save plot as PDF. If None, displays plot instead
    
    Returns:
        DataFrame with insertion sequences and amplicon counts
    """

    def calc_distance_between_flanks(query_positions: list) -> int:
        ordered = sorted(query_positions)
        return ordered[2]-ordered[1]

    results = pd.read_csv(blast_result_file, sep='\t')

    if len(results) == 0:
        print(f"No hits in {blast_result_file}.")
        return None

    # Pivot table: qseqid as rows, suffixed columns by left/right
    pivoted = results.pivot(index='qseqid', columns='sseqid')
    pivoted.columns = [f'{col[0]}_{col[1]}' for col in pivoted.columns]
    pivoted = pivoted.reset_index()

    pivoted[['first_boundary', 'second_boundary']] = pivoted.apply(
        lambda x: pd.Series(get_insert_boundaries(
            [
                x['qstart_left_flank'],
                x['qend_left_flank'],
                x['qstart_right_flank'],
                x['qend_right_flank']
            ]
            )
        ), axis=1
    )

    pivoted['dist_between_flanks'] = pivoted.apply(
        lambda x: calc_distance_between_flanks(
            [
                x['qstart_left_flank'],
                x['qend_left_flank'],
                x['qstart_right_flank'],
                x['qend_right_flank']
            ]
        ), axis=1)

    # Extract insertion sequence
    def get_insertion(row):
        start = row['first_boundary'] - 1
        end = row['second_boundary']
        return row['sequence_left_flank'][start:end]
    
    pivoted['insertion_seq'] = pivoted.apply(get_insertion, axis=1)

    # Calculate number of amplicons
    pivoted['n_amplicons'] = (
        pivoted['dist_between_flanks'] / expect_amplicon_len
    )

    # Plot distribution of n_amplicons
    category_series = categorize_insertion_multiples(pivoted['n_amplicons'])
    category_counts = category_series.value_counts()
    category_counts.index = category_counts.index.map(str)
    category_counts = category_counts.sort_index()
    # Split and sort indices
    int_indices = [i for i in category_counts.index if i.isdigit()]
    str_indices = [i for i in category_counts.index if not i.isdigit()]
    # Sort numerically and alphabetically
    sorted_indices = sorted(int_indices, key=int) + sorted(str_indices)
    # Reindex with sorted order
    category_counts = category_counts.reindex(sorted_indices)

    fig, ax = plt.subplots()
    plot = ax.bar(category_counts.index, category_counts.values)
    ax.set_ylabel('count')
    ax.bar_label(plot, label_type='center')
    ax.set_title(os.path.basename(blast_result_file))

    if out_pdf:
        with PdfPages(out_pdf) as pdf:
            plt.show()
            pdf.savefig(fig)
            plt.close(fig)
    else:
        plt.show()

    return pivoted
