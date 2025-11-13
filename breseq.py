import subprocess
import os
import pandas as pd


def run_breseq(reference, output_dir, *fastq_files, polymorphism_prediction=True, fold_coverage=None):
    """
    Run breseq on given fastq files against a reference genome.
    
    Args:
        reference: Path to reference genome
        output_dir: Directory for breseq output
        *fastq_files: One or more FASTQ files. Can be strings or lists.
        polymorphism_prediction: Whether to enable polymorphism prediction
        fold_coverage: Coverage limit for reads
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Flatten fastq_files in case any elements are lists
    flattened_files = []
    for item in fastq_files:
        if isinstance(item, list):
            flattened_files.extend(item)
        else:
            flattened_files.append(item)
            
    cmd = [
        'breseq',
        '-r', reference,
    ] + flattened_files + [
        '--genbank-field-for-seq-id', 'ACCESSION',
        '--output', output_dir,
        '-j', '10'
    ]
    if polymorphism_prediction:
        cmd.append('--polymorphism-prediction')
    if fold_coverage:
        cmd += ['--limit-fold-coverage', str(fold_coverage)]

    print(f"Running breseq with command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=output_dir)

    # return output_dir

    # Add #=TITLE sample_name to genome diff file
    gd_file = os.path.join(output_dir, 'data', 'output.gd')
    # Read the existing content
    with open(gd_file, 'r') as f:
        content = f.readlines()
    # Add the title line at the beginning
    title_line = f'#=TITLE\t{os.path.basename(output_dir)}\n'
    content.insert(0, title_line)
    # Write back the modified content
    with open(gd_file, 'w') as f:
        f.writelines(content)

    return gd_file
    

def run_bam2cov(sample_dir, region):
    bam2cov_dir = os.path.join(sample_dir, 'BAM2COV')
    os.makedirs(bam2cov_dir, exist_ok=True)
    outfile = os.path.join(bam2cov_dir, region.replace(':', '_').replace('-', '_'))
    bam2cov_cmd = [
        'breseq', 'BAM2COV',
        '-b', os.path.join(sample_dir, 'data', 'reference.bam'),  # Default BAM output location
        '-o', outfile,
        '-r', region,
        '-t'  # table of coverage instead of plot
    ]

    subprocess.run(bam2cov_cmd, check=True, cwd=sample_dir)

    return outfile+'.tab'


def parse_region_average_cov_from_file(filepath):
    """
    Extract the 'region_average_cov' value from a coverage tab file.
    """
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#\tregion_average_cov"):
                # Split by tabs and take the last field
                return float(line.strip().split("\t")[-1])


def compare_gdiff(reference, out_file, gdiffs, format='HTML'):
    compare_cmd = [
        'gdtools', 'ANNOTATE',
        '-r', reference,
        '-o', out_file,
        '-f', format,
    ] + gdiffs

    subprocess.run(compare_cmd, check=True)

    return out_file


def count_mutations(path_to_gdiff, count_pols=True, output_file=None):

    if output_file is None:
        output_file = os.path.join(os.path.dirname(path_to_gdiff), 'count.csv')

    count_cmd = [
        'gdtools', 'COUNT',
        '-r', os.path.join(os.path.dirname(path_to_gdiff), 'reference.fasta'),
        '-o', output_file,
        path_to_gdiff

    ]

    if count_pols:
        count_cmd.append('-p')

    subprocess.run(count_cmd, check=True)

    count = pd.read_csv(output_file)['total'].iloc[0]

    return count


def apply_mutations_to_ref(reference, input_gd, ouput_ref):
    apply_cmd = [
        'gdtools', 'APPLY',
        '-o', ouput_ref,
        '-f', 'GENBANK',
        '-r', reference,
        input_gd
    ]

    subprocess.run(apply_cmd, check=True)

    return ouput_ref