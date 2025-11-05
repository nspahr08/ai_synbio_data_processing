import os
import subprocess
import tempfile
from Bio import SeqIO
import shutil
from binfo_utils import convert_gbk_to_gff3


def map_reads(fwd_reads, rvs_reads, reference, output_dir, keep_index=False, index_dir=None, sample_name=None):
    temp_dir = tempfile.mkdtemp()
    if keep_index:
        bt2_index = index_reference(reference, index_dir=index_dir)
    else:
        bt2_index = index_reference(reference, index_dir=temp_dir)
    sam_file = map_to_indexed_ref(fwd_reads, rvs_reads, bt2_index, output_dir, sample_name=sample_name)
    bam_file = sort_and_index_bam(sam_file)
    if temp_dir and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    return bam_file


def index_reference(reference: str, index_dir):
    print(f"Building bowtie2 index in {index_dir} from reference {reference}")
    # Determine file type
    ref_ext = os.path.splitext(reference)[1].lower()
    fasta_path = reference
    if ref_ext in ['.gb', '.gbk', '.genbank']:
        fasta_path = os.path.join(index_dir, 'reference.fasta')
        SeqIO.convert(reference, 'genbank', fasta_path, 'fasta')
    index = os.path.join(index_dir, 'reference')
    subprocess.run([
        'bowtie2-build', fasta_path, index
    ], check=True)
    return index


def map_to_indexed_ref(
        fwd_reads, rvs_reads, bowtie2_index, output_dir, sample_name='output'):
    """
    Map reads to an indexed reference using bowtie2.
    
    Args:
        fwd_reads: Forward reads file path(s), single or list
        rvs_reads: Reverse reads file path(s), single or list
        bowtie2_index: Path to bowtie2 index (prefix only)
        output_dir: Directory for output files
        sample_name: Prefix for output files
    """
    # Convert inputs to lists if they're strings
    if isinstance(fwd_reads, str):
        fwd_reads = [fwd_reads]
    if isinstance(rvs_reads, str):
        rvs_reads = [rvs_reads]
        
    # Ensure equal number of forward and reverse read files
    if len(fwd_reads) != len(rvs_reads):
        raise ValueError("Number of forward and reverse read files must match")
        
    print(f"Mapping {fwd_reads} and {rvs_reads} to reference {bowtie2_index}")
    output_sam = os.path.join(output_dir, sample_name + '.sam')
    log_file = os.path.join(output_dir, sample_name + '_bowtie2.log')
    
    # Build command with multiple input files
    cmd = [
        'bowtie2',
        '-x', bowtie2_index,
        '-1', ','.join(fwd_reads),  # bowtie2 accepts comma-separated lists
        '-2', ','.join(rvs_reads),
        '-S', output_sam
    ]
    
    with open(log_file, 'w') as log_fh:
        subprocess.run(cmd, check=True, stderr=log_fh)
    return output_sam


def count_mapped_reads(bam_file):
    cmd = [
        'samtools', 'view',
        '-c',  # count reads in bam file
        '-F', '4',  # filter out reads with 0x0004 bit flag set, i.e. unmapped reads
        bam_file
    ]
    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,  # captures both stdout and stderr
        text=True  # converts output to string instead of bytes
    )
    # stdout contains the count as a string
    return int(result.stdout.strip())


def parse_bowtie2_log(log_file):
    """
    Parse bowtie2 log file to extract the overall alignment rate.
    
    Args:
        log_file: Path to the bowtie2 log file
        
    Returns:
        float: The overall alignment rate as a percentage, or None if not found
    """
    try:
        with open(log_file, 'r') as f:
            for line in f:
                if "overall alignment rate" in line:
                    # Extract percentage from alignment rate line
                    return float(line.split('%')[0].strip())
    except (FileNotFoundError, ValueError, IndexError):
        return None
    return None


def sort_and_index_bam(sam_file):
    print(f"Sorting and indexing BAM file from {sam_file}")
    output_bam = os.path.splitext(sam_file)[0] + '.sorted.bam'
    # Convert SAM to BAM, sort
    cmd_sort = [
        'samtools', 'sort', '-o', output_bam, sam_file
    ]
    subprocess.run(cmd_sort, check=True)
    # Index BAM
    cmd_index = [
        'samtools', 'index', output_bam
    ]
    subprocess.run(cmd_index, check=True)
    return output_bam


def run_fadu(bam_file, reference_file, output_dir, fadu_folder=None):
    
    temp_dir = None

    if fadu_folder is None:
        fadu_folder = '/Users/nataschaspahr/code/FADU/'

    reference_file_ext = os.path.splitext(reference_file)[1].lower()

    if reference_file_ext in ['.gb', '.gbk', '.genbank']:
        temp_dir = tempfile.mkdtemp()
        gff3 = convert_gbk_to_gff3(reference_file, temp_dir+'/reference.gff3')

    elif reference_file_ext not in ['.gff3', '.gff']:
        msg = ("Reference file must be in GenBank (.gb, .gbk) "
               "or GFF3 (.gff3) format.")
        raise ValueError(msg)

    else:
        gff3 = reference_file

    fadu_cmd = [
        'julia',
        f'--project={fadu_folder}/fadu_pkgs',
        f'{fadu_folder}/fadu.jl',
        '-M',  # remove multimapping reads
        '-g', gff3,
        '-b', bam_file,
        '-o', output_dir,
        '-f', 'gene',
        '-a', 'ID'  # feature name tag in gff3 to use
    ]
    
    print(f"Running FADU with command: {' '.join(fadu_cmd)}")
    subprocess.run(fadu_cmd, check=True)

    if temp_dir is not None:
        shutil.rmtree(temp_dir)

    return output_dir


def map_and_feature_count(
        fwd_reads,
        rvs_reads,
        reference,
        output_dir,
        keep_index=False,
        index_dir=None,
        sample_name=None,
        fadu_folder=None):
    """
    Maps reads to a reference and performs feature counting using FADU.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_file = map_reads(
        fwd_reads,
        rvs_reads,
        reference,
        output_dir,
        keep_index=keep_index,
        index_dir=index_dir,
        sample_name=sample_name
    )
    
    output_dir = run_fadu(
        bam_file,
        reference,
        output_dir,
        fadu_folder=fadu_folder
    )

    return output_dir
