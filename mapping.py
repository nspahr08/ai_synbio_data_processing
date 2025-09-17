import os
import subprocess
import tempfile
from Bio import SeqIO
import shutil


def map_reads(fwd_reads, rvs_reads, reference, output_dir, keep_index=False, index_dir=None):
    temp_dir = tempfile.mkdtemp()
    if keep_index:
        bt2_index = index_reference(reference, index_dir=index_dir)
    else:
        bt2_index = index_reference(reference, index_dir=temp_dir)
    sam_file = map_to_indexed_ref(fwd_reads, rvs_reads, bt2_index, output_dir)
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


def map_to_indexed_ref(fwd_reads, rvs_reads, bowtie2_index, output_dir):
    print(f"Mapping {fwd_reads} and {rvs_reads} to reference {bowtie2_index}")
    output_sam = os.path.join(output_dir, 'output.sam')
    # bowtie2 expects index prefix without file extension
    cmd = [
        'bowtie2',
        '-x', bowtie2_index,
        '-1', fwd_reads,
        '-2', rvs_reads,
        '-S', output_sam
    ]
    subprocess.run(cmd, check=True)
    return output_sam


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
