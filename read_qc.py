import os
import subprocess


def run_fastqc(path_to_file):
    output_dir = os.path.dirname(path_to_file)
    fastqc_cmd = [
        'fastqc', path_to_file,
    ]

    print(f"Running FastQC on {path_to_file}")
    subprocess.run(fastqc_cmd, check=True, cwd=output_dir)


def run_filtlong(input_fastq, output_dir, min_length=1000, keep_percent=90):
    outfile = os.path.basename(input_fastq).replace('.fastq.gz', '_filtered.fastq.gz')
    outfile_path = os.path.join(output_dir, outfile)

    filtlong_cmd = [
        'filtlong',
        '--min_length', str(min_length),
        '--keep_percent', str(keep_percent),
        input_fastq
    ]

    print(f"Running filtlong on {input_fastq}, output will be saved to {output_dir}")

    with open(outfile_path, 'wb') as out_f:
        filtlong_proc = subprocess.Popen(filtlong_cmd, stdout=subprocess.PIPE, cwd=output_dir)
        gzip_proc = subprocess.Popen(['gzip'], stdin=filtlong_proc.stdout, stdout=out_f)
        filtlong_proc.stdout.close()  # Allow filtlong_proc to receive a SIGPIPE if gzip_proc exits.
        gzip_proc.communicate()

    return outfile_path


def run_nanocomp(files, sample_names, output_dir, threads=10):

    nanocomp_cmd = [
        'NanoComp',
        '-t', str(threads),
        '--fastq'
        ] + files + ['--names'] + sample_names

    subprocess.run(nanocomp_cmd, check=True, cwd=output_dir)
