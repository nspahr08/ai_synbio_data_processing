import os
import subprocess


def run_fastqc(path_to_file, threads=10):
    output_dir = os.path.dirname(path_to_file)
    fastqc_cmd = [
        'fastqc', path_to_file,
        '-t', str(threads)
    ]

    print(f"Running FastQC on {path_to_file}")
    subprocess.run(fastqc_cmd, check=True, cwd=output_dir)


def run_fastp(path_to_fwd, path_to_rev, path_to_fwd_out, path_to_rev_out, threads=10):
    # Make sure that base file name for R1 and R2 are the same
    # otherwise, use R1 base name, but print a warning
    output_name_fwd = path_to_fwd_out.split('_R1')[0]
    output_name_rev = path_to_rev_out.split('_R2')[0]
    if output_name_fwd != output_name_rev:
        print("Warning: R1 and R2 file names do not match. Using R1 base name for output.")
    output_dir = os.path.dirname(output_name_fwd)
    
    fastp_cmd = [
        'fastp',
        '-i', path_to_fwd,
        '-I', path_to_rev,
        '-o', path_to_fwd_out,
        '-O', path_to_rev_out,
        '-w', str(threads),
        '-j', os.path.join(output_dir, output_name_fwd + '_fastp.json'),
        '-h', os.path.join(output_dir, output_name_fwd + '_fastp.html')
    ]

    print(f"Running fastp on {path_to_fwd} and {path_to_rev}")
    subprocess.run(fastp_cmd, check=True, cwd=output_dir)


def run_filtlong(input_fastq, output_dir, min_length=1000, keep_percent=90):
    outfile = os.path.basename(input_fastq).replace('.fastq.gz', '_filtered.fastq.gz')
    outfile_path = os.path.join(output_dir, outfile)
    log_file = os.path.join(output_dir, os.path.basename(input_fastq).replace('.fastq.gz', '_filtlong.log'))

    filtlong_cmd = [
        'filtlong',
        '--min_length', str(min_length),
        '--keep_percent', str(keep_percent),
        input_fastq
    ]

    print(f"Running filtlong on {input_fastq}, output will be saved to {output_dir}")

    with open(outfile_path, 'wb') as out_f, open(log_file, 'w') as log_f:
        filtlong_proc = subprocess.Popen(
            filtlong_cmd,
            stdout=subprocess.PIPE,
            stderr=log_f,
            cwd=output_dir
            )
        gzip_proc = subprocess.Popen(['gzip'], stdin=filtlong_proc.stdout, stdout=out_f)
        filtlong_proc.stdout.close()  # Allow filtlong_proc to receive a SIGPIPE if gzip_proc exits.
        gzip_proc.communicate()

    return outfile_path


def run_nanocomp(files, sample_names, output_dir, threads=10):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    nanocomp_cmd = [
        'NanoComp',
        '-t', str(threads),
        '--fastq'
        ] + files + ['--names'] + sample_names

    subprocess.run(nanocomp_cmd, check=True, cwd=output_dir)


def run_multiqc(folder, output_dir=None):
    multiqc_cmd = [
        'multiqc',
        folder,
        '--force'
    ]
    output_dir = folder if output_dir is None else output_dir
    subprocess.run(multiqc_cmd, cwd=output_dir)

    return os.path.join(folder, 'multiqc_report.html')