"""
Python module for managing breseq command-line parameters and execution.

This module provides two main classes:
- Breseq_params: Manages breseq parameter settings
- Breseq: Manages breseq execution for samples
"""

import os
import re
import json
import hashlib
import subprocess
import shlex
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
import csv
import zlib
import base64


def _get_ref_genomes_path() -> str:
    """Get REF_GENOMES path from environment variable or use default."""
    return os.environ.get('REF_GENOMES', '/storage/synbio/ai_synbio_data/reference_data/genomes')


def _convert_param_name_to_python(name: str) -> str:
    """Convert breseq parameter name to Python attribute name.
    
    Examples:
        --polymorphism-prediction -> polymorphism_prediction
        -r -> reference
        --read-min-length -> read_min_length
    """
    # Remove leading dashes
    name = name.lstrip('-')
    # Replace hyphens with underscores
    name = name.replace('-', '_')
    return name


def _convert_python_name_to_param(name: str) -> str:
    """Convert Python attribute name back to breseq parameter format.
    
    Examples:
        polymorphism_prediction -> --polymorphism-prediction
        reference -> -r
        read_min_length -> --read-min-length
    """
    # Special case for short flags
    if name == 'reference':
        return '-r'
    elif name == 'name':
        return '-n'
    elif name == 'num_processors':
        return '-j'
    elif name == 'output':
        return '-o'
    elif name == 'polymorphism_prediction':
        return '-p'
    elif name == 'nanopore':
        return '-x'
    elif name == 'contig_reference':
        return '-c'
    elif name == 'junction_only_reference':
        return '-s'
    elif name == 'targeted_sequencing':
        return '-t'
    elif name == 'header_genome_diff':
        return '-g'
    elif name == 'keep_intermediates':
        return '-k'
    else:
        # Convert underscores to hyphens and add double dash
        return '--' + name.replace('_', '-')


class Breseq_params:
    """Manages breseq command-line parameters.
    
    Required parameters:
        - reference: Reference genome file from REF_GENOMES folder
        - num_processors (int): Number of processors
        - polymorphism_prediction (bool): Whether to use polymorphism mode
    
    All other parameters are optional with defaults from breseq documentation.
    """
    
    # Required parameters
    REQUIRED_PARAMS = ['reference', 'num_processors', 'polymorphism_prediction']
    
    def __init__(self, reference: str,
                 num_processors: int,
                 polymorphism_prediction: bool,
                 **kwargs):
        """Initialize Breseq_params.
        
        Can be created from existing Breseq object or built by hand.
        If built by hand, all required parameters must be provided.
        
        Args:
            reference: Reference genome filename (must be in REF_GENOMES)
            num_processors: Number of processors to use
            polymorphism_prediction: Whether to use polymorphism mode (True/False)
            **kwargs: Any other breseq parameters
        """
        self._ref_genomes_path = _get_ref_genomes_path()

        # Set required parameters
        self.reference = reference
        self.num_processors = num_processors
        self.polymorphism_prediction = polymorphism_prediction
        
        # Set optional parameters with defaults
        self._set_defaults()
        
        # Set any additional parameters from kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        # Validate if all required params are set
        self._validate_required()
    
    def _set_defaults(self):
        """Set default values for all optional parameters."""
        # Basic options
        self.name = None  # -n, DEFAULT=<none>
        
        # Read file options
        self.limit_fold_coverage = 0  # -l, DEFAULT=OFF
        self.aligned_sam = False  # --aligned-sam
        self.read_min_length = 18  # --read-min-length, DEFAULT=18
        self.read_max_same_base_fraction = 0.9  # --read-max-same-base-fraction, DEFAULT=0.9
        self.read_max_N_fraction = 0.5  # --read-max-N-fraction, DEFAULT=0.5
        self.long_read_trigger_length = 1000  # --long-read-trigger-length, DEFAULT=1000
        self.long_read_split_length = 200  # --long-read-split-length, DEFAULT=200
        self.long_read_distribute_remainder = False  # --long-read-distribute-remainder
        self.genbank_field_for_seq_id = None  # --genbank-field-for-seq-id, DEFAULT=AUTOMATIC
        
        # Reference file options
        self.contig_reference = None  # -c, DEFAULT=0
        self.junction_only_reference = []  # -s, can be multiple, DEFAULT=0
        self.targeted_sequencing = False  # -t
        self.user_evidence_gd = None  # --user-evidence-gd
        
        # Read alignment options
        self.minimum_mapping_quality = 0  # -m, DEFAULT=0
        self.base_quality_cutoff = 3  # -b, DEFAULT=3
        self.quality_score_trim = False  # --quality-score-trim
        self.require_match_length = 0  # --require-match-length, DEFAULT=0
        self.require_match_fraction = 0.9  # --require-match-fraction, DEFAULT=0.9
        self.maximum_read_mismatches = None  # --maximum-read-mismatches, DEFAULT=OFF
        
        # Bowtie2 options
        self.bowtie2_scoring = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals"  # DEFAULT
        self.bowtie2_stage1 = "--local -i S,1,0.25 --score-min L,1,0.9 -k 2000"  # DEFAULT
        self.bowtie2_stage2 = "--local -i S,1,0.25 --score-min L,6,0.2 -k 2000"  # DEFAULT
        self.bowtie2_junction = "--local -i S,1,0.25 --score-min L,1,0.70 -k 2000"  # DEFAULT
        
        # Junction evidence options
        self.no_junction_prediction = False  # --no-junction-prediction
        self.junction_indel_split_length = 3  # --junction-indel-split-length, DEFAULT=3
        self.junction_alignment_pair_limit = 100000  # --junction-alignment-pair-limit, DEFAULT=100000
        self.junction_minimum_candidates = 100  # --junction-minimum-candidates, DEFAULT=100
        self.junction_maximum_candidates = 5000  # --junction-maximum-candidates, DEFAULT=5000
        self.junction_candidate_length_factor = 0.1  # --junction-candidate-length-factor, DEFAULT=0.1
        self.junction_minimum_candidate_pos_hash_score = 2  # DEFAULT=2
        self.junction_score_cutoff = 3.0  # --junction-score-cutoff, DEFAULT=3.0
        self.junction_minimum_pos_hash_score = None  # DEFAULT varies by mode
        self.junction_minimum_side_match = None  # DEFAULT varies by mode
        self.junction_minimum_pr_no_read_start_per_position = 0.1  # DEFAULT=0.1
        self.junction_allow_suboptimal_matches = False  # --junction-allow-suboptimal-matches
        
        # Missing coverage options
        self.deletion_coverage_seed_cutoff = 0  # DEFAULT=0
        self.deletion_coverage_propagation_cutoff = 0  # DEFAULT=0
        self.call_mutations_overlapping_MC = False  # --call-mutations-overlapping-MC
        
        # Consensus RA evidence options
        self.consensus_score_cutoff = 10  # DEFAULT=10
        self.consensus_frequency_cutoff = 0.8  # DEFAULT=0.8
        self.consensus_minimum_variant_coverage = 0  # DEFAULT=0
        self.consensus_minimum_total_coverage = 0  # DEFAULT=0
        self.consensus_minimum_variant_coverage_each_strand = 0  # DEFAULT=0
        self.consensus_minimum_total_coverage_each_strand = 0  # DEFAULT=0
        self.consensus_reject_indel_homopolymer_length = 0  # DEFAULT=OFF
        self.consensus_reject_surrounding_homopolymer_length = 0  # DEFAULT=OFF
        
        # Polymorphism RA evidence options
        self.polymorphism_score_cutoff = None  # DEFAULT varies by mode
        self.polymorphism_frequency_cutoff = None  # DEFAULT varies by mode
        self.polymorphism_minimum_variant_coverage = 0  # DEFAULT=0
        self.polymorphism_minimum_total_coverage = 0  # DEFAULT=0
        self.polymorphism_minimum_variant_coverage_each_strand = None  # DEFAULT varies by mode
        self.polymorphism_minimum_total_coverage_each_strand = 0  # DEFAULT=0
        self.polymorphism_bias_cutoff = 0  # DEFAULT=OFF
        self.polymorphism_no_indels = False  # --polymorphism-no-indels
        self.polymorphism_reject_indel_homopolymer_length = None  # DEFAULT varies by mode
        self.polymorphism_reject_surrounding_homopolymer_length = None  # DEFAULT varies by mode
        
        # Output options
        self.max_displayed_reads = 100  # DEFAULT=100
        self.brief_html_output = False  # --brief-html-output
        self.header_genome_diff = None  # -g
        self.no_javascript = False  # --no-javascript
        
        # Pipeline control options
        self.skip_RA_MC_prediction = False  # --skip-RA-MC-prediction
        self.skip_JC_prediction = False  # --skip-JC-prediction
        self.skip_MC_prediction = False  # --skip-MC-prediction
        
        # Debugging options
        self.keep_intermediates = False  # -k
        self.per_position_file = False  # --per-position-file
        self.junction_debug = False  # --junction-debug
        
        # Experimental options
        self.cnv = False  # --cnv
        self.cnv_tile_size = 500  # DEFAULT=500
        self.cnv_ignore_redundant = False  # --cnv-ignore-redundant
        
        # Nanopore preset (affects multiple parameters)
        self.nanopore = False  # -x
    
    def _validate_required(self):
        """Validate that all required parameters are set."""
        missing = []
        for param in self.REQUIRED_PARAMS:
            if not hasattr(self, param) or getattr(self, param) is None:
                missing.append(param)
        
        if missing:
            raise ValueError(f"Missing required parameters: {', '.join(missing)}")
    
    def _validate_reference(self, ref_file: str):
        """Validate that reference file exist in REF_GENOMES folder."""
        ref_genomes_path = Path(self._ref_genomes_path)
        if not ref_genomes_path.exists():
            raise FileNotFoundError(f"REF_GENOMES directory does not exist: {ref_genomes_path}")
        
        ref_path = ref_genomes_path / ref_file
        if not ref_path.exists():
            raise FileNotFoundError(
                f"Reference file '{ref_file}' not found in REF_GENOMES directory: {ref_genomes_path}"
            )

    def __setattr__(self, name: str, value):
        """Intercept setting of `reference` to validate existence.
        """
        if name == 'reference' and value is not None:
            self._validate_reference(value)
        super().__setattr__(name, value)
    
    def _get_all_params_dict(self) -> Dict[str, Any]:
        """Get all parameters as a dictionary for hashing."""
        params = {}
        for key in dir(self):
            # Skip private/internal attributes
            if key.startswith('_'):
                continue

            # Skip class-level properties/descriptors (e.g., @property "version_name")
            cls_attr = getattr(self.__class__, key, None)
            if isinstance(cls_attr, property):
                continue

            # Safely get the attribute value; some descriptors may raise on access
            try:
                value = getattr(self, key)
            except Exception:
                # If accessing the attribute raises, skip it
                continue

            # Skip callables (methods)
            if callable(value):
                continue

            # Skip None values for optional params (they use defaults)
            if value is not None:
                params[key] = value
        return params
    
    def _generate_version_hash(self) -> str:
        """Generate a hash from all parameter values for version naming."""
        params_dict = self._get_all_params_dict()
        # Sort by key for consistent hashing
        sorted_params = sorted(params_dict.items())
        # Convert to string representation
        params_str = json.dumps(sorted_params, sort_keys=True)
        # Generate hash
        hash_obj = hashlib.sha256(params_str.encode())
        # Return first 10 characters of hex digest
        return hash_obj.hexdigest()[:10]
    
    @property
    def version_name(self) -> str:
        """Get version name based on parameter hash."""
        # Use a reversible, self-contained encoding of the params so the
        # version name can be decoded back to the original params.
        params_dict = self._get_all_params_dict()
        return _params_to_encoded_version_name(params_dict)


def _params_to_encoded_version_name(params: Dict[str, Any], prefix: str = 'breseq_params_') -> str:
    """Encode params dict into a compact, filesystem-safe version name.

    Method: deterministic JSON -> zlib.compress -> urlsafe base64 (no padding).
    The resulting string is safe for most filesystems (URL-safe base64 uses
    only letters, digits, '-', '_').
    """
    # Deterministic JSON representation
    j = json.dumps(params, sort_keys=True, separators=(',', ':')).encode('utf-8')
    # Compress to reduce length
    compressed = zlib.compress(j, level=6)
    # URL-safe base64 encode and strip padding
    b64 = base64.urlsafe_b64encode(compressed).decode('ascii').rstrip('=')
    return prefix + b64


def _encoded_version_name_to_params(version_name: str, prefix: str = 'breseq_params_') -> Dict[str, Any]:
    """Decode a version_name produced by _params_to_encoded_version_name back to params dict.

    Raises ValueError if the prefix is wrong or decoding fails.
    """
    if not version_name.startswith(prefix):
        raise ValueError(f"version_name does not start with expected prefix '{prefix}'")
    b64 = version_name[len(prefix):]
    # Restore padding for base64
    padding = '=' * ((4 - len(b64) % 4) % 4)
    try:
        compressed = base64.urlsafe_b64decode(b64 + padding)
        j = zlib.decompress(compressed)
        return json.loads(j.decode('utf-8'))
    except Exception as e:
        raise ValueError(f"Failed to decode version_name: {e}") from e
    
    def to_command_args(self) -> List[str]:
        """Convert parameters to breseq command-line arguments list.
        
        Returns:
            List of command-line arguments (without 'breseq' command itself)
        """
        args = []
        
        # Required parameters
        # -r reference (single)
        if getattr(self, 'reference', None):
            args.extend(['-r', os.path.join(self._ref_genomes_path, self.reference)])
        
        # -j num_processors
        args.extend(['-j', str(self.num_processors)])
        
        # -p polymorphism_prediction (flag)
        if self.polymorphism_prediction:
            args.append('-p')
        
        # Optional parameters (only include if not default)
        # -n name
        if self.name is not None and self.name:
            args.extend(['-n', self.name])
        
        # Read file options
        if self.limit_fold_coverage > 0:
            args.extend(['-l', str(self.limit_fold_coverage)])
        
        if self.aligned_sam:
            args.append('--aligned-sam')
        
        if self.read_min_length != 18:
            args.extend(['--read-min-length', str(self.read_min_length)])
        
        if self.read_max_same_base_fraction != 0.9:
            args.extend(['--read-max-same-base-fraction', str(self.read_max_same_base_fraction)])
        
        if self.read_max_N_fraction != 0.5:
            args.extend(['--read-max-N-fraction', str(self.read_max_N_fraction)])
        
        if self.long_read_trigger_length != 1000:
            args.extend(['--long-read-trigger-length', str(self.long_read_trigger_length)])
        
        if self.long_read_split_length != 200:
            args.extend(['--long-read-split-length', str(self.long_read_split_length)])
        
        if self.long_read_distribute_remainder:
            args.append('--long-read-distribute-remainder')
        
        if self.genbank_field_for_seq_id:
            args.extend(['--genbank-field-for-seq-id', self.genbank_field_for_seq_id])
        
        # Reference file options
        if self.contig_reference:
            args.extend(['-c', self.contig_reference])
        
        if self.junction_only_reference:
            for ref in self.junction_only_reference:
                args.extend(['-s', ref])
        
        if self.targeted_sequencing:
            args.append('-t')
        
        if self.user_evidence_gd:
            args.extend(['--user-evidence-gd', self.user_evidence_gd])
        
        # Read alignment options
        if self.minimum_mapping_quality != 0:
            args.extend(['-m', str(self.minimum_mapping_quality)])
        
        if self.base_quality_cutoff != 3:
            args.extend(['-b', str(self.base_quality_cutoff)])
        
        if self.quality_score_trim:
            args.append('--quality-score-trim')
        
        if self.require_match_length != 0:
            args.extend(['--require-match-length', str(self.require_match_length)])
        
        if self.require_match_fraction != 0.9:
            args.extend(['--require-match-fraction', str(self.require_match_fraction)])
        
        if self.maximum_read_mismatches is not None:
            args.extend(['--maximum-read-mismatches', str(self.maximum_read_mismatches)])
        
        # Bowtie2 options (only if different from defaults)
        default_bowtie2_scoring = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals"
        if self.bowtie2_scoring != default_bowtie2_scoring:
            args.extend(['--bowtie2-scoring', self.bowtie2_scoring])
        
        default_bowtie2_stage1 = "--local -i S,1,0.25 --score-min L,1,0.9 -k 2000"
        if self.bowtie2_stage1 != default_bowtie2_stage1:
            args.extend(['--bowtie2-stage1', self.bowtie2_stage1])
        
        default_bowtie2_stage2 = "--local -i S,1,0.25 --score-min L,6,0.2 -k 2000"
        if self.bowtie2_stage2 != default_bowtie2_stage2:
            args.extend(['--bowtie2-stage2', self.bowtie2_stage2])
        
        default_bowtie2_junction = "--local -i S,1,0.25 --score-min L,1,0.70 -k 2000"
        if self.bowtie2_junction != default_bowtie2_junction:
            args.extend(['--bowtie2-junction', self.bowtie2_junction])
        
        # Junction evidence options
        if self.no_junction_prediction:
            args.append('--no-junction-prediction')
        
        if self.junction_indel_split_length != 3:
            args.extend(['--junction-indel-split-length', str(self.junction_indel_split_length)])
        
        if self.junction_alignment_pair_limit != 100000:
            args.extend(['--junction-alignment-pair-limit', str(self.junction_alignment_pair_limit)])
        
        if self.junction_minimum_candidates != 100:
            args.extend(['--junction-minimum-candidates', str(self.junction_minimum_candidates)])
        
        if self.junction_maximum_candidates != 5000:
            args.extend(['--junction-maximum-candidates', str(self.junction_maximum_candidates)])
        
        if self.junction_candidate_length_factor != 0.1:
            args.extend(['--junction-candidate-length-factor', str(self.junction_candidate_length_factor)])
        
        if self.junction_minimum_candidate_pos_hash_score != 2:
            args.extend(['--junction-minimum-candidate-pos-hash-score', str(self.junction_minimum_candidate_pos_hash_score)])
        
        if self.junction_score_cutoff != 3.0:
            args.extend(['--junction-score-cutoff', str(self.junction_score_cutoff)])
        
        if self.junction_minimum_pos_hash_score is not None:
            args.extend(['--junction-minimum-pos-hash-score', str(self.junction_minimum_pos_hash_score)])
        
        if self.junction_minimum_side_match is not None:
            args.extend(['--junction-minimum-side-match', str(self.junction_minimum_side_match)])
        
        if self.junction_minimum_pr_no_read_start_per_position != 0.1:
            args.extend(['--junction-minimum-pr-no-read-start-per-position', str(self.junction_minimum_pr_no_read_start_per_position)])
        
        if self.junction_allow_suboptimal_matches:
            args.append('--junction-allow-suboptimal-matches')
        
        # Missing coverage options
        if self.deletion_coverage_seed_cutoff != 0:
            args.extend(['--deletion-coverage-seed-cutoff', str(self.deletion_coverage_seed_cutoff)])
        
        if self.deletion_coverage_propagation_cutoff != 0:
            args.extend(['--deletion-coverage-propagation-cutoff', str(self.deletion_coverage_propagation_cutoff)])
        
        if self.call_mutations_overlapping_MC:
            args.append('--call-mutations-overlapping-MC')
        
        # Consensus RA evidence options
        if self.consensus_score_cutoff != 10:
            args.extend(['--consensus-score-cutoff', str(self.consensus_score_cutoff)])
        
        if self.consensus_frequency_cutoff != 0.8:
            args.extend(['--consensus-frequency-cutoff', str(self.consensus_frequency_cutoff)])
        
        if self.consensus_minimum_variant_coverage != 0:
            args.extend(['--consensus-minimum-variant-coverage', str(self.consensus_minimum_variant_coverage)])
        
        if self.consensus_minimum_total_coverage != 0:
            args.extend(['--consensus-minimum-total-coverage', str(self.consensus_minimum_total_coverage)])
        
        if self.consensus_minimum_variant_coverage_each_strand != 0:
            args.extend(['--consensus-minimum-variant-coverage-each-strand', str(self.consensus_minimum_variant_coverage_each_strand)])
        
        if self.consensus_minimum_total_coverage_each_strand != 0:
            args.extend(['--consensus-minimum-total-coverage-each-strand', str(self.consensus_minimum_total_coverage_each_strand)])
        
        if self.consensus_reject_indel_homopolymer_length != 0:
            args.extend(['--consensus-reject-indel-homopolymer-length', str(self.consensus_reject_indel_homopolymer_length)])
        
        if self.consensus_reject_surrounding_homopolymer_length != 0:
            args.extend(['--consensus-reject-surrounding-homopolymer-length', str(self.consensus_reject_surrounding_homopolymer_length)])
        
        # Polymorphism RA evidence options
        if self.polymorphism_score_cutoff is not None:
            args.extend(['--polymorphism-score-cutoff', str(self.polymorphism_score_cutoff)])
        
        if self.polymorphism_frequency_cutoff is not None:
            args.extend(['--polymorphism-frequency-cutoff', str(self.polymorphism_frequency_cutoff)])
        
        if self.polymorphism_minimum_variant_coverage != 0:
            args.extend(['--polymorphism-minimum-variant-coverage', str(self.polymorphism_minimum_variant_coverage)])
        
        if self.polymorphism_minimum_total_coverage != 0:
            args.extend(['--polymorphism-minimum-total-coverage', str(self.polymorphism_minimum_total_coverage)])
        
        if self.polymorphism_minimum_variant_coverage_each_strand is not None:
            args.extend(['--polymorphism-minimum-variant-coverage-each-strand', str(self.polymorphism_minimum_variant_coverage_each_strand)])
        
        if self.polymorphism_minimum_total_coverage_each_strand != 0:
            args.extend(['--polymorphism-minimum-total-coverage-each-strand', str(self.polymorphism_minimum_total_coverage_each_strand)])
        
        if self.polymorphism_bias_cutoff != 0:
            args.extend(['--polymorphism-bias-cutoff', str(self.polymorphism_bias_cutoff)])
        
        if self.polymorphism_no_indels:
            args.append('--polymorphism-no-indels')
        
        if self.polymorphism_reject_indel_homopolymer_length is not None:
            args.extend(['--polymorphism-reject-indel-homopolymer-length', str(self.polymorphism_reject_indel_homopolymer_length)])
        
        if self.polymorphism_reject_surrounding_homopolymer_length is not None:
            args.extend(['--polymorphism-reject-surrounding-homopolymer-length', str(self.polymorphism_reject_surrounding_homopolymer_length)])
        
        # Output options
        if self.max_displayed_reads != 100:
            args.extend(['--max-displayed-reads', str(self.max_displayed_reads)])
        
        if self.brief_html_output:
            args.append('--brief-html-output')
        
        if self.header_genome_diff:
            args.extend(['-g', self.header_genome_diff])
        
        if self.no_javascript:
            args.append('--no-javascript')
        
        # Pipeline control options
        if self.skip_RA_MC_prediction:
            args.append('--skip-RA-MC-prediction')
        
        if self.skip_JC_prediction:
            args.append('--skip-JC-prediction')
        
        if self.skip_MC_prediction:
            args.append('--skip-MC-prediction')
        
        # Debugging options
        if self.keep_intermediates:
            args.append('-k')
        
        if self.per_position_file:
            args.append('--per-position-file')
        
        if self.junction_debug:
            args.append('--junction-debug')
        
        # Experimental options
        if self.cnv:
            args.append('--cnv')
        
        if self.cnv_tile_size != 500:
            args.extend(['--cnv-tile-size', str(self.cnv_tile_size)])
        
        if self.cnv_ignore_redundant:
            args.append('--cnv-ignore-redundant')
        
        # Nanopore preset
        if self.nanopore:
            args.append('-x')
        
        return args


class Breseq:
    """Manages breseq execution for a sample.
    
    Can be created from an existing run or for a new run.
    """
    
    def __init__(self, sample_path: Optional[str] = None, params: Optional[Breseq_params] = None):
        """Initialize Breseq for a new run.
        
        Args:
            sample_path: Path to sample sequence data directory
            params: Breseq_params object with all required parameters
        """
        if sample_path is None or params is None:
            raise ValueError("For new runs, both sample_path and params must be provided. "
                           "Use Breseq.from_existing() to load from existing run.")
        
        self.sample_path = Path(sample_path)
        if not self.sample_path.exists():
            raise FileNotFoundError(f"Sample path does not exist: {self.sample_path}")
        
        self.params = params
        self.params._validate_required()

        self.reference_path = Path(_get_ref_genomes_path()) / self.params.reference
        
        # Compute output folder
        self.output_folder = self.sample_path / 'breseq' / self.params.version_name
        self.exists = self.output_folder.exists() and (self.output_folder / 'output' / 'output.done').exists()
        # Per-instance cache of region average coverage values
        # Key: region string (e.g., "NC_005966:1-1000") -> float or None
        self.region_average_cov = {}
    
    @classmethod
    def from_existing(cls, output_folder: Union[str, Path]) -> 'Breseq':
        """Create Breseq object from existing breseq run.
        
        Args:
            output_folder: Path to breseq output folder (e.g., breseq_params_v1)
            
        Returns:
            Breseq object with parameters loaded from log.txt and summary.json
        """
        output_folder = Path(output_folder)
        
        # Find log.txt (could be in output/ subfolder or directly in output_folder)
        log_path = output_folder / 'output' / 'log.txt'
        if not log_path.exists():
            log_path = output_folder / 'log.txt'
        if not log_path.exists():
            raise FileNotFoundError(f"Could not find log.txt in {output_folder}")
        
        # Find summary.json
        summary_path = output_folder / 'output' / 'summary.json'
        if not summary_path.exists():
            summary_path = output_folder / 'summary.json'
        
        # Parse parameters
        params = cls._parse_params_from_log(log_path)
        
        # If summary.json exists, use it to fill in missing parameters
        if summary_path.exists():
            cls._parse_params_from_summary(summary_path, params)
        
        # Infer sample_path from log.txt or output folder structure
        sample_path = cls._infer_sample_path(log_path, output_folder)
        
        # Create Breseq object
        breseq = cls.__new__(cls)
        breseq.sample_path = Path(sample_path)
        breseq.params = params
        # Keep a copy of the reference on the Breseq object for easy
        # access.
        breseq.reference = breseq.params.reference
        breseq.reference_path = Path(_get_ref_genomes_path()) / breseq.reference
        breseq.output_folder = output_folder

        breseq.exists = True
        # Initialize per-instance cache for region coverage
        breseq.region_average_cov = {}

        return breseq
    
    @staticmethod
    def _parse_params_from_log(log_path: Path) -> Breseq_params:
        """Parse breseq parameters from log.txt file.
        
        Args:
            log_path: Path to log.txt file
            
        Returns:
            Breseq_params object with parsed parameters
        """
        with open(log_path, 'r') as f:
            content = f.read()
        
        # Find the breseq command line
        # Look for line starting with "breseq"
        lines = content.split('\n')
        command_line = None
        for line in lines:
            if line.strip().startswith('breseq'):
                command_line = line.strip()
                break
        
        if not command_line:
            raise ValueError(f"Could not find breseq command in {log_path}")
        
        # Parse command line arguments
        # Use shlex to properly handle quoted arguments
        try:
            parts = shlex.split(command_line)
        except ValueError:
            # Fallback: simple split if shlex fails
            parts = command_line.split()
        
        # Remove 'breseq' command
        if parts[0] == 'breseq':
            parts = parts[1:]
        
        # Parse arguments into a dictionary
        args_dict = {}
        ref_genomes_path = _get_ref_genomes_path()
        i = 0
        while i < len(parts):
            arg = parts[i]
            
            # Handle flags (no value)
            if arg in ['-p', '--polymorphism-prediction']:
                args_dict['polymorphism_prediction'] = True
                i += 1
            elif arg in ['-t', '--targeted-sequencing']:
                args_dict['targeted_sequencing'] = True
                i += 1
            elif arg in ['-x', '--nanopore']:
                args_dict['nanopore'] = True
                i += 1
            elif arg == '--aligned-sam':
                args_dict['aligned_sam'] = True
                i += 1
            elif arg == '--long-read-distribute-remainder':
                args_dict['long_read_distribute_remainder'] = True
                i += 1
            elif arg == '--quality-score-trim':
                args_dict['quality_score_trim'] = True
                i += 1
            elif arg == '--no-junction-prediction':
                args_dict['no_junction_prediction'] = True
                i += 1
            elif arg == '--junction-allow-suboptimal-matches':
                args_dict['junction_allow_suboptimal_matches'] = True
                i += 1
            elif arg == '--call-mutations-overlapping-MC':
                args_dict['call_mutations_overlapping_MC'] = True
                i += 1
            elif arg == '--polymorphism-no-indels':
                args_dict['polymorphism_no_indels'] = True
                i += 1
            elif arg == '--brief-html-output':
                args_dict['brief_html_output'] = True
                i += 1
            elif arg == '--no-javascript':
                args_dict['no_javascript'] = True
                i += 1
            elif arg in ['--skip-RA-MC-prediction', '--skip-RA-MC-prediction']:
                args_dict['skip_RA_MC_prediction'] = True
                i += 1
            elif arg == '--skip-JC-prediction':
                args_dict['skip_JC_prediction'] = True
                i += 1
            elif arg == '--skip-MC-prediction':
                args_dict['skip_MC_prediction'] = True
                i += 1
            elif arg in ['-k', '--keep-intermediates']:
                args_dict['keep_intermediates'] = True
                i += 1
            elif arg == '--per-position-file':
                args_dict['per_position_file'] = True
                i += 1
            elif arg == '--junction-debug':
                args_dict['junction_debug'] = True
                i += 1
            elif arg == '--cnv':
                args_dict['cnv'] = True
                i += 1
            elif arg == '--cnv-ignore-redundant':
                args_dict['cnv_ignore_redundant'] = True
                i += 1
            # Handle parameters with values
            elif arg in ['-r', '--reference']:
                # Only a single reference is supported. If multiple -r are
                # present in the log, raise an error.
                if i + 1 < len(parts):
                    ref_path = parts[i + 1]
                    # Extract filename from full path
                    if ref_genomes_path in ref_path:
                        ref_file = os.path.relpath(ref_path, ref_genomes_path)
                    else:
                        ref_file = os.path.basename(ref_path)
                    if 'reference' in args_dict:
                        raise ValueError('Multiple -r/--reference entries found; only one reference is supported')
                    args_dict['reference'] = ref_file
                    i += 2
                else:
                    i += 1
            elif arg in ['-n', '--name']:
                if i + 1 < len(parts):
                    args_dict['name'] = parts[i + 1]
                    i += 2
                else:
                    i += 1
            elif arg in ['-j', '--num-processors']:
                if i + 1 < len(parts):
                    args_dict['num_processors'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg in ['-o', '--output']:
                # Skip output path (not stored in params)
                if i + 1 < len(parts):
                    i += 2
                else:
                    i += 1
            elif arg in ['-c', '--contig-reference']:
                if i + 1 < len(parts):
                    args_dict['contig_reference'] = parts[i + 1]
                    i += 2
                else:
                    i += 1
            elif arg in ['-s', '--junction-only-reference']:
                if 'junction_only_reference' not in args_dict:
                    args_dict['junction_only_reference'] = []
                if i + 1 < len(parts):
                    args_dict['junction_only_reference'].append(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--user-evidence-gd':
                if i + 1 < len(parts):
                    args_dict['user_evidence_gd'] = parts[i + 1]
                    i += 2
                else:
                    i += 1
            elif arg in ['-l', '--limit-fold-coverage']:
                if i + 1 < len(parts):
                    args_dict['limit_fold_coverage'] = float(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--read-min-length':
                if i + 1 < len(parts):
                    args_dict['read_min_length'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--read-max-same-base-fraction':
                if i + 1 < len(parts):
                    args_dict['read_max_same_base_fraction'] = float(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--read-max-N-fraction':
                if i + 1 < len(parts):
                    args_dict['read_max_N_fraction'] = float(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--long-read-trigger-length':
                if i + 1 < len(parts):
                    args_dict['long_read_trigger_length'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--long-read-split-length':
                if i + 1 < len(parts):
                    args_dict['long_read_split_length'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--genbank-field-for-seq-id':
                if i + 1 < len(parts):
                    args_dict['genbank_field_for_seq_id'] = parts[i + 1]
                    i += 2
                else:
                    i += 1
            elif arg in ['-m', '--minimum-mapping-quality']:
                if i + 1 < len(parts):
                    args_dict['minimum_mapping_quality'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg in ['-b', '--base-quality-cutoff']:
                if i + 1 < len(parts):
                    args_dict['base_quality_cutoff'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--require-match-length':
                if i + 1 < len(parts):
                    args_dict['require_match_length'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--require-match-fraction':
                if i + 1 < len(parts):
                    args_dict['require_match_fraction'] = float(parts[i + 1])
                    i += 2
                else:
                    i += 1
            elif arg == '--maximum-read-mismatches':
                if i + 1 < len(parts):
                    args_dict['maximum_read_mismatches'] = int(parts[i + 1])
                    i += 2
                else:
                    i += 1
            # Skip FASTQ file arguments (they're not parameters)
            elif arg.endswith('.fastq') or arg.endswith('.fastq.gz') or arg.endswith('.fq') or arg.endswith('.fq.gz'):
                i += 1
            else:
                # Try to handle other arguments generically
                # Check if it's a known parameter format
                if arg.startswith('--'):
                    param_name = _convert_param_name_to_python(arg)
                    if i + 1 < len(parts) and not parts[i + 1].startswith('-'):
                        # Has a value
                        try:
                            # Try to convert to appropriate type
                            value = parts[i + 1]
                            # Try int
                            try:
                                value = int(value)
                            except ValueError:
                                # Try float
                                try:
                                    value = float(value)
                                except ValueError:
                                    pass  # Keep as string
                            args_dict[param_name] = value
                            i += 2
                        except:
                            i += 1
                    else:
                        # Boolean flag
                        args_dict[param_name] = True
                        i += 1
                else:
                    i += 1
        
        # Ensure polymorphism_prediction is set (default to False if not found)
        if 'polymorphism_prediction' not in args_dict:
            args_dict['polymorphism_prediction'] = False
        
        # Create Breseq_params object
        return Breseq_params(**args_dict)
    
    @staticmethod
    def _parse_params_from_summary(summary_path: Path, params: Breseq_params):
        """Parse additional parameters from summary.json and update params.
        
        Args:
            summary_path: Path to summary.json file
            params: Breseq_params object to update
        """
        with open(summary_path, 'r') as f:
            summary = json.load(f)
        
        if 'options' not in summary:
            return
        
        options = summary['options']
        
        # Map summary.json structure to parameters
        # workflow section
        if 'workflow' in options:
            wf = options['workflow']
            if 'num_processors' in wf and not hasattr(params, 'num_processors'):
                params.num_processors = wf['num_processors']
            if 'genbank_field_for_seq_id' in wf:
                params.genbank_field_for_seq_id = wf['genbank_field_for_seq_id']
            if 'polymorphism_prediction' in wf:
                params.polymorphism_prediction = wf.get('polymorphism_prediction', False)
        
        # mutation_identification section
        if 'mutation_identification' in options:
            mi = options['mutation_identification']
            if 'polymorphism_prediction' in mi:
                params.polymorphism_prediction = mi['polymorphism_prediction']
            if 'base_quality_cutoff' in mi:
                params.base_quality_cutoff = mi['base_quality_cutoff']
            if 'mutation_log10_e_value_cutoff' in mi:
                params.consensus_score_cutoff = mi['mutation_log10_e_value_cutoff']
            if 'consensus_frequency_cutoff' in mi:
                params.consensus_frequency_cutoff = mi['consensus_frequency_cutoff']
            if 'consensus_minimum_variant_coverage' in mi:
                params.consensus_minimum_variant_coverage = mi['consensus_minimum_variant_coverage']
            if 'consensus_minimum_total_coverage' in mi:
                params.consensus_minimum_total_coverage = mi['consensus_minimum_total_coverage']
            if 'consensus_minimum_variant_coverage_each_strand' in mi:
                params.consensus_minimum_variant_coverage_each_strand = mi['consensus_minimum_variant_coverage_each_strand']
            if 'consensus_minimum_total_coverage_each_strand' in mi:
                params.consensus_minimum_total_coverage_each_strand = mi['consensus_minimum_total_coverage_each_strand']
            if 'consensus_reject_indel_homopolymer_length' in mi:
                params.consensus_reject_indel_homopolymer_length = mi['consensus_reject_indel_homopolymer_length']
            if 'consensus_reject_surrounding_homopolymer_length' in mi:
                params.consensus_reject_surrounding_homopolymer_length = mi['consensus_reject_surrounding_homopolymer_length']
            if 'polymorphism_log10_e_value_cutoff' in mi:
                params.polymorphism_score_cutoff = mi['polymorphism_log10_e_value_cutoff']
            if 'polymorphism_frequency_cutoff' in mi:
                params.polymorphism_frequency_cutoff = mi['polymorphism_frequency_cutoff']
            if 'polymorphism_minimum_variant_coverage' in mi:
                params.polymorphism_minimum_variant_coverage = mi['polymorphism_minimum_variant_coverage']
            if 'polymorphism_minimum_total_coverage' in mi:
                params.polymorphism_minimum_total_coverage = mi['polymorphism_minimum_total_coverage']
            if 'polymorphism_minimum_variant_coverage_each_strand' in mi:
                params.polymorphism_minimum_variant_coverage_each_strand = mi['polymorphism_minimum_variant_coverage_each_strand']
            if 'polymorphism_minimum_total_coverage_each_strand' in mi:
                params.polymorphism_minimum_total_coverage_each_strand = mi['polymorphism_minimum_total_coverage_each_strand']
            if 'polymorphism_bias_p_value_cutoff' in mi:
                params.polymorphism_bias_cutoff = mi['polymorphism_bias_p_value_cutoff']
            if 'no_indel_polymorphisms' in mi:
                params.polymorphism_no_indels = mi['no_indel_polymorphisms']
            if 'polymorphism_reject_indel_homopolymer_length' in mi:
                params.polymorphism_reject_indel_homopolymer_length = mi['polymorphism_reject_indel_homopolymer_length']
            if 'polymorphism_reject_surrounding_homopolymer_length' in mi:
                params.polymorphism_reject_surrounding_homopolymer_length = mi['polymorphism_reject_surrounding_homopolymer_length']
            if 'deletion_coverage_seed_cutoff' in mi:
                params.deletion_coverage_seed_cutoff = mi['deletion_coverage_seed_cutoff']
            if 'deletion_coverage_propagation_cutoff' in mi:
                params.deletion_coverage_propagation_cutoff = mi['deletion_coverage_propagation_cutoff']
            if 'targeted_sequencing' in mi:
                params.targeted_sequencing = mi['targeted_sequencing']
            if 'quality_score_trim' in mi:
                params.quality_score_trim = bool(mi['quality_score_trim'])
        
        # read_alignment section
        if 'read_alignment' in options:
            ra = options['read_alignment']
            if 'bowtie2_scoring' in ra:
                params.bowtie2_scoring = ra['bowtie2_scoring']
            if 'bowtie2_stage1' in ra:
                params.bowtie2_stage1 = ra['bowtie2_stage1']
            if 'bowtie2_stage2' in ra:
                params.bowtie2_stage2 = ra['bowtie2_stage2']
            if 'bowtie2_junction' in ra:
                params.bowtie2_junction = ra['bowtie2_junction']
            if 'minimum_mapping_quality' in ra:
                params.minimum_mapping_quality = ra['minimum_mapping_quality']
            if 'require_match_fraction' in ra:
                params.require_match_fraction = ra['require_match_fraction']
            if 'require_match_length' in ra:
                params.require_match_length = ra['require_match_length']
            if 'maximum_read_mismatches' in ra and ra['maximum_read_mismatches'] >= 0:
                params.maximum_read_mismatches = ra['maximum_read_mismatches']
        
        # read_file section
        if 'read_file' in options:
            rf = options['read_file']
            if 'read_file_read_length_min' in rf:
                params.read_min_length = rf['read_file_read_length_min']
            if 'read_file_max_same_base_fraction' in rf:
                params.read_max_same_base_fraction = rf['read_file_max_same_base_fraction']
            if 'read_file_max_N_fraction' in rf:
                params.read_max_N_fraction = rf['read_file_max_N_fraction']
            if 'read_file_long_read_trigger_length' in rf:
                params.long_read_trigger_length = rf['read_file_long_read_trigger_length']
            if 'read_file_long_read_split_length' in rf:
                params.long_read_split_length = rf['read_file_long_read_split_length']
            if 'read_file_long_read_distribute_remainder' in rf:
                params.long_read_distribute_remainder = rf['read_file_long_read_distribute_remainder']
            if 'read_file_coverage_fold_limit' in rf and rf['read_file_coverage_fold_limit'] > 0:
                params.limit_fold_coverage = rf['read_file_coverage_fold_limit']
            if 'aligned_sam_mode' in rf:
                params.aligned_sam = rf['aligned_sam_mode']
        
        # candidate_junction section
        if 'candidate_junction' in options:
            cj = options['candidate_junction']
            if 'preprocess_junction_min_indel_split_length' in cj:
                params.junction_indel_split_length = cj['preprocess_junction_min_indel_split_length']
            if 'maximum_junction_sequence_passed_alignment_pairs_to_consider' in cj:
                params.junction_alignment_pair_limit = cj['maximum_junction_sequence_passed_alignment_pairs_to_consider']
            if 'minimum_candidate_junctions' in cj:
                params.junction_minimum_candidates = cj['minimum_candidate_junctions']
            if 'maximum_candidate_junctions' in cj:
                params.junction_maximum_candidates = cj['maximum_candidate_junctions']
            if 'maximum_candidate_junction_length_factor' in cj:
                params.junction_candidate_length_factor = cj['maximum_candidate_junction_length_factor']
            if 'minimum_candidate_junction_pos_hash_score' in cj:
                params.junction_minimum_candidate_pos_hash_score = cj['minimum_candidate_junction_pos_hash_score']
        
        # alignment_resolution section
        if 'alignment_resolution' in options:
            ar = options['alignment_resolution']
            if 'junction_pos_hash_neg_log10_p_value_cutoff' in ar:
                params.junction_score_cutoff = ar['junction_pos_hash_neg_log10_p_value_cutoff']
            if 'minimum_alignment_resolution_pos_hash_score' in ar:
                params.junction_minimum_pos_hash_score = ar['minimum_alignment_resolution_pos_hash_score']
            if 'junction_minimum_side_match' in ar:
                params.junction_minimum_side_match = ar['junction_minimum_side_match']
            if 'minimum_pr_no_read_start_per_position' in ar:
                params.junction_minimum_pr_no_read_start_per_position = ar['minimum_pr_no_read_start_per_position']
            if 'junction_allow_suboptimal_matches' in ar:
                params.junction_allow_suboptimal_matches = ar['junction_allow_suboptimal_matches']
        
        # output section
        if 'output' in options:
            out = options['output']
            if 'max_displayed_reads' in out:
                params.max_displayed_reads = out['max_displayed_reads']
            if 'header_genome_diff_file_name' in out and out['header_genome_diff_file_name']:
                params.header_genome_diff = out['header_genome_diff_file_name']
            if 'no_javascript' in out:
                params.no_javascript = out['no_javascript']
        
        # Add more mappings as needed
    
    @staticmethod
    def _infer_sample_path(log_path: Path, output_folder: Path) -> str:
        """Infer sample path from log.txt or output folder structure.
        
        Args:
            log_path: Path to log.txt
            output_folder: Path to breseq output folder
            
        Returns:
            Inferred sample path
        """
        # Try to extract from log.txt first
        with open(log_path, 'r') as f:
            content = f.read()
        
        # Look for FASTQ file paths in the command
        fastq_pattern = r'([^\s]+\.(fastq|fq)(\.gz)?)'
        matches = re.findall(fastq_pattern, content)
        if matches:
            # Get the first FASTQ file path
            first_fastq = matches[0][0]
            fastq_path = Path(first_fastq)
            # Assume sample path is the parent of the directory containing FASTQ files
            # Common structure: .../trimmed/file.fastq.gz or .../received/file.fastq.gz
            if fastq_path.parent.name in ['trimmed', 'received']:
                sample_path = fastq_path.parent.parent
            else:
                sample_path = fastq_path.parent
            return str(sample_path)
        
        # Fallback: infer from output folder structure
        # output_folder should be: sample_path/breseq/breseq_params_vX
        if output_folder.parent.name == 'breseq':
            return str(output_folder.parent.parent)
        
        # Last resort: return parent of output_folder
        return str(output_folder.parent)
    
    def run(self, overwrite: bool = False):
        """Run breseq with the configured parameters.
        
        Args:
            overwrite: If True, overwrite existing results. If False, raise error if exists.
        
        Raises:
            FileExistsError: If output exists and overwrite=False
            RuntimeError: If breseq command fails
        """
        # Check if output already exists
        if self.exists and not overwrite:
            raise FileExistsError(
                f"Breseq output already exists at {self.output_folder}. "
                "Set overwrite=True to overwrite."
            )
        
        # Create output directory
        self.output_folder.mkdir(parents=True, exist_ok=True)
        
        # Find FASTQ files in sample_path
        fastq_files = []
        for pattern in ['*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz']:
            fastq_files.extend(list(self.sample_path.glob(pattern)))
        
        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in {self.sample_path}")
        
        # Sort FASTQ files for consistent ordering
        fastq_files = sorted(fastq_files)
        
        # Build command
        cmd = ['breseq'] + self.params.to_command_args()
        cmd.extend([str(f) for f in fastq_files])
        cmd.extend(['-o', str(self.output_folder)])
        
        # Run breseq
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            self.exists = True

            # Add #=TITLE sample_name to genome diff file
            gd_file = os.path.join(str(self.output_folder), 'data', 'output.gd')
            # Read the existing content
            with open(gd_file, 'r') as f:
                content = f.readlines()
            # Add the title line at the beginning
            title_line = f'#=TITLE\t{os.path.basename(output_dir)}\n'
            content.insert(0, title_line)
            # Write back the modified content
            with open(gd_file, 'w') as f:
                f.writelines(content)

        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Breseq command failed with return code {e.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Error output: {e.stderr}"
            ) from e

    def run_bam2cov(self, region: str) -> str:
        """Run breseq BAM2COV for this Breseq run and a region.

        Requires that the Breseq run exists (``self.exists`` is True). The
        function will use the BAM produced by the run at
        ``<sample_dir>/data/reference.bam``.
        """
        if not getattr(self, 'exists', False):
            raise RuntimeError(
                "Breseq output not found for this object. Run '\n"
                "breseq_obj.run() to produce results before running BAM2COV."
            )

        # Use the run's output folder for BAM2COV outputs
        bam2cov_dir = Path(self.output_folder) / 'BAM2COV'
        bam2cov_dir.mkdir(parents=True, exist_ok=True)

        # Determine BAM input path from the run output
        bam_path = Path(self.output_folder) / 'data' / 'reference.bam'
        if not bam_path.exists():
            raise FileNotFoundError(
                f"BAM file not found for this run: {bam_path}. "
                "Ensure the Breseq run produced the BAM file."
            )

        # Determine fasta input path from the run output
        fasta_path = Path(self.output_folder) / 'data' / 'reference.fasta'
        if not bam_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found for this run: {fasta_path}. "
                "Ensure there is a reference.fasta file in the data directory."
            )

        # Build output base name and path
        outfile_base = region.replace(':', '_').replace('-', '_')
        outfile = bam2cov_dir / outfile_base

        bam2cov_cmd = [
            'breseq', 'BAM2COV',
            '-b', str(bam_path),
            '-o', str(outfile),
            '-r', region,
            '-f', str(fasta_path),
            '-t',
        ]

        try:
            subprocess.run(
                bam2cov_cmd,
                check=True,
                cwd=str(self.sample_path),
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"BAM2COV command failed (rc={e.returncode}).\n"
                f"Command: {' '.join(bam2cov_cmd)}\n"
                f"stderr: {e.stderr}"
            ) from e

        return str(outfile) + '.tab'

    def get_region_average_coverage(
        self,
        region: str,
        run_if_missing: bool = False,
    ) -> Optional[float]:
        """Return average coverage for a region and cache it on this object.

        Behavior:
          - If the coverage .tab file for ``region`` exists under
            ``<output_folder>/BAM2COV``, the file is parsed and the value
            returned.
          - If not present and ``run_if_missing`` is True, ``run_bam2cov`` is
            invoked to produce the file, which is then parsed.

        The parsed value is stored in ``self.region_average_cov[region]``.
        """
        # Return cached value if present
        if region in self.region_average_cov:
            return self.region_average_cov[region]

        bam2cov_dir = Path(self.output_folder) / 'BAM2COV'
        outfile_base = region.replace(':', '_').replace('-', '_')
        tab_path = bam2cov_dir / (outfile_base + '.tab')

        def _parse_tab(path: Union[str, Path]) -> Optional[float]:
            p = Path(path)
            if not p.exists():
                return None
            try:
                with open(p, 'r') as fh:
                    for line in fh:
                        if line.startswith('#') and 'region_average_cov' in line:
                            parts = line.strip().split('\t')
                            if parts:
                                try:
                                    return float(parts[-1])
                                except ValueError:
                                    return None
            except Exception:
                return None
            return None

        if tab_path.exists():
            val = _parse_tab(tab_path)
            self.region_average_cov[region] = val
            return val

        if not run_if_missing:
            self.region_average_cov[region] = None
            return None

        # Run BAM2COV to create the .tab file, then parse
        tabfile = self.run_bam2cov(region)
        val = _parse_tab(tabfile)
        self.region_average_cov[region] = val
        return val

    def count_mutations(
        self,
        output_csv: Optional[Union[str, Path]] = None,
        detailed_output: Optional[Union[str, Path]] = None,
        verbose: bool = False,
    ) -> Dict[str, Optional[int]]:
        """Run `gdtools COUNT` on the run's GenomeDiff and parse totals.

        Writes CSV outputs under the run folder when paths are not provided.

        Args:
            output_csv: Path for the consensus COUNT CSV (optional).
            detailed_output: Path for the polymorphism COUNT CSV (optional).
            verbose: If True, run gdtools with -v.

        Returns:
            Dict with keys 'consensus' and 'polymorphisms' mapping to counts
            (or None if not available).
        """
        # Only consider the run's main genome-diff file: data/output.gd
        gd_path = Path(self.output_folder) / 'data' / 'output.gd'
        if not gd_path.exists():
            raise FileNotFoundError(f"Could not find data/output.gd in {self.output_folder}")

        def _run_count(out_path: Path, extra_args: Optional[List[str]] = None):
            cmd = ['gdtools', 'COUNT']
            cmd.extend(['-r', str(self.reference_path)])
            if verbose:
                cmd.append('-v')
            if extra_args:
                cmd.extend(extra_args)
            cmd.extend(['-o', str(out_path)])
            cmd.append(str(gd_path))
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"gdtools COUNT failed (rc={e.returncode}).\n"
                    f"cmd: {' '.join(cmd)}\n"
                    f"stderr: {e.stderr}"
                ) from e

        def _parse_count_csv(path: Path) -> Optional[int]:
            if not path.exists():
                return None
            with open(path, 'r') as fh:
                reader = csv.DictReader(fh)
                if not reader.fieldnames:
                    return None
                key = None
                for h in reader.fieldnames:
                    if h and h.lower() == 'total':
                        key = h
                        break
                if key is None:
                    return None
                total = 0
                for row in reader:
                    v = row.get(key)
                    if not v:
                        continue
                    try:
                        total += int(float(v))
                    except Exception:
                        continue
                return total

        out_cons = Path(output_csv) if output_csv else Path(self.output_folder) / 'gdtools_count.csv'
        _run_count(out_cons)
        cons_total = _parse_count_csv(out_cons)
        self.consensus_mutation_count = cons_total
        self.gdtools_count_csv = str(out_cons)

        poly_total: Optional[int] = None
        if getattr(self.params, 'polymorphism_prediction', False):
            out_poly = Path(detailed_output) if detailed_output else Path(self.output_folder) / 'gdtools_count_polymorphisms.csv'
            _run_count(out_poly, extra_args=['-p'])
            poly_total = _parse_count_csv(out_poly) - cons_total
            self.polymorphism_mutation_count = poly_total
            self.gdtools_count_polymorphisms_csv = str(out_poly)
        else:
            self.polymorphism_mutation_count = None

        return {'consensus': cons_total, 'polymorphisms': poly_total}

    def count_reads(self) -> Dict[str, Optional[int]]:
        """Parse input and mapped read counts from data/output.gd.

        Looks for header lines in the GenomeDiff file like:
        #=INPUT-READS\t<integer>
        #=MAPPED-READS\t<integer>

        Stores the parsed values on the Breseq object as
        ``self.input_read_count`` and ``self.mapped_read_count`` and also
        returns them in a dict: {'input': int|None, 'mapped': int|None}.

        Raises:
            FileNotFoundError: if data/output.gd does not exist for this run.
        """
        gd_path = Path(self.output_folder) / 'data' / 'output.gd'
        if not gd_path.exists():
            raise FileNotFoundError(f"Could not find data/output.gd in {self.output_folder}")

        input_reads: Optional[int] = None
        mapped_reads: Optional[int] = None

        # Read header-style lines that start with '#=' and parse key/value
        try:
            with open(gd_path, 'r') as fh:
                for line in fh:
                    if not line.startswith('#='):
                        # Once we pass headers, stop parsing for performance
                        # (GenomeDiff puts headers at the top)
                        break
                    # Remove leading '#=' then split on whitespace or tab
                    rest = line[2:].strip()
                    if not rest:
                        continue
                    parts = re.split(r"\s+", rest, maxsplit=1)
                    if len(parts) < 2:
                        continue
                    key = parts[0].upper()
                    val = parts[1].strip()
                    # Some keys may include underscores; normalize
                    key = key.replace('_', '-')
                    # Try to parse integer value
                    try:
                        num = int(float(val))
                    except Exception:
                        continue

                    if key == 'INPUT-READS' and input_reads is None:
                        input_reads = num
                    elif key == 'MAPPED-READS' and mapped_reads is None:
                        mapped_reads = num
        except Exception:
            # If anything goes wrong reading, leave values as None
            input_reads = input_reads
            mapped_reads = mapped_reads

        # Attach to object for convenience
        self.input_read_count = input_reads
        self.mapped_read_count = mapped_reads

        return {'input': input_reads, 'mapped': mapped_reads}

    def apply_mutations(
        self,
        output_path: Union[str, Path],
        format: str = 'GENBANK',
        seq_ids: Optional[List[str]] = None,
        polymorphism_mode: Optional[bool] = None,
        applied_gd: Optional[Union[str, Path]] = None,
        verbose: bool = False,
    ) -> str:
        """Apply mutations from data/output.gd to the reference using gdtools APPLY.

        By default this writes a GenBank file with the applied reference.

        Args:
            output_path: where to write the applied reference. If None, a
                sensible default under the run folder is used.
            format: Output format: one of 'GENBANK', 'FASTA', 'GFF3'.
            seq_ids: Optional list of sequence IDs to keep (passed with -s).
            polymorphism_mode: If True, include polymorphic mutations. If
                None, defaults to the run's polymorphism_prediction param.
            applied_gd: Optional path to write an updated GenomeDiff file
                (gdtools --applied-gd).
            verbose: If True, pass -v to gdtools.

        Returns:
            The path to the applied reference file (string).
        """
        # Locate the run's GenomeDiff
        gd_path = Path(self.output_folder) / 'data' / 'output.gd'
        if not gd_path.exists():
            raise FileNotFoundError(f"Could not find data/output.gd in {self.output_folder}")

        # Normalize format
        fmt = (format or 'GENBANK').upper()
        if fmt not in ('GENBANK', 'FASTA', 'GFF3'):
            raise ValueError("format must be one of 'GENBANK', 'FASTA', or 'GFF3'")

        # output_path is mandatory; normalize to Path
        out_path = Path(output_path)

        # Determine polymorphism mode default
        if polymorphism_mode is None:
            polymorphism_mode = bool(getattr(self.params, 'polymorphism_prediction', False))

        # Build command
        cmd = ['gdtools', 'APPLY']
        cmd.extend(['-r', str(self.reference_path)])
        cmd.extend(['-o', str(out_path)])
        cmd.extend(['-f', fmt])
        if seq_ids:
            for sid in seq_ids:
                cmd.extend(['-s', sid])
        if polymorphism_mode:
            cmd.append('-p')
        if applied_gd:
            cmd.extend(['--applied-gd', str(applied_gd)])
        if verbose:
            cmd.append('-v')
        cmd.append(str(gd_path))

        # Ensure output folder exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"gdtools APPLY failed (rc={e.returncode}).\n"
                f"cmd: {' '.join(cmd)}\n"
                f"stderr: {e.stderr}"
            ) from e

        # Record attributes on the object
        self.applied_reference = str(out_path)
        if applied_gd:
            self.applied_gd_path = str(applied_gd)

        return str(out_path)

