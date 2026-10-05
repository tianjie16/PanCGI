set -euo pipefail

usage() {
    printf '%s\n' \
      'Usage: run_pancgi.sh --prepared DIR --gfa FILE --hal FILE --genomes FILE --paths FILE --out-dir DIR [options]' \
      '  --hal-runtime docker|native  Explicit HAL execution mode (default: docker)' \
      '  --docker-image IMAGE        Required only for docker; already installed image' \
      '  --docker-bin PATH           Docker executable (default: docker)' \
      '  --python PATH               Python executable (default: python3)' \
      '  --hal-only-exclusions FILE  Reviewed HAL-only contigs with no GFA path and no CGI/SV input' \
      '  --threads N                 Graph, feature and analysis workers (default: 1)' \
      '  --hal-threads N             Concurrent HAL workers (default: 1)' \
      '  --hal-cpus N                CPUs per Docker HAL worker (default: 1)' \
      '  --hal-memory SIZE           Memory per Docker HAL worker, e.g. 1g (default: 1g)' \
      '  --hal-stats PATH            Explicit halStats executable (default: halStats)' \
      '  --hal2fasta PATH            Explicit hal2fasta executable (default: hal2fasta)' \
      '  --hal-liftover PATH          Explicit halLiftover executable (default: halLiftover)' \
      '  --scratch DIR               Existing temporary directory; output work remains under out-dir' \
      '  --help                      Show this help' \
      'Scientific profile: identity=0.80; length ratio=0.80; all-pairs clique; Parasail 32-bit; multi-scale graph genotyping.' \
      'Native execution inherits limits from its scheduler/container; threads alone are not a memory limit.'
}

GFA=
PREPARED=
HAL=
GENOMES=
PATHS=
HAL_ONLY_EXCLUSIONS=
OUTDIR=
PYTHON_BIN=python3
THREADS=1
HAL_THREADS=1
DOCKER_BIN=docker
DOCKER_IMAGE=
HAL_STATS=halStats
HAL2FASTA=hal2fasta
HAL_LIFTOVER=halLiftover
HAL_RUNTIME=docker
HAL_CPUS=1
HAL_MEMORY=1g
SCRATCH=

while [[ $# -gt 0 ]]; do
    if [[ "$1" != --help && "$1" != -h && $# -lt 2 ]]; then
        printf 'Missing value for option: %s\n' "$1" >&2
        exit 2
    fi
    case "$1" in
        --gfa) GFA=$2; shift 2 ;;
        --prepared) PREPARED=$2; shift 2 ;;
        --hal) HAL=$2; shift 2 ;;
        --genomes) GENOMES=$2; shift 2 ;;
        --paths) PATHS=$2; shift 2 ;;
        --hal-only-exclusions) HAL_ONLY_EXCLUSIONS=$2; shift 2 ;;
        --out-dir) OUTDIR=$2; shift 2 ;;
        --python) PYTHON_BIN=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --hal-threads) HAL_THREADS=$2; shift 2 ;;
        --docker-bin) DOCKER_BIN=$2; shift 2 ;;
        --docker-image) DOCKER_IMAGE=$2; shift 2 ;;
        --hal-stats) HAL_STATS=$2; shift 2 ;;
        --hal2fasta) HAL2FASTA=$2; shift 2 ;;
        --hal-liftover) HAL_LIFTOVER=$2; shift 2 ;;
        --hal-runtime) HAL_RUNTIME=$2; shift 2 ;;
        --hal-cpus) HAL_CPUS=$2; shift 2 ;;
        --hal-memory) HAL_MEMORY=$2; shift 2 ;;
        --scratch) SCRATCH=$2; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

for value in PREPARED GFA HAL GENOMES PATHS OUTDIR; do
    if [[ -z "${!value}" ]]; then
        printf 'Missing required option for %s\n' "$value" >&2
        usage >&2
        exit 2
    fi
done

for path in "$GFA" "$HAL" "$GENOMES" "$PATHS"; do
    if [[ ! -s "$path" ]]; then
        printf 'Input file is missing or empty: %s\n' "$path" >&2
        exit 2
    fi
done

for value in THREADS HAL_THREADS HAL_CPUS; do
    if [[ ! "${!value}" =~ ^[1-9][0-9]*$ ]]; then
        printf '%s must be a positive integer: %s\n' "$value" "${!value}" >&2
        exit 2
    fi
done

case "$HAL_RUNTIME" in
    docker)
        [[ -n "$DOCKER_IMAGE" ]] || { printf 'Docker HAL requires --docker-image\n' >&2; exit 2; }
        DOCKER_IMAGE=$("$DOCKER_BIN" image inspect "$DOCKER_IMAGE" --format '{{.Id}}')
        [[ "$DOCKER_IMAGE" =~ ^sha256:[a-f0-9]{64}$ ]] || exit 2
        [[ "$HAL_MEMORY" =~ ^[1-9][0-9]*[mMgG]$ ]] || { printf 'Invalid HAL memory limit\n' >&2; exit 2; }
        ;;
    native)
        [[ -z "$DOCKER_IMAGE" ]] || { printf 'Native HAL cannot use --docker-image\n' >&2; exit 2; }
        for program in "$HAL_STATS" "$HAL2FASTA" "$HAL_LIFTOVER"; do
            command -v "$program" >/dev/null || { printf 'Missing HAL executable: %s\n' "$program" >&2; exit 2; }
        done
        ;;
    *) printf 'Invalid --hal-runtime: %s\n' "$HAL_RUNTIME" >&2; exit 2 ;;
esac
export PANCGI_HAL_RUNTIME=$HAL_RUNTIME PANCGI_HAL_CPUS=$HAL_CPUS PANCGI_HAL_MEMORY=$HAL_MEMORY
if [[ -n "$SCRATCH" ]]; then
    [[ -d "$SCRATCH" && -w "$SCRATCH" ]] || { printf 'Scratch directory must exist and be writable\n' >&2; exit 2; }
    export TMPDIR=$SCRATCH
fi

if [[ -e "$OUTDIR" ]]; then
    printf 'Output path already exists: %s\n' "$OUTDIR" >&2
    exit 2
fi

PACKAGE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PIPELINE=${PACKAGE_DIR}/cpgi_nr_prod.py
MAPPING=${PACKAGE_DIR}/pancgi_mapping.py
INPUT_PREP=${PACKAGE_DIR}/pancgi_inputs.py
DEPENDENCY_CHECK=${PACKAGE_DIR}/scripts/check_runtime_dependencies.py
PREPARED_INPUT_VALIDATOR=${PACKAGE_DIR}/scripts/validate_prepared_inputs.py
HAL_PIPELINE=${PACKAGE_DIR}/pancgi_hal.py
WORKDIR=${OUTDIR}/work
MAPPINGDIR=${WORKDIR}/mapping
VALIDATED=${MAPPINGDIR}/validated
INPUTDIR=${WORKDIR}/inputs
PATHBED_DIR=${WORKDIR}/pathbed
HALDIR=${WORKDIR}/hal_liftover
RESULTDIR=${WORKDIR}/results_internal
CATALOG=${INPUTDIR}/sample_catalog.tsv
GFA_INVENTORY=${MAPPINGDIR}/gfa_paths.tsv
HAL_INVENTORY=${MAPPINGDIR}/hal_sequences.tsv

mkdir -p "$WORKDIR" "$RESULTDIR" "$MAPPINGDIR"
export PANCGI_PARAMETERS_DIR=$WORKDIR/parameters
export PYTHONHASHSEED=0
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

run_step() {
    local name=$1
    shift
    "$PYTHON_BIN" "$PACKAGE_DIR/scripts/run_stage.py" "$WORKDIR/stages" "$name" "$@"
}

run_step 00_validate_dependencies \
    "$PYTHON_BIN" "$DEPENDENCY_CHECK"

run_step 01_validate_preparation \
    "$PYTHON_BIN" "$PACKAGE_DIR/pancgi_preparation.py" validate \
    --prepared "$PREPARED" --gfa "$GFA" --hal "$HAL" --out-dir "$MAPPINGDIR"

exclusion_args=()
if [[ -n "$HAL_ONLY_EXCLUSIONS" ]]; then
    exclusion_args=(--hal-only-exclusions "$HAL_ONLY_EXCLUSIONS")
fi

run_step 03_validate_inputs \
    "$PYTHON_BIN" "$PACKAGE_DIR/pancgi_contract.py" \
    --gfa "$GFA" --hal "$HAL" \
    --genomes "$GENOMES" \
    --paths "$PATHS" \
    --gfa-inventory "$GFA_INVENTORY" \
    --hal-inventory "$HAL_INVENTORY" \
    --out-dir "$VALIDATED" "${exclusion_args[@]}"

run_step 04_unfold_graph \
    "$PYTHON_BIN" "$PIPELINE" unfold-graph \
    --gfa "$GFA" \
    --contigs "$VALIDATED/contigs.validated.tsv" \
    --gfa-inventory "$GFA_INVENTORY" \
    --expected-lengths "$VALIDATED/path_lengths.json" \
    --out-dir "$PATHBED_DIR" \
    --threads "$THREADS" \
    --compression gzip

run_step 05_prepare_inputs \
    "$PYTHON_BIN" "$INPUT_PREP" \
    --genomes "$VALIDATED/genomes.validated.tsv" \
    --contigs "$VALIDATED/contigs.validated.tsv" \
    --out-dir "$INPUTDIR" \
    --pathbed-dir "$PATHBED_DIR"

run_step 05_validate_prepared_inputs \
    "$PYTHON_BIN" "$PREPARED_INPUT_VALIDATOR" \
    --catalog "$CATALOG" \
    --genomes "$VALIDATED/genomes.validated.tsv" \
    --contigs "$VALIDATED/contigs.validated.tsv" \
    --pathbed-dir "$PATHBED_DIR"

run_step 06_make_cpgi_fasta \
    "$PYTHON_BIN" "$PIPELINE" make-cpgi-fasta \
    --catalog "$CATALOG" \
    --contigs "$VALIDATED/contigs.validated.tsv" \
    --hal "$HAL" \
    --log-dir "$WORKDIR/hal_fasta_logs" \
    --threads "$HAL_THREADS" \
    --docker-bin "$DOCKER_BIN" \
    --docker-image "$DOCKER_IMAGE" \
    --hal2fasta "$HAL2FASTA" \
    --missing-report "$WORKDIR/cpgi_fasta_missing.tsv"

run_step 07_hal_liftover \
    "$PYTHON_BIN" "$HAL_PIPELINE" \
    --catalog "$CATALOG" \
    --contigs "$VALIDATED/contigs.validated.tsv" \
    --hal "$HAL" \
    --out-dir "$HALDIR" \
    --threads "$HAL_THREADS" \
    --docker-bin "$DOCKER_BIN" \
    --docker-image "$DOCKER_IMAGE" \
    --hal-stats "$HAL_STATS" \
    --hal-liftover "$HAL_LIFTOVER"

run_step 08_build_features \
    "$PYTHON_BIN" "$PIPELINE" build-features-prod \
    --catalog "$CATALOG" \
    --out "$WORKDIR/features.jsonl.gz" \
    --excluded "$WORKDIR/features.excluded.tsv.gz" \
    --hal-psl-dir "$HALDIR/psl" \
    --threads "$THREADS" \
    --min-graph-cov 0.95 \
    --flank-bp 1000 \
    --flank-max-steps 32 \
    --hal-min-coverage 0.5 \
    --hal-min-identity 0.0 \
    --sv-contig-coordinate-base 0

run_step 08_index_features \
    "$PYTHON_BIN" "$PACKAGE_DIR/pancgi_features.py" \
    --features "$WORKDIR/features.jsonl.gz" --out "$WORKDIR/features.sqlite" \
    --anchor-k 3 --max-mid-anchors 8

run_step 09_cluster_loci \
    "$PYTHON_BIN" "$PIPELINE" cluster-loci-prod \
    --features "$WORKDIR/features.sqlite" \
    --out-locus "$RESULTDIR/locus_seed.tsv.gz" \
    --out-members "$RESULTDIR/locus_seed_members.tsv.gz"

run_step 10_polish_loci \
    "$PYTHON_BIN" "$PIPELINE" polish-loci \
    --features "$WORKDIR/features.sqlite" \
    --in-locus "$RESULTDIR/locus_seed.tsv.gz" \
    --in-members "$RESULTDIR/locus_seed_members.tsv.gz" \
    --out-locus "$RESULTDIR/locus_polished.tsv.gz" \
    --out-members "$RESULTDIR/locus_polished_members.tsv.gz" \
    --threads "$THREADS" \
    --feature-store-db "$WORKDIR/features.step04_polish.sqlite"

run_step 11_anchor_loci \
    "$PYTHON_BIN" "$PIPELINE" anchor-loci-primary \
    --features "$WORKDIR/features.sqlite" \
    --in-locus "$RESULTDIR/locus_polished.tsv.gz" \
    --in-members "$RESULTDIR/locus_polished_members.tsv.gz" \
    --out-locus "$RESULTDIR/locus_anchored.tsv.gz" \
    --out-members "$RESULTDIR/locus_anchored_members.tsv.gz" \
    --nonref-site-window-bp 100 \
    --hal-nonref-min-coverage 0.5 \
    --threads "$THREADS" \
    --parallel-backend process \
    --mp-start-method fork

run_step 12_cluster_alleles \
    "$PYTHON_BIN" "$PIPELINE" cluster-alleles-prod \
    --features "$WORKDIR/features.sqlite" \
    --locus-catalog "$RESULTDIR/locus_anchored.tsv.gz" \
    --locus-members "$RESULTDIR/locus_anchored_members.tsv.gz" \
    --out-allele "$RESULTDIR/allele_catalog.tsv.gz" \
    --out-allele-members "$RESULTDIR/allele_members.tsv.gz" \
    --identity 0.80 \
    --min-len-ratio 0.80 \
    --seq-backend parasail \
    --very-long-threshold 100000 \
    --very-long-backend external \
    --very-long-external-template "\"$PYTHON_BIN\" \"$PACKAGE_DIR/wfa_longalign_wrapper.py\" --seq1 \"{seq1}\" --seq2 \"{seq2}\" --log \"$RESULTDIR/wfa_external_calls.tsv\"" \
    --allele-rep-strategy asm_seq_medoid \
    --allele-cluster-mode allpairs_clique \
    --threads "$THREADS" \
    --parallel-backend process \
    --mp-start-method fork

run_step 13_strict_genotype \
    "$PYTHON_BIN" "$PIPELINE" strict-genotype-prod \
    --features "$WORKDIR/features.sqlite" \
    --locus-catalog "$RESULTDIR/locus_anchored.tsv.gz" \
    --locus-members "$RESULTDIR/locus_anchored_members.tsv.gz" \
    --allele-catalog "$RESULTDIR/allele_catalog.tsv.gz" \
    --allele-members "$RESULTDIR/allele_members.tsv.gz" \
    --catalog "$CATALOG" \
    --out-prefix "$RESULTDIR/pancgi" \
    --threads "$THREADS"

run_step 15_publish_results \
    "$PYTHON_BIN" "$PACKAGE_DIR/pancgi_results.py" \
    --internal-results "$RESULTDIR" --validated "$VALIDATED" \
    --inputs "$INPUTDIR" --out-dir "$OUTDIR/results"

printf '[%s] PanCGI completed: %s\n' "$(date '+%F %T')" "$OUTDIR"
