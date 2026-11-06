#!/usr/bin/env bash
# Profile a DMRG++ input with Valgrind Massif and render a text report.
#
# Usage:
#   ./massif.sh [--dmrg PATH] [--output-dir DIR] [--] INPUT [DMRG_ARGUMENT...]
#
# DMRG_EXECUTABLE and MASSIF_OUTPUT_DIR provide the same overrides as the
# corresponding options.  The executable defaults to a configured build tree.

set -euo pipefail

usage()
{
    cat <<'EOF'
Usage: massif.sh [--dmrg PATH] [--output-dir DIR] [--] INPUT [DMRG_ARGUMENT...]

Profile INPUT with DMRG++ under Valgrind Massif.  Additional arguments are
passed to dmrg after "-f INPUT".  Results are written to a new subdirectory of
MASSIF_OUTPUT_DIR (default: massif-results) and its path is printed on success.

Options:
  --dmrg PATH          DMRG++ executable (or set DMRG_EXECUTABLE)
  --output-dir DIR     Parent directory for results (or set MASSIF_OUTPUT_DIR)
  -h, --help           Show this help text
EOF
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
dmrg_executable=${DMRG_EXECUTABLE:-}
output_dir=${MASSIF_OUTPUT_DIR:-"${script_dir}/massif-results"}

while (( $# > 0 )); do
    case $1 in
    --dmrg)
        (( $# >= 2 )) || { echo "error: --dmrg requires a path" >&2; exit 2; }
        dmrg_executable=$2
        shift 2
        ;;
    --output-dir)
        (( $# >= 2 )) || { echo "error: --output-dir requires a path" >&2; exit 2; }
        output_dir=$2
        shift 2
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    --)
        shift
        break
        ;;
    -*)
        echo "error: unknown option: $1" >&2
        usage >&2
        exit 2
        ;;
    *)
        break
        ;;
    esac
done

(( $# > 0 )) || { echo "error: an input file is required" >&2; usage >&2; exit 2; }
input_file=$1
shift

[[ -f $input_file ]] || { echo "error: input file not found: $input_file" >&2; exit 2; }
command -v valgrind >/dev/null || { echo "error: valgrind is not installed" >&2; exit 127; }
command -v ms_print >/dev/null || { echo "error: ms_print is not installed" >&2; exit 127; }

if [[ -z $dmrg_executable ]]; then
    for candidate in \
        "${script_dir}/build-callgrind/dmrg/dmrg" \
        "${script_dir}/build/dmrg/dmrg" \
        "${script_dir}/dmrg/dmrg"; do
        if [[ -x $candidate ]]; then
            dmrg_executable=$candidate
            break
        fi
    done
fi

[[ -n $dmrg_executable && -x $dmrg_executable ]] || {
    echo "error: cannot execute DMRG++ binary: ${dmrg_executable:-not found}" >&2
    echo "       configure one with --dmrg PATH or DMRG_EXECUTABLE" >&2
    exit 2
}

input_name=$(basename -- "$input_file")
input_stem=${input_name%.*}
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
result_dir="${output_dir}/${input_stem}-${timestamp}-$$"
massif_output="${result_dir}/massif.out"
report="${result_dir}/massif.txt"
valgrind_log="${result_dir}/valgrind.log"
command_file="${result_dir}/command.txt"

mkdir -p -- "$result_dir"
printf '%q ' valgrind --tool=massif --time-unit=ms --detailed-freq=1 \
    --max-snapshots=200 --massif-out-file="$massif_output" --log-file="$valgrind_log" \
    "$dmrg_executable" -f "$input_file" "$@" > "$command_file"
printf '\n' >> "$command_file"

valgrind --tool=massif \
    --time-unit=ms \
    --detailed-freq=1 \
    --max-snapshots=200 \
    --massif-out-file="$massif_output" \
    --log-file="$valgrind_log" \
    -- "$dmrg_executable" -f "$input_file" "$@"

ms_print "$massif_output" > "$report"

printf 'Massif results: %s\n' "$result_dir"
printf '  Report:       %s\n' "$report"
printf '  Raw profile:  %s\n' "$massif_output"
printf '  Valgrind log: %s\n' "$valgrind_log"
