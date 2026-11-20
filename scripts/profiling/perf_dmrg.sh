#!/bin/bash

set -euo pipefail

if (( $# == 0 )); then
    echo "Usage: $0 <command> [arguments...]" >&2
    exit 2
fi

perf record -e cycles:u -g --call-graph dwarf,16384 -- "$@"
perf report
