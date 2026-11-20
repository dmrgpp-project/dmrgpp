#!/bin/bash

set -euo pipefail

computeProcsForMake()
{
    local nprocs

    if [[ $(uname -s) == Darwin ]]; then
        nprocs=$(sysctl -n hw.ncpu)
    else
        nprocs=$(nproc)
    fi

    if (( nprocs > 3 )); then
        echo $((nprocs - 2))
    else
        echo 1
    fi
}

computeProcsForCallgrind()
{
    local nprocs

    if [[ $(uname -s) == Darwin ]]; then
        nprocs=$(sysctl -n hw.ncpu)
    else
        nprocs=$(nproc)
    fi

    if (( nprocs > 4 )); then
        echo 4
    else
        echo "${nprocs}"
    fi
}

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

case ${1:-} in
"")
    build_dir="${source_dir}/build"
    mkdir -p "${build_dir}"
    cmake -S "${source_dir}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release
    cmake --build "${build_dir}" -j "$(computeProcsForMake)"
    ;;
callgrind)
    (( $# == 1 )) || { echo "Usage: $0 [callgrind]" >&2; exit 2; }
    build_dir="${source_dir}/build-callgrind"
    cmake -S "${source_dir}" -B "${build_dir}" \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DBUILD_TESTING=OFF
    cmake --build "${build_dir}" --target dmrg -j "$(computeProcsForCallgrind)"
    ;;
*)
    echo "Usage: $0 [callgrind]" >&2
    exit 2
    ;;
esac
