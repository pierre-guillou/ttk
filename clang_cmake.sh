#!/usr/bin/env bash

function conf_builddir_libcxx {
    local builddir=$1
    rm -rf $builddir
    mkdir $builddir
    CXX=clang++ \
       CC=clang \
       LD=ld.lld \
       CXXFLAGS="-nostdinc++ -I/usr/include/c++/v1" \
       LDFLAGS="-stdlib=libc++" \
       cmake \
       -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE \
       -B $builddir
}

function build_macOS {
    local builddir='build_macOS'
    conf_builddir_libcxx $builddir
    cmake --build $builddir
}

function check_macOS {
    local builddir='build_macOS'
    conf_builddir_libcxx $builddir
    jq -r '.[] | select(.file | contains("'$builddir'") | not) | .file' \
       $builddir/compile_commands.json | xargs -n1 -P 4 clang-check -p $builddir
}

function conf_builddir {
    local builddir=$1
    rm -rf $builddir
    mkdir $builddir
    CXX=clang++ \
       CC=clang \
       cmake \
       -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE \
       -B $builddir
}

function build_clang {
    local builddir='build_clang'
    conf_builddir $builddir
    cmake --build $builddir
}

function check_clang {
    local builddir='build_clang'
    conf_builddir $builddir
    jq -r '.[] | select(.file | contains("'$builddir'") | not) | .file' \
       $builddir/compile_commands.json | xargs -n1 -P 4 clang-check -p $builddir
}

function check_tidy {
    local builddir='build_clang'
    conf_builddir $builddir
    jq -r '.[] | select(.file | contains("'$builddir'") | not) | .file' \
       $builddir/compile_commands.json \
        | xargs -n1 -P 4 clang-tidy -p $builddir -header-filter="ttk/core"
}

function check_sa {
    local builddir='build_clang'
    conf_builddir $builddir
    jq -r '.[] | select(.file | contains("'$builddir'") | not) | .file' \
       $builddir/compile_commands.json \
        | xargs -n1 -P 4 \
                clang-check \
                -p $builddir -analyze -extra-arg -Xclang -extra-arg -analyzer-output=text
}

function scan_build {
    local builddir='build_scan'
    rm -rf $builddir
    mkdir $builddir
    scan-build cmake \
               -DCMAKE_BUILD_TYPE=Debug \
               -B $builddir
    scan-build make -C $builddir
}

$1
