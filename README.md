# zmine
A GPU miner for the Zcash cryptocurrency.

This project is licensed under the TGPPL. See LICENSE for more information.

This project is currently under development. It is not intended for end users before the 28th (zcashminers.org contest close). Please report any issues you find on this repository.

## Building

    mkdir libs
    cd libs

[Download and build libsodium](https://download.libsodium.org/doc/installation/index.html) in the `libs/ directory.

    cd ..
    make

## Running

    export LD_LIBRARY_PATH=`pwd`/libs/libsodium-1.0.11/src/libsodium/.libs/:/usr/local/cuda-7.5/lib64
    ./a.out

## TODO

* Eliminate ij_buf in favor of indices per hash.
  * Merge two sorted lists and each xor iteration.
  * Iterate (N-1) to check for distinct indices.
* Integrate Zcash testing API.
* Change lookup table to use same method as final functions.

[![Build Status](https://travis-ci.org/douggard/zmine.svg?branch=master)](https://travis-ci.org/douggard/zmine)
