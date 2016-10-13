# zmine
A GPU miner for the Zcash cryptocurrency.

## Building

    mkdir libs
    cd libs

[Download and build libsodium](https://download.libsodium.org/doc/installation/index.html) in the `libs/ directory.

    cd ..
    make

## Running

    export LD_LIBRARY_PATH=`pwd`/libs/libsodium-1.0.11/src/libsodium/.libs/:/usr/local/cuda-7.5/lib64
    ./a.out
