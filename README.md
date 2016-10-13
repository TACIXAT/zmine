# zmine
A GPU miner for the Zcash cryptocurrency.

## Building

Download and build libsodium in libs/ director.

    make

## Running

    export LD_LIBRARY_PATH=`pwd`/libs/libsodium-1.0.11/src/libsodium/.libs/:/usr/local/cuda-7.5/lib64
    ./a.out
