# zmine
A GPU miner for the Zcash cryptocurrency.

This project is licensed under the TGPPL. See LICENSE for more information.

This project is currently under development. If you are not a developer planning on contributing back to this repository, please wait until after the zcashminers.org contest closes (10/28) to ask me for help building or running the code.

## Building

    mkdir libs
    cd libs

[Download and build libsodium](https://download.libsodium.org/doc/installation/index.html) in the `libs/ directory.

    cd ..
    make

## Running

    export LD_LIBRARY_PATH=`pwd`/libs/libsodium-1.0.11/src/libsodium/.libs/:/usr/local/cuda-7.5/lib64
    ./a.out
