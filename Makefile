all: program

BASE_PATH = $(shell pwd)

INCLUDES = -I$(BASE_PATH)/libs/libsodium-1.0.11/src/libsodium/include/ -I/usr/local/cuda/include/
NACL_LIBBIES = -L$(BASE_PATH)/libs/libsodium-1.0.11/src/libsodium/.libs/ -lsodium
NV_LIBBIES = -L/usr/local/cuda/lib64 -lcuda -lcudart
LIBBIES = $(NV_LIBBIES) $(NACL_LIBBIES)
OPTS = #--std=c99

program: clean csolv
	g++ -g $(OPTS) main.c cblake2b.o cutil.o cxor.o cfinal.o $(INCLUDES) $(LIBBIES)

csolv:
	nvcc -lineinfo -arch=compute_52 -code=sm_52 $(INCLUDES) -c cblake2b.cu -c cutil.cu -c cxor.cu -c cfinal.cu 

clean: 
	rm -rf *.o a.out
