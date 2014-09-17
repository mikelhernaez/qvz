# Makefile for building C programs to do encoding, decoding, and clustering

SRC=well.c codebook.c main.c util.c lines.c quantizer.c pmf.c distortion.c qv_stream.c qv_compressor.c arith.c os_stream.c
CLUSTER_SRC=cluster.c util.c

OBJ=$(SRC:.c=.o)
CLUSTER_OBJ=$(CLUSTER_SRC:.c=.o)

CC=gcc
RM=rm -f

CLFAGS=-O3 -Wall -fomit-frame-pointer -I../include -DLINUX -ggdb
LDFLAGS=-lc -lm -lrt

%.o : %.c
	$(CC) $(CLFAGS) -c $<

all : qvz

cluster : $(CLUSTER_OBJ)
	$(CC) $(LDFLAGS) $(CLUSTER_OBJ) -o cluster

qvz : $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o qvz

clean :
	$(RM) *.o qvz
