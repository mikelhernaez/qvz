# Makefile for building C programs to do encoding, decoding, and clustering

SRC=well.c codebook.c main.c util.c lines.c quantizer.c pmf.c distortion.c qv_stream.c qv_compressor.c arith.c os_stream.c
CLUSTER_SRC=cluster.c util.c

OBJ=$(SRC:.c=.o)
CLUSTER_OBJ=$(CLUSTER_SRC:.c=.o)

CC=gcc
RM=rm -f

CFLAGS=-O3 -Wall -I../include -DLINUX
LDFLAGS=-lc -lm -lrt

%.o : %.c
	$(CC) $(CFLAGS) -c $<

all : qvz

cluster : $(CLUSTER_OBJ)
	$(CC) $(CLUSTER_OBJ) -o cluster $(LDFLAGS)

qvz : $(OBJ)
	$(CC)  $(OBJ) -o qvz $(LDFLAGS)

debug : CFLAGS += -DDEBUG -ggdb
debug : qvz

clean :
	$(RM) *.o qvz
