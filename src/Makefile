# Makefile for building C programs to do encoding, decoding, and clustering

ENC_SRC=well.c codebook.c main.c util.c
DEC_SRC=well.c codebook.c decoder.c util.c
CLUSTER_SRC=cluster.c util.c
TEST_SRC=codebook.c distortion.c lines.c pmf.c quantizer.c test.c util.c well.c

ENC_OBJ=$(ENC_SRC:.c=.o)
DEC_OBJ=$(DEC_SRC:.c=.o)
CLUSTER_OBJ=$(CLUSTER_SRC:.c=.o)
TEST_OBJ=$(TEST_SRC:.c=.o)

CC=gcc
RM=rm -f

CLFAGS=-O3 -Wall -fomit-frame-pointer -I../include -DLINUX
LDFLAGS=-lc -lm

%.o : %.c
	$(CC) $(CLFAGS) -c $<

all : encoder decoder cluster test

test : $(TEST_OBJ)
	$(CC) $(LDFLAGS) $(TEST_OBJ) -o test

cluster : $(CLUSTER_OBJ)
	$(CC) $(LDFLAGS) $(CLUSTER_OBJ) -o cluster

encoder : $(ENC_OBJ)
	$(CC) $(LDFLAGS) $(ENC_OBJ) -o encoder

decoder : $(DEC_OBJ)
	$(CC) $(LDFLAGS) $(DEC_OBJ) -o decoder

clean :
	$(RM) *.o encoder decoder cluster test
