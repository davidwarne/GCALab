#!/bin/make

CC = gcc
OPTS = -O2
#OPTS = -g
SRC =  test.c GCA.c ../libMesh/mesh.c ../libMesh/vectorMath.c
OBJS = $(SRC:.c=.o)
INC = -I../libMesh/ -I./
BIN = Test
LIBS = -lm -lpthread 
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) *.o
