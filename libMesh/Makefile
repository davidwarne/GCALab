#!/bin/make

CC = gcc
OPTS = -O2 -DNO_THREADS
#OPTS = -g
SRC =  tester.c mesh.c vectorMath.c
OBJS = $(SRC:.c=.o)
INC = 
BIN = Testmesh
LIBS = -lm 
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS)
