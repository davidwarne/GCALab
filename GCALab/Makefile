#!/bin/make

CC = gcc
OPTS = -O2 -DNO_THREADS
#OPTS = -g
SRC =  GCALab.c GCALab_fio.c ../libGCA/GCA.c ../libMesh/mesh.c ../libMesh/vectorMath.c
OBJS = $(SRC:.c=.o)
INC = -I../libMesh/ -I../libGCA/ -I./
BIN = GCALab
LIBS = -lm -lpthread
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS)
