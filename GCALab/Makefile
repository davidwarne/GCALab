#!/bin/make

CC = gcc
OPTS = -O2 -DNO_THREADS -DWITH_GRAPHICS
#OPTS = -g -DWITH_GRAPHICS -DNO_THREADS
SRC =  GCALab.c GCALab_fio.c ../libGCA/GCA.c ../libMesh/mesh.c ../libMesh/vectorMath.c ../libBitMap/BitMapFile.c ../libBitMap/BitMapWriter.c
OBJS = $(SRC:.c=.o)
INC = -I../libMesh/ -I../libGCA/ -I ../libBitMap -I./ 
BIN = GCALab
LIBS = -lm -lpthread -lg -lglut -lGLU
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS)
