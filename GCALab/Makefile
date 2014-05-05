#!/bin/make

CC = gcc
OPTS = -O2 -DWITH_GRAPHICS
#OPTS = -g -DWITH_GRAPHICS 
#OPTS = -g
SRC =  GCALab.c GCALab_fio.c 
OBJS = $(SRC:.c=.o)
INC = -I../libMesh/ -I../libGCA/ -I ../libBitMap -I./ 
BIN = GCALab
LIBS = -lm -lpthread -lg -lglut -lGL -lGLU -L../libBitMap -lbitmap -L../libMesh -lmesh -L../libGCA -lGCA
#PROFILE = -g -pg


.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS)
