#!/bin/make

NAME = bitmap
TARGET = lib$(NAME)
STATIC = $(TARGET).a
SHARED = $(TARGET).so

# compiler options
CC = gcc
OPTS = -O2 -fPIC -DNO_THREADS
#OPTS = -g

#archive options
AR = ar
AROPTS = -rcvs

SRC =  BitMapFile.c BitMapReader.c BitMapWriter.c 
OBJS = $(SRC:.c=.o)
LIBS = -lm  
#PROFILE = -g -pg


.SUFFIXES: .c .o

all:
	make $(STATIC)
	make $(SHARED)

$(STATIC): $(OBJS)
	$(AR) $(AROPTS) $@ $?

$(SHARED): $(OBJS)
	$(CC) -shared -o $@ $? $(LIBS)

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ 

clean:
	set nonomatch; rm -f $(OBJS) $(SHARED) $(STATIC)
