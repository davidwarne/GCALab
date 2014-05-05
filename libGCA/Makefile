#!/bin/make

NAME = GCA
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

SRC =  GCA.c
OBJS = $(SRC:.c=.o)
TESTSRC = test.c
INC = -I ../libMesh 
BIN = Test
LIBS = -lm -L ../libMesh -lmesh 
#PROFILE = -g -pg


.SUFFIXES: .c .o

all:
	make $(STATIC)
	make $(SHARED)
	make $(BIN)

$(STATIC): $(OBJS)
	$(AR) $(AROPTS) $@ $?

$(SHARED): $(OBJS)
	$(CC) -shared -o $@ $? $(LIBS)

.c.o:
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) $(TESTSRC)  -o $(BIN) $(LIBS) $(INC) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS) $(SHARED) $(STATIC)
