#!/bin/make

NAME = mesh
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

SRC =  mesh.c vectorMath.c
OBJS = $(SRC:.c=.o)
TESTSRC = tester.c
INC = 
BIN = Testmesh
LIBS = -lm 
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
	$(CC) $(OPTS) $(PROFILE)  $(OBJS) $(TESTSRC)  -o $(BIN) $(LIBS) 
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) $(OBJS) $(SHARED) $(STATIC)
