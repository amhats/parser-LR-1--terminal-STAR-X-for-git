#COMPILE WITH STATIC LIBRARY
CC = gcc
CFLAGS = -lbpasparser -lgmp -lm 
LIBRARY = ../lib/
TARGET = test
INCLUDE = ../inlucde/bpas_parser.h
TEST_INCLUDE = /home/user/maple2017/extern/include/
TEST_LIBS = -lmaplec -lmaple -lhf 

all:$(TARGET)

$(TARGET): parser_test.c  
	$(CC) -I$(INCLUDE) -I$(TEST_INCLUDE) -o $@ $^ $(CFLAGS) $(TEST_LIBS)

clean:
	rm *.o test