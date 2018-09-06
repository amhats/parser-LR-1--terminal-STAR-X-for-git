#COMPILE WITH STATIC LIBRARY
CC = gcc
CFLAGS = -lbpasparser -lgmp -lm 
# LIBRARY = ../lib/
TARGET = myparser
INCLUDE = ../inlucde/bpas_parser.h
TEST_INCLUDE = /home/user/maple2017/extern/include/
TEST_LIBS = -lmaplec -lmaple -lhf 

all:$(TARGET) test

$(TARGET): myparser.c  
	$(CC) -I$(INCLUDE) -I$(TEST_INCLUDE) -o $@ $^ $(CFLAGS) $(TEST_LIBS)

test:
	./myparser
clean:
	rm *.o myparser 




















# #COMPILE WITHOUT STATIC LIBRARY
# CC = gcc
# CFLAGS = -lgmp -lm 
# DFLAGS = -g
# TARGET = myparser
# INCLUDE = ../inlucde/bpas_parser.h
# SOURCE = ../src/bpas_parser.c ../SMQP_Support-AA.c ../parser_errno.c ../parser_helper.c ../parser_lex.yy.c ../parser_grammar.tab.c
# TEST_INCLUDE = /home/user/maple2017/extern/include/
# TEST_LIBS = -lmaplec -lmaple -lhf 
# all:$(TARGET)

# $(TARGET): myparser.c $(SOURCE) 
# 	$(CC) -I$(INCLUDE) -I$(TEST_INCLUDE) -o $@ $^ $(CFLAGS) $(DFLAGS) $(TEST_LIBS)



# clean:
# 	rm *.o myparser 
















# #STATIC LINKED LIBRARY NOT WORKING
# CC = gcc
# CFLAGS = -lgmp -lm 
# DFLAGS = -g
# TARGET = myparser
# LIBRARY = ../lib/
# INCLUDE = ../include/
# # SOURCE = ../src/bpas_parser.c ../SMQP_Support-AA.c ../parser_errno.c ../parser_helper.c ../parser_lex.yy.c ../parser_grammar.tab.c
# TEST_INCLUDE = /home/user/maple2017/extern/include/
# TEST_LIBS = -lmaplec -lmaple -lhf 
# all:$(TARGET)

# $(TARGET): myparser.c #$(SOURCE) 
# 	$(CC) -I$(INC) -I$(INCLUDE) -I$(TEST_INCLUDE) -o $@ $^ -L$(LIBRARY) $(CFLAGS) $(DFLAGS) $(TEST_LIBS) -lbpasparser -lparser



# clean:
# 	rm *.o myparser 



















