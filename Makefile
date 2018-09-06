#STATIC LIBRARY - SIMPLE
CC = gcc
BISON = bison 
FLEX = flex 
DFLAGS = 
STATICLIB = lib/libbpasparser.a
BINARY = parser
# TEST_INCLUDE = /home/user/maple2017/extern/include/
# TEST_LIBS = -lmaplec -lmaple -lhf
# OBJ = SMQP_Support-AA.o parser_errno.o parser_helper.o parser_print.o parser_lex.yy.o parser_grammar.tab.o bpas_parser.o
# INCLUDE = include/bpas_parser.h SMQP_Support-AA.h SMZP_Support.h parser_errno.h parser_helper.h parser_print.h parser_grammar.tab.h 
# SOURCE = src/bpas_parser.c SMQP_Support-AA.c SMZP_Support.c parser_errno.c parser_helper.c parser_print.c parser_lex.yy.c parser_grammar.tab.c
# OBJ = parser_main.o SMQP_Support-AA.o parser_errno.o parser_helper.o parser_print.o parser_lex.yy.o parser_grammar.tab.o bpas_parser.o
INCLUDE = SMQP_Support-AA.h SMZP_Support.h parser_errno.h parser_helper.h parser_print.h parser_grammar.tab.h 
SOURCE = parser_main.c SMQP_Support-AA.c SMZP_Support.c parser_errno.c parser_helper.c parser_print.c parser_lex.yy.c parser_grammar.tab.c 
BISON_HEADER_SRC = parser_grammar.tab.h parser_grammar.tab.c
FLEX_SRC = parser_lex.yy.c

SUBDIR = test

all: $(BISON_HEADER_SRC) $(FLEX_SRC) $(BINARY) 

debug: DFLAGS += -DPARSER_DEBUG -g
debug: $(BISON_HEADER_SRC) $(FLEX_SRC) $(BINARY) 

$(FLEX_SRC): parser_scanner.l
	$(FLEX) -o $@ $<

$(BISON_HEADER_SRC): parser_grammar.y
	$(BISON) -vd $<

# $(STATICLIB): $(OBJ)
# 	ar -cvqs $@ $^

$(BINARY): $(SOURCE) 
	$(CC) -o $@ $^ $(DFLAGS) -lgmp -lreadline -lhistory

$(OBJ) :  $(SOURCE)
	$(CC) -c $^ 



# $(OBJ) :  $(SOURCE)
# 	$(CC) -I$(TEST_INCLUDE) -c $^ $(TEST_LIBS)

.PHONY:

test:
	cd tests  && $(MAKE) && ./test
example:
	cd examples && $(MAKE)

clean:
	rm *.o *.tab.o *.yy.o *.gch *.a *.so *.output lib/*.a *.tab.* *.yy.* examples/myparser tests/test parser

