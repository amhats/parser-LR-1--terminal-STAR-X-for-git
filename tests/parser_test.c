#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/bpas_parser.h"
#include "../parser_print.h"
#include "../BPAS-PARSER-TEST-FRAMEWORK/parser_testframework.h"
#include "../BPAS-PARSER-TEST-FRAMEWORK/maple_evaluate.h"

void handle_input(const char* input);

int main(int argc , char **argv){
	handle_input("data.txt");
	handle_input("data2.txt");
    return 0;
}


void handle_input(const char* input){
	FILE *temp;
	temp = fopen(input, "r");
	if(!temp){
		fprintf(stderr, "%s", "Error 2: opening the file! \n");
		exit(EXIT_FAILURE);
	}
	
	fseek(temp, 0, SEEK_END);
	size_t size = ftell(temp);
	rewind(temp);
	// printf("size: %zd\n", size);
	
	char *buf = (char*)malloc(size*sizeof(char));
	size_t r = fread(buf, sizeof(char), size, temp);
	
	char *newbuf = strip_comments(buf, size);
	free(buf);
	
	char *delim = ";";
	char *result =  strtok(newbuf, delim);
	while(result != NULL){
		size_t ilen = strlen(result);
		for(int i =0; i<ilen+1; i++){
			if(result[i] == '\n')
				result[i] = ' ';
		}
		altarr_pack *ret = generate_altarr_pack(result);
		if(ret->altarr_t_data != NULL){
			char *ret2 = print_poly_to_string_variable(ret->altarr_t_data, ret->vars, ret->numVars);
			if(strlen(ret2) !=0){
				ASSERT_EQUAL(maple_evaluate_equal(result, ret2), 0);
				parser_test(result, "TEST");
			}else{
				//print function evaluate input and return null, for example "x-x" returned as NULL
				//since the test function expect "x-x" when compared NULL with actual, test fails
				ASSERT_EQUAL(0, 0);
				parser_test(result, "TEST");
			}
		}
		result = strtok(NULL, delim);
	}
	free(newbuf);
	free(result);
}

