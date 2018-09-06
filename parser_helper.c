#include "parser_helper.h"


powervar *create_power_var(char *var, degrees_t exp){
    powervar *temp_powervar = (powervar*)malloc(sizeof(powervar));
    if(temp_powervar == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
	temp_powervar->var = (char*)malloc(strlen(var)*sizeof(char)+1);
    if(temp_powervar->var == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
	memcpy(temp_powervar->var, var, strlen(var)+1);
	temp_powervar->exp = exp;

	return temp_powervar;
}


term *create_term(int arraySize){
	term* temp_term = (term*)malloc(sizeof(term));
    if(temp_term == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    temp_term->exp = calloc(arraySize, sizeof(degrees_t));
    if(temp_term->exp == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    mpq_init(temp_term->coef);
    mpq_set_str(temp_term->coef, "1", 10);
	return temp_term;
}

term_z *create_term_z(int arraySize){
	term_z* temp_term = (term_z*)malloc(sizeof(term_z));
    if(temp_term == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    temp_term->exp = calloc(arraySize, sizeof(degrees_t));
    if(temp_term->exp == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }
    mpz_init(temp_term->coef);
    mpz_set_str(temp_term->coef, "1", 10);
	return temp_term;
}

degrees_t* generate_long_int_exponents(int arraySize ){
    degrees_t* temp_exponents = calloc(arraySize, sizeof(unsigned long int));
    if(!temp_exponents){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
    }

    return temp_exponents;
}

void free2DArray(char** arr, int size){
    if(arr != NULL){
        for(int i=0; i<size; i++)
            free(arr[i]);
        free(arr);
    }
}

void fill_term_exponent(char **orderedVarArray, degrees_t *exponentToFill, char* currentVar, degrees_t currentExp, int *arraySize){
    int pos = -1;
    //The for loop is a bottle neck, unfortunately it is required to reset the array because
    //of a wierd bug in the exponentToFill array. he problem is, instead of keeping exponentToFill
    //filled with zero, the program fills the certain index with random number if thier length is
    //greater than 3. This index (if they are greater than 3) will be filled with the value 1 or 33
    //instead of 0.
    // for(int i=0; i<numVariables; i++){
    //     exponentToFill[i] = 0;
    // }
    for(int i=0; i<(*arraySize); i++){
        if(strcmp(currentVar, orderedVarArray[i])==0){ //unsafe function
            pos = i;
            break;
        }
    }
    if(pos >= 0){
        if(exponentToFill[pos] > 0){//added may 4 , when input [x] x*x only returned x
            exponentToFill[pos] += currentExp;
        }else{
            exponentToFill[pos] = currentExp;
        }
    }else{
        // TODO: printf("new varaiable is discovered\n"); 
        //if the variable used unspecified, mostly the code resides here 
        //to enumirate through the array to update the exponents
        // printf("variable is not found ::)\n");
        // This condition is useful if the data structure linked-list
	}
}

int check_if_it_exists(char** vars, char* v, int numvar){
    if(vars){
        for(int i=0; i<numvar; ++i){
            if(strcmp(vars[i], v)==0){
                return 1;
            }
        }
        return 0;
    }
    return 0;
}


char** push_back_dynamic(char **varArray, int *numvar, char* var){
    if(!check_if_it_exists(varArray, var, *numvar)){
        (*numvar)++;
        if(!varArray){
            varArray = (char**)malloc((*numvar)*sizeof(char*));
        }else{
            varArray = realloc(varArray, (*numvar)*sizeof(char*));
        }
        varArray[(*numvar)-1] = (char*)malloc(strlen(var)+1*sizeof(char));
        memcpy(varArray[(*numvar)-1], var, strlen(var)+1);
    }

    return varArray;
}


void deep_degrees_cpy(degrees_t *dst, degrees_t *src, int size){
    for(int i=0; i<size; i++){
        dst[i] = src[i];
    }
}


void clear_term(term *tempTerm, int numVar){
	for(int i=0; i<numVar; i++){
        tempTerm->exp[i] = 0;
    }
}

void clear_term_z(term_z *tempTerm, int numVar){
	for(int i=0; i<numVar; i++){
        tempTerm->exp[i] = 0;
    }
}


void free_term(term* temp){
	free(temp->exp);
    mpq_clear(temp->coef);
	free(temp);
}

void free_term_z(term_z* temp){
	free(temp->exp);
    mpz_clear(temp->coef);
	free(temp);
}


static void addTermSafe_AA(AltArr_t* aa, degrees_t d, const ratNum_t coef){
	if (AA_SIZE(aa) >= aa->alloc) {
		aa->alloc *= 2;
		aa->elems = (AAElem_t*) realloc(aa->elems, sizeof(AAElem_t)*aa->alloc);
	}
	mpq_init(aa->elems[AA_SIZE(aa)].coef);
	mpq_set(aa->elems[AA_SIZE(aa)].coef, coef);
	aa->elems[AA_SIZE(aa)].degs = d;
	++(AA_SIZE(aa));
}

static void addTermSafe_AAZ(AltArrZ_t* aa, degrees_t d, const mpz_t coef){
	if (AA_SIZE(aa) >= aa->alloc) {
		aa->alloc *= 2;
		aa->elems = (AAZElem_t*) realloc(aa->elems, sizeof(AAZElem_t)*aa->alloc);
	}
	mpz_init(aa->elems[AA_SIZE(aa)].coef);
	mpz_set(aa->elems[AA_SIZE(aa)].coef, coef);
	aa->elems[AA_SIZE(aa)].degs = d;
	++(AA_SIZE(aa));
}


void add_packed_degree_term_to_smqp_aa(AltArr_t *aa, degrees_t *deg, ratNum_t coef, int numvar){
    
    if(numvar == 0){ //will dump error message in the console complaining EXP and NVARS are empty
        numvar = 1;
    }
    int *size = getExpOffsetArray(numvar);
    degrees_t d = 0;

    for(int i = 0; i < numvar; ++i){
        d |= (deg[i] << size[i]);
    }

    addTermSafe_AA(aa, d, coef);
    free(size);
}

void add_packed_degree_term_to_smzp_aaz(AltArrZ_t *aa, degrees_t *deg, mpz_t coef, int numvar){
    
    if(numvar == 0){ //will dump error message in the console complaining EXP and NVARS are empty
        numvar = 1;
    }
    int *size = getExpOffsetArray(numvar);
    degrees_t d = 0;

    for(int i = 0; i < numvar; ++i){
        d |= (deg[i] << size[i]);
    }

    addTermSafe_AAZ(aa, d, coef);
    free(size);
}


degrees_t* unpacked_degs(degrees_t packed_degs, int numvars){
    if(numvars == 0){ //will dump error message in the console complaining EXP and NVARS are empty
        numvars = 1;
    }
    degrees_t* __restrict__ oldmasks = getExpMaskArray(numvars);
    int* __restrict__ oldsizes = getExpOffsetArray(numvars);
    
    degrees_t* ret = (degrees_t*)calloc(numvars, sizeof(degrees_t));
    if(!ret){
        fprintf(stderr, "%s\n", "error callocing");
        return NULL;
    }
    for(int j=0; j<numvars; j++){
        ret[j] = GET_NTH_EXP(packed_degs, oldmasks[j], oldsizes[j]);
    } 
    
    free(oldmasks);
    free(oldsizes);
    return ret;
}

char* strip_comments(const char* buf, size_t size){
    char *newbuf = (char*)malloc((size)*sizeof(char));
    int i =0;
    int j = 0;
    while(i<size-1){
        if(buf[i]=='/'&&buf[i+1]=='*'){
            i = i+2;
            while(i<size){
                if(buf[i]=='*' && buf[i+1]=='/'){
                    break;
                }
                i++;
            }                
            i = i+2;            
        }
        newbuf[j] = buf[i];
        j++;
        i++;
    }
    
    return newbuf;
}