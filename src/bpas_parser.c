#include "../include/bpas_parser.h"
#include "../parser_errno.h"

extern gmp_type g_coef_type;
extern void *altarr_data;
extern char** g_variables;
extern int g_num_variables;
extern void yy_scan_string(char*);
 

static int total_var_size(char** vars, const int numVars){
    int size = 0;
    for(int i=0; i<numVars; i++){
        size += strlen(vars[i]);
    }
    return size;
}

static char *create_var_list(char **vars, const int numVars){
    int size = total_var_size(vars, numVars);
    int t_size = size + (numVars-1) + 3;
    char *var_list = (char*)calloc(t_size, sizeof(char));
    
    if(!var_list){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    
    strncat(var_list, "[", strlen("["));
    for(int i=0; i<numVars; i++){
        strncat(var_list, vars[i], strlen(vars[i]));
        if(i != (numVars-1)){
            strncat(var_list, ",", strlen(","));
        }
    }
    strncat(var_list, "]", strlen("]"));
    strncat(var_list, "\0", strlen("\0"));
    return var_list;
}

char **create_dynamic_str_array(int arrSize, ...){
    char **temp_str_arr = (char**)malloc(sizeof(char*));
    if(!temp_str_arr){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    va_list args;
    va_start(args, arrSize);
    for(int i=0; i<arrSize; i++){
        char *temp_var = strdup(va_arg(args, char*));
        if(!temp_var){
            parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
            exit(PARSER_NOALLOC);
        }
        temp_str_arr[i] = (char*)malloc(sizeof(char)*strlen(temp_var)+1);
        strncpy(temp_str_arr[i], temp_var, strlen(temp_var)+1);
    }
    va_end(args);
    
    return temp_str_arr;
}

AltArr_t* generate_altarr(const char* poly_str){
    g_coef_type = GMP_MPQ;

    char *temp_polystr = strdup(poly_str);
    yy_scan_string(temp_polystr);
    yyparse();
    free(temp_polystr);
    return (AltArr_t*)altarr_data;
}

AltArrZ_t* generate_altarr_z(const char* poly_str){
    g_coef_type = GMP_MPZ;

    char *temp_polystr = strdup(poly_str);
    yy_scan_string(temp_polystr);
    yyparse();
    free(temp_polystr);
    return (AltArrZ_t*)altarr_data;
}

altarr_pack* generate_altarr_pack(const char* poly_str){
    g_coef_type = GMP_MPQ;

    char *temp_polystr = strdup(poly_str);
    altarr_pack *pack = (altarr_pack*)malloc(sizeof(altarr_pack));
    yy_scan_string(temp_polystr);
    yyparse();
    pack->altarr_pack_type.altarr_t_data = (AltArr_t*)altarr_data;
    pack->vars = g_variables;
    pack->numVars = g_num_variables;
    free(temp_polystr);
    return pack;
}

altarr_pack* generate_altarr_pack_z(const char* poly_str){
    g_coef_type = GMP_MPZ;

    char *temp_polystr = strdup(poly_str);
    altarr_pack *pack = (altarr_pack*)malloc(sizeof(altarr_pack));
    yy_scan_string(temp_polystr);
    yyparse();
    pack->altarr_pack_type.altarrz_t_data = (AltArrZ_t*)altarr_data;
    pack->vars = g_variables;
    pack->numVars = g_num_variables;
    free(temp_polystr);
    return pack;
}

AltArr_t* generate_altarr_var_defined(const char* poly_str, char** variables, int num_var){
    char *temp_var_list = create_var_list(variables, num_var);
    char *temp_poly = (char*)malloc(sizeof(char)*(strlen(temp_var_list) + strlen(poly_str)) + 1);
    if(!temp_poly){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }

    strncpy(temp_poly, temp_var_list, strlen(temp_var_list));
    strncpy(temp_poly, poly_str, strlen(poly_str)+1);

    yy_scan_string(temp_poly);
    yyparse();
    free(temp_var_list);
    free(temp_poly);
    
    return (AltArr_t*)altarr_data;
}









/* AltArr_t* generate_altarr_var_defined(const char* poly_str, const char** variables, int num_var){
    int var_size = 0;
    for(int i=0; i<num_var; i++){
        var_size += sizeof(variables[i]);
    }
    char *temp_var = (char*)calloc((num_var+var_size+2), sizeof(char));
    if(temp_var == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }
    strncat(temp_var, "[", sizeof("["));
    for(int i=0; i<num_var; i++){
        strncat(temp_var, variables[i], sizeof(variables[i]));
        if(i != (num_var-1)){
            strncat(temp_var, ",", sizeof(","));
        }
    }
    strncat(temp_var, "]", sizeof("]"));

    char *temp_str = (char*)malloc((sizeof(temp_var)+sizeof(poly_str)+1)*sizeof(char));
    if(temp_str == NULL){
        parser_error_with_errno_reason(PARSER_NOALLOC, __FILE__, __func__, __LINE__);
        exit(PARSER_NOALLOC);
    }

    char* temp_polystr = strdup(poly_str);
    strncat(temp_str, temp_var, sizeof(temp_var));
    strncat(temp_str, temp_polystr, sizeof(poly_str));
    strncat(temp_str, "\0", sizeof("\0"));

    yy_scan_string(temp_str);
    yyparse();

    free(temp_polystr);
    free(temp_var);
    free(temp_str); 
    
    return altarr_data;
} */