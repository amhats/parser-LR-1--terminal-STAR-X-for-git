#ifndef PRINT_C
#define PRINT_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "SMQP_Support-AA.h"
#include "parser_type.h"
#include "parser_helper.h"
#include "parser_errno.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    # define __BEGIN_DECLS extern "C" {
    # define __END_DECLS }
#else
    # define __BEGIN_DECLS /* empty */
    # define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/**
 * @brief AltArr_t* array is converted to a string polynomial containing the variables
 *  and retuned as string or character array.
 * 
 * @param altarr - AltArr_t* or polynomial
 * @param var - variables
 * @param numVar - number of variables
 * @return char* - string 
 */
char* print_poly_to_string_variable(AltArr_t *altarr, char** var, int numVar );
char* print_poly_to_string_variable_z(AltArrZ_t *altarr, char** var, int numVar );

/**
 * @brief AltArr_t* array is converted to a string polynomial containing the variables
 *  and printed to file.
 * 
 * @param l - AltArr_t* or polynomial to print to file
 * @param var - set of variables in the polynomial
 * @param numVar - number of variables
 * @param filename - filename to print to
 * @param insertVarsAtTheBeginning - flag to add variable name at the benginning of the polynomials(eg if set it add [x,y,z]..).
 */
void print_poly_to_file(AltArr_t *l, char** var, int numVar, char *filename, int insertVarsAtTheBeginning);
void print_poly_to_file_z(AltArrZ_t *l, char** var, int numVar, char *filename, int insertVarsAtTheBeginning);

/**
 * @brief print AltArr_t* or polynomail with the variables.
 * 
 * @param l - AltArr* or polynomial to print
 * @param var - set of variables
 * @param numVar - number of variables
 */
void print_poly_to_terminal(AltArr_t *l, char** var, int numVar );
void print_poly_to_terminal_z(AltArrZ_t *l, char** var, int numVar );

/**
 * @brief print AltArr_t* or polynomail without the variables. Directly copied from
 * SMPQ_Support-AA, to avoid compiler warning adpoted a new function name.
 * 
 * @param aa - AltArr* or polynomial to print.
 */
void print_naked_AltArr_t_poly(AltArr_t* aa);
void print_naked_AltArr_t_poly_z(AltArrZ_t* aa);

/**
 * @brief Simply print the AltArr_t* with the variables, directly copied from 
 * SMPQ_Support-AA, to avoid compiler warning adpoted a new function name.
 * 
 * @param aa - array of polynomial
 * @param vars 
 * @param numvars 
 * @param msg 
 */
void print_naked_AltArr_t_poly2(AltArr_t* aa, char **vars, int numvars, char* msg);
void print_naked_AltArr_t_poly2_z(AltArrZ_t* aa, char **vars, int numvars, char* msg);

/**
 * @brief Print a polynomial term or type term (as defined in the type.h).
 * 
 * @param t - term as defined in type.h
 * @param vars - variable arrays
 * @param numvars - number of variables
 * @param message - include message to accompany the term being printed
 */
void print_term(term* t, char **vars, int numvars, char* message);
void print_term_z(term_z* t, char **vars, int numvars, char* message);

__END_DECLS

#endif