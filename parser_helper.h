#ifndef HELPER_C
#define HELPER_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "SMQP_Support-AA.h"
#include "SMZP_Support.h"
#include "parser_type.h"
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
 * @brief Check if variable exist in the array 
 * 
 * @param vars - arrays of variables
 * @param v - variable to search for
 * @param numvar - number of variables
 * @return int - 1 or 0(does not exist)
 */
int check_if_it_exists(char** vars, char* v, int numvar);

/**
 * @brief Unpack a packed array debugging and printing purpose.
 * 
 * @param packed_degs 
 * @param numvars 
 * @return degrees_t* 
 */
degrees_t* unpacked_degs(degrees_t packed_degs, int numvars);

/**
 * Removes any string surrounded by / followed by *
 * followed by the string followed by * and / from a file or
 * a string input.
 * 
 * @param buf buffer containing the string.
 * @param size size of the buffer
 * @return buffer of string with out the comment
 */
char* strip_comments(const char* buf, size_t size);

/**
 * @brief free string of arrays from dynamically located memory.
 * 
 * @param arr 
 * @param size 
 */
void free2DArray(char** arr, int size);

/**
 * @brief Set the term type to NULL or Zero (struct defined in type.h).
 * 
 * @param tempTerm 
 * @param numVar 
 */
void clear_term(term *tempTerm, int numVar);
void clear_term_z(term_z *tempTerm, int numVar);

/**
 * @brief Free or release the term type (struct defined in type.h).
 * 
 * @param temp 
 */
void free_term(term* temp);
void free_term_z(term_z* temp);

/**
 * @brief Create a power variable type (struct defined in type.h)
 * 
 * @param var 
 * @param exp 
 * @return powervar* 
 */
powervar *create_power_var(char *var, degrees_t exp);
/**
 * @brief Create a term type (struct defined in type.h)
 * Refere type.h
 * 
 * @param arraySize 
 * @return term* 
 */
term* create_term(int arraySize);
term_z* create_term_z(int arraySize);

/**
 * @brief Generate give size of long int array.
 * 
 * @param arraySize 
 * @return degrees_t* 
 */
degrees_t* generate_long_int_exponents(int arraySize);

/**
 * @brief The function takes a globle term type, for the duration of scanning the production
 * term as defined in parser_grammar.y, it fills the global term with currently scanned variable
 * at its proper order. If the variable already exists it adds the current expoenent variable to 
 * the existing variable.
 * 
 * @param orderedVarArray  - arrays of variables ordered, 
 * @param exponentToFill  - expoenent to fill to orderedVarArray
 * @param currentVar - current variable scanned
 * @param currentExp - current exponent scanned
 * @param arraySize - array size of the variables
 */
void fill_term_exponent(char **orderedVarArray, degrees_t *exponentToFill, char* currentVar, degrees_t currentExp, int *arraySize);
/**
 * @brief create arrays of string dynamically. The arrays size is adjusted as a new
 *  item is introduced. 
 *      Note: This function requires varArray to be set to NULL if it is empty.
 *      If an empty varArray variable is not set to NULL, segmentation falut is eminent.
 * 
 * @param varArray 
 * @param numvar 
 * @param var 
 * @return char** 
 */
char** push_back_dynamic(char **varArray, int *numvar, char* var);
/**
 * @brief Deep copy source array to destination.
 * 
 * @param dst 
 * @param src 
 * @param size 
 */
void deep_degrees_cpy(degrees_t *dst, degrees_t *src, int size);

/**
 * @brief pack the raw degree's of a term to 2-words or 4-words before adding to
 * AltArr_t* .
 * 
 * @param aa 
 * @param deg 
 * @param coef 
 * @param numvar 
 */
void add_packed_degree_term_to_smqp_aa(AltArr_t *aa, degrees_t *deg, ratNum_t coef, int numvar);
void add_packed_degree_term_to_smzp_aaz(AltArrZ_t *aa, degrees_t *deg, mpz_t coef, int numvar);

//The following function is not for SMQP_Support-AA, currently not defined in SMQP_Support-AA
/**
 * @brief Safely add term to AltArr_t* , the function is not defined in the dependency library
 * SMQP_Support-AA currently.
 * 
 * @param aa - AltArr_t* type to add to.
 * @param d - degree's of the term to be added
 * @param coef - Coefficient of the term
 */
static void addTermSafe_AA(AltArr_t* aa, degrees_t d, const ratNum_t coef);
static void addTermSafe_AAZ(AltArrZ_t* aa, degrees_t d, const mpz_t codef);

__END_DECLS

#endif