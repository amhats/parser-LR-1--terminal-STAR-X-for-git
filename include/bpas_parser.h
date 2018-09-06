#ifndef BPAS_PARSER_H
#define BPAS_PARSER_H

#define _GNU_SOURCE
#define __STDC_WANT_LIB_EXT2__ 1
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "../parser_helper.h"
#include "../parser_type.h"
#include "../parser_grammar.tab.h"
#include "../SMQP_Support-AA.h"
#include "../SMZP_Support.h"

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

char **create_dynamic_str_array(int arrSize, ...);
AltArr_t* generate_altarr(const char* poly_str);
AltArrZ_t* generate_altarr_z(const char* poly_str);
AltArr_t* generate_altarr_var_defined(const char* poly_str,  char** variables, int num_var);
altarr_pack* generate_altarr_pack(const char* poly_str);
altarr_pack* generate_altarr_pack_z(const char* poly_str);

__END_DECLS

#endif

