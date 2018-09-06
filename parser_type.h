#ifndef PARSER_TYPE_H
#define PARSER_TYPE_H

#include "SMQP_Support-AA.h"
#include "SMZP_Support.h"

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

typedef struct _term{
    ratNum_t coef;
	degrees_t *exp;
}term;

typedef struct _term_z{
	mpz_t coef;
	degrees_t *exp;
}term_z;

typedef struct _powervar{
	char *var;
	degrees_t exp;
}powervar;

typedef struct _altarr_pack{
	union{
		AltArr_t* altarr_t_data;
		AltArrZ_t* altarrz_t_data;
	}altarr_pack_type;
	char** vars;
	int numVars;
}altarr_pack;

typedef enum _gmp_type{
	GMP_MPQ = 0,
	GMP_MPZ
}gmp_type;

__END_DECLS

#endif
