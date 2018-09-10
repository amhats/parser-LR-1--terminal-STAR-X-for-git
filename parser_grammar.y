%{
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser_type.h"
#include "parser_helper.h"
#include "parser_debug_helper.h"
#include "parser_errno.h"
#include "parser_print.h"

#define DEFAULT_AA_SIZE 10

//coefficient default type (enum typedef gmp_type)
gmp_type g_coef_type = GMP_MPQ;
//Argument option variables
int arg_test = 0;
char* arg_filename = NULL;
char* current_input = NULL; //grab current input for testing purpose


extern int yylex();
void yyerror(char*);

void *altarr_data = NULL;

int is_negative_var = 0;
int is_new_var = 0; // new variable encountered on var unspecified poly input
int is_all_var_defined = 0; // user defined all the variables in the poly

char** g_variables = NULL;
int g_num_variables = 0;
%}

%union{
	char *string_type;
	int integer_value_type;
	powervar *powervar_type;
	void *term_type;
	void *AltArr_t_type;
};

%token <string_type> RATNUM VAR NUM 
%token MULTIPLY PLUS MINUS L_BRACE R_BRACE LS_BRACE RS_BRACE DIVIDE POWER COMMA UMINUS


%left MINUS PLUS
%left MULTIPLY DIVIDE
%left L_BRACE R_BRACE
%right POWER
%left UMINUS


%type <string_type> variable coef exponent polynomial
%type <powervar_type> powerVariable
%type <string_type> other
%type <term_type>  term
%type <AltArr_t_type> poly

%start polynomial

%%
polynomial		: polynomial poly 								{
																	// altarr_data = $2;
																	// return PARSER_SUCCESS;	
																	if(g_coef_type == GMP_MPQ){
																		print_poly_to_terminal($2, g_variables, g_num_variables);	
																		print_naked_AltArr_t_poly($2);
																	}else{
																		print_poly_to_terminal_z($2, g_variables, g_num_variables);	
																	}	

																	free2DArray(g_variables, g_num_variables);
																	g_variables = NULL;
																	g_num_variables = 0;
																}
				| polynomial other 								{}
                |                                       		{}
				;
poly			:MINUS poly %prec UMINUS						{
																	if(g_coef_type == GMP_MPQ){
																		if(is_new_var){
																			expandNumVars_AA($2, g_num_variables);
																		}	
																		negatePolynomial_AA($2);
																		mergeSortPolynomial_AA($2);
																		$$ = $2;
																	}else if(g_coef_type == GMP_MPZ){
																		if(is_new_var){
																			expandNumVars_AAZ($2, g_num_variables);
																		}	
																		negatePolynomial_AAZ($2);
																		mergeSortPolynomial_AAZ($2);
																		$$ = $2;
																	}
																	
																	#ifdef PARSER_DEBUG
																		printf_magenta("poly: $2: -poly");
																		if(g_coef_type == GMP_MPQ)
																			print_poly_to_terminal($2, g_variables, g_num_variables);
																		else if(g_coef_type == GMP_MPQ)
																			print_poly_to_terminal_z($2, g_variables, g_num_variables);
																		printf_magenta("poly: $$: -poly");
																		if(g_coef_type == GMP_MPQ)
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		else if(g_coef_type == GMP_MPZ)
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																	#endif
																}
				|poly MINUS poly 								{
																	if(g_coef_type == GMP_MPQ){
																		if(is_new_var){
																			expandNumVars_AA($1, g_num_variables);
																		}
																		AltArr_t *temp_aa_result = makePolynomial_AA(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AA($1);
																		printf("DELETE ME $1>>>>>> ");
																		print_poly_to_terminal($1, g_variables, g_num_variables);
																		negatePolynomial_AA($3);
																		mergeSortPolynomial_AA($3);
																		printf("DELETE ME $3>>>>>> ");
																		print_poly_to_terminal($3, g_variables, g_num_variables);
																		temp_aa_result = addPolynomials_AA($1, $3, g_num_variables);
																		printf("DELETE ME temp_aa_result >>>>>> ");
																		print_naked_AltArr_t_poly2(temp_aa_result, g_variables, g_num_variables, "***");
																		print_naked_AltArr_t_poly(temp_aa_result);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$: poly-poly");
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AA($1);
																		freePolynomial_AA($3);
																	}else if(g_coef_type == GMP_MPZ){
																		if(is_new_var){
																			expandNumVars_AAZ($1, g_num_variables);
																		}
																		AltArrZ_t *temp_aa_result = makePolynomial_AAZ(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AAZ($1);
																		negatePolynomial_AAZ($3);
																		mergeSortPolynomial_AAZ($3);
																		temp_aa_result = addPolynomials_AAZ($1, $3, g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal_z($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal_z($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$: poly-poly");
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AAZ($1);
																		freePolynomial_AAZ($3);
																	}
																}
				|poly PLUS poly									{
																	if(g_coef_type == GMP_MPQ){
																		if(is_new_var){
																			expandNumVars_AA($1, g_num_variables);
																		}
																		AltArr_t *temp_aa_result = makePolynomial_AA(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AA($1);
																		mergeSortPolynomial_AA($3);
																		temp_aa_result = addPolynomials_AA($1, $3, g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$: poly+poly");
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AA($1);
																		freePolynomial_AA($3);
																	} else if(g_coef_type == GMP_MPZ){
																		if(is_new_var){
																			expandNumVars_AAZ($1, g_num_variables);
																		}
																		AltArrZ_t *temp_aa_result = makePolynomial_AAZ(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AAZ($1);
																		mergeSortPolynomial_AAZ($3);
																		temp_aa_result = addPolynomials_AAZ($1, $3, g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal_z($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal_z($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$: poly+poly");
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AAZ($1);
																		freePolynomial_AAZ($3);
																	}
																}
				| poly MULTIPLY poly							{
																	if(g_coef_type == GMP_MPQ){
																		if(is_new_var){
																			expandNumVars_AA($1, g_num_variables);
																		}
																		AltArr_t *temp_aa_result = makePolynomial_AA(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AA($1);
																		mergeSortPolynomial_AA($3);
																		temp_aa_result = multiplyPolynomials_AA($1, $3, g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$ :poly * poly");
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AA($1);
																		freePolynomial_AA($3);
																	}else if(g_coef_type == GMP_MPZ){
																		if(is_new_var){
																			expandNumVars_AAZ($1, g_num_variables);
																		}
																		AltArrZ_t *temp_aa_result = makePolynomial_AAZ(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AAZ($1);
																		mergeSortPolynomial_AAZ($3);
																		temp_aa_result = multiplyPolynomials_AAZ($1, $3, g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: $1: poly");
																			print_poly_to_terminal_z($1, g_variables, g_num_variables);
																			printf_magenta("poly: $3: poly");
																			print_poly_to_terminal_z($3, g_variables, g_num_variables);
																			printf_magenta("poly: $$ :poly * poly");
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AAZ($1);
																		freePolynomial_AAZ($3);
																	}
																}
				| L_BRACE poly R_BRACE							{
																	$$ = $2;

																	#ifdef PARSER_DEBUG
																		printf_magenta("poly: (poly)");
																		if(g_coef_type == GMP_MPQ)
																			print_poly_to_terminal($2, g_variables, g_num_variables);
																		else if(g_coef_type == GMP_MPZ)
																			print_poly_to_terminal_z($2, g_variables, g_num_variables);
																		printf_magenta("poly: $$: (poly)");
																		if(g_coef_type == GMP_MPQ)
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		else if(g_coef_type == GMP_MPZ)
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																	#endif
																}
				|L_BRACE poly R_BRACE POWER exponent			{
																	if(g_coef_type == GMP_MPQ){
																		AltArr_t *temp_aa_result = makePolynomial_AA(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AA($2);
																		temp_aa_result = exponentiatePoly_AA($2, atoi($5), g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: poly^");
																			print_poly_to_terminal($2, g_variables, g_num_variables);
																			printf("NUM: %s\n", $5);
																			printf_magenta("poly: $$: poly^NUM");
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AA($2);
																	}else if(g_coef_type == GMP_MPZ){
																		AltArrZ_t *temp_aa_result = makePolynomial_AAZ(DEFAULT_AA_SIZE, g_num_variables);
																		mergeSortPolynomial_AAZ($2);
																		temp_aa_result = exponentiatePoly_AAZ($2, atoi($5), g_num_variables);
																		$$ = temp_aa_result;

																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: poly^");
																			print_poly_to_terminal_z($2, g_variables, g_num_variables);
																			printf("NUM: %s\n", $5);
																			printf_magenta("poly: $$: poly^NUM");
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																		#endif
																		freePolynomial_AAZ($2);
																	}
																}
				|term											{
																	if(g_coef_type == GMP_MPQ){
																		AltArr_t *temp_aa = makePolynomial_AA(DEFAULT_AA_SIZE, g_num_variables);
																		term *temp_term = $1;
																		add_packed_degree_term_to_smqp_aa(temp_aa, temp_term->exp, temp_term->coef, g_num_variables);
																		$$ = temp_aa;
																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: term");
																			print_term($1, g_variables, g_num_variables, "");
																			printf_magenta("poly: $$: term");
																			print_naked_AltArr_t_poly($$);
																			print_poly_to_terminal($$, g_variables, g_num_variables);
																		#endif
																		free_term(temp_term);
																	}else if(g_coef_type == GMP_MPZ){
																		AltArrZ_t *temp_aa = makePolynomial_AAZ(DEFAULT_AA_SIZE, g_num_variables);	
																		term_z *temp_term_z = $1;
																		add_packed_degree_term_to_smzp_aaz(temp_aa, temp_term_z->exp, temp_term_z->coef, g_num_variables);
																		$$ = temp_aa;
																		#ifdef PARSER_DEBUG
																			printf_magenta("poly: term");
																			print_term_z($1, g_variables, g_num_variables, "");
																			printf_magenta("poly: $$: term");
																			print_naked_AltArr_t_poly_z($$);
																			print_poly_to_terminal_z($$, g_variables, g_num_variables);
																		#endif
																		free_term_z(temp_term_z);
																	}
																}
				;
term			: coef											{
																	if(g_coef_type == GMP_MPQ){
																		term *local_term = create_term(g_num_variables);
																		mpq_t temp_coef;
																		mpq_init(temp_coef);
																		mpq_set_str(temp_coef, $1, 10);
																		mpq_mul(local_term->coef, local_term->coef, temp_coef);
																		$$ = local_term;
																	}else if(g_coef_type == GMP_MPZ){
																		term_z *local_term_z = create_term_z(g_num_variables);
																		mpz_t temp_coef;
																		mpz_init(temp_coef);
																		mpz_set_str(temp_coef, $1, 10);
																		mpz_mul(local_term_z->coef, local_term_z->coef, temp_coef);
																		$$ = local_term_z;
																	}
																	
																	#ifdef PARSER_DEBUG
																		printf_blue("term: coef ");
																		printf("%s\n", $1);
																		print_term($$, g_variables, g_num_variables, "term: $$: coef");
																	#endif
																	free($1);
																}
				| powerVariable									{
																	if(g_coef_type == GMP_MPQ){
																		term *local_term = create_term(g_num_variables);
																		if(is_negative_var){
																			mpq_set_str(local_term->coef, "-1", 10);
																			is_negative_var = 0;
																		}else{ //un-necessary, coef already 1 ahen term is created
																			mpq_set_str(local_term->coef, "1", 10);
																		}
																		fill_term_exponent(g_variables, local_term->exp, $1->var, $1->exp, &g_num_variables);
																		$$ = local_term;
																	}else if(g_coef_type == GMP_MPZ){
																		term_z *local_term_z = create_term_z(g_num_variables);
																		if(is_negative_var){
																			mpz_set_str(local_term_z->coef, "-1", 10);
																			is_negative_var = 0;
																		}else{ //un-necessary, coef already 1 ahen term is created
																			mpz_set_str(local_term_z->coef, "1", 10);
																		}
																		fill_term_exponent(g_variables, local_term_z->exp, $1->var, $1->exp, &g_num_variables);
																		$$ = local_term_z;
																	}

																	#ifdef PARSER_DEBUG
																		printf_blue("term: powerVaribale ");
																		printf("%s^%llu\n", $1->var, $1->exp);
																		if(g_coef_type == GMP_MPQ)
																			print_term($$, g_variables, g_num_variables, "term: $$: powerVariable");
																		else
																			print_term_z($$, g_variables, g_num_variables, "term: $$: powerVariable");
																	#endif
																	free($1->var);
																	free($1);
																}
				;				
powerVariable	: variable										{
																	is_new_var = 0;
																	if(!is_all_var_defined && !check_if_it_exists(g_variables, $1, g_num_variables)){
																		g_variables = push_back_dynamic(g_variables, &g_num_variables, $1);
																		is_new_var = 1;
																	}
																	powervar *temp_power_var = create_power_var($1, 1);
																	$$ = temp_power_var;
																	#ifdef PARSER_DEBUG
																		printf_green("powervariable: variable ");
																		printf("%s\n", $1);
																		printf_green("powervariable: $$: variable ");
																		printf("%s^%llu\n", $$->var, $$->exp);
																	#endif
																	free($1);
																}
				| variable POWER exponent						{
																	is_new_var = 0;
																	if(!is_all_var_defined && !check_if_it_exists(g_variables, $1, g_num_variables)){
																		g_variables = push_back_dynamic(g_variables, &g_num_variables, $1);
																		is_new_var = 1;
																	}
																	powervar *temp_power_var = create_power_var($1, strtoul($3, NULL, 10));
																	$$ = temp_power_var;
																	#ifdef PARSER_DEBUG
																		printf_green("powervariable: variable^exponent ");
																		printf("%s^%s\n", $1, $3);
																		printf_green("powervariable: $$: variable^exponent ");
																		printf("%s^%llu\n", $$->var, $$->exp);
																	#endif
																	free($1);
																}
				;
variable 		: VAR											{
																	$$ = strdup($1);
																	#ifdef PARSER_DEBUG
																		printf_blue("variable: VAR ");
																		printf(" %s\n", $1);
																	#endif
																	free($1);
																}
				;
exponent		: NUM											{
																	$$ = strdup($1);
																	#ifdef PARSER_DEBUG
																		printf_blue("exponent: NUM ");
																		printf(" %s\n", $1);
																	#endif
																	free($1);
																}
				;
coef			: NUM											{
																	$$ = strdup($1);
																	#ifdef PARSER_DEBUG
																		printf_blue("coef: NUM ");
																		printf(" %s\n", $1);
																	#endif
																	free($1);
																}
				| RATNUM										{
																	$$ = strdup($1);
																	#ifdef PARSER_DEBUG
																		printf_blue("coef: RATNUM ");
																		printf(" %s\n", $1);
																	#endif	
																	free($1);
																}
				;
other 			: variable 										{ 
																	g_variables = push_back_dynamic(g_variables, &g_num_variables, $1);
																}
				|  LS_BRACE other 								{}
				|  other  RS_BRACE								{ is_all_var_defined = 1;}
				|  other COMMA variable							{
																	g_variables = push_back_dynamic(g_variables, &g_num_variables, $3);
																	free($3);
																}
				;

%%

void yyerror(char *msg){
	return;
}

