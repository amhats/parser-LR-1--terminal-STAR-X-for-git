#include <stdio.h>
#include <stdlib.h>
#include "../include/bpas_parser.h"
#include "../parser_helper.h"
#include "../parser_print.h"
#include "../BPAS-PARSER-TEST-FRAMEWORK/parser_testframework.h"
#include "../BPAS-PARSER-TEST-FRAMEWORK/maple_evaluate.h"

int main(int argc, char **argv){
    char *poly = "-(x^6*y - x*y^6)^6";
    AltArr_t *result = generate_altarr(poly);
    print_naked_AltArr_t_poly(result);
    freePolynomial_AA(result);

    char *poly_z = "-(x^6*y - x*y^6)^6";
    // AltArrZ_t *result_z = generate_altarr_z(poly_z);
    // print_naked_AltArr_t_poly_z(result_z);
    // freePolynomial_AAZ(result_z);


    char *poly2 = "(x+y+z^2-m)^2";
    int numvars = 4;
    char **vars = create_dynamic_str_array(numvars, "x", "y", "z", "m");
    AltArr_t *result2 = generate_altarr_var_defined(poly2, vars, numvars);
    print_naked_AltArr_t_poly(result2);
    print_poly_to_terminal(result2, vars, numvars);

    // char *poly3 = "-(x^6*y - x*y^6)^6";
    char *poly3 = "-x-y";
    int numvars3 = 2;
    char **vars3 = create_dynamic_str_array(numvars3, "x", "y");
    AltArr_t *result3 = generate_altarr(poly3);
    print_naked_AltArr_t_poly(result3);
    print_poly_to_terminal(result3, vars3, numvars3);
    char *ret = print_poly_to_string_variable(result3, vars3, numvars3);
    ASSERT_EQUAL(maple_evaluate_equal(poly3, ret), 0);
    parser_test(poly3, "TEST");

    // const char *poly4 = "((x*y+z)*(y*z+x)+z*x+y)*((y*z+x)*(x*z+y)+x*y+z)+(x*z+y)*(x*y+z)+y*z+x";
    const char *poly4 = "9*z^44*z^5*x*x^4 + 12*x*y*z*z*w*w*w";
    int numvars4 = 4;
    char **vars4 = create_dynamic_str_array(numvars4, "x", "y", "z", "w");
    AltArr_t *result4 = generate_altarr(poly4);
    print_poly_to_terminal(result4, vars4, numvars3);
    char *ret2 = print_poly_to_string_variable(result4, vars4, numvars4);
    ASSERT_EQUAL(maple_evaluate_equal(poly4, ret2), 0);
    parser_test(poly4, "TEST");

    return 0;
}