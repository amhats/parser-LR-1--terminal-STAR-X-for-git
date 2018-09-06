#include <stdio.h>
#include "parser_testframework.h"
#include "maple_evaluate.h"



int main(int argc, char** argv) {
    ASSERT_EQUAL(1, 0);
    parser_test("1 and 0", NULL);
    ASSERT_STRING_EQUAL("JOKE", "JOKE");
    parser_test("JOKE", NULL);
    ASSERT_STRING_EQUAL("JOKE", "joke");
    parser_test("JOKE and joke", NULL);
    ASSERT_EQUAL(maple_evaluate_equal("[x]-x", "-x"), 0);
    parser_test("-x", NULL);
    ASSERT_EQUAL(maple_evaluate_equal("[x,y](x+y)^2", "x^2+2*x*y+y^2"), 0);
    parser_test("(x+y)^2", NULL);
    ASSERT_STRING_EQUAL(maple_evaluate("[x,y,z,t](x+y+z+2*t)^2"), "4*t^2+4*t*x+4*t*y+4*t*z+x^2+2*x*y+2*x*z+y^2+2*y*z+z^2");
    parser_test("(x+y+z+2*t)^2", NULL);
    ASSERT_STRING_EQUAL(maple_evaluate("[x,y,z,t](1-x-y-z+t)^2"), "t^2-2*t*x-2*t*y-2*t*z+x^2+2*x*y+2*x*z+y^2+2*y*z+z^2+2*t-2*x-2*y-2*z+1");
    parser_test("(1-x-y-z+t)^2", NULL);
    parser_report();
    return (EXIT_SUCCESS);
}
