struct RationalPolynomial*number_field_element_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*base, struct Integer*exponent,
    struct RationalPolynomial*generator_minimal_polynomial);

#include "integer/rational/rational_polynomial/number_field_element/number_field_polynomial/number_field_polynomial.h"