struct Number*number_copy(struct Stack*output_stack, struct Number*a);
struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value);
int8_t number_formal_compare(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct Number*a, struct Number*b);
struct RationalPolynomial*number_annulling_polynomial_to_minimal_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*annulling_polynomial, struct Number*a);
size_t number_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a);

#include "number/number_add/number_add.h"
#include "number/number_divide/number_divide.h"
#include "number/number_estimate/number_estimate.h"
#include "number/number_exponentiate/number_exponentiate.h"
#include "number/number_multiply/number_multiply.h"
#include "number/roots_of_unity/roots_of_unity.h"