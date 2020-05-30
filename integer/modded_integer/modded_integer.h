struct Integer*integer_get_residue(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic);
struct Integer*modded_integer_negate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend,
    struct Integer*characteristic);
struct Integer*modded_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic);
struct Integer*modded_integer_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);

#include "integer/modded_integer/modded_polynomial/modded_polynomial.h"