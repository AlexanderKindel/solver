struct Float float_copy(struct Stack*output_stack, struct Float*a);
struct Float float_reduce(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*significand, struct Integer*exponent);
bool float_equals(struct Float*a, struct Float*b);
struct Float float_get_magnitude(struct Stack*output_stack, struct Float*a);
struct Float float_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct Float float_negate(struct Stack*output_stack, struct Float*a);
struct Float float_subtract(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*minuend, struct Float*subtrahend);
struct Float float_multiply(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct Float float_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float base, struct Integer*exponent);
int8_t float_get_sign(struct Float*a);
int8_t float_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Float*a, struct Float*b);
struct Float float_get_min(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct Float float_get_max(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct FloatInterval float_estimate_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float*a, struct Rational*interval_size,
    struct Integer*index);
struct Rational float_to_rational(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float*a);

#include "integer/float/float_interval/float_interval.h"