struct Rational rational_copy(struct Stack*output_stack, struct Rational*a);
struct Rational rational_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*numerator, struct Integer*denominator);
bool rational_equals(struct Rational*a, struct Rational*b);
struct Rational rational_get_magnitude(struct Stack*output_stack, struct Rational*a);
struct Rational rational_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Rational*a, struct Rational*b);
struct Rational rational_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b);
struct Rational rational_negate(struct Stack*output_stack, struct Rational*a);
struct Rational rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*minuend, struct Rational*subtrahend);
struct Rational rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational rational_unreduced_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational rational_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b);
struct Rational rational_get_reciprocal(struct Stack*output_stack, struct Rational*a);
struct Rational rational_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Rational*divisor);
struct Rational rational_integer_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Integer*divisor);
struct Rational rational_double(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a);
struct Rational rational_halve(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a);
struct Rational rational_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational base, struct Integer*exponent);
int8_t rational_get_sign(struct Rational*a);
int8_t rational_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Rational*a, struct Rational*b);
struct Rational rational_get_min(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational rational_get_max(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational rational_get_argument(struct Stack*output_stack, struct Rational*a);
struct RationalInterval rational_estimate_cosine(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);
struct RationalInterval rational_estimate_sine(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);
struct RationalInterval rational_estimate_arctangent(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);
struct RationalInterval rational_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*y, struct Rational*x,
    struct Rational interval_size);
struct FloatInterval rational_get_float_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);

#include "integer/rational/gaussian_rational/gaussian_rational.h"
#include "integer/rational/matrix/matrix.h"
#include "integer/rational/rational_interval/rational_interval.h"
#include "integer/rational/rational_polynomial/rational_polynomial.h"