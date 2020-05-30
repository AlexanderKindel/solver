struct GaussianRational gaussian_rational_copy(struct Stack*output_stack,
    struct GaussianRational*a);
bool gaussian_rational_equals(struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational gaussian_rational_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational gaussian_rational_negate(struct Stack*output_stack,
    struct GaussianRational*a);
struct GaussianRational gaussian_rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*minuend,
    struct GaussianRational*subtrahend);
struct GaussianRational gaussian_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational gaussian_rational_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct Rational*b);

#include "integer/rational/gaussian_rational/gaussian_rational_polynomial/gaussian_rational_polynomial.h"