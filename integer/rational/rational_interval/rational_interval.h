struct RationalInterval rational_interval_copy(struct Stack*output_stack,
    struct RationalInterval*a);
struct Rational rational_interval_get_max_magnitude(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a);
struct RationalInterval rational_interval_add(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*a, struct RationalInterval*b);
struct RationalInterval rational_interval_negate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a);
struct RationalInterval rational_interval_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct RationalInterval*b);
void rational_interval_expand_bounds(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalInterval*a,
    struct Rational*bound_candidate);
struct Rational rational_interval_get_factor_interval_size(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*factor_a,
    struct RationalInterval*factor_b, struct Rational*product_interval_size);
struct RationalInterval rational_interval_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size);
struct FloatInterval rational_interval_to_float_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*a,
    struct Rational*bound_interval_size);

#include "integer/rational/rational_interval/pi/pi.h"