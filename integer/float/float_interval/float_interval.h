struct FloatInterval float_interval_copy(struct Stack*output_stack, struct FloatInterval*a);
struct Float float_interval_get_max_magnitude(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a);
struct FloatInterval float_interval_add(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b);
struct FloatInterval float_interval_negate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a);
struct FloatInterval float_interval_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*minuend,
    struct FloatInterval*subtrahend);
struct FloatInterval float_interval_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b);
bool float_intervals_are_disjoint(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct FloatInterval*a, struct FloatInterval*b);
bool region_a_contains_b(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b);
struct RationalInterval float_interval_to_rational_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*a);

#include "integer/float/float_interval/region/region.h"