#include "declarations.h"

bool regions_are_disjoint(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b)
{
    return float_intervals_are_disjoint(local_stack_a, local_stack_b, &a->real_interval,
        &b->real_interval) ||
        float_intervals_are_disjoint(local_stack_a, local_stack_b, &a->imaginary_interval,
            &b->imaginary_interval);
}

bool region_a_contains_b(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b)
{
    return float_interval_a_contains_b(local_stack_a, local_stack_b, &a->real_interval,
        &b->real_interval) &&
        float_interval_a_contains_b(local_stack_a, local_stack_b, &a->imaginary_interval,
            &b->imaginary_interval);
}