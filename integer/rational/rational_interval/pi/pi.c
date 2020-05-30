#include "declarations.h"

void pi_refine_interval(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Rational*pi_interval_size)
{
    struct Rational sum = rational_add(local_stack, output_stack,
        &(struct Rational) { &one, integer_add(local_stack, pi_eight_k, INT(5, 1)) },
        &(struct Rational) { &one, integer_add(local_stack, pi_eight_k, INT(6, 1)) });
    struct Rational term = rational_reduce(local_stack, output_stack, INT(2, 1),
        integer_add(local_stack, pi_eight_k, INT(4, 1)));
    sum = rational_add(local_stack, output_stack, &term, &sum);
    sum = rational_subtract(local_stack, output_stack,
        &(struct Rational) { INT(4, 1), integer_add(local_stack, pi_eight_k, &one) }, &sum);
    struct Rational quotient =
        rational_integer_divide(local_stack, output_stack, &sum, pi_sixteen_to_the_k);
    pi.min = rational_add(output_stack, local_stack, &pi.min, &quotient);
    *pi_interval_size = rational_multiply(output_stack, local_stack, pi_interval_size,
        &(struct Rational) { &one, INT(16, 1) });
    pi_eight_k = integer_add(output_stack, pi_eight_k, INT(8, 1));
    pi_sixteen_to_the_k =
        integer_multiply(output_stack, local_stack, pi_sixteen_to_the_k, INT(16, 1));
    local_stack->cursor = local_stack->start;
}

void pi_set_stacks(struct Stack**out_old_stack, struct Stack**out_new_stack)
{
    if (pi.min.denominator < pi_stack_b.start)
    {
        *out_old_stack = &pi_stack_a;
        *out_new_stack = &pi_stack_b;
    }
    else
    {
        *out_old_stack = &pi_stack_b;
        *out_new_stack = &pi_stack_a;
    }
}

void pi_estimate(struct Rational*interval_size)
{
    struct Stack*old_stack;
    struct Stack*new_stack;
    pi_set_stacks(&old_stack, &new_stack);
    struct Rational pi_interval_size = rational_subtract(new_stack, old_stack, &pi.max, &pi.min);
    if (rational_compare(&pi_stack_a, &pi_stack_b, &pi_interval_size, interval_size) > 0)
    {
        do
        {
            pi_refine_interval(new_stack, old_stack, &pi_interval_size);
            SWAP(old_stack, new_stack, struct Stack*);
        } while (rational_compare(&pi_stack_a, &pi_stack_b, &pi_interval_size, interval_size) > 0);
        pi.max = rational_add(old_stack, new_stack, &pi.min, &pi_interval_size);
    }
    else
    {
        new_stack->cursor = new_stack->start;
    }
}

void pi_shrink_interval_to_one_side_of_value(struct Rational*value)
{
    if (rational_compare(&pi_stack_a, &pi_stack_b, &pi.min, value) <= 0 &&
        rational_compare(&pi_stack_a, &pi_stack_b, value, &pi.max) <= 0)
    {
        struct Stack*old_stack;
        struct Stack*new_stack;
        pi_set_stacks(&old_stack, &new_stack);
        struct Rational pi_interval_size =
            rational_subtract(old_stack, new_stack, &pi.max, &pi.min);
        do
        {
            pi_refine_interval(new_stack, old_stack, &pi_interval_size);
            pi.max = rational_add(new_stack, old_stack, &pi.min, &pi_interval_size);
            SWAP(old_stack, new_stack, struct Stack*);
        } while (rational_compare(&pi_stack_a, &pi_stack_b, &pi.min, value) <= 0 &&
            rational_compare(&pi_stack_a, &pi_stack_b, value, &pi.max) <= 0);
    }
}