#include "declarations.h"

struct Rational*rational_place_value(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = &rational_one;
    while (rational_compare(output_stack, local_stack, out, a) > 0)
    {
        out->denominator = integer_doubled(local_stack, out->denominator);
    }
    out->denominator = integer_doubled(output_stack, out->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval*number_refine_rational_estimate_interval(
    rational_estimate_getter get_rational_estimate, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalInterval*estimate,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval*out;
    if (rational_compare(output_stack, local_stack, interval_size,
        rational_subtract(local_stack, output_stack, estimate->max, estimate->min)) < 0)
    {
        out = get_rational_estimate(output_stack, local_stack, a, interval_size);
    }
    else
    {
        out = rational_interval_copy(output_stack, estimate);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval*number_halve_rational_estimate_interval(
    rational_estimate_getter get_rational_estimate, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalInterval*estimate)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval*out = get_rational_estimate(output_stack, local_stack, a,
        rational_half(local_stack, output_stack,
            rational_subtract(local_stack, output_stack, estimate->max, estimate->min)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_refine_float_estimate_interval(float_estimate_getter get_float_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct FloatInterval*estimate, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*out;
    if (rational_compare(output_stack, local_stack, interval_size,
        float_to_rational(local_stack, output_stack,
            float_subtract(local_stack, output_stack, estimate->max, estimate->min))) < 0)
    {
        out = get_float_estimate(output_stack, local_stack, a, interval_size);
    }
    else
    {
        out = float_interval_copy(output_stack, estimate);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_halve_float_estimate_interval(float_estimate_getter get_float_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct FloatInterval*estimate)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*current_interval_size =
        float_subtract(local_stack, output_stack, estimate->max, estimate->min);
    struct FloatInterval* out = get_float_estimate(output_stack, local_stack, a,
        float_to_rational(local_stack, output_stack,
            &(struct Float){current_interval_size->significand,
                integer_add(local_stack, current_interval_size->significand, &INT(1, -))}));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_float_estimate_from_rational(
    rational_estimate_getter get_rational_estimate, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*rational_estimate_interval_size =
        rational_place_value(local_stack, output_stack, interval_size);
    struct FloatInterval*out = rational_interval_to_float_interval(output_stack, local_stack,
        get_rational_estimate(local_stack, output_stack, a, rational_estimate_interval_size),
        rational_estimate_interval_size);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval*number_argument_cosine_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval*argument_estimate =
        number_rational_argument_estimate(local_stack, output_stack, a, argument_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_min_cosine,
        argument_estimate->min, argument_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_max_cosine,
        argument_estimate->max, argument_interval_size);
    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
    out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
        argument_min_cosine.max, argument_max_cosine.max));
    pi_shrink_interval_to_one_side_of_value(argument_estimate->min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate->max);
    if (rational_compare(output_stack, local_stack, argument_estimate->min, pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, pi.max, argument_estimate->max) <= 0)
    {
        out->min->numerator = integer_initialize(output_stack, 1, -1);
        out->min->denominator = &one;
    }
    else
    {
        out->min = rational_copy(output_stack, rational_max(output_stack, local_stack,
            argument_min_cosine.min, argument_max_cosine.min));
        out->min = rational_copy(output_stack, out->min);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval*number_argument_sine_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval*argument_estimate = number_rational_argument_estimate(local_stack,
        output_stack, a, argument_interval_size);
    struct RationalInterval argument_min_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_min_sine, argument_estimate->min,
        argument_interval_size);
    struct RationalInterval argument_max_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_max_sine, argument_estimate->max,
        argument_interval_size);
    struct RationalInterval argument_estimate_multiple =
    { rational_doubled(local_stack, output_stack, argument_estimate->min),
    rational_doubled(local_stack, output_stack, argument_estimate->max) };
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.max);
    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
    if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min, pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, pi.max, argument_estimate_multiple.max) <= 0)
    {
        out->max = &rational_one;
        out->min = rational_copy(output_stack,
            rational_max(output_stack, local_stack, argument_min_sine.min, argument_max_sine.min));
    }
    else
    {
        argument_estimate_multiple.min = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.min, &INT(3, +));
        argument_estimate_multiple.max = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.max, &INT(3, +));
        pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.max);
        if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
            pi.min) <= 0 && rational_compare(output_stack, local_stack, pi.max,
                argument_estimate_multiple.max) <= 0)
        {
            out->min->numerator = integer_initialize(output_stack, 1, -1);
            out->min->denominator = &one;
            out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
        else
        {
            out->min = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.min, argument_max_sine.min));
            out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval*number_rectangular_part_from_polar_form(
    struct RationalInterval*(trig_function)(struct Stack*, struct Stack*, struct Number*,
        struct Rational*),
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval*magnitude_estimate =
        number_rational_magnitude_estimate(local_stack, output_stack, a, &rational_one);
    struct Rational*factor_interval_size = rational_divide(local_stack, output_stack, interval_size,
        rational_integer_add(local_stack, output_stack, magnitude_estimate->max, &INT(2, +)));
    magnitude_estimate =
        number_refine_rational_estimate_interval(number_rational_magnitude_estimate, local_stack,
            output_stack, a, magnitude_estimate, factor_interval_size);
    struct RationalInterval*trig_value =
        trig_function(local_stack, output_stack, a, factor_interval_size);
    struct RationalInterval*out =
        rational_interval_multiply(output_stack, local_stack, trig_value, magnitude_estimate);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_estimate_part_sum(float_estimate_getter get_element_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*element_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size,
            integer_from_size_t(local_stack, a->element_count));
    struct FloatInterval*out = &(struct FloatInterval) { &float_zero, &float_zero };
    struct FloatInterval*element_estimate;
    for (size_t i = a->element_count; i-- > 1;)
    {
        out = float_interval_add(local_stack, output_stack, out,
            get_element_estimate(local_stack, output_stack, a->elements[i],
                element_estimate_interval_size));
    }
    out = float_interval_add(output_stack, local_stack, out,
        get_element_estimate(local_stack, output_stack, a->elements[0],
            element_estimate_interval_size));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_float_real_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval*out =
            number_float_magnitude_estimate(local_stack, output_stack, a, interval_size);
        if (a->value->numerator->sign < 0)
        {
            out = float_interval_negative(output_stack, local_stack, out);
        }
        else
        {
            out = float_interval_copy(output_stack, out);
        }
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '^':
        return number_float_estimate_from_rational(number_rational_real_part_estimate, output_stack,
            local_stack, a, interval_size);
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RectangularEstimate rectangular_estimate;
        number_rectangular_estimate(local_stack, output_stack, &rectangular_estimate, a,
            interval_size);
        struct FloatInterval*out =
            float_interval_copy(output_stack, rectangular_estimate.real_part_estimate);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
        return number_estimate_part_sum(number_float_real_part_estimate, output_stack, local_stack,
            a, interval_size);
    default:
        crash("Number operation not recognized.");
    }
}

struct RationalInterval*number_rational_real_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        return number_rectangular_part_from_polar_form(number_argument_cosine_estimate,
            output_stack, local_stack, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval*out = float_interval_to_rational_interval(output_stack, local_stack,
            number_float_real_part_estimate(local_stack, output_stack, a, interval_size));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
}

struct FloatInterval*number_float_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
    {
        struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
        out->min = &float_zero;
        out->max = &float_zero;
        return out;
    }
    case '^':
        return number_float_estimate_from_rational(number_rational_imaginary_part_estimate,
            output_stack, local_stack, a, interval_size);
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RectangularEstimate rectangular_estimate;
        number_rectangular_estimate(local_stack, output_stack, &rectangular_estimate, a,
            interval_size);
        struct FloatInterval*out =
            float_interval_copy(output_stack, rectangular_estimate.imaginary_part_estimate);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
        return number_estimate_part_sum(number_float_imaginary_part_estimate, output_stack,
            local_stack, a, interval_size);
    default:
        crash("Number operation not recognized.");
    }
}

struct RationalInterval*number_rational_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        return number_rectangular_part_from_polar_form(number_argument_sine_estimate, output_stack,
            local_stack, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval*out = float_interval_to_rational_interval(output_stack, local_stack,
            number_float_imaginary_part_estimate(local_stack, output_stack, a, interval_size));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
}

void number_rectangular_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RectangularEstimate*out, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '*')
    {
        struct RectangularEstimate*factor_estimates =
            ARRAY_ALLOCATE(local_stack, a->element_count, struct RectangularEstimate);
        struct Float*interval_scale = &float_one;
        struct Float*term_to_subtract = &float_one;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            number_rectangular_estimate(local_stack, output_stack, factor_estimates + i,
                a->elements[i], &rational_one);
            struct Float*bound_max = float_max(output_stack, local_stack,
                float_interval_max_magnitude(local_stack, output_stack,
                    factor_estimates[i].real_part_estimate),
                float_interval_max_magnitude(local_stack, output_stack,
                    factor_estimates[i].imaginary_part_estimate));
            interval_scale = float_multiply(local_stack, output_stack, interval_scale,
                float_add(local_stack, output_stack, bound_max, &float_one));
            term_to_subtract =
                float_multiply(local_stack, output_stack, term_to_subtract, bound_max);
        }
        interval_size = rational_divide(local_stack, output_stack, interval_size,
            float_to_rational(local_stack, output_stack,
                float_subtract(local_stack, output_stack, interval_scale, term_to_subtract)));
        out->real_part_estimate = ALLOCATE(output_stack, struct FloatInterval);
        out->real_part_estimate->min = &float_zero;
        out->real_part_estimate->max = &float_zero;
        out->imaginary_part_estimate = ALLOCATE(output_stack, struct FloatInterval);
        out->imaginary_part_estimate->min = &float_zero;
        out->imaginary_part_estimate->max = &float_zero;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            factor_estimates[i].real_part_estimate =
                number_refine_float_estimate_interval(number_float_real_part_estimate, local_stack,
                    output_stack, a->elements[i], factor_estimates[i].real_part_estimate,
                    interval_size);
            factor_estimates[i].imaginary_part_estimate =
                number_refine_float_estimate_interval(number_float_imaginary_part_estimate,
                    local_stack, output_stack, a->elements[i],
                    factor_estimates[i].imaginary_part_estimate, interval_size);
            struct FloatInterval*new_real_part = float_interval_subtract(local_stack, output_stack,
                float_interval_multiply(local_stack, output_stack, out->real_part_estimate,
                    factor_estimates[i].real_part_estimate),
                float_interval_multiply(local_stack, output_stack, out->imaginary_part_estimate,
                    factor_estimates[i].imaginary_part_estimate));
            out->imaginary_part_estimate = float_interval_add(local_stack, output_stack,
                float_interval_multiply(local_stack, output_stack, out->real_part_estimate,
                    factor_estimates[i].imaginary_part_estimate),
                float_interval_multiply(local_stack, output_stack, out->imaginary_part_estimate,
                    factor_estimates[i].real_part_estimate));
            out->real_part_estimate = new_real_part;
        }
        out->real_part_estimate = float_interval_copy(output_stack, out->real_part_estimate);
        out->imaginary_part_estimate =
            float_interval_copy(output_stack, out->imaginary_part_estimate);
    }
    else
    {
        out->real_part_estimate =
            number_float_real_part_estimate(output_stack, local_stack, a, interval_size);
        out->imaginary_part_estimate =
            number_float_imaginary_part_estimate(output_stack, local_stack, a, interval_size);
    }
}

void number_rectangular_estimate_halve_dimensions(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct Number*a)
{
    out->real_part_estimate = number_halve_float_estimate_interval(number_float_real_part_estimate,
        output_stack, local_stack, a, out->real_part_estimate);
    out->imaginary_part_estimate =
        number_halve_float_estimate_interval(number_float_imaginary_part_estimate, output_stack,
            local_stack, a, out->imaginary_part_estimate);
}

struct FloatInterval*number_rectangular_part_estimate_for_magnitude_estimate(
    struct EstimateGetters*part_estimate_getters, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*out = part_estimate_getters->fl(output_stack, local_stack, a,
        rational_divide(local_stack, output_stack, interval_size,
            rational_integer_add(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    rational_interval_max_magnitude(local_stack, output_stack, 
                        part_estimate_getters->rational(local_stack, output_stack, a,
                            &rational_one)), &INT(4, +)), &INT(2, +))));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_magnitude_estimate_bound_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Float*real_part_bound,
    struct Float*imaginary_part_bound, struct Rational*interval_size)
{
    float_estimate_root(output_stack, local_stack, out_min, out_max,
        float_add(local_stack, output_stack,
            float_multiply(local_stack, output_stack, real_part_bound, real_part_bound),
            float_multiply(local_stack, output_stack, imaginary_part_bound, imaginary_part_bound)),
        interval_size, &INT(2, +));
}

struct FloatInterval*number_float_magnitude_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    switch (a->operation)
    {
    case 'r':
        return rational_float_estimate(output_stack, local_stack, a->value, interval_size);
    case '^':
    {
        struct FloatInterval*radicand_magnitude_estimate = 
            number_float_magnitude_estimate(local_stack, output_stack,
                a->radicand, rational_integer_divide(local_stack, output_stack,
                    rational_exponentiate(local_stack, output_stack, interval_size, a->index),
                    &INT(3, +)));
        struct Rational*bound_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct Float*unused_estimate_bound;
        struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
        float_estimate_root(local_stack, output_stack, &out->min, &unused_estimate_bound,
            radicand_magnitude_estimate->min, bound_interval_size, a->index);
        float_estimate_root(local_stack, output_stack, &unused_estimate_bound, &out->max,
            radicand_magnitude_estimate->max, bound_interval_size, a->index);
        out->min = float_copy(output_stack, out->min);
        out->max = float_copy(output_stack, out->max);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        struct FloatInterval**magnitude_estimates =
            ARRAY_ALLOCATE(local_stack, a->element_count, struct FloatInterval*);
        struct Float*interval_scale = &float_one;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            magnitude_estimates[i] = number_float_magnitude_estimate(local_stack, output_stack,
                a->elements[i], &rational_one);
            interval_scale = float_multiply(local_stack, output_stack, interval_scale,
                float_add(local_stack, output_stack, &float_one,
                    float_interval_max_magnitude(local_stack, output_stack,
                        magnitude_estimates[i])));
        }
        interval_size = rational_divide(local_stack, output_stack, interval_size,
            float_to_rational(local_stack, output_stack, interval_scale));
        struct FloatInterval*product = &(struct FloatInterval){ &float_one, &float_one };
        for (size_t i = 0; i < a->element_count; ++i)
        {
            magnitude_estimates[i] =
                number_refine_float_estimate_interval(number_float_argument_estimate, local_stack,
                    output_stack, a->elements[i], magnitude_estimates[i], interval_size);
            product = float_interval_multiply(local_stack, output_stack, product,
                magnitude_estimates[i]);
        }
        struct FloatInterval*out = float_interval_copy(output_stack, product);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct FloatInterval*real_part_estimate =
            number_rectangular_part_estimate_for_magnitude_estimate(&real_estimate_getters,
                local_stack, output_stack, a, interval_size);
        struct FloatInterval*imaginary_part_estimate =
            number_rectangular_part_estimate_for_magnitude_estimate(&imaginary_estimate_getters,
                local_stack, output_stack, a, interval_size);
        struct Float*unused_estimate_bound;
        struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
        number_magnitude_estimate_bound_interval(local_stack, output_stack, &out->min,
            &unused_estimate_bound, real_part_estimate->min, imaginary_part_estimate->min,
            interval_size);
        number_magnitude_estimate_bound_interval(local_stack, output_stack, &unused_estimate_bound,
            &out->max, real_part_estimate->max, imaginary_part_estimate->max, interval_size);
        out->min = float_copy(output_stack, out->min);
        out->max = float_copy(output_stack, out->max);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct RationalInterval*number_rational_magnitude_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval*out = float_interval_to_rational_interval(output_stack, local_stack,
        number_float_magnitude_estimate(local_stack, output_stack, a, interval_size));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*number_float_argument_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '*')
    {
        return number_estimate_part_sum(number_float_argument_estimate, output_stack, local_stack,
            a, interval_size);
    }
    return number_float_estimate_from_rational(number_rational_argument_estimate, output_stack,
        local_stack, a, interval_size);
}

struct RationalInterval*number_rational_argument_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
    {
        struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
        if (a->value->numerator->sign < 0)
        {
            pi_estimate(interval_size);
            out->min = pi.min;
            out->max = pi.max;
        }
        else
        {
            out->min = &rational_zero;
            out->max = &rational_zero;
        }
        return out;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval*radicand_rational_argument_estimate =
            number_rational_argument_estimate(local_stack, output_stack, a->radicand,
                rational_integer_multiply(local_stack, output_stack, interval_size, a->index));
        struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
        out->min = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate->min, a->index);
        out->max = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate->max, a->index);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval*out = float_interval_to_rational_interval(output_stack, local_stack,
            number_float_argument_estimate(local_stack, output_stack, a, interval_size));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval*real_part_estimate =
            number_rational_real_part_estimate(local_stack, output_stack, a, &rational_one);
        struct RationalInterval*imaginary_part_estimate =
            number_rational_imaginary_part_estimate(local_stack, output_stack, a, &rational_one);
        if (imaginary_part_estimate->min->numerator->sign !=
            imaginary_part_estimate->max->numerator->sign)
        {
            while (true)
            {
                struct Rational*imaginary_max_magnitude =
                    rational_interval_max_magnitude(local_stack, output_stack,
                        imaginary_part_estimate);
                if (rational_polynomial_root_count_in_rectangle(local_stack, output_stack,
                    a->minimal_polynomial, real_part_estimate,
                    &(struct RationalInterval){rational_negative(local_stack,
                        imaginary_max_magnitude), imaginary_max_magnitude}) == 1)
                {
                    while (real_part_estimate->min->numerator->sign ==
                        -real_part_estimate->max->numerator->sign)
                    {
                        real_part_estimate = number_halve_rational_estimate_interval(
                            number_rational_real_part_estimate, local_stack, output_stack, a,
                            real_part_estimate);
                    }
                    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
                    if (real_part_estimate->max->numerator->sign > 0)
                    {
                        out->min = &rational_zero;
                        out->max = &rational_zero;
                    }
                    else
                    {
                        pi_estimate(interval_size);
                        out->min = pi.min;
                        out->max = pi.max;
                    }
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                imaginary_part_estimate =
                    number_halve_rational_estimate_interval(number_rational_imaginary_part_estimate,
                        local_stack, output_stack, a, imaginary_part_estimate);
                if (imaginary_part_estimate->min->numerator->sign ==
                    imaginary_part_estimate->max->numerator->sign)
                {
                    break;
                }
                real_part_estimate =
                    number_halve_rational_estimate_interval(number_rational_real_part_estimate,
                        local_stack, output_stack, a, real_part_estimate);
            }
        }
        interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct Rational*imaginary_part_interval_size =
            rational_interval_factor_interval_size(local_stack, output_stack,
                imaginary_part_estimate, &(struct RationalInterval){
                    rational_reciprocal(local_stack, real_part_estimate->max),
                    rational_reciprocal(local_stack, real_part_estimate->min)}, interval_size);
        imaginary_part_estimate =
            number_refine_rational_estimate_interval(number_rational_imaginary_part_estimate,
                local_stack, output_stack, a, imaginary_part_estimate,
                imaginary_part_interval_size);
        struct Rational*real_part_interval_size = rational_multiply(local_stack, output_stack,
            real_part_estimate->min, imaginary_part_interval_size);
        real_part_interval_size = rational_divide(local_stack, output_stack,
            rational_multiply(local_stack, output_stack, real_part_estimate->min,
                real_part_interval_size),
            rational_subtract(local_stack, output_stack, &rational_one, real_part_interval_size));
        real_part_estimate =
            number_refine_rational_estimate_interval(number_rational_real_part_estimate,
                local_stack, output_stack, a, real_part_estimate, real_part_interval_size);
        struct RationalInterval*out = rational_interval_estimate_atan2(output_stack, local_stack,
            imaginary_part_estimate, real_part_estimate, interval_size);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}