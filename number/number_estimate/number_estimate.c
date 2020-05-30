#include "declarations.h"

typedef struct RationalInterval(*RationalEstimateGetter)(struct Stack*restrict,
    struct Stack*restrict, struct Number*, struct Rational*);

typedef struct FloatInterval(*FloatEstimateGetter)(struct Stack*restrict, struct Stack*restrict,
    struct Number*, struct Rational*);

struct EstimateGetters
{
    RationalEstimateGetter rational;
    FloatEstimateGetter fl;
};

struct EstimateGetters real_estimate_getters =
{ number_get_rational_real_part_estimate, number_get_float_real_part_estimate };
struct EstimateGetters imaginary_estimate_getters =
{ number_get_rational_imaginary_part_estimate, number_get_float_imaginary_part_estimate };
struct EstimateGetters magnitude_estimate_getters =
{ number_get_rational_magnitude_estimate, number_get_float_magnitude_estimate };

#define ARGUMENT_ESTIMATE_IS_RADIANS(interval) interval.max.denominator
#define ARGUMENT_ESTIMATE_MAKE_REVOLUTIONS(interval) interval.max.denominator = 0

struct RationalInterval number_refine_rational_estimate_interval(
    RationalEstimateGetter get_rational_estimate, struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct RationalInterval*estimate,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out;
    struct Rational current_interval_size =
        rational_subtract(local_stack, output_stack, &estimate->max, &estimate->min);
    if (rational_compare(output_stack, local_stack, interval_size, &current_interval_size) < 0)
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

struct RationalInterval number_halve_rational_estimate_interval(
    RationalEstimateGetter get_rational_estimate, struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct RationalInterval*estimate)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational current_interval_size =
        rational_subtract(local_stack, output_stack, &estimate->max, &estimate->min);
    struct Rational halved_interval_size =
        rational_halve(local_stack, output_stack, &current_interval_size);
    struct RationalInterval out =
        get_rational_estimate(output_stack, local_stack, a, &halved_interval_size);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval number_refine_float_estimate_interval(FloatEstimateGetter get_float_estimate,
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct FloatInterval*estimate, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float float_current_interval_size =
        float_subtract(local_stack, output_stack, &estimate->max, &estimate->min);
    struct Rational current_interval_size =
        float_to_rational(local_stack, output_stack, &float_current_interval_size);
    struct FloatInterval out;
    if (rational_compare(output_stack, local_stack, interval_size, &current_interval_size) < 0)
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

struct FloatInterval number_get_float_estimate_from_rational_estimate(
    RationalEstimateGetter get_rational_estimate, struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational rational_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, INT(3, 1));
    struct RationalInterval rational_estimate =
        get_rational_estimate(local_stack, output_stack, a, &rational_estimate_interval_size);
    struct FloatInterval out = rational_interval_to_float_interval(output_stack, local_stack,
        &rational_estimate, &rational_estimate_interval_size);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval revolutions_to_radians(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    if (!a->numerator->value_count)
    {
        return (struct RationalInterval) { rational_zero, rational_zero };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational doubled_a = rational_double(local_stack, output_stack, a);
    struct Rational pi_estimate_interval_size =
        rational_divide(local_stack, output_stack, interval_size, &doubled_a);
    pi_estimate(&pi_estimate_interval_size);
    struct RationalInterval out =
    { rational_multiply(output_stack, local_stack, &doubled_a, &pi.min),
        rational_multiply(output_stack, local_stack, &doubled_a, &pi.max) };
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval number_estimate_cosine_of_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, INT(3, 1));
    struct RationalInterval argument_estimate =
        number_estimate_argument(local_stack, output_stack, a, &argument_interval_size);
    struct RationalInterval out;
    if (!ARGUMENT_ESTIMATE_IS_RADIANS(argument_estimate))
    {
        if (!argument_estimate.min.numerator->value_count)
        {
            local_stack->cursor = local_stack_savepoint;
            return (struct RationalInterval) { rational_one, rational_one };
        }
        if (integer_equals(argument_estimate.min.numerator, &one))
        {
            if (integer_equals(argument_estimate.min.denominator, INT(2, 1)))
            {
                local_stack->cursor = local_stack_savepoint;
                out.min.numerator = integer_initialize(output_stack, 1, -1);
                out.min.denominator = &one;
                out.max = out.min;
                return out;
            }
            if (integer_equals(argument_estimate.min.denominator, INT(4, 1)))
            {
                local_stack->cursor = local_stack_savepoint;
                return (struct RationalInterval) { rational_zero, rational_zero };
            }
        }
        else if (integer_equals(argument_estimate.min.numerator, INT(3, 1)) &&
            integer_equals(argument_estimate.min.denominator, INT(4, 1)))
        {
            local_stack->cursor = local_stack_savepoint;
            return (struct RationalInterval) { rational_zero, rational_zero };
        }
        argument_estimate = revolutions_to_radians(local_stack, output_stack,
            &argument_estimate.min, &argument_interval_size);
    }
    struct RationalInterval argument_min_cosine = rational_estimate_cosine(local_stack,
        output_stack, &argument_estimate.min, &argument_interval_size);
    struct RationalInterval argument_max_cosine = rational_estimate_cosine(local_stack,
        output_stack, &argument_estimate.max, &argument_interval_size);
    out.max = rational_get_max(output_stack, local_stack, &argument_min_cosine.max,
        &argument_max_cosine.max);
    pi_shrink_interval_to_one_side_of_value(&argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(&argument_estimate.max);
    if (rational_compare(output_stack, local_stack, &argument_estimate.min, &pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, &pi.max, &argument_estimate.max) <= 0)
    {
        out.min.numerator = integer_initialize(output_stack, 1, -1);
        out.min.denominator = &one;
    }
    else
    {
        out.min = rational_get_min(output_stack, local_stack, &argument_min_cosine.min,
            &argument_max_cosine.min);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval number_estimate_sine_of_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, INT(3, 1));
    struct RationalInterval argument_estimate =
        number_estimate_argument(local_stack, output_stack, a, &argument_interval_size);
    struct RationalInterval out;
    if (!ARGUMENT_ESTIMATE_IS_RADIANS(argument_estimate))
    {
        if (!argument_estimate.min.numerator->value_count)
        {
            local_stack->cursor = local_stack_savepoint;
            return (struct RationalInterval) { rational_zero, rational_zero };
        }
        if (integer_equals(argument_estimate.min.numerator, &one))
        {
            if (integer_equals(argument_estimate.min.denominator, INT(2, 1)))
            {
                local_stack->cursor = local_stack_savepoint;
                return (struct RationalInterval) { rational_zero, rational_zero };
            }
            if (integer_equals(argument_estimate.min.denominator, INT(4, 1)))
            {
                local_stack->cursor = local_stack_savepoint;
                return (struct RationalInterval) { rational_one, rational_one };
            }
        }
        else if (integer_equals(argument_estimate.min.numerator, INT(3, 1)) &&
            integer_equals(argument_estimate.min.denominator, INT(4, 1)))
        {
            local_stack->cursor = local_stack_savepoint;
            out.min.numerator = integer_initialize(output_stack, 1, -1);
            out.min.denominator = &one;
            out.max = out.min;
            return out;
        }
        argument_estimate = revolutions_to_radians(local_stack, output_stack,
            &argument_estimate.min, &argument_interval_size);
    }
    struct RationalInterval argument_min_sine = rational_estimate_sine(local_stack, output_stack,
        &argument_estimate.min, &argument_interval_size);
    struct RationalInterval argument_max_sine = rational_estimate_sine(local_stack, output_stack,
        &argument_estimate.max, &argument_interval_size);
    struct RationalInterval argument_estimate_multiple =
    { rational_double(local_stack, output_stack, &argument_estimate.min),
    rational_double(local_stack, output_stack, &argument_estimate.max) };
    pi_shrink_interval_to_one_side_of_value(&argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(&argument_estimate_multiple.max);
    if (rational_compare(output_stack, local_stack, &argument_estimate_multiple.min,
        &pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, &pi.max, &argument_estimate_multiple.max) <= 0)
    {
        out.max = rational_one;
        out.min = rational_get_min(output_stack, local_stack, &argument_min_sine.min,
            &argument_max_sine.min);
    }
    else
    {
        argument_estimate_multiple.min = rational_integer_divide(local_stack, output_stack,
            &argument_estimate_multiple.min, INT(3, 1));
        argument_estimate_multiple.max = rational_integer_divide(local_stack, output_stack,
            &argument_estimate_multiple.max, INT(3, 1));
        pi_shrink_interval_to_one_side_of_value(&argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(&argument_estimate_multiple.max);
        if (rational_compare(output_stack, local_stack, &argument_estimate_multiple.min,
            &pi.min) <= 0 &&
            rational_compare(output_stack, local_stack, &pi.max,
                &argument_estimate_multiple.max) <= 0)
        {
            out.min.numerator = integer_initialize(output_stack, 1, -1);
            out.min.denominator = &one;
            out.max = rational_get_max(output_stack, local_stack, &argument_min_sine.max,
                &argument_max_sine.max);
        }
        else
        {
            out.min = rational_get_min(output_stack, local_stack, &argument_min_sine.min,
                &argument_max_sine.min);
            out.max = rational_get_max(output_stack, local_stack, &argument_min_sine.max,
                &argument_max_sine.max);
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval number_polar_form_to_rectangular_part(
    struct RationalInterval(trig_function)(struct Stack*restrict, struct Stack*restrict,
        struct Number*, struct Rational*),
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval magnitude_estimate =
        number_get_rational_magnitude_estimate(local_stack, output_stack, a, &rational_one);
    struct Rational sum =
        rational_integer_add(local_stack, output_stack, &magnitude_estimate.max, INT(2, 1));
    struct Rational factor_interval_size =
        rational_divide(local_stack, output_stack, interval_size, &sum);
    magnitude_estimate =
        number_refine_rational_estimate_interval(number_get_rational_magnitude_estimate,
            local_stack, output_stack, a, &magnitude_estimate, &factor_interval_size);
    struct RationalInterval trig_estimate =
        trig_function(local_stack, output_stack, a, &factor_interval_size);
    struct RationalInterval out =
        rational_interval_multiply(output_stack, local_stack, &trig_estimate, &magnitude_estimate);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval number_estimate_part_sum(FloatEstimateGetter get_element_estimate,
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational element_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size,
            size_t_to_integer(local_stack, a->element_count));
    struct FloatInterval out = { float_zero, float_zero };
    struct FloatInterval element_estimate;
    for (size_t i = a->element_count; i-- > 1;)
    {
        element_estimate = get_element_estimate(local_stack, output_stack, a->elements[i],
            &element_estimate_interval_size);
        out = float_interval_add(local_stack, output_stack, &out, &element_estimate);
    }
    element_estimate = get_element_estimate(local_stack, output_stack, a->elements[0],
        &element_estimate_interval_size);
    out = float_interval_add(output_stack, local_stack, &out, &element_estimate);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval number_get_float_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        return rational_get_float_estimate(output_stack, local_stack, &a->value, interval_size);
    case '^':
        return number_get_float_estimate_from_rational_estimate(
            number_get_rational_real_part_estimate, output_stack, local_stack, a, interval_size);
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Region region_estimate =
            number_get_region_estimate(local_stack, output_stack, a, interval_size);
        struct FloatInterval out =
            float_interval_copy(output_stack, &region_estimate.real_interval);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
        return number_estimate_part_sum(number_get_float_real_part_estimate, output_stack,
            local_stack, a, interval_size);
    default:
        crash("Number operation not recognized.");
    }
}

struct RationalInterval number_get_rational_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        return number_polar_form_to_rectangular_part(number_estimate_cosine_of_argument,
            output_stack, local_stack, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval float_estimate =
            number_get_float_real_part_estimate(local_stack, output_stack, a, interval_size);
        struct RationalInterval out =
            float_interval_to_rational_interval(output_stack, local_stack, &float_estimate);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
}

struct FloatInterval number_get_float_imaginary_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        return (struct FloatInterval) { float_zero, float_zero };
    case '^':
        return number_get_float_estimate_from_rational_estimate(
            number_get_rational_imaginary_part_estimate, output_stack, local_stack, a,
            interval_size);
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Region region_estimate =
            number_get_region_estimate(local_stack, output_stack, a, interval_size);
        struct FloatInterval out =
            float_interval_copy(output_stack, &region_estimate.imaginary_interval);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
        return number_estimate_part_sum(number_get_float_imaginary_part_estimate, output_stack,
            local_stack, a, interval_size);
    default:
        crash("Number operation not recognized.");
    }
}

struct RationalInterval number_get_rational_imaginary_part_estimate(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        return number_polar_form_to_rectangular_part(number_estimate_sine_of_argument, output_stack,
            local_stack, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval float_estimate =
            number_get_float_imaginary_part_estimate(local_stack, output_stack, a, interval_size);
        struct RationalInterval out =
            float_interval_to_rational_interval(output_stack, local_stack, &float_estimate);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
}

struct Region number_get_region_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '*')
    {
        struct Region*factor_estimates =
            ARRAY_ALLOCATE(local_stack, a->element_count, struct Region);
        struct Float interval_scale = float_one;
        struct Float term_to_subtract = float_one;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            factor_estimates[i] = number_get_region_estimate(local_stack, output_stack,
                a->elements[i], &rational_one);
            struct Float real_max_magnitude = float_interval_get_max_magnitude(local_stack,
                output_stack, &factor_estimates[i].real_interval);
            struct Float imaginary_max_magnitude = float_interval_get_max_magnitude(local_stack,
                output_stack, &factor_estimates[i].imaginary_interval);
            struct Float bound_max = float_get_max(output_stack, local_stack, &real_max_magnitude,
                &imaginary_max_magnitude);
            struct Float sum = float_add(local_stack, output_stack, &bound_max, &float_one);
            interval_scale = float_multiply(local_stack, output_stack, &interval_scale, &sum);
            term_to_subtract =
                float_multiply(local_stack, output_stack, &term_to_subtract, &bound_max);
        }
        struct Float float_difference =
            float_subtract(local_stack, output_stack, &interval_scale, &term_to_subtract);
        struct Rational difference =
            float_to_rational(local_stack, output_stack, &float_difference);
        struct Rational refined_interval_size =
            rational_divide(local_stack, output_stack, interval_size, &difference);
        struct Region out = { { float_one, float_one }, { float_zero, float_zero } };
        for (size_t i = 0; i < a->element_count; ++i)
        {
            factor_estimates[i].real_interval =
                number_refine_float_estimate_interval(number_get_float_real_part_estimate,
                    local_stack, output_stack, a->elements[i], &factor_estimates[i].real_interval,
                    &refined_interval_size);
            factor_estimates[i].imaginary_interval =
                number_refine_float_estimate_interval(number_get_float_imaginary_part_estimate,
                    local_stack, output_stack, a->elements[i],
                    &factor_estimates[i].imaginary_interval, &refined_interval_size);
            struct FloatInterval product_a = float_interval_multiply(local_stack, output_stack,
                &out.real_interval, &factor_estimates[i].real_interval);
            struct FloatInterval product_b = float_interval_multiply(local_stack, output_stack,
                &out.imaginary_interval, &factor_estimates[i].imaginary_interval);
            struct FloatInterval new_real_part =
                float_interval_subtract(local_stack, output_stack, &product_a, &product_b);
            product_a = float_interval_multiply(local_stack, output_stack, &out.real_interval,
                &factor_estimates[i].imaginary_interval);
            product_b = float_interval_multiply(local_stack, output_stack, &out.imaginary_interval,
                &factor_estimates[i].real_interval);
            out.imaginary_interval =
                float_interval_add(local_stack, output_stack, &product_a, &product_b);
            out.real_interval = new_real_part;
        }
        out.real_interval = float_interval_copy(output_stack, &out.real_interval);
        out.imaginary_interval = float_interval_copy(output_stack, &out.imaginary_interval);
        return out;
    }
    else
    {
        return (struct Region) {
            number_get_float_real_part_estimate(output_stack, local_stack, a, interval_size),
            number_get_float_imaginary_part_estimate(output_stack, local_stack, a, interval_size) };
    }
}

void number_halve_region_estimate_dimensions(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Region*estimate, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float real_interval_size = float_subtract(local_stack, output_stack,
        &estimate->real_interval.max, &estimate->real_interval.min);
    struct Float imaginary_interval_size = float_subtract(local_stack, output_stack,
        &estimate->imaginary_interval.max, &estimate->imaginary_interval.min);
    struct Float interval_size;
    if (!real_interval_size.significand->value_count)
    {
        interval_size = imaginary_interval_size;
    }
    else if (!imaginary_interval_size.significand->value_count)
    {
        interval_size = real_interval_size;
    }
    else
    {
        interval_size = float_get_min(local_stack, output_stack, &real_interval_size,
            &imaginary_interval_size);
    }
    struct Rational rational_interval_size =
        float_to_rational(local_stack, output_stack, &interval_size);
    rational_interval_size = rational_halve(local_stack, output_stack, &rational_interval_size);
    *estimate = number_get_region_estimate(output_stack, local_stack, a, &rational_interval_size);
    local_stack->cursor = local_stack_savepoint;
}

struct FloatInterval number_get_rectangular_part_estimate_for_magnitude_estimate(
    struct EstimateGetters*part_estimate_getters, struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval a_estimate =
        part_estimate_getters->rational(local_stack, output_stack, a, &rational_one);
    struct Rational a_estimate_max_magnitude =
        rational_interval_get_max_magnitude(local_stack, output_stack, &a_estimate);
    struct Rational multiple =
        rational_integer_multiply(local_stack, output_stack, &a_estimate_max_magnitude, INT(4, 1));
    struct Rational sum = rational_integer_add(local_stack, output_stack, &multiple, INT(2, 1));
    struct Rational quotient = rational_divide(local_stack, output_stack, interval_size, &sum);
    struct FloatInterval out = part_estimate_getters->fl(output_stack, local_stack, a, &quotient);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_get_magnitude_estimate_bound_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float*out_min, struct Float*out_max,
    struct Float*real_part_bound, struct Float*imaginary_part_bound, struct Rational*interval_size)
{
    struct Float product_a =
        float_multiply(local_stack, output_stack, real_part_bound, real_part_bound);
    struct Float product_b =
        float_multiply(local_stack, output_stack, imaginary_part_bound, imaginary_part_bound);
    struct Float sum = float_add(local_stack, output_stack, &product_a, &product_b);
    struct FloatInterval root =
        float_estimate_root(output_stack, local_stack, &sum, interval_size, INT(2, 1));
    *out_min = root.min;
    *out_max = root.max;
}

struct FloatInterval number_get_float_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    switch (a->operation)
    {
    case 'r':
    {
        struct Rational a_magnitude = rational_get_magnitude(local_stack, &a->value);
        struct FloatInterval out =
            rational_get_float_estimate(output_stack, local_stack, &a_magnitude, interval_size);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '^':
    {
        struct Rational power =
            rational_exponentiate(local_stack, output_stack, *interval_size, a->index);
        struct Rational quotient =
            rational_integer_divide(local_stack, output_stack, &power, INT(3, 1));
        struct FloatInterval radicand_magnitude_estimate = 
            number_get_float_magnitude_estimate(local_stack, output_stack, a->radicand, &quotient);
        struct Rational refined_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, INT(3, 1));
        struct FloatInterval out = float_estimate_root(local_stack, output_stack, 
            &radicand_magnitude_estimate.min, &refined_interval_size, a->index);
        out.min = float_copy(output_stack, &out.min);
        struct FloatInterval out_max = float_estimate_root(local_stack, output_stack,
            &radicand_magnitude_estimate.max, &refined_interval_size, a->index);
        out.max = float_copy(output_stack, &out_max.max);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        struct FloatInterval*magnitude_estimates =
            ARRAY_ALLOCATE(local_stack, a->element_count, struct FloatInterval);
        struct Float interval_scale = float_one;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            magnitude_estimates[i] = number_get_float_magnitude_estimate(local_stack, output_stack,
                a->elements[i], &rational_one);
            struct Float estimate_max_magnitude = float_interval_get_max_magnitude(local_stack,
                output_stack, magnitude_estimates + i);
            struct Float sum =
                float_add(local_stack, output_stack, &float_one, &estimate_max_magnitude);
            interval_scale = float_multiply(local_stack, output_stack, &interval_scale, &sum);
        }
        struct Rational rational_interval_scale =
            float_to_rational(local_stack, output_stack, &interval_scale);
        struct Rational refined_interval_size =
            rational_divide(local_stack, output_stack, interval_size, &rational_interval_scale);
        struct FloatInterval out = { float_one, float_one };
        for (size_t i = 0; i < a->element_count; ++i)
        {
            magnitude_estimates[i] =
                number_refine_float_estimate_interval(number_get_float_magnitude_estimate,
                    local_stack, output_stack, a->elements[i], magnitude_estimates + i,
                    &refined_interval_size);
            out = float_interval_multiply(local_stack, output_stack, &out, magnitude_estimates + i);
        }
        out = float_interval_copy(output_stack, &out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        struct Rational refined_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, INT(3, 1));
        struct FloatInterval real_part_estimate =
            number_get_rectangular_part_estimate_for_magnitude_estimate(&real_estimate_getters,
                local_stack, output_stack, a, &refined_interval_size);
        struct FloatInterval imaginary_part_estimate =
            number_get_rectangular_part_estimate_for_magnitude_estimate(&imaginary_estimate_getters,
                local_stack, output_stack, a, &refined_interval_size);
        struct Float unused_estimate_bound;
        struct FloatInterval out;
        number_get_magnitude_estimate_bound_interval(local_stack, output_stack, &out.min,
            &unused_estimate_bound, &real_part_estimate.min, &imaginary_part_estimate.min,
            &refined_interval_size);
        number_get_magnitude_estimate_bound_interval(local_stack, output_stack,
            &unused_estimate_bound, &out.max, &real_part_estimate.max,
                &imaginary_part_estimate.max, &refined_interval_size);
        out = float_interval_copy(output_stack, &out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    crash("Number operation not recognized.");
}

struct RationalInterval number_get_rational_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval magnitude_estimate =
        number_get_float_magnitude_estimate(local_stack, output_stack, a, interval_size);
    struct RationalInterval out =
        float_interval_to_rational_interval(output_stack, local_stack, &magnitude_estimate);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval real_number_estimate_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*real_part_estimate, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out = *real_part_estimate;
    while (out.min.numerator->sign == -out.max.numerator->sign)
    {
        out = number_halve_rational_estimate_interval(number_get_rational_real_part_estimate,
            local_stack, output_stack, a, &out);
    }
    out.min = rational_get_argument(output_stack, &out.min);
    ARGUMENT_ESTIMATE_MAKE_REVOLUTIONS(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//If the max field of the return value is 0, the min field is the exact value of the argument in
//revolutions. Otherwise, the return value is an estimate of the argument in radians.
struct RationalInterval number_estimate_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
    {
        return (struct RationalInterval) { rational_get_argument(output_stack, &a->value), { 0 } };
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational product =
            rational_integer_multiply(local_stack, output_stack, interval_size, a->index);
        struct RationalInterval radicand_argument_estimate =
            number_estimate_argument(local_stack, output_stack, a->radicand, &product);
        struct RationalInterval out;
        out.min = rational_integer_divide(output_stack, local_stack,
            &radicand_argument_estimate.min, a->index);
        if (ARGUMENT_ESTIMATE_IS_RADIANS(radicand_argument_estimate))
        {
            out.max = rational_integer_divide(output_stack, local_stack,
                &radicand_argument_estimate.max, a->index);
        }
        else
        {
            ARGUMENT_ESTIMATE_MAKE_REVOLUTIONS(out);
        }
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational factor_argument_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size,
                size_t_to_integer(local_stack, a->element_count));
        struct Rational revolutions = rational_zero;
        struct RationalInterval radians = { rational_zero, rational_zero };
        for (size_t i = 0; i < a->element_count; ++i)
        {
            struct RationalInterval factor_argument = number_estimate_argument(local_stack,
                output_stack, a->elements[i], &factor_argument_interval_size);
            if (ARGUMENT_ESTIMATE_IS_RADIANS(factor_argument))
            {
                radians =
                    rational_interval_add(local_stack, output_stack, &radians, &factor_argument);
            }
            else
            {
                revolutions =
                    rational_add(local_stack, output_stack, &revolutions, &factor_argument.min);
            }
        }
        if (radians.max.numerator->value_count)
        {
            revolutions = rational_double(local_stack, output_stack, &revolutions);
            struct Rational pi_estimate_interval_size = rational_divide(local_stack, output_stack,
                &factor_argument_interval_size, &revolutions);
            pi_estimate(&pi_estimate_interval_size);
            struct RationalInterval out = revolutions_to_radians(local_stack, output_stack,
                &revolutions, &factor_argument_interval_size);
            out = rational_interval_add(output_stack, local_stack, &radians, &out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct RationalInterval out;
        out.min = rational_copy(output_stack, &revolutions);
        ARGUMENT_ESTIMATE_MAKE_REVOLUTIONS(out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval real_part_estimate =
            number_get_rational_real_part_estimate(local_stack, output_stack, a, &rational_one);
        struct RationalInterval imaginary_part_estimate =
            number_get_rational_imaginary_part_estimate(local_stack, output_stack, a,
                &rational_one);
        if (imaginary_part_estimate.min.numerator->sign !=
            imaginary_part_estimate.max.numerator->sign)
        {
            while (true)
            {
                struct Rational imaginary_max_magnitude =
                    rational_interval_get_max_magnitude(local_stack, output_stack,
                        &imaginary_part_estimate);
                if (rational_polynomial_count_roots_in_rectangle(local_stack, output_stack,
                    a->minimal_polynomial, &real_part_estimate,
                    &(struct RationalInterval) { rational_negate(local_stack,
                        &imaginary_max_magnitude), imaginary_max_magnitude }) == 1)
                {
                    struct RationalInterval out = real_number_estimate_argument(output_stack,
                        local_stack, &real_part_estimate, a);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                imaginary_part_estimate = number_halve_rational_estimate_interval(
                    number_get_rational_imaginary_part_estimate, local_stack, output_stack, a,
                        &imaginary_part_estimate);
                if (imaginary_part_estimate.min.numerator->sign ==
                    imaginary_part_estimate.max.numerator->sign)
                {
                    break;
                }
                real_part_estimate =
                    number_halve_rational_estimate_interval(number_get_rational_real_part_estimate,
                        local_stack, output_stack, a, &real_part_estimate);
            }
        }
        if (!imaginary_part_estimate.min.numerator->sign &&
            !imaginary_part_estimate.max.numerator->sign)
        {
            struct RationalInterval out =
                real_number_estimate_argument(output_stack, local_stack, &real_part_estimate, a);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Rational refined_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, INT(3, +1));
        struct Rational imaginary_part_interval_size =
            rational_interval_get_factor_interval_size(local_stack, output_stack,
                &imaginary_part_estimate, &(struct RationalInterval) {
                rational_get_reciprocal(local_stack, &real_part_estimate.max),
                rational_get_reciprocal(local_stack, &real_part_estimate.min) },
                &refined_interval_size);
        imaginary_part_estimate =
            number_refine_rational_estimate_interval(number_get_rational_imaginary_part_estimate,
                local_stack, output_stack, a, &imaginary_part_estimate,
                &imaginary_part_interval_size);
        struct Rational real_part_interval_size = rational_multiply(local_stack, output_stack,
            &real_part_estimate.min, &imaginary_part_interval_size);
        struct Rational product = rational_multiply(local_stack, output_stack,
            &real_part_estimate.min, &real_part_interval_size);
        struct Rational difference =
            rational_subtract(local_stack, output_stack, &rational_one, &real_part_interval_size);
        real_part_interval_size = rational_divide(local_stack, output_stack, &product, &difference);
        real_part_estimate =
            number_refine_rational_estimate_interval(number_get_rational_real_part_estimate,
                local_stack, output_stack, a, &real_part_estimate, &real_part_interval_size);
        struct RationalInterval out = rational_interval_estimate_atan2(output_stack, local_stack,
            &imaginary_part_estimate, &real_part_estimate, &refined_interval_size);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    crash("Number operation not recognized.");
}

struct RationalInterval number_estimate_argument_in_radians(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out =
        number_estimate_argument(local_stack, output_stack, a, interval_size);
    if (ARGUMENT_ESTIMATE_IS_RADIANS(out))
    {
        out = rational_interval_copy(output_stack, &out);
    }
    else
    {
        out = revolutions_to_radians(output_stack, local_stack, &out.min, interval_size);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}