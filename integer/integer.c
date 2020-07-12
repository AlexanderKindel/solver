#include "declarations.h"

size_t integer_get_size(size_t value_count)
{
    return sizeof(struct Integer) + value_count * sizeof(uint32_t);
}

struct Integer*integer_allocate(struct Stack*output_stack, size_t value_count)
{
    return stack_slot_allocate(output_stack, integer_get_size(value_count),
        _Alignof(struct Integer));
}

struct Integer*integer_copy(struct Stack*output_stack, struct Integer*a)
{
    size_t size = integer_get_size(a->value_count);
    struct Integer*out = stack_slot_allocate(output_stack, size, _Alignof(struct Integer));
    memcpy(out, a, size);
    return out;
}

struct Integer*integer_initialize(struct Stack*output_stack, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_allocate(output_stack, 1);
    out->value_count = 1;
    out->sign = sign;
    out->value[0] = value;
    return out;
}

struct Integer*char_to_integer(struct Stack*output_stack, char a)
{
    uint32_t uint_a = a - '0';
    if (uint_a)
    {
        struct Integer*out = integer_allocate(output_stack, 1);
        out->value_count = 1;
        out->sign = 1;
        out->value[0] = uint_a;
        return out;
    }
    else
    {
        return &zero;
    }
}

void integer_trim_leading_zeroes(struct Integer*a)
{
    for (size_t i = a->value_count; i > 0; --i)
    {
        if (a->value[i - 1] != 0)
        {
            a->value_count = i;
            return;
        }
    }
    a->value_count = 0;
    a->sign = 0;
}

struct Integer*size_t_to_integer(struct Stack*output_stack, size_t a)
{
    size_t value_count = sizeof(size_t) / sizeof(uint32_t);
    struct Integer*out = integer_allocate(output_stack, value_count);
    out->value_count = value_count;
    out->sign = 1;
    *(size_t*)out->value = a;
    integer_trim_leading_zeroes(out);
    output_stack->cursor = out->value + out->value_count;
    return out;
}

size_t integer_to_size_t(struct Integer*a)
{
    ASSERT(a->sign >= 0, "integer_to_size_t called on a negative Integer.");
    switch (a->value_count)
    {
    case 0:
        return 0;
    case 1:
        return a->value[0];
    case 2:
        return *(size_t*)a->value;
    default:
        crash("integer_to_size_t called on an Integer too large to fit in a size_t.");
    }
}

bool integer_equals(struct Integer*a, struct Integer*b)
{
    if (a->sign != b->sign)
    {
        return false;
    }
    if (a->value_count != b->value_count)
    {
        return false;
    }
    for (int i = 0; i < a->value_count; ++i)
    {
        if (a->value[i] != b->value[i])
        {
            return false;
        }
    }
    return true;
}

struct Integer*integer_get_magnitude(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    if (out->sign < 0)
    {
        out->sign = 1;
    }
    return out;
}

void calculate_sum_values(struct Integer*short_integer, struct Integer*sum)
{
    for (int i = 0; i < short_integer->value_count; ++i)
    {
        uint64_t remainder = short_integer->value[i];
        for (int j = i; j < sum->value_count; ++j)
        {
            uint64_t sum_value = sum->value[j] + remainder;
            sum->value[j] = (uint32_t)sum_value;
            remainder = (sum_value & 0xffffffff00000000) >> 32;
            if (remainder == 0)
            {
                break;
            }
        }
    }
}

void integer_get_twos_complement(struct Integer*a)
{
    for (size_t i = 0; i < a->value_count; ++i)
    {
        a->value[i] = ~a->value[i];
    }
    for (size_t i = 0; i < a->value_count; ++i)
    {
        uint32_t power = 1;
        for (size_t j = 0; j < 32; ++j)
        {
            a->value[i] ^= power;
            if ((a->value[i] & power) != 0)
            {
                return;
            }
            power = power << 1;
        }
    }
}

//Assumes that, at a->values, there is enough allocated but unused memory for a->value_count to
//increase to max(a->value_count, b->value_count) + 1.
void integer_add_to_a_in_place(struct Integer*restrict a, struct Integer*restrict b)
{
    if (b->sign == 0)
    {
        return;
    }
    if (a->sign == 0)
    {
        memcpy(a, b, integer_get_size(b->value_count));
        return;
    }
    if (a->value_count < b->value_count)
    {
        memset(&a->value[a->value_count], 0,
            (b->value_count - a->value_count + 1) * sizeof(uint32_t));
        a->value_count = b->value_count + 1;
    }
    else
    {
        a->value[a->value_count] = 0;
        a->value_count += 1;
    }
    if (a->sign == b->sign)
    {
        calculate_sum_values(b, a);
    }
    else
    {
        a->sign = b->sign;
        integer_get_twos_complement(a);
        calculate_sum_values(b, a);
        if (a->value[a->value_count - 1] != 0)
        {
            integer_get_twos_complement(a);
            a->sign *= -1;
        }
    }
    integer_trim_leading_zeroes(a);
}

struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b)
{
    struct Integer*out = integer_copy(output_stack, a);
    stack_slot_allocate(output_stack, (b->value_count + 1) * sizeof(uint32_t),
        _Alignof(struct Integer));
    integer_add_to_a_in_place(out, b);
    output_stack->cursor = &out->value[out->value_count];
    return out;
}

struct Integer*integer_negate(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    out->sign = -out->sign;
    return out;
}

struct Integer*integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out =
        integer_add(output_stack, minuend, integer_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = stack_slot_allocate(output_stack,
        integer_get_size(a->value_count + b->value_count), _Alignof(struct Integer));
    out->value_count = 0;
    out->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)a->value[i] * b->value[j];
            size_t shift = i + j;
            struct Integer*integer_component = integer_allocate(local_stack, shift + 2);
            integer_component->value_count = shift;
            integer_component->sign = 0;
            memset(&integer_component->value, 0, shift * sizeof(uint32_t));
            if (product_component)
            {
                integer_component->value[shift] = (uint32_t)product_component;
                integer_component->value_count += 1;
                integer_component->sign = 1;
                uint32_t high_bytes = (product_component & 0xffffffff00000000) >> 32;
                if (high_bytes)
                {
                    integer_component->value[shift + 1] = high_bytes;
                    integer_component->value_count += 1;
                }
            }
            integer_add_to_a_in_place(out, integer_component);
        }
    }
    out->sign = a->sign * b->sign;
    output_stack->cursor = &out->value[out->value_count];
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void integer_upshift(struct Integer*a, uint8_t shift)
{
    for (size_t i = a->value_count; i-- > 1;)
    {
        a->value[i] = a->value[i] << shift | a->value[i - 1] >> (32 - shift);
    }
    a->value[0] = a->value[0] << shift;
}

void integer_downshift(struct Integer*a, uint8_t shift)
{
    for (size_t i = 0; i < a->value_count - 1; ++i)
    {
        a->value[i] = a->value[i] >> shift | a->value[i + 1] << (32 - shift);
    }
    a->value[a->value_count - 1] = a->value[a->value_count - 1] >> shift;
}

uint32_t integer_get_leading_digit_place(struct Integer*a)
{
    if (a->value_count)
    {
        return get_leading_digit_place(a->value[a->value_count - 1]);
    }
    return 0;
}

struct IntegerDivision integer_euclidean_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*dividend, struct Integer*divisor)
{
    uint32_t dividend_leading_digit_place = integer_get_leading_digit_place(dividend);
    uint32_t divisor_leading_digit_place = integer_get_leading_digit_place(divisor);
    if (dividend->value_count < divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
            dividend_leading_digit_place < divisor_leading_digit_place))
    {
        return (struct IntegerDivision) { &zero, integer_copy(output_stack, dividend) };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*divisor_magnitude = integer_allocate(local_stack, dividend->value_count);
    divisor_magnitude->value_count = dividend->value_count;
    divisor_magnitude->sign = 1;
    size_t quotient_value_index = dividend->value_count - divisor->value_count;
    memset(&divisor_magnitude->value, 0, quotient_value_index * sizeof(uint32_t));
    memcpy(&divisor_magnitude->value[quotient_value_index], &divisor->value,
        divisor->value_count * sizeof(uint32_t));
    int shift = dividend_leading_digit_place - divisor_leading_digit_place;
    uint32_t quotient_digit;
    if (shift > 0)
    {
        integer_upshift(divisor_magnitude, shift);
        quotient_digit = 1 << shift;
    }
    else if (shift < 0)
    {
        shift *= -1;
        integer_downshift(divisor_magnitude, shift);
        quotient_digit = 1 << (32 - shift);
        quotient_value_index -= 1;
    }
    else
    {
        quotient_digit = 1;
    }
    struct IntegerDivision out = { integer_allocate(output_stack, dividend->value_count),
        integer_get_magnitude(local_stack, dividend) };
    out.quotient->value_count = dividend->value_count;
    out.quotient->sign = dividend->sign * divisor->sign;
    memset(&out.quotient->value, 0, out.quotient->value_count * sizeof(uint32_t));
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer*difference =
                integer_subtract(local_stack, output_stack, out.remainder, divisor_magnitude);
            if (difference->sign >= 0)
            {
                out.quotient->value[quotient_value_index] |= quotient_digit;
                out.remainder = difference;
            }
            if (quotient_digit == 1)
            {
                if (quotient_value_index == 0)
                {
                    goto break_both_loops;
                }
                quotient_digit = 0x80000000;
                quotient_value_index -= 1;
            }
            else
            {
                quotient_digit = quotient_digit >> 1;
            }
            integer_halve_in_place(divisor_magnitude);
        }
    }
break_both_loops:
    integer_trim_leading_zeroes(out.quotient);
    output_stack->cursor = out.quotient->value + out.quotient->value_count;
    out.remainder = integer_copy(output_stack, out.remainder);
    if (out.remainder->sign != 0)
    {
        out.remainder->sign = dividend->sign;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*integer_double(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out;
    if (a->value[a->value_count - 1] & 0x80000000)
    {
        out = integer_allocate(output_stack, a->value_count + 1);
        out->value_count = a->value_count + 1;
        out->value[a->value_count] = 0;
    }
    else
    {
        out = integer_allocate(output_stack, a->value_count);
        out->value_count = a->value_count;
    }
    out->sign = a->sign;
    memcpy(&out->value, &a->value, a->value_count * sizeof(uint32_t));
    integer_upshift(out, 1);
    return out;
}

void integer_halve_in_place(struct Integer*a)
{
    integer_downshift(a, 1);
    integer_trim_leading_zeroes(a);
}

struct Integer*integer_halve(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    integer_halve_in_place(out);
    return out;
}

struct Integer*integer_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return &one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = &one;
    while (true)
    {
        if (exponent->value[0] & 1)
        {
            out = integer_multiply(local_stack, output_stack, out, base);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = integer_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = integer_multiply(local_stack, output_stack, base, base);
    }
}

//Leaves a significant amount of excess allocations on output_stack.
struct ExtendedGCDInfo integer_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo out = { b, &zero, &one, &zero, &one };
    while (!integer_equals(a, &zero))
    {
        struct IntegerDivision division =
            integer_euclidean_divide(output_stack, local_stack, out.gcd, a);
        struct Integer*m = integer_subtract(output_stack, local_stack, out.a_coefficient,
            integer_multiply(local_stack, output_stack, out.b_over_gcd, division.quotient));
        struct Integer*n = integer_subtract(output_stack, local_stack, out.b_coefficient,
            integer_multiply(local_stack, output_stack, out.a_over_gcd, division.quotient));
        out.gcd = a;
        a = division.remainder;
        out.a_coefficient = out.b_over_gcd;
        out.b_coefficient = out.a_over_gcd;
        out.b_over_gcd = m;
        out.a_over_gcd = n;
    }
    out.b_over_gcd = integer_negate(output_stack, out.b_over_gcd);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*integer_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!integer_equals(b, &zero))
    {
        struct Integer*c = b;
        b = integer_euclidean_divide(local_stack, output_stack, a, b).remainder;
        a = c;
    }
    a = integer_copy(output_stack, a);
    local_stack->cursor = local_stack_savepoint;
    return a;
}

int8_t integer_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Integer*a, struct Integer*b)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    int8_t out = integer_subtract(local_stack_a, local_stack_b, a, b)->sign;
    local_stack_a->cursor = local_stack_a_savepoint;
    return out;
}

struct Integer*integer_get_lcm(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = integer_multiply(output_stack, local_stack,
        integer_euclidean_divide(local_stack, output_stack, a,
            integer_get_gcd(local_stack, output_stack, a, b)).quotient,
        b);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Rounded up to the nearest integer.
struct Integer*integer_take_square_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval float_square_root = float_estimate_root(local_stack, output_stack,
        &(struct Float) { a, &zero }, &rational_one, INT(2, 1));
    struct Rational rational_square_root =
        float_to_rational(local_stack, output_stack, &float_square_root.max);
    struct Integer*out;
    if (integer_equals(rational_square_root.denominator, &one))
    {
        out = integer_copy(output_stack, rational_square_root.numerator);
    }
    else
    {
        out = integer_add(output_stack, &one,
            integer_euclidean_divide(local_stack, output_stack, rational_square_root.numerator,
                rational_square_root.denominator).quotient);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*get_next_prime(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b)
{
    if (next_prime < prime_stack.cursor)
    {
        struct Integer*out = *next_prime;
        ++next_prime;
        return out;
    }
    void*local_stack_a_savepoint = local_stack_a->cursor;
    *next_prime =
        integer_add(local_stack_a, *(next_prime - 1), *(struct Integer**)prime_stack.start);
    while (true)
    {
        struct Integer**potential_factor = prime_stack.start;
        while (true)
        {
            if ((integer_euclidean_divide(local_stack_a, local_stack_b, *next_prime,
                *potential_factor).remainder)->value_count == 0)
            {
                break;
            }
            ++potential_factor;
            if (potential_factor == next_prime)
            {
                *next_prime = integer_copy(&permanent_stack, *next_prime);
                struct Integer*out = *next_prime;
                ++next_prime;
                local_stack_a->cursor = local_stack_a_savepoint;
                return out;
            }
        }
        *next_prime = integer_add(local_stack_a, *next_prime, *(struct Integer**)prime_stack.start);
    }
}

size_t size_t_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    size_t**out, size_t a)
{
    void*local_stack_savepoint = local_stack->cursor;
    *out = ALLOCATE(output_stack, size_t);
    size_t factor_count = 0;
    next_prime = prime_stack.start;
    size_t square_root = integer_to_size_t(integer_take_square_root(local_stack, output_stack,
        size_t_to_integer(local_stack, a)));
    while (true)
    {
        size_t prime = integer_to_size_t(get_next_prime(output_stack, local_stack));
        if (prime > square_root)
        {
            break;
        }
        if (a % prime == 0)
        {
            (*out)[factor_count] = prime;
            factor_count += 1;
            ALLOCATE(output_stack, size_t);
            do
            {
                a = a / prime;
            } while (a % prime == 0);
        }
    }
    if (a != 1)
    {
        (*out)[factor_count] = a;
    }
    else
    {
        output_stack->cursor = *out + factor_count;
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t integer_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Factor**out, struct Integer*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    *out = align_cursor(output_stack, _Alignof(struct Factor));
    size_t factor_count = 0;
    next_prime = prime_stack.start;
    struct Integer*square_root = integer_take_square_root(local_stack, output_stack, a);
    while (true)
    {
        struct Integer*prime = get_next_prime(output_stack, local_stack);
        if (integer_compare(output_stack, local_stack, prime, square_root) > 0)
        {
            break;
        }
        struct IntegerDivision division =
            integer_euclidean_divide(local_stack, output_stack, a, prime);
        if (division.remainder->value_count == 0)
        {
            unaligned_stack_slot_allocate(output_stack, sizeof(struct Factor));
            struct Factor*factor = *out + factor_count;
            factor->value = prime;
            factor->multiplicity = &zero;
            factor_count += 1;
            do
            {
                factor->multiplicity = integer_add(local_stack, factor->multiplicity, &one);
                a = division.quotient;
                division = integer_euclidean_divide(local_stack, output_stack, a, factor->value);
            } while (division.remainder->value_count == 0);
        }
    }
    if (!integer_equals(a, &one))
    {
        unaligned_stack_slot_allocate(output_stack, sizeof(struct Factor));
        (*out)[factor_count].value = a;
        (*out)[factor_count].multiplicity = &one;
        factor_count += 1;
    }
    for (size_t i = 0; i < factor_count; ++i)
    {
        (*out)[i].value = integer_copy(output_stack, (*out)[i].value);
        (*out)[i].multiplicity = integer_copy(output_stack, (*out)[i].multiplicity);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t integer_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*a)
{
    if (a->sign == 0)
    {
        *(char*)ALLOCATE(output_stack, char) = '0';
        return 1;
    }
    void*local_stack_savepoint = local_stack->cursor;
    char*buffer_start = stack_slot_allocate(output_stack, 10 * a->value_count + 1, 1);
    char*next_char = output_stack->cursor;
    struct Integer*quotient = a;
    struct Integer*ten = INT(10, 1);
    while (quotient->sign != 0)
    {
        struct IntegerDivision division =
            integer_euclidean_divide(local_stack, output_stack, quotient, ten);
        quotient = division.quotient;
        next_char -= 1;
        if (division.remainder->sign != 0)
        {
            *next_char = division.remainder->value[0] + '0';
        }
        else
        {
            *next_char = '0';
        }
    }
    size_t char_count = (char*)output_stack->cursor - next_char;
    if (a->sign < 0)
    {
        *buffer_start = '-';
        memcpy(buffer_start + 1, next_char, char_count);
        char_count += 1;
    }
    else
    {
        memcpy(buffer_start, next_char, char_count);
    }
    output_stack->cursor = buffer_start + char_count;
    local_stack->cursor = local_stack_savepoint;
    return char_count;
}

#include "integer/float/float.c"
#include "integer/integer_polynomial/integer_polynomial.c"
#include "integer/modded_integer/modded_integer.c"
#include "integer/rational/rational.c"
