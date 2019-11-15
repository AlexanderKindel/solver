#include "declarations.h"

size_t integer_size(size_t value_count)
{
    return sizeof(struct Integer) + value_count * sizeof(uint32_t);
}

struct Integer*integer_allocate(struct Stack*output_stack, size_t value_count)
{
    return stack_slot_allocate(output_stack, integer_size(value_count), _Alignof(struct Integer));
}

struct Integer*integer_copy(struct Stack*output_stack, struct Integer*a)
{
    size_t size = integer_size(a->value_count);
    struct Integer*out = stack_slot_allocate(output_stack, size, _Alignof(struct Integer));
    memcpy(out, a, size);
    return out;
}

struct Integer*integer_initialize(struct Stack*stack, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_allocate(stack, 1);
    out->value_count = 1;
    out->sign = sign;
    out->value[0] = value;
    return out;
}

struct Integer*integer_from_char(struct Stack*output_stack, char value)
{
    uint32_t uint_value = value - '0';
    if (uint_value)
    {
        struct Integer*out = integer_allocate(output_stack, 1);
        out->value_count = 1;
        out->sign = 1;
        out->value[0] = uint_value;
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

struct Integer*integer_from_size_t(struct Stack*output_stack, size_t value)
{
    size_t value_count = sizeof(size_t) / sizeof(uint32_t);
    struct Integer*out = integer_allocate(output_stack, value_count);
    out->value_count = value_count;
    out->sign = 1;
    *(size_t*)out->value = value;
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

bool integer_equals_zero(struct Integer*a)
{
    return a->value_count == 0;
}

struct Integer*integer_magnitude(void*output_arena, struct Integer*a)
{
    struct Integer*out = integer_copy(output_arena, a);
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

void twos_complement(struct Integer*a)
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
void integer_add_to_a_in_place(struct Integer*a, struct Integer*b)
{
    if (b->sign == 0)
    {
        return;
    }
    if (a->sign == 0)
    {
        memcpy(a, b, integer_size(b->value_count));
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
        twos_complement(a);
        calculate_sum_values(b, a);
        if (a->value[a->value_count - 1] != 0)
        {
            twos_complement(a);
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

struct Integer*integer_negative(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    out->sign = -out->sign;
    return out;
}

struct Integer*integer_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*minuend, struct Integer*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out =
        integer_add(output_stack, minuend, integer_negative(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = stack_slot_allocate(output_stack,
        integer_size(a->value_count + b->value_count), _Alignof(struct Integer));
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
            if (product_component > 0)
            {
                integer_component->value[shift] = (uint32_t)product_component;
                integer_component->value_count += 1;
                integer_component->sign = 1;
                uint32_t high_bytes = (product_component & 0xffffffff00000000) >> 32;
                if (high_bytes > 0)
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

size_t leading_digit_place(struct Integer*a)
{
    if (!a->value_count)
    {
        return 0;
    }
    uint32_t divisor_leading_digit = 0x80000000;
    uint32_t last_value = a->value[a->value_count - 1];
    size_t i = 31;
    while (true)
    {
        if ((last_value & divisor_leading_digit) != 0)
        {
            return i;
        }
        divisor_leading_digit = divisor_leading_digit >> 1;
        --i;
    }
}

void integer_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerDivision*out, struct Integer*dividend, struct Integer*divisor)
{
    size_t dividend_leading_digit_place = leading_digit_place(dividend);
    size_t divisor_leading_digit_place = leading_digit_place(divisor);
    if (dividend->value_count > divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
            dividend_leading_digit_place >= divisor_leading_digit_place))
    {
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
        out->quotient = integer_allocate(output_stack, dividend->value_count);
        out->quotient->value_count = dividend->value_count;
        out->quotient->sign = dividend->sign * divisor->sign;
        memset(&out->quotient->value, 0, out->quotient->value_count * sizeof(uint32_t));
        out->remainder = integer_magnitude(local_stack, dividend);
        while (true)
        {
            for (int i = 32; i > 0; --i)
            {
                struct Integer*difference =
                    integer_subtract(local_stack, output_stack, out->remainder, divisor_magnitude);
                if (difference->sign >= 0)
                {
                    out->quotient->value[quotient_value_index] |= quotient_digit;
                    out->remainder = difference;
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
                integer_halve(divisor_magnitude);
            }
        }
        break_both_loops:
        integer_trim_leading_zeroes(out->quotient);
        output_stack->cursor = &out->quotient->value[out->quotient->value_count];
        out->remainder = integer_copy(output_stack, out->remainder);
        if (out->remainder->sign != 0)
        {
            out->remainder->sign = dividend->sign;
        }
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        out->quotient = &zero;
        out->remainder = integer_copy(output_stack, dividend);
    }
}

struct Integer*integer_euclidean_quotient(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*dividend, struct Integer*divisor)
{
    struct IntegerDivision division;
    integer_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return division.quotient;
}

struct Integer*integer_euclidean_remainder(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*dividend, struct Integer*divisor)
{
    struct IntegerDivision division;
    integer_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return division.remainder;
}

struct Integer*integer_doubled(struct Stack*output_stack, struct Integer*a)
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

void integer_halve(struct Integer*a)
{
    integer_downshift(a, 1);
    integer_trim_leading_zeroes(a);
}

struct Integer*integer_half(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    integer_halve(out);
    return out;
}

int8_t integer_compare(struct Stack*stack_a, struct Stack*stack_b, struct Integer*a,
    struct Integer*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = integer_subtract(stack_a, stack_b, a, b)->sign;
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Integer*integer_lcm(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a,
    struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = integer_multiply(output_stack, local_stack,
        integer_euclidean_quotient(local_stack, output_stack, a,
            integer_gcd(local_stack, output_stack, a, b)), b);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Rounded up to the nearest integer.
struct Integer*integer_square_root(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*float_square_root_min;
    struct Float*float_square_root_max;
    float_estimate_root(local_stack, output_stack, &float_square_root_min, &float_square_root_max,
        &(struct Float){ a, &zero }, &rational_one, INT(2, 1));
    struct Rational*rational_square_root =
        float_to_rational(local_stack, output_stack, float_square_root_max);
    struct Integer*out;
    if (integer_equals(rational_square_root->denominator, &one))
    {
        out = integer_copy(output_stack, rational_square_root->numerator);
    }
    else
    {
        out = integer_add(output_stack, &one,
            integer_euclidean_quotient(local_stack, output_stack, rational_square_root->numerator,
                rational_square_root->denominator));
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*get_prime(struct Stack*stack_a, struct Stack*stack_b, size_t index)
{
    if (primes[index])
    {
        return primes[index];
    }
    void*stack_a_savepoint = stack_a->cursor;
    primes[index] = integer_add(stack_a, primes[index - 1], primes[0]);
    while (true)
    {
        size_t potential_factor_index = 1;
        while (true)
        {
            if ((integer_euclidean_remainder(stack_a, stack_b, primes[index],
                primes[potential_factor_index]))->value_count == 0)
            {
                break;
            }
            ++potential_factor_index;
            if (potential_factor_index == index)
            {
                stack_a->cursor = stack_a_savepoint;
                primes[index] = integer_copy(&permanent_stack, primes[index]);
                return primes[index];
            }
        }
        primes[index] = integer_add(stack_a, primes[index], primes[0]);
    }
}

size_t size_t_factor(struct Stack*output_stack, struct Stack*local_stack, size_t**out, size_t a)
{
    void*local_stack_savepoint = local_stack->cursor;
    *out = ALLOCATE(output_stack, size_t);
    size_t factor_count = 0;
    size_t prime_index = 0;
    size_t prime = 2;
    size_t square_root = integer_to_size_t(integer_square_root(local_stack, output_stack,
        integer_from_size_t(local_stack, a)));
    while (prime <= square_root)
    {
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
        ++prime_index;
        prime = integer_to_size_t(get_prime(output_stack, local_stack, prime_index));
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

size_t integer_factor(struct Stack*output_stack, struct Stack*local_stack, struct Factor**out,
    struct Integer*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    *out = array_start(output_stack, _Alignof(struct Factor));
    size_t factor_count = 0;
    size_t prime_index = 0;
    struct Integer*prime = primes[0];
    struct Integer*square_root = integer_square_root(local_stack, output_stack, a);
    while (integer_compare(output_stack, local_stack, prime, square_root) <= 0)
    {
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division, a, prime);
        if (division.remainder->value_count == 0)
        {
            extend_array(output_stack, sizeof(struct Factor));
            struct Factor*factor = *out + factor_count;
            factor->value = prime;
            factor->multiplicity = &zero;
            factor_count += 1;
            do
            {
                factor->multiplicity = integer_add(local_stack, factor->multiplicity, &one);
                a = division.quotient;
                integer_euclidean_divide(local_stack, output_stack, &division, a, factor->value);
            } while (division.remainder->value_count == 0);
        }
        ++prime_index;
        prime = get_prime(output_stack, local_stack, prime_index);
    }
    if (!integer_equals(a, &one))
    {
        extend_array(output_stack, sizeof(struct Factor));
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

size_t integer_string(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a)
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
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division, quotient, ten);
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
