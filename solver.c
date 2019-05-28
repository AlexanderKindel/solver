#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#ifdef _DEBUG
#define ABORT(message) printf(message); abort();
#else
#define ABORT(message)
#endif

DWORD page_size;
size_t arena_size = 8192;

size_t next_page_boundary(size_t address)
{
    size_t distance_from_page_start = address % page_size;
    if (distance_from_page_start)
    {
        return address + page_size - distance_from_page_start;
    }
    return address;
}

void*stack_slot_new(void**stack_cursor, size_t slot_size)
{
    void*slot = *stack_cursor;
    size_t boundary = next_page_boundary(*stack_cursor);
    *stack_cursor = (size_t)*stack_cursor + slot_size;
    if (*stack_cursor > boundary)
    {
        VirtualAlloc(boundary, (size_t)*stack_cursor - boundary, MEM_COMMIT, PAGE_READWRITE);
    }
    return slot;
}

void rewind_stack_cursor(void**cursor, void*cursor_target)
{
    size_t boundary = next_page_boundary(cursor_target);
    if (*cursor > boundary)
    {
        VirtualFree(boundary, (size_t)*cursor - boundary, MEM_DECOMMIT);
    }
    *cursor = cursor_target;
}

struct Integer
{
    size_t value_count;
    int8_t sign;
    uint32_t value;//The first element of an array of length value_count.
};

struct Integer zero = { 0, 0, 0 };
#define INT(value, sign) (struct Integer){ 1, sign, value }

size_t integer_size(size_t value_count)
{
    return sizeof(struct Integer) + (value_count - 1) * sizeof(uint32_t);
}

struct Integer*stack_integer_slot_new(void**stack_cursor, size_t value_count)
{
    return stack_slot_new(stack_cursor, integer_size(value_count));
}

struct IntegerSlot
{
    struct IntegerSlot*next_slot;
    struct Integer integer;
};

struct IntegerSlotPage
{
    struct IntegerSlotPage*next_page;
    struct IntegerSlot memory;
};

struct IntegerPool
{
    struct IntegerSlotPage*first_page;
    struct IntegerSlot*cursor;
    struct IntegerSlot*free_list;
};

struct IntegerPool*permanent_integer_pool;

void integer_pool_new(struct IntegerPool*out)
{
    out->first_page = VirtualAlloc(0, page_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    out->cursor = &out->first_page->memory;
}

struct IntegerPool*get_value_count_pool(struct IntegerPool*pool, size_t value_count)
{
    if (value_count == 1)
    {
        return pool;
    }
    struct IntegerPool*next_pool = pool + 1;
    if (!next_pool->first_page)
    {
        integer_pool_new(next_pool);
    }
    return get_value_count_pool(next_pool, value_count - 1);
}

struct Integer*pool_integer_slot_new(struct IntegerPool*pool, size_t value_count)
{
    if (!value_count)
    {
        return &zero;
    }
    struct IntegerPool*value_count_pool = get_value_count_pool(pool, value_count);
    struct IntegerSlot*slot = pool->free_list;
    if (slot)
    {
        pool->free_list = slot->next_slot;
    }
    else
    {
        slot = pool->cursor;
        pool->cursor = (size_t)(pool->cursor + 1) + (value_count - 1) * sizeof(uint32_t);
        size_t boundary = next_page_boundary(slot);
        if (pool->cursor > boundary)
        {
            struct IntegerSlotPage*new_page =
                VirtualAlloc(0, page_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
            pool->cursor = &new_page->memory;
            ((struct IntegerSlotPage*)(boundary - page_size))->next_page = new_page;
        }
    }
    return &slot->integer;
}

void pool_integer_free(struct IntegerPool*pool, struct Integer*a)
{
    if (a->value_count)
    {
        struct IntegerPool*value_count_pool = get_value_count_pool(pool, a->value_count);
        struct IntegerSlot*slot = (size_t)a - sizeof(struct IntegerSlot*);
        slot->next_slot = value_count_pool->free_list;
        value_count_pool->free_list = slot;
    }
}

struct Integer*integer_new(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    uint32_t value, int8_t sign)
{
    struct Integer*out = integer_slot_allocator(memory, 1);
    if (value == 0)
    {
        out->value_count = 0;
        out->sign = 0;
    }
    else
    {
        out->value_count = 1;
        out->value = value;
        out->sign = sign;
    }
    return out;
}

struct Integer*integer_copy(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    struct Integer*a)
{
    struct Integer*copy = integer_slot_allocator(memory, a->value_count);
    memcpy(copy, a, integer_size(a->value_count));
    return copy;
}

void integer_move(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    struct Integer**a)
{
    *a = integer_copy(integer_slot_allocator, memory, *a);
}

void integer_move_from_pool(struct IntegerPool*pool, void**stack_cursor, struct Integer**a)
{
    struct Integer*pool_copy = *a;
    integer_move(stack_integer_slot_new, stack_cursor, a);
    pool_integer_free(pool, pool_copy);
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
        if ((&a->value)[i] != (&b->value)[i])
        {
            return false;
        }
    }
    return true;
}

bool integer_equals_one(struct Integer*a)
{
    return a->value_count == 1 && a->sign > 0 && a->value == 1;
}

void calculate_sum_values(struct Integer*short_integer, struct Integer*sum)
{
    for (int i = 0; i < short_integer->value_count; ++i)
    {
        uint64_t remainder = (&short_integer->value)[i];
        for (int j = i; j < sum->value_count; ++j)
        {
            uint64_t sum_value = (&sum->value)[j] + remainder;
            (&sum->value)[j] = (uint32_t)sum_value;
            remainder = (sum_value & 0xffffffff00000000) >> 32;
            if (remainder == 0)
            {
                break;
            }
        }
    }
}

void twos_complement(struct Integer*integer)
{
    for (int i = 0; i < integer->value_count; ++i)
    {
        (&integer->value)[i] = ~(&integer->value)[i];
    }
    for (int i = 0; i < integer->value_count; ++i)
    {
        uint32_t power = 1;
        for (int j = 0; j < 32; ++j)
        {
            (&integer->value)[i] ^= power;
            if (((&integer->value)[i] & power) != 0)
            {
                return;
            }
            power = power << 1;
        }
    }
}

void trim_leading_zeroes(struct Integer*a)
{
    for (int i = a->value_count - 1; i >= 0; --i)
    {
        if ((&a->value)[i] != 0)
        {
            a->value_count = i + 1;
            return;
        }
    }
    a->value_count = 0;
    a->sign = 0;
}

struct Integer*integer_add(void**stack_cursor, struct Integer*a, struct Integer*b)
{
    if (a->sign == 0)
    {
        return integer_copy(stack_integer_slot_new, stack_cursor, b);
    }
    if (b->sign == 0)
    {
        return integer_copy(stack_integer_slot_new, stack_cursor, a);
    }
    struct Integer*long_integer = a;
    struct Integer*short_integer = b;
    if (a->value_count < b->value_count)
    {
        long_integer = b;
        short_integer = a;
    }
    struct Integer*sum = stack_integer_slot_new(stack_cursor, long_integer->value_count + 1);
    sum->value_count = long_integer->value_count + 1;
    sum->sign = short_integer->sign;
    memcpy(&sum->value, &long_integer->value, long_integer->value_count * sizeof(uint32_t));
    (&sum->value)[long_integer->value_count] = 0;
    if (short_integer->sign == long_integer->sign)
    {
        calculate_sum_values(short_integer, sum);
    }
    else
    {
        twos_complement(sum);
        calculate_sum_values(short_integer, sum);
        if ((&sum->value)[long_integer->value_count] != 0)
        {
            twos_complement(sum);
            sum->sign *= -1;
        }
    }
    trim_leading_zeroes(sum);
    rewind_stack_cursor(stack_cursor, &sum->value + sum->value_count);
    return sum;
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
        memset(&a->value + a->value_count, 0,
            (b->value_count - a->value_count + 1) * sizeof(uint32_t));
        a->value_count = b->value_count + 1;
    }
    else
    {
        a->value_count += 1;
    }
    if (a->sign == b->sign)
    {
        calculate_sum_values(b, a);
    }
    else
    {
        twos_complement(a);
        calculate_sum_values(b, a);
        if ((&a->value)[a->value_count - 1] != 0)
        {
            twos_complement(a);
            a->sign *= -1;
        }
    }
    trim_leading_zeroes(a);
}

struct Integer*integer_subtract(void**output_stack_cursor, void*scratch_stack_cursor,
    struct Integer*minuend, struct Integer*subtrahend)
{
    struct Integer*negative =
        integer_copy(stack_integer_slot_new, &scratch_stack_cursor, subtrahend);
    negative->sign *= -1;
    return integer_add(output_stack_cursor, minuend, negative);
}

int8_t integer_compare(void*stack_a_cursor, void*stack_b_cursor, struct Integer*a, struct Integer*b)
{
    return integer_subtract(&stack_a_cursor, stack_b_cursor, a, b)->sign;
}

struct Integer*integer_multiply(void**output_stack_cursor, void*scratch_stack_cursor,
    struct Integer*a, struct Integer*b)
{
    struct Integer*product =
        stack_slot_new(output_stack_cursor, integer_size(a->value_count + b->value_count));
    product->value_count = 0;
    product->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)(&a->value)[i] * (&b->value)[j];
            size_t shift = i + j;
            struct Integer*integer_component =
                stack_slot_new(&scratch_stack_cursor, integer_size(shift + 2));
            integer_component->value_count = shift;
            integer_component->sign = 0;
            memset(&integer_component->value, 0, shift * sizeof(uint32_t));
            if (product_component > 0)
            {
                (&integer_component->value)[shift] = (uint32_t)product_component;
                integer_component->value_count += 1;
                integer_component->sign = 1;
                uint32_t high_bytes = (product_component & 0xffffffff00000000) >> 32;
                if (high_bytes > 0)
                {
                    (&integer_component->value)[shift + 1] = high_bytes;
                    integer_component->value_count += 1;
                }
            }
            integer_add_to_a_in_place(product, integer_component);
        }
    }
    product->sign = a->sign * b->sign;
    rewind_stack_cursor(output_stack_cursor, &product->value + product->value_count);
    return product;
}

struct Division
{
    struct Integer*quotient;
    struct Integer*remainder;
};

void integer_upshift(struct Integer*a, uint8_t shift)
{
    for (int i = a->value_count - 1; i > 0; --i)
    {
        (&a->value)[i] = (&a->value)[i] << shift | (&a->value)[i - 1] >> 32 - shift;
    }
    a->value = a->value << shift;
}

void integer_downshift(struct Integer*a, uint8_t shift)
{
    for (size_t i = 0; i < a->value_count - 1; ++i)
    {
        (&a->value)[i] = (&a->value)[i] >> shift | (&a->value)[i + 1] << 32 - shift;
    }
    (&a->value)[a->value_count - 1] = (&a->value)[a->value_count - 1] >> shift;
}

//Can't be used on pool Integers because it can change a->value_count, and pool_integer_free relies
//on an Integer's value_count being the same as when it was allocated. Use integer_half instead.
void integer_halve(struct Integer*a)
{  
    integer_downshift(a, 1);
    trim_leading_zeroes(a);
}

void calculate_division_values(void**output_stack_cursor, void*scratch_stack_cursor,
    struct Integer*positive_divisor, struct Division*division, size_t quotient_value_index,
    uint32_t quotient_digit)
{
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer*difference = integer_subtract(output_stack_cursor,
                scratch_stack_cursor, division->remainder, positive_divisor);
            if (difference->sign >= 0)
            {
                (&division->quotient->value)[quotient_value_index] |= quotient_digit;
                division->remainder = difference;
            }
            if (quotient_digit == 1)
            {
                if (quotient_value_index == 0)
                {
                    return;
                }
                quotient_digit = 0x80000000;
                quotient_value_index -= 1;
            }
            else
            {
                quotient_digit = quotient_digit >> 1;
            }
            integer_halve(positive_divisor);
        }
    }
}

int leading_digit_place(struct Integer*a)
{
    if (!a->value_count)
    {
        return 0;
    }
    uint32_t divisor_leading_digit = 0x80000000;
    uint32_t last_value = (&a->value)[a->value_count - 1];
    for (int i = 31; i >= 0; --i)
    {
        if ((last_value & divisor_leading_digit) != 0)
        {
            return i;
        }
        divisor_leading_digit = divisor_leading_digit >> 1;
    }
}

void integer_euclidean_divide(void**output_stack_cursor, void*scratch_stack_cursor,
    struct Division*out, struct Integer*dividend, struct Integer*divisor)
{
    int8_t quotient_sign = dividend->sign * divisor->sign;
    int dividend_leading_digit_place = leading_digit_place(dividend);
    int divisor_leading_digit_place = leading_digit_place(divisor);
    if (dividend->value_count > divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
            dividend_leading_digit_place >= divisor_leading_digit_place))
    {
        struct Integer*positive_divisor =
            stack_integer_slot_new(&scratch_stack_cursor, dividend->value_count);
        positive_divisor->value_count = dividend->value_count;
        positive_divisor->sign = 1;
        size_t quotient_value_index = dividend->value_count - divisor->value_count;
        memset(&positive_divisor->value, 0, quotient_value_index * sizeof(uint32_t));
        memcpy(&positive_divisor->value + quotient_value_index, &divisor->value,
            divisor->value_count * sizeof(uint32_t));
        int shift = dividend_leading_digit_place - divisor_leading_digit_place;
        uint32_t quotient_digit;
        if (shift > 0)
        {
            integer_upshift(positive_divisor, shift);
            quotient_digit = 1 << shift;
        }
        else if (shift < 0)
        {
            shift *= -1;
            integer_downshift(positive_divisor, shift);
            quotient_digit = 1 << 32 - shift;
            quotient_value_index -= 1;
        }
        else
        {
            quotient_digit = 1;
        }
        int8_t dividend_sign = dividend->sign;
        dividend->sign = 1;
        out->quotient = stack_integer_slot_new(output_stack_cursor, dividend->value_count);
        out->quotient->value_count = dividend->value_count;
        out->quotient->sign = quotient_sign;
        memset(&out->quotient->value, 0, out->quotient->value_count * sizeof(uint32_t));
        out->remainder = dividend;
        calculate_division_values(&scratch_stack_cursor, *output_stack_cursor, positive_divisor,
            out, quotient_value_index, quotient_digit);
        trim_leading_zeroes(out->quotient);
        *(uint32_t**)output_stack_cursor += out->quotient->value_count;
        integer_move(stack_integer_slot_new, output_stack_cursor, &out->remainder);
        if (out->remainder->sign != 0)
        {
            out->remainder->sign = dividend->sign;
        }
        dividend->sign = dividend_sign;
    }
    else
    {
        out->quotient = &zero;
        out->remainder = integer_copy(stack_integer_slot_new, output_stack_cursor, dividend);
    }
}

struct Integer*integer_doubled(void**stack_cursor, struct Integer*a)
{
    struct Integer*out;
    if ((&a->value)[a->value_count - 1] & 0x80000000)
    {
        out = stack_integer_slot_new(stack_cursor, a->value_count + 1);
        out->value_count = a->value_count + 1;
        (&out->value)[a->value_count] = 0;
    }
    else
    {
        out = stack_integer_slot_new(stack_cursor, a->value_count);
        out->value_count = a->value_count;
    }
    out->sign = a->sign;
    memcpy(&out->value, &a->value, a->value_count * sizeof(uint32_t));
    integer_upshift(out, 1);
    return out;
}

struct Integer*integer_half(void**stack_cursor, struct Integer*a)
{
    struct Integer*out = integer_copy(stack_integer_slot_new, stack_cursor, a);
    integer_halve(out);
    return out;
}

struct Integer*integer_exponentiate(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Integer*base,
    struct Integer*exponent)
{
    struct Integer*exponentiation =
        integer_new(stack_integer_slot_new, &scratch_stack_cursor, 1, 1);
    struct Integer*base_to_a_power_of_two = base;
    while (exponent->sign > 0)
    {
        if (exponent->value & 1)
        {
            exponentiation = integer_multiply(&scratch_stack_cursor, output_stack_cursor,
                exponentiation, base_to_a_power_of_two);
        }
        base_to_a_power_of_two = integer_multiply(&scratch_stack_cursor, output_stack_cursor,
            base_to_a_power_of_two, base_to_a_power_of_two);
        exponent = integer_half(&scratch_stack_cursor, exponent);
    }
    return integer_copy(integer_slot_allocator, memory, exponentiation);
}

void integer_string(void**output_stack_cursor, void*scratch_stack_cursor, struct Integer*a)
{
    if (a->sign == 0)
    {
        *(char*)stack_slot_new(output_stack_cursor, 1) = '0';
        return;
    }
    if (a->sign < 0)
    {
        *(char*)stack_slot_new(output_stack_cursor, 1) = '-';
    }
    char*buffer_start = stack_slot_new(output_stack_cursor, 10 * a->value_count + 1);
    char*next_char = *output_stack_cursor;
    struct Integer*quotient = a;
    struct Integer power = { 1, 1, 10 };
    while (quotient->sign != 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, *output_stack_cursor, &division, quotient,
            &power);
        quotient = division.quotient;
        next_char -= 1;
        if (division.remainder->sign != 0)
        {
            *next_char = division.remainder->value + '0';
        }
        else
        {
            *next_char = '0';
        }
    }
    size_t char_count = *(char**)output_stack_cursor - next_char;
    memcpy(buffer_start, next_char, char_count);
    rewind_stack_cursor(output_stack_cursor, buffer_start + char_count);
}

struct Integer*get_gcd(void**output_stack_cursor, void*scratch_stack_cursor, struct Integer*a,
    struct Integer*b)
{
    while (b->value_count > 0)
    {
        struct Integer*c = b;
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, *output_stack_cursor, &division, a, b);
        b = division.remainder;
        a = c;
    }
    return integer_copy(stack_integer_slot_new, output_stack_cursor, a);
}

struct ExtendedGCDInfo
{
    struct Integer*gcd;
    struct Integer*a_coefficient;
    struct Integer*b_coefficient;
    struct Integer*a_over_gcd;
    struct Integer*b_over_gcd;
};

void extended_gcd(void**output_stack_cursor, void*scratch_stack_cursor, struct ExtendedGCDInfo*out,
    struct Integer*a, struct Integer*b)
{
    out->a_coefficient = &zero;
    out->b_coefficient = integer_new(stack_integer_slot_new, &scratch_stack_cursor, 1, 1);
    out->a_over_gcd = &zero;
    out->b_over_gcd = out->b_coefficient;
    while (a->value_count)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, *output_stack_cursor, &division, b, a);
        struct Integer*m = integer_subtract(&scratch_stack_cursor, *output_stack_cursor,
            out->a_coefficient, integer_multiply(&scratch_stack_cursor, *output_stack_cursor,
                out->b_over_gcd, division.quotient));
        struct Integer*n =
            integer_subtract(&scratch_stack_cursor, *output_stack_cursor, out->b_coefficient,
                integer_multiply(&scratch_stack_cursor, *output_stack_cursor, out->a_over_gcd,
                    division.quotient));
        b = a;
        a = division.remainder;
        out->a_coefficient = out->b_over_gcd;
        out->b_coefficient = out->a_over_gcd;
        out->b_over_gcd = m;
        out->a_over_gcd = n;
    }
    out->gcd = integer_copy(stack_integer_slot_new, output_stack_cursor, b);
    integer_move(stack_integer_slot_new, output_stack_cursor, &out->a_coefficient);
    integer_move(stack_integer_slot_new, output_stack_cursor, &out->b_coefficient);
    integer_move(stack_integer_slot_new, output_stack_cursor, &out->a_over_gcd);
    integer_move(stack_integer_slot_new, output_stack_cursor, &out->b_over_gcd);
    out->b_over_gcd->sign *= -1;
}

struct Integer*primes;

struct Integer*next_prime_slot(struct Integer*prime)
{
    return (size_t)prime + integer_size(prime->value_count);
}

struct Integer*next_prime(void*stack_a_cursor, void*stack_b_cursor, struct Integer*prime)
{
    struct Integer*next_slot = next_prime_slot(prime);
    size_t boundary = next_page_boundary(prime);
    if (next_slot == boundary)
    {
        VirtualAlloc(boundary, page_size, MEM_COMMIT, PAGE_READWRITE);
    }
    if (next_slot->value_count)
    {
        return next_slot;
    }
    memcpy(next_slot, prime, integer_size(prime->value_count));
    integer_add_to_a_in_place(next_slot, primes);
    while (true)
    {
        struct Integer*potential_factor = primes;
        while (true)
        {
            struct Division division;
            integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division, next_slot,
                potential_factor);
            if (division.remainder->value_count == 0)
            {
                break;
            }
            potential_factor = next_prime_slot(potential_factor);
            if (potential_factor == next_slot)
            {
                return next_slot;
            }
        }
        integer_add_to_a_in_place(next_slot, primes);
    }
}

struct Rational
{
    struct Integer*numerator;
    struct Integer*denominator;
};

struct RationalInterval
{
    struct Rational min;
    struct Rational max;
};

bool rational_new(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*out,
    struct Integer*numerator, struct Integer*denominator)
{
    if (denominator->sign == 0)
    {
        printf("Tried to divide by 0.");
        return false;
    }
    if (numerator->sign == 0)
    {
        out->numerator = &zero;
        out->denominator = integer_new(integer_slot_allocator, memory, 1, 1);
        return true;
    }
    struct ExtendedGCDInfo gcd_info;
    extended_gcd(&scratch_stack_cursor, output_stack_cursor, &gcd_info, numerator, denominator);
    out->numerator = gcd_info.a_over_gcd;
    out->denominator = gcd_info.b_over_gcd;
    if (out->denominator->sign < 0)
    {
        out->numerator->sign *= -1;
        out->denominator->sign *= -1;
    }
    integer_move(integer_slot_allocator, memory, &out->numerator);
    integer_move(integer_slot_allocator, memory, &out->denominator);
    return true;
}

void rational_free(struct IntegerPool*integer_pool, struct Rational*a)
{
    pool_integer_free(integer_pool, a->numerator);
    pool_integer_free(integer_pool, a->denominator);
}

void rational_move(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    struct Rational*a)
{
    integer_move(integer_slot_allocator, memory, &a->numerator);
    integer_move(integer_slot_allocator, memory, &a->denominator);
}

void rational_move_from_pool(struct IntegerPool*pool, void**stack_cursor, struct Rational*a)
{
    integer_move_from_pool(pool, stack_cursor, &a->numerator);
    integer_move_from_pool(pool, stack_cursor, &a->denominator);
}

bool rational_out_add(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, uint32_t*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Rational*b)
{
    return rational_new(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        out,
        integer_add(&scratch_stack_cursor,
            integer_multiply(&scratch_stack_cursor,
                output_stack_cursor, a->numerator, b->denominator),
            integer_multiply(&scratch_stack_cursor, output_stack_cursor, b->numerator,
                a->denominator)),
        integer_multiply(&scratch_stack_cursor, output_stack_cursor, a->denominator,
            b->denominator));
}

struct Rational rational_add(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*a, struct Rational*b)
{
    struct Rational sum;
    rational_out_add(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        &sum, a, b);
    return sum;
}

bool rational_out_subtract(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, uint32_t*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Rational*b)
{
    return rational_new(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        out,
        integer_subtract(&scratch_stack_cursor, output_stack_cursor,
            integer_multiply(&scratch_stack_cursor,
                output_stack_cursor, a->numerator, b->denominator),
            integer_multiply(&scratch_stack_cursor, output_stack_cursor, b->numerator,
                a->denominator)),
        integer_multiply(&scratch_stack_cursor, output_stack_cursor, a->denominator,
            b->denominator));
}

struct Rational rational_subtract(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*stack_a_cursor, void*stack_b_cursor, struct Rational*a, struct Rational*b)
{
    struct Rational difference;
    rational_out_subtract(integer_slot_allocator, memory, stack_a_cursor, stack_b_cursor,
        &difference, a, b);
    return difference;
}

bool rational_out_multiply(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Rational*b)
{
    return rational_new(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        out,
        integer_multiply(&scratch_stack_cursor, output_stack_cursor, a->numerator, b->numerator),
        integer_multiply(&scratch_stack_cursor, output_stack_cursor, a->denominator,
            b->denominator));
}

struct Rational rational_multiply(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*stack_a_cursor, void*stack_b_cursor, struct Rational*a, struct Rational*b)
{
    struct Rational product;
    rational_out_multiply(integer_slot_allocator, memory, stack_a_cursor, stack_b_cursor, &product,
        a, b);
    return product;
}

bool rational_integer_out_multiply(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Integer*b)
{
    return rational_new(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        out, integer_multiply(&scratch_stack_cursor, output_stack_cursor, a->numerator, b),
        a->denominator);
}

bool rational_integer_out_divide(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*out,
    struct Rational*dividend, struct Integer*divisor)
{
    return rational_new(integer_slot_allocator, memory, output_stack_cursor, scratch_stack_cursor,
        out, dividend->numerator, integer_multiply(&scratch_stack_cursor, output_stack_cursor,
            dividend->denominator, divisor));
}

int8_t rational_compare(void*stack_a_cursor, void*stack_b_cursor, struct Rational*a,
    struct Rational*b)
{
    return integer_compare(stack_a_cursor, stack_b_cursor,
        integer_multiply(&stack_a_cursor, stack_b_cursor, a->numerator, b->denominator),
        integer_multiply(&stack_a_cursor, stack_b_cursor, a->denominator, b->numerator));
}

struct Rational pi_estimate_min;
struct Rational pi_interval_size;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

void pi_refine_interval(void*stack_a_cursor, void*stack_b_cursor, struct Rational*interval_size)
{
    int8_t interval_size_comparison = rational_compare(stack_a_cursor, stack_b_cursor,
        &pi_interval_size, interval_size);
    if (interval_size_comparison > 0)
    {
        struct Rational sixteenth = { &INT(1, 1), &INT(16, 1) };
        rational_move_from_pool(permanent_integer_pool, &stack_a_cursor, &pi_estimate_min);
        rational_move_from_pool(permanent_integer_pool, &stack_a_cursor,
            &pi_interval_size);
        integer_move_from_pool(permanent_integer_pool, &stack_a_cursor, &pi_sixteen_to_the_k);
        integer_move_from_pool(permanent_integer_pool, &stack_a_cursor, &pi_eight_k);
        while (interval_size_comparison > 0)
        {
            struct Rational term_a =
                { &INT(4, 1), integer_add(&stack_a_cursor, pi_eight_k, &INT(1, 1)) };
            struct Rational term_b;
            rational_new(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor, stack_b_cursor,
                &term_b, &INT(2, 1), integer_add(&stack_a_cursor, pi_eight_k, &INT(4, 1)));
            struct Rational sum;
            rational_out_subtract(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &sum, &term_a, &term_b);
            struct Rational term_c =
                { &INT(1, 1), integer_add(&stack_a_cursor, pi_eight_k, &INT(5, 1)) };
            sum = rational_subtract(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &sum, &term_c);
            struct Rational term_d =
                { &INT(1, 1), integer_add(&stack_a_cursor, pi_eight_k, &INT(6, 1)) };
            sum = rational_subtract(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &sum, &term_d);
            struct Rational term;
            rational_integer_out_divide(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &term, &sum, pi_sixteen_to_the_k);
            pi_estimate_min = rational_add(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &pi_estimate_min, &term);
            pi_interval_size = rational_multiply(stack_integer_slot_new, &stack_a_cursor,
                stack_a_cursor, stack_b_cursor, &pi_interval_size, &sixteenth);
            pi_eight_k = integer_add(&stack_a_cursor, pi_eight_k, &INT(8, 1));
            pi_sixteen_to_the_k =
                integer_multiply(&stack_a_cursor, stack_b_cursor, pi_sixteen_to_the_k, &INT(16, 1));
            interval_size_comparison =
                rational_compare(stack_a_cursor, stack_b_cursor, &pi_interval_size, interval_size);
        }
        rational_move(pool_integer_slot_new, permanent_integer_pool, &pi_estimate_min);
        rational_move(pool_integer_slot_new, permanent_integer_pool, &pi_interval_size);
        integer_move(pool_integer_slot_new, permanent_integer_pool, &pi_sixteen_to_the_k);
        integer_move(pool_integer_slot_new, permanent_integer_pool, &pi_eight_k);
    }
}

struct Float
{
    struct Integer*significand;
    struct Integer*exponent;
};

void float_new(struct Integer*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Float*out,
    struct Integer*significand, struct Integer*exponent)
{
    while (!(significand->value & 1))
    {
        significand = integer_half(&scratch_stack_cursor, significand);
        exponent = integer_add(&scratch_stack_cursor, exponent, &INT(1, 1));
    }
    integer_move(integer_slot_allocator, memory, &significand);
    integer_move(integer_slot_allocator, memory, &exponent);
    out->significand = significand;
    out->exponent = exponent;
}

void float_move(struct IntegerPool*integer_pool, struct Float*a)
{
    integer_move(pool_integer_slot_new, integer_pool, &a->significand);
    integer_move(pool_integer_slot_new, integer_pool, &a->exponent);
}

void float_move_from_pool(struct IntegerPool*integer_pool, void**stack_cursor, struct Float*a)
{
    integer_move_from_pool(integer_pool, stack_cursor, &a->significand);
    integer_move_from_pool(integer_pool, stack_cursor, &a->exponent);
}

void float_out_add(void**output_stack_cursor, void*scratch_stack_cursor, struct Float*out,
    struct Float*a, struct Float*b)
{
    struct Integer*exponent_difference =
        integer_subtract(&scratch_stack_cursor, *output_stack_cursor, a->exponent, b->exponent);
    if (exponent_difference->sign > 0)
    {
        struct Integer two = { 1, 1, 2 };
        out->significand = integer_add(output_stack_cursor, b->significand,
            integer_multiply(&scratch_stack_cursor, *output_stack_cursor, a->significand,
                integer_exponentiate(stack_integer_slot_new, &scratch_stack_cursor,
                    scratch_stack_cursor, *output_stack_cursor, &two, exponent_difference)));
        out->exponent = b->exponent;
    }
    else if (exponent_difference->sign < 0)
    {
        float_out_add(output_stack_cursor, scratch_stack_cursor, out, a, b);
    }
    else
    {
        float_new(stack_integer_slot_new, output_stack_cursor, *output_stack_cursor,
            scratch_stack_cursor, out,
            integer_add(&scratch_stack_cursor, a->significand, b->significand), a->exponent);
    }
}

void float_out_subtract(void**output_stack_cursor, void*scratch_stack_cursor, struct Float*out,
    struct Float*minuend, struct Float*subtrahend)
{
    struct Float negative = { integer_copy(stack_integer_slot_new, &scratch_stack_cursor,
        subtrahend->significand), subtrahend->exponent };
    negative.significand->sign *= -1;
    float_out_add(output_stack_cursor, scratch_stack_cursor, out, minuend, &subtrahend);
}

void float_to_rational(void**output_stack_cursor, void*scratch_stack_cursor, struct Rational*out,
    struct Float*a)
{
    struct Integer two = { 1, 1, 2 };
    if (a->exponent->value_count < 0)
    {
        out->numerator =
            integer_new(stack_integer_slot_new, output_stack_cursor, 1, a->significand->sign);
        struct Integer*positive_exponent =
            integer_copy(stack_integer_slot_new, &scratch_stack_cursor, a->exponent);
        positive_exponent->sign = 1;
        out->denominator = integer_multiply(output_stack_cursor, scratch_stack_cursor,
            integer_exponentiate(stack_integer_slot_new, &scratch_stack_cursor,
                scratch_stack_cursor, *output_stack_cursor, &two, positive_exponent),
            a->significand);
        out->denominator->sign = 1;
    }
    else
    {
        out->numerator = integer_multiply(output_stack_cursor, scratch_stack_cursor,
            integer_exponentiate(stack_integer_slot_new, &scratch_stack_cursor,
                scratch_stack_cursor, *output_stack_cursor, &two, a->exponent), a->significand);
        out->denominator = integer_new(stack_integer_slot_new, output_stack_cursor, 1, 1);
    }
}

void float_free(struct IntegerPool*integer_pool, struct Float*a)
{
    pool_integer_free(integer_pool, a->significand);
    pool_integer_free(integer_pool, a->exponent);
}

struct FloatInterval
{
    struct Float min;
    struct Float max;
};

void float_interval_move(struct IntegerPool*integer_pool, struct FloatInterval*a)
{
    if (a->min.significand)
    {
        float_move(integer_pool, &a->min);
        float_move(integer_pool, &a->max);
    }
}

struct Number
{
    union
    {
        struct
        {//When operation == 'r'.
            struct Rational value;
            struct Integer*magnitude_estimate_denominator;
            struct Integer*magnitude_estimate_remainder;
        };
        struct
        {//When operation != 'r'.
            struct Number*left;
            struct Number*right;
            struct FloatInterval imaginary_part_estimate;
        };
    };
    struct FloatInterval real_part_estimate;
    struct FloatInterval argument_estimate;
    struct FloatInterval magnitude_estimate;
    struct Number*previous;
    struct Number*next;
    char operation;
};

struct Number*number_slot_new(struct Number**pool_cursor, struct Number**free_list)
{
    struct Number*out = *free_list;
    if (out)
    {
        *free_list = out->next;
    }
    else
    {
        out = stack_slot_new(pool_cursor, sizeof(struct Number));
    }
    memset(out, 0, sizeof(struct Number));
    return out;
}

struct Number*number_copy(struct Number**pool_cursor, struct Number**free_list,
    struct IntegerPool*integer_pool, struct Number*a)
{
    struct Number*out = number_slot_new(pool_cursor, free_list);
    memcpy(out, a, sizeof(struct Number));
    if (a->operation == 'r')
    {
        rational_move(pool_integer_slot_new, integer_pool, &out->value);
    }
    else
    {
        out->left = number_copy(pool_cursor, free_list, integer_pool, a->left);
        out->right = number_copy(pool_cursor, free_list, integer_pool, a->right);
        float_interval_move(integer_pool, &a->imaginary_part_estimate);
    }
    float_interval_move(integer_pool, &a->real_part_estimate);
    float_interval_move(integer_pool, &a->argument_estimate);
    float_interval_move(integer_pool, &a->magnitude_estimate);
    return out;
}

void number_root_free(struct Number**free_list, struct IntegerPool*integer_pool, struct Number*a)
{
    if (a->operation == 'r')
    {
        rational_free(integer_pool, &a->value);
    }
    a->next = *free_list;
    *free_list = a;
    if (a->argument_estimate.min.significand)
    {
        float_free(integer_pool, &a->argument_estimate.min);
        float_free(integer_pool, &a->argument_estimate.max);
    }
}

void number_free(struct Number**free_list, struct IntegerPool*integer_pool, struct Number*a)
{
    if (a->operation == 'r')
    {
        rational_free(integer_pool, &a->value);
    }
    else
    {
        number_free(free_list, integer_pool, a->left);
        number_free(free_list, integer_pool, a->right);
    }
    a->next = *free_list;
    *free_list = a;
    if (a->argument_estimate.min.significand)
    {
        float_free(integer_pool, &a->argument_estimate.min);
        float_free(integer_pool, &a->argument_estimate.max);
    }
}

bool get_input(struct Number**number_pool_cursor, struct Number**number_free_list,
    struct IntegerPool*integer_pool, void*stack_a_cursor, void*stack_b_cursor)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    struct Number*previous = 0;
    while (true)
    {
        struct Number*number = number_slot_new(number_pool_cursor, number_free_list);
        number->previous = previous;
        number->next = *number_pool_cursor;
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            number->value.numerator =
                integer_new(stack_integer_slot_new, &stack_a_cursor, next_char - '0', 1);
            struct Integer ten = { 1, 1, 10 };
            next_char = getchar();
            while (isdigit(next_char))
            {
                number->value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
                    number->value.numerator, &ten);
                struct Integer digit = { 1, 1, next_char - '0' };
                number->value.numerator =
                    integer_add(&stack_a_cursor, number->value.numerator, &digit);
                next_char = getchar();
            }
            integer_move(pool_integer_slot_new, integer_pool, &number->value.numerator);
        }
        else
        {
            switch (next_char)
            {
            case '+':
                break;
            case '-':
                break;
            case '*':
                break;
            case '/':
                break;
            case '^':
                break;
            case '(':
                break;
            case ')':
                break;
            default:
                printf("\"%c\"%s", next_char, " is an invalid character.");
                while (getchar() != '\n')
                {}
                return false;
            }
            number->operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            number->next = 0;
            return true;
        }
        previous = number;
    }
}

bool is_numeric(struct Number*a)
{
    if (!a)
    {
        return false;
    }
    return a->operation == 'r' || a->left;
}

bool is_unparsed_operation(struct Number*a, char operation)
{
    return a->operation == operation && !a->left;
}

bool parse_binary_operation(struct Number*a, char operation)
{
    if (!is_unparsed_operation(a, operation))
    {
        return true;
    }
    if (!a->previous || !is_numeric(a->previous))
    {
        printf("%c missing left operand.", operation);
        return false;
    }
    if (a->next && is_numeric(a->next))
    {
        a->left = a->previous;
        a->right = a->next;
        a->previous = a->previous->previous;
        if (a->previous)
        {
            a->previous->next = a;
        }
        a->next = a->next->next;
        if (a->next)
        {
            a->next->previous = a;
        }
    }
    else
    {
        printf("%c missing right operand.", operation);
        return false;
    }
    return true;
}

void rewind_to_first(struct Number**a)
{
    while ((*a)->previous)
    {
        *a = (*a)->previous;
    }
}

bool parse_binary_operation_pair(struct Number**number, char operation_a, char operation_b)
{
    rewind_to_first(number);
    while (true)
    {
        if (!parse_binary_operation(*number, operation_a))
        {
            return false;
        }
        else if (!parse_binary_operation(*number, operation_b))
        {
            return false;
        }
        if (!(*number)->next)
        {
            return true;
        }
        *number = (*number)->next;
    }
}

struct Number*parse_input(struct Number**number_pool_cursor, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Number*input)
{
    if (!input)
    {
        printf("Empty expression.");
        return 0;
    }
    while (true)
    {
        if (input->operation == ')')
        {
            printf("Unmatched ).");
            return 0;
        }
        if (input->operation == '(')
        {
            int unmatched_paren_count = 1;
            struct Number*nested_number = input;
            while (unmatched_paren_count > 0)
            {
                nested_number = nested_number->next;
                if (!nested_number)
                {
                    printf("Unmatched (.");
                    return 0;
                }
                if (nested_number->operation == '(')
                {
                    unmatched_paren_count += 1;
                }
                else if (nested_number->operation == ')')
                {
                    unmatched_paren_count -= 1;
                }
            }
            struct Number*previous = input->previous;
            struct Number*next = nested_number->next;
            input->next->previous = 0;
            nested_number->previous->next = 0;
            struct Number*parsed_nested_expression =
                parse_input(number_pool_cursor, number_free_list, integer_pool, input->next);
            if (parsed_nested_expression)
            {
                parsed_nested_expression->previous = previous;
                if (previous)
                {
                    previous->next = parsed_nested_expression;
                }
                parsed_nested_expression->next = next;
                if (next)
                {
                    next->previous = parsed_nested_expression;
                }
                number_root_free(number_free_list, integer_pool, input);
                number_root_free(number_free_list, integer_pool, nested_number);
                input = parsed_nested_expression;
            }
            else
            {
                return 0;
            }
        }
        if (!input->next)
        {
            break;
        }
        input = input->next;
    }
    rewind_to_first(&input);
    while (true)
    {
        struct Number*next_number = input->next;
        if (is_unparsed_operation(input, '+') && !is_numeric(input->previous))
        {
            if (!next_number || (!is_numeric(next_number) &&
                !is_unparsed_operation(next_number, '+') &&
                !is_unparsed_operation(next_number, '-')))
            {
                printf("+ missing right operand.");
                return 0;
            }
            next_number->previous = input->previous;
            if (input->previous)
            {
                input->previous->next = next_number;
            }
            number_root_free(number_free_list, integer_pool, input);
        }
        if (!next_number)
        {
            break;
        }
        input = next_number;
    }
    rewind_to_first(&input);
    while (true)
    {
        if (!parse_binary_operation(input, '^'))
        {
            return 0;
        }
        if (!input->next)
        {
            break;
        }
        input = input->next;
    }
    rewind_to_first(&input);
    while (true)
    {
        if (is_unparsed_operation(input, '-') && !is_numeric(input->previous))
        {
            if (!input->next)
            {
                printf("- missing right operand.");
                return 0;
            }
            if (is_numeric(input->next))
            {
                input->operation = 'r';
                input->value.numerator = integer_new(pool_integer_slot_new, integer_pool, 1, -1);
                input->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
                struct Number*times = number_slot_new(number_pool_cursor, number_free_list);
                times->operation = '*';
                times->left = input;
                times->right = input->next;
                times->next = input->next->next;
                if (input->previous)
                {
                    input->previous->next = times;
                }
                input = times;
                if (input->next)
                {
                    input->next->previous = input;
                }
                else
                {
                    break;
                }
                input = input->next;
            }
            else if (is_unparsed_operation(input->next, '-'))
            {
                struct Number*new_next = input->next->next;
                if (!new_next)
                {
                    printf("- missing right operand.");
                    return 0;
                }
                new_next->previous = input->previous;
                if (input->previous)
                {
                    input->previous->next = new_next;
                }
                number_root_free(number_free_list, integer_pool, input->next);
                number_root_free(number_free_list, integer_pool, input);
                input = new_next;
            }
            else
            {
                printf("- missing left operand.");
                return 0;
            }
        }
        else
        {
            if (!input->next)
            {
                break;
            }
            input = input->next;
        }
    }
    if (!parse_binary_operation_pair(&input, '*', '/'))
    {
        return 0;
    }
    if (!parse_binary_operation_pair(&input, '+', '-'))
    {
        return 0;
    }
    return input;
}

bool estimate_needs_refinement(void*stack_a_cursor, void*stack_b_cursor,
    struct FloatInterval*estimate, struct Rational*interval_size)
{
    if (estimate->min.significand)
    {
        struct Float current_interval_size;
        float_out_subtract(&stack_a_cursor, stack_b_cursor, &current_interval_size, &estimate->max,
            &estimate->min);
        struct Rational rational_current_interval_size;
        float_to_rational(&stack_a_cursor, stack_b_cursor, &rational_current_interval_size,
            &current_interval_size);
        if (rational_compare(stack_a_cursor, stack_b_cursor, interval_size,
            &rational_current_interval_size) >= 0)
        {
            return false;
        }
    }
    return true;
}

void rational_continue_float_estimate(void**output_stack_cursor, void*scratch_stack_cursor,
    struct IntegerPool*integer_pool, struct FloatInterval*out, struct Integer**estimate_denominator,
    struct Integer**estimate_remainder, struct Rational*a, struct Rational*interval_size)
{
    while (integer_compare(scratch_stack_cursor, *output_stack_cursor, interval_size->denominator,
        integer_multiply(&scratch_stack_cursor, *output_stack_cursor, *estimate_denominator,
            interval_size->numerator)) > 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, *output_stack_cursor, &division,
            integer_doubled(&scratch_stack_cursor, *estimate_remainder), a->denominator);
        out->min.significand = integer_add(&scratch_stack_cursor, division.quotient,
            integer_doubled(&scratch_stack_cursor, out->min.significand));
        *estimate_denominator = integer_doubled(&scratch_stack_cursor, *estimate_denominator);
        *estimate_remainder = division.remainder;
        out->min.exponent = integer_add(&scratch_stack_cursor, out->min.exponent, &INT(1, -1));
        if ((*estimate_remainder)->value_count == 0)
        {
            out->max = out->min;
            return;
        }
    }
    float_new(pool_integer_slot_new, integer_pool, scratch_stack_cursor, *output_stack_cursor,
        &out->min, out->min.significand, out->min.exponent);
    float_new(pool_integer_slot_new, integer_pool, scratch_stack_cursor, *output_stack_cursor,
        &out->max, integer_add(&scratch_stack_cursor, out->min.significand, &INT(1, 1)),
        out->min.exponent);
    integer_move(stack_integer_slot_new, output_stack_cursor, estimate_denominator);
    integer_move(stack_integer_slot_new, output_stack_cursor, estimate_remainder);
}

void rational_float_estimate(void*stack_a_cursor, void*stack_b_cursor,
    struct IntegerPool*integer_pool, struct FloatInterval*out, struct Rational*a,
    struct Rational*interval_size)
{
    struct Integer*numerator_magnitude =
        integer_copy(stack_integer_slot_new, &stack_a_cursor, a->numerator);
    numerator_magnitude->sign *= numerator_magnitude->sign;
    struct Division division;
    integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division, numerator_magnitude,
        a->denominator);
    struct Integer*estimate_remainder = division.remainder;
    struct Integer*estimate_denominator = 
        integer_new(stack_integer_slot_new, &stack_a_cursor, 1, 1);
    out->min.significand = division.quotient;
    out->min.exponent = &zero;
    rational_continue_float_estimate(&stack_a_cursor, stack_b_cursor, integer_pool, out,
        &estimate_denominator, &estimate_remainder, a, interval_size);
}

struct FloatInterval*float_magnitude_estimate(void*stack_a_cursor, void*stack_b_cursor,
    struct IntegerPool*integer_pool, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(stack_a_cursor, stack_b_cursor, &a->magnitude_estimate,
        interval_size))
    {
        return &a->magnitude_estimate;
    }
    if (a->operation == 'r')
    {
        if (!a->magnitude_estimate.min.significand)
        {
            struct Integer*numerator_magnitude =
                integer_copy(stack_integer_slot_new, &stack_a_cursor, a->value.numerator);
            numerator_magnitude->sign = 1;
            struct Division division;
            integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                numerator_magnitude, a->value.denominator);
            a->magnitude_estimate_denominator =
                integer_new(stack_integer_slot_new, &stack_a_cursor, 1, 1);
            a->magnitude_estimate.min.significand = division.quotient;
            a->magnitude_estimate.min.exponent = &zero;
            a->magnitude_estimate_remainder = division.remainder;
        }
        else
        {
            float_move_from_pool(integer_pool, &stack_a_cursor, &a->magnitude_estimate.min);
            float_move_from_pool(integer_pool, &stack_a_cursor, &a->magnitude_estimate.max);
            integer_move_from_pool(integer_pool, &stack_a_cursor,
                &a->magnitude_estimate_denominator);
            integer_move_from_pool(integer_pool, &stack_a_cursor, &a->magnitude_estimate_remainder);
        }
        rational_continue_float_estimate(stack_a_cursor, stack_b_cursor, integer_pool,
            &a->magnitude_estimate, a, a->magnitude_estimate_denominator,
            a->magnitude_estimate_remainder, interval_size);
        integer_move(pool_integer_slot_new, integer_pool, &a->magnitude_estimate_denominator);
        integer_move(pool_integer_slot_new, integer_pool, &a->magnitude_estimate_remainder);
        return &a->magnitude_estimate;
    }
    float_free(integer_pool, &a->magnitude_estimate.min);
    float_free(integer_pool, &a->magnitude_estimate.max);
    ABORT("float_magnitude_estimate case not yet implemented.");
}

void rational_argument_estimate(void**output_stack_cursor, void*scratch_stack_cursor,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->sign < 0)
        {
            pi_refine_interval(*output_stack_cursor, scratch_stack_cursor, interval_size);
            out->min = pi_estimate_min;
            rational_out_add(stack_integer_slot_new, output_stack_cursor, *output_stack_cursor,
                scratch_stack_cursor, &out->max, &pi_estimate_min, &pi_interval_size);
        }
        else
        {
            out->min.numerator = &zero;
            out->min.denominator = integer_new(stack_integer_slot_new, output_stack_cursor, 1, 1);
            out->max.numerator = &zero;
            out->max.denominator = out->min.denominator;
        }
        return;
    case '^':
    {
        struct Rational radicand_argument_interval_size;
        rational_integer_out_multiply(stack_integer_slot_new, &scratch_stack_cursor,
            scratch_stack_cursor, *output_stack_cursor, &radicand_argument_interval_size,
            interval_size, a->right->value.denominator);
        struct RationalInterval radicand_rational_argument_estimate;
        rational_argument_estimate(&scratch_stack_cursor, *output_stack_cursor, integer_pool,
            &radicand_rational_argument_estimate, a->left, &radicand_argument_interval_size);
        rational_integer_out_divide(stack_integer_slot_new, output_stack_cursor,
            *output_stack_cursor, scratch_stack_cursor, &out->min,
            &radicand_rational_argument_estimate.min, a->right->value.denominator);
        rational_integer_out_divide(stack_integer_slot_new, output_stack_cursor,
            *output_stack_cursor, scratch_stack_cursor, &out->max,
            &radicand_rational_argument_estimate.max, a->right->value.denominator);
        return;
    }
    default:
        ABORT("rational_argument_estimate case not yet implemented.");
    }
}

struct FloatInterval*float_argument_estimate(void*stack_a_cursor, void*stack_b_cursor,
    struct IntegerPool*integer_pool, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(stack_a_cursor, stack_b_cursor, &a->argument_estimate,
        interval_size))
    {
        return &a->magnitude_estimate;
    }
    if (a->operation == '*')
    {
        struct Integer two = { 1, 1, 2 };
        struct Rational factor_estimate_interval_size;
        rational_integer_out_divide(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
            stack_b_cursor, &factor_estimate_interval_size, interval_size, &two);
        struct FloatInterval*left_factor_argument_estimate = float_argument_estimate(stack_a_cursor,
            stack_a_cursor, integer_pool, a->left, &factor_estimate_interval_size);
        struct FloatInterval*right_factor_argument_estimate =
            float_argument_estimate(stack_a_cursor, stack_b_cursor, integer_pool, a->right,
                &factor_estimate_interval_size);
        float_out_add(&stack_a_cursor, stack_b_cursor, &a->argument_estimate.min,
            &left_factor_argument_estimate->min, &right_factor_argument_estimate->min);
        float_out_add(&stack_a_cursor, stack_b_cursor, &a->argument_estimate.max,
            &left_factor_argument_estimate->max, &right_factor_argument_estimate->max);
        return &a->magnitude_estimate;
    }
    struct Rational rational_estimate_interval_size = { &INT(1, 1), &INT(1, 1) };
    while (rational_compare(stack_a_cursor, stack_b_cursor, &rational_estimate_interval_size,
        interval_size) > 0)
    {
        rational_estimate_interval_size.denominator =
            integer_doubled(&stack_a_cursor, rational_estimate_interval_size.denominator);
    }
    rational_estimate_interval_size.denominator =
        integer_doubled(&stack_a_cursor, rational_estimate_interval_size.denominator);
    struct RationalInterval rational_estimate;
    rational_argument_estimate(&stack_a_cursor, stack_b_cursor, integer_pool, &rational_estimate, a,
        &rational_estimate_interval_size);
    rational_float_estimate(stack_a_cursor, stack_b_cursor, integer_pool, &a->argument_estimate,
        &rational_estimate.min, &rational_estimate_interval_size);
    struct FloatInterval max_estimate;
    rational_float_estimate(stack_a_cursor, stack_b_cursor, integer_pool, &max_estimate,
        &rational_estimate.max, &rational_estimate_interval_size);
    a->argument_estimate.max = max_estimate.max;
    return &a->argument_estimate;
}

struct Number*number_add(struct Number**number_pool_cursor, struct Number**number_free_list,
    struct IntegerPool*integer_pool, void*stack_a_cursor, void*stack_b_cursor, struct Number*a,
    struct Number*b)
{
    if (a->operation == 'r')
    {
        if (a->value.numerator->value_count == 0)
        {
            number_free(number_free_list, integer_pool, a);
            return b;
        }
        if (b->operation == 'r')
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = 'r';
            if (!rational_out_add(pool_integer_slot_new, integer_pool, stack_a_cursor,
                stack_b_cursor, &out->value, &a->value, &b->value))
            {
                return 0;
            }
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            return out;
        }
    }
    else if (a->operation == '^')
    {
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    else if (a->operation == '*')
    {
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    else if (a->operation == '+')
    {
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '+':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    return number_add(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, b, a);
}

struct Number*number_multiply(struct Number**number_pool_cursor, struct Number**number_free_list,
    struct IntegerPool*integer_pool, void*stack_a_cursor, void*stack_b_cursor, struct Number*a,
    struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            return out;
        }
        if (integer_equals_one(a->value.numerator) && integer_equals_one(a->value.denominator))
        {
            number_free(number_free_list, integer_pool, a);
            return b;
        }
        if (b->operation == 'r')
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = 'r';
            if (!rational_out_multiply(pool_integer_slot_new, integer_pool, stack_a_cursor,
                stack_b_cursor, &out->value, &a->value, &b->value))
            {
                return 0;
            }
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            return out;
        }
        break;
    case '^':
        switch (b->operation)
        {
        case 'r':
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '*';
            out->left = b;
            out->right = a;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return out;
        }
        default:
            break;
        }
    case '*':
        switch (b->operation)
        {
        case 'r':
        {
            if (a->left->operation == 'r')
            {
                struct NumberSlot*coefficient = number_multiply(number_pool_cursor,
                    number_free_list, integer_pool, stack_a_cursor, stack_b_cursor, b, a->left);
                if (!coefficient)
                {
                    return 0;
                }
                number_free(number_free_list, integer_pool, b);
                struct NumberSlot*out = number_multiply(number_pool_cursor, number_free_list,
                    integer_pool, stack_a_cursor, stack_b_cursor, coefficient, a->right);
                number_root_free(number_free_list, integer_pool, a);
                return out;
            }
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return true;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            struct Number*out = number_multiply(number_pool_cursor, number_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, a->right, b);
            if (!out)
            {
                return 0;
            }
            out = number_multiply(number_pool_cursor, number_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, a->left, out);
            number_root_free(number_free_list, integer_pool, a);
            return out;
        }
        default:
            break;
        }
    case '+':
    {
        struct Number*b_copy = number_copy(number_pool_cursor, number_free_list, integer_pool, b);
        if (b->operation == 'r')
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            out->operation = '+';
            out->left = number_slot_new(number_pool_cursor, number_free_list);
            out->left = number_multiply(number_pool_cursor, number_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, a->left, b);
            if (!out->left)
            {
                return 0;
            }
            out->right = number_multiply(number_pool_cursor, number_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, a->right, b_copy);
            if (!out->right)
            {
                return 0;
            }
            number_root_free(number_free_list, integer_pool, a);
            return out;
        }
        struct Number*left = number_multiply(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a->left, b);
        if (!left)
        {
            return 0;
        }
        struct Number*right = number_multiply(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a->right, b_copy);
        if (!right)
        {
            return 0;
        }
        number_root_free(number_free_list, integer_pool, a);
        return number_add(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor, left, right);
    }
    }
    return number_multiply(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, b, a);
}

bool number_reciprocal(struct Number*a)
{
    if (a->operation == 'r')
    {
        if (!a->value.numerator->value_count)
        {
            printf("Tried to divide by 0.");
            return false;
        }
        struct Integer*old_numerator = a->value.numerator;
        a->value.numerator = a->value.denominator;
        a->value.denominator = old_numerator;
        return true;
    }
    else
    {
        ABORT("number_reciprocal called on Number type for which it is not yet implemented.");
    }
}

struct Factor
{
    struct Integer*value;
    struct Integer*multiplicity;
};

struct Number*number_exponentiate(struct Number**number_pool_cursor,
    struct Number**number_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*base, struct Number*exponent)
{
    if (exponent->operation != 'r')
    {
        printf("The input expression contains an exponentiation whose exponent is not both real "
            "and rational; this program doesn't handle transcendental numbers.");
        return 0;
    }
    if (exponent->value.numerator->sign < 0)
    {
        if (!number_reciprocal(base))
        {
            return 0;
        }
        exponent->value.numerator->sign = 1;
        return number_exponentiate(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, base, exponent);
    }
    if (base->operation == 'r')
    {
        if (base->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
            number_free(number_free_list, integer_pool, base);
            number_free(number_free_list, integer_pool, exponent);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            return true;
        }
        if (integer_equals_one(base->value.denominator))
        {
            if (integer_equals_one(exponent->value.denominator))
            {
                struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
                out->operation = 'r';
                out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
                out->value.numerator = integer_exponentiate(pool_integer_slot_new, integer_pool,
                    stack_a_cursor, stack_b_cursor, base->value.numerator,
                    exponent->value.numerator);
                number_free(number_free_list, integer_pool, base);
                number_free(number_free_list, integer_pool, exponent);
                return out;
            }
            struct Integer*radicand = integer_exponentiate(stack_integer_slot_new, &stack_a_cursor,
                stack_a_cursor, stack_b_cursor, base->value.numerator, exponent->value.numerator);
            number_free(number_free_list, integer_pool, base);
            int8_t radicand_sign = radicand->sign;
            radicand->sign = 1;
            size_t factor_count = 0;
            struct Factor*factors = stack_b_cursor;
            struct Integer*prime = primes;
            while (integer_compare(stack_a_cursor, stack_b_cursor, prime, radicand) <= 0)
            {
                struct Division division;
                integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division, radicand,
                    prime);
                if (division.remainder->value_count == 0)
                {
                    factor_count += 1;
                    struct Factor*factor = stack_slot_new(&stack_b_cursor, sizeof(struct Factor));
                    factor->value = prime;
                    factor->multiplicity = &zero;
                    while (division.remainder->value_count == 0)
                    {
                        factor->multiplicity =
                            integer_add(&stack_a_cursor, factor->multiplicity, &INT(1, 1));
                        radicand = division.quotient;
                        integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                            radicand, factor->value);
                    }
                }
                prime = next_prime(stack_a_cursor, stack_b_cursor, prime);
            }
            struct Integer*coefficient = integer_new(stack_integer_slot_new, &stack_a_cursor, 1, 1);
            struct Integer*multiplicity_gcd = exponent->value.denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Division multiplicity_reduction;
                integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->value.denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(&stack_a_cursor, stack_b_cursor, coefficient,
                    integer_exponentiate(stack_integer_slot_new, &stack_a_cursor, stack_a_cursor,
                        stack_b_cursor, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = get_gcd(&stack_a_cursor, stack_b_cursor, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            if (radicand_sign > 0)
            {
                struct Division reduced_degree;
                integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &reduced_degree,
                    exponent->value.denominator, multiplicity_gcd);
                if (integer_equals_one(reduced_degree.quotient))
                {
                    number_free(number_free_list, integer_pool, exponent);
                    struct Number*out = number_slot_new(number_pool_cursor, number_free_list);
                    out->operation = 'r';
                    out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move(pool_integer_slot_new, integer_pool, &out->value.numerator);
                    return out;
                }
                pool_integer_free(integer_pool, exponent->value.denominator);
                exponent->value.denominator = reduced_degree.quotient;
                integer_move(pool_integer_slot_new, integer_pool, &exponent->value.denominator);
            }
            struct Number*number_coefficient =
                number_slot_new(number_pool_cursor, number_free_list);
            number_coefficient->operation = 'r';
            number_coefficient->value.numerator = coefficient;
            integer_move(pool_integer_slot_new, integer_pool, &number_coefficient->value.numerator);
            number_coefficient->value.denominator =
                integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            struct Number*number_radicand = number_slot_new(number_pool_cursor, number_free_list);
            number_radicand->operation = 'r';
            number_radicand->value.denominator =
                integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            number_radicand->value.numerator =
                integer_new(stack_integer_slot_new, &stack_a_cursor, 1, radicand_sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Division reduced_multiplicity;
                integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &reduced_multiplicity,
                    factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(stack_integer_slot_new,
                    &stack_a_cursor, stack_a_cursor, stack_b_cursor, factors[factor_index].value,
                    reduced_multiplicity.quotient);
                number_radicand->value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
                    number_radicand->value.numerator, exponentiation);
            }
            integer_move(pool_integer_slot_new, integer_pool, &number_radicand->value.numerator);
            struct Number*surd = number_slot_new(number_pool_cursor, number_free_list);
            surd->operation = '^';
            surd->left = number_radicand;
            surd->right = exponent;
            return number_multiply(number_pool_cursor, number_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, number_coefficient, surd);
        }
        struct Number*new_denominator = number_slot_new(number_pool_cursor, number_free_list);
        new_denominator->operation = 'r';
        new_denominator->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        new_denominator->value.numerator = integer_exponentiate(pool_integer_slot_new, integer_pool,
            stack_a_cursor, stack_b_cursor, base->value.denominator, exponent->value.numerator);
        struct Number*new_numerator_base = number_slot_new(number_pool_cursor, number_free_list);
        new_numerator_base->operation = 'r';
        new_numerator_base->value.denominator =
            integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        new_numerator_base->value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
            base->value.numerator, integer_exponentiate(stack_integer_slot_new, &stack_a_cursor,
                stack_a_cursor, stack_b_cursor, base->value.denominator,
                integer_subtract(&stack_a_cursor, stack_b_cursor, exponent->value.denominator,
                    &INT(1, 1))));
        integer_move(pool_integer_slot_new, integer_pool, &new_numerator_base->value.numerator);
        struct Number*new_numerator = number_exponentiate(number_pool_cursor, number_free_list,
            integer_pool, stack_a_cursor, stack_b_cursor, new_numerator_base, exponent);
        number_free(number_free_list, integer_pool, new_numerator_base);
        number_reciprocal(new_denominator);
        struct Number*out = number_multiply(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, new_numerator, new_denominator);
        number_free(number_free_list, integer_pool, new_numerator);
        number_free(number_free_list, integer_pool, new_denominator);
        return out;
    }
    return true;
}

struct Number*number_evaluate_root(struct Number**number_pool_cursor,
    struct Number**number_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*a)
{
    switch (a->operation)
    {
    case '+':
        return number_add(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor, a->left, a->right);
    case '-':
    {
        struct Number*negative = number_slot_new(number_pool_cursor, number_free_list);
        negative->operation = 'r';
        negative->value.numerator = integer_new(pool_integer_slot_new, integer_pool, 1, -1);
        negative->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        struct Number*product = number_multiply(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, negative, a->right);
        if (!product)
        {
            return 0;
        }
        return number_add(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor, a->left, product);
    }
    case '*':
        return number_multiply(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor, a->left, a->right);
    case '/':
        number_reciprocal(a->right);
        return number_multiply(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor, a->left, a->right);
    case '^':
        return number_exponentiate(number_pool_cursor, number_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a->left, a->right);
    }
}

struct Number*number_evaluate(struct Number**number_pool_cursor, struct Number**number_free_list,
    struct IntegerPool*integer_pool, uint32_t*stack_a_cursor, uint32_t*stack_b_cursor,
    struct Number*a)
{
    if (a->operation == 'r')
    {
        return a;
    }
    a->left = number_evaluate(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, a->left);
    if (!a->left)
    {
        return 0;
    }
    a->right = number_evaluate(number_pool_cursor, number_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, a->right);
    if (!a->right)
    {
        return 0;
    }
    return number_evaluate_root(number_pool_cursor, number_free_list, integer_pool,
        stack_a_cursor, stack_b_cursor, a);
}

void print_number(void*stack_a_cursor, void*stack_b_cursor, struct Number*number)
{
    if (number->operation == 'r')
    {
        char*string = stack_a_cursor;
        integer_string(&stack_a_cursor, stack_b_cursor, number->value.numerator);
        if (!integer_equals_one(number->value.denominator))
        {
            *(char*)stack_slot_new(&stack_a_cursor, 1) = '/';
            integer_string(&stack_a_cursor, stack_b_cursor, number->value.denominator);
        }
        *(char*)stack_slot_new(&stack_a_cursor, 1) = 0;
        printf("%s", string);
    }
    else
    {
        printf("(");
        print_number(stack_a_cursor, stack_b_cursor, number->left);
        printf(")");
        printf("%c", number->operation);
        printf("(");
        print_number(stack_a_cursor, stack_b_cursor, number->right);
        printf(")");
    }
}

void init()
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwPageSize;
    primes = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    VirtualAlloc(primes, page_size, MEM_COMMIT, PAGE_READWRITE);
    struct Integer*prime = primes;
    integer_new(stack_integer_slot_new, &prime, 2, 1);
    integer_new(stack_integer_slot_new, &prime, 3, 1);
    permanent_integer_pool = VirtualAlloc(0, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    integer_pool_new(permanent_integer_pool);
    pi_estimate_min.numerator = integer_new(pool_integer_slot_new, permanent_integer_pool, 47, 1);
    pi_estimate_min.denominator = integer_new(pool_integer_slot_new, permanent_integer_pool, 15, 1);
    pi_interval_size.numerator =
        integer_new(pool_integer_slot_new, permanent_integer_pool, 1696, 1);
    pi_interval_size.denominator =
        integer_new(pool_integer_slot_new, permanent_integer_pool, 12285, 1);
    pi_sixteen_to_the_k = integer_new(pool_integer_slot_new, permanent_integer_pool, 16, 1);
    pi_eight_k = integer_new(pool_integer_slot_new, permanent_integer_pool, 8, 1);
}

int main()
{
    init();
    struct Number*number_pool = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    struct Number*number_pool_cursor = number_pool;
    struct Number*number_free_list = 0;
    void*stack_a = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_a_cursor = stack_a;
    void*stack_b = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_b_cursor = stack_b;
    while (true)
    {
        struct IntegerPool*integer_pool =
            VirtualAlloc(0, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
        integer_pool_new(integer_pool);
        if (get_input(&number_pool_cursor, &number_free_list, integer_pool, stack_a_cursor,
            stack_b_cursor))
        {
            struct Number*number = number_pool;
            while (number->next)
            {
                if ((number->operation == 'r' && number->next->operation == '(') ||
                    (number->operation == ')' &&
                    (number->next->operation == 'r' || number->next->operation == '(')))
                {
                    struct Number*times = number_slot_new(&number_pool_cursor, &number_free_list);
                    times->operation = '*';
                    times->previous = number;
                    times->next = number->next;
                    number->next->previous = times;
                    number->next = times;
                }
                number = number->next;
            }
            struct Number*input =
                parse_input(&number_pool_cursor, &number_free_list, integer_pool, number_pool);
            if (input)
            {
                struct Number*evaluation = number_evaluate(&number_pool_cursor,
                    &number_free_list, integer_pool, stack_a_cursor, stack_b_cursor, input);
                if (evaluation)
                {
                    printf("=\n");
                    print_number(stack_a_cursor, stack_b_cursor, evaluation);
                }
            }
        }
        printf("\n\n");
        rewind_stack_cursor(&number_pool_cursor, number_pool);
        number_free_list = 0;
        rewind_stack_cursor(&stack_a_cursor, stack_a);
        rewind_stack_cursor(&stack_b_cursor, stack_b);
        while (integer_pool->first_page)
        {
            struct IntegerPool*next_value_count_pool = integer_pool + 1;
            struct IntegerSlotPage*page = integer_pool->first_page;
            while (page)
            {
                struct IntegerSlotPage*next_page = page->next_page;
                VirtualFree(page, 0, MEM_RELEASE);
                page = next_page;
            }
            integer_pool = next_value_count_pool;
        }
    }
    return 0;
}