#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#ifdef _DEBUG
#define ASSERT(condition, message) if(!(condition)){printf(message); *(int*)0 = 0;}
#else
#define ASSERT(condition, message)
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

uint32_t*int_stack_slot_new(void**stack_cursor, size_t value_count)
{
    return stack_slot_new(stack_cursor, value_count * sizeof(uint32_t));
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

struct SizedIntPage
{
    struct SizedIntPage*next_page;
    uint32_t memory;
};

struct IntPool
{
    uint32_t*cursor;
    uint32_t*free_list;
    struct IntPool*next_value_count_pool;
    struct SizedIntPage page;
};

struct IntPool*sized_int_pool_new()
{
    struct IntPool*pool = VirtualAlloc(0, page_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    pool->cursor = &pool->page.memory;
    return pool;
}

struct IntPool*get_value_count_pool(struct IntPool*int_pool, size_t value_count)
{
    if (value_count == 1)
    {
        return int_pool;
    }
    if (!int_pool->next_value_count_pool)
    {
        int_pool->next_value_count_pool = sized_int_pool_new();
    }
    return get_value_count_pool(int_pool->next_value_count_pool, value_count - 1);
}

uint32_t**int_pool_slot_next(uint32_t*slot, size_t value_count)
{
    return slot + value_count;
}

uint32_t*int_pool_slot_new(struct IntPool*int_pool, size_t value_count)
{
    if (!value_count)
    {
        return 0;
    }
    struct IntPool*pool = get_value_count_pool(int_pool, value_count);
    struct uint32_t*allocation = pool->free_list;
    if (allocation)
    {
        pool->free_list = *int_pool_slot_next(allocation, value_count);
    }
    else
    {
        allocation = pool->cursor;
        size_t boundary = next_page_boundary(allocation);
        pool->cursor = (size_t)(pool->cursor + value_count) + sizeof(uint32_t*);
        if (pool->cursor > boundary)
        {
            struct SizedIntPage*new_page =
                VirtualAlloc(0, page_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
            pool->cursor = &new_page->memory;
            size_t previous_page_address = boundary - page_size;
            if (previous_page_address == int_pool)
            {
                int_pool->page.next_page = new_page;
            }
            else
            {
                ((struct SizedIntPage*)previous_page_address)->next_page = new_page;
            }
        }
        *int_pool_slot_next(allocation, value_count) = 0;
    }
    return allocation;
}

struct Integer
{
    size_t value_count;
    uint32_t*values;
    int8_t sign;
};

void integer_new(uint32_t*(allocator)(void*, size_t), void*memory, struct Integer*out,
    int32_t value)
{
    if (value == 0)
    {
        memset(out, 0, sizeof(struct Integer));
        return;
    }
    out->value_count = 1;
    out->values = allocator(memory, 1);
    if (value > 0)
    {
        out->values[0] = value;
        out->sign = 1;
    }
    else
    {
        out->values[0] = -value;
        out->sign = -1;
    }
}

void int_pool_slot_free(struct IntPool*int_pool, struct Integer*a)
{
    struct IntPool*pool = get_value_count_pool(int_pool, a->value_count);
    *int_pool_slot_next(a->values, a->value_count) = pool->free_list;
    pool->free_list = a->values;
}

void integer_copy(void**stack_cursor, struct Integer*out, struct Integer*a)
{
    ASSERT(a != out, "integer_copy a parameter aliased with out parameter.");
    out->value_count = a->value_count;
    out->values = int_stack_slot_new(stack_cursor, a->value_count);
    out->sign = a->sign;
    memcpy(out->values, a->values, a->value_count * sizeof(uint32_t));
}

void integer_copy_values_to_allocator(uint32_t*(allocator)(void*, size_t), void*memory,
    struct Integer*a)
{
    uint32_t*slot = allocator(memory, a->value_count);
    memcpy(slot, a->values, a->value_count * sizeof(uint32_t));
    a->values = slot;
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
        if (a->values[i] != b->values[i])
        {
            return false;
        }
    }
    return true;
}

bool integer_equals_one(struct Integer*a)
{
    return a->value_count == 1 && a->values[0] == 1 && a->sign > 0;
}

void calculate_sum_values(struct Integer*short_integer, struct Integer*sum)
{
    for (int i = 0; i < short_integer->value_count; ++i)
    {
        uint64_t remainder = short_integer->values[i];
        for (int j = i; j < sum->value_count; ++j)
        {
            uint64_t sum_value = sum->values[j] + remainder;
            sum->values[j] = (uint32_t)sum_value;
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
        integer->values[i] = ~integer->values[i];
    }
    for (int i = 0; i < integer->value_count; ++i)
    {
        uint32_t power = 1;
        for (int j = 0; j < 32; ++j)
        {
            integer->values[i] ^= power;
            if ((integer->values[i] & power) != 0)
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
        if (a->values[i] != 0)
        {
            a->value_count = i + 1;
            return;
        }
    }
    a->value_count = 0;
    a->sign = 0;
}

void integer_out_add(void**stack_cursor, struct Integer*out, struct Integer*a, struct Integer*b)
{
    ASSERT(a != out && b != out, "integer_out_add operand aliased with out parameter.");
    if (a->sign == 0)
    {
        integer_copy(stack_cursor, out, b);
        return;
    }
    if (b->sign == 0)
    {
        integer_copy(stack_cursor, out, a);
        return;
    }
    struct Integer*long_integer = a;
    struct Integer*short_integer = b;
    if (a->value_count < b->value_count)
    {
        long_integer = b;
        short_integer = a;
    }
    out->sign = short_integer->sign;
    out->value_count = long_integer->value_count + 1;
    out->values = int_stack_slot_new(stack_cursor, out->value_count);
    memcpy(out->values, long_integer->values, long_integer->value_count * sizeof(uint32_t));
    out->values[long_integer->value_count] = 0;
    if (short_integer->sign == long_integer->sign)
    {
        calculate_sum_values(short_integer, out);
    }
    else
    {
        twos_complement(out);
        calculate_sum_values(short_integer, out);
        if (out->values[long_integer->value_count] != 0)
        {
            twos_complement(out);
            out->sign *= -1;
        }
    }
    trim_leading_zeroes(out);
    rewind_stack_cursor(stack_cursor, out->values + out->value_count);
}

struct Integer integer_add(void**stack_cursor, struct Integer*a, struct Integer*b)
{
    struct Integer sum;
    integer_out_add(stack_cursor, &sum, a, b);
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
        a->value_count = b->value_count;
        a->sign = b->sign;
        memcpy(a->values, b->values, b->value_count * sizeof(uint32_t));
        return;
    }
    if (a->value_count < b->value_count)
    {
        memset(a->values + a->value_count, 0,
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
        if (a->values[a->value_count - 1] != 0)
        {
            twos_complement(a);
            a->sign *= -1;
        }
    }
    trim_leading_zeroes(a);
}

void integer_out_subtract(void**stack_cursor, struct Integer*out, struct Integer*minuend,
    struct Integer*subtrahend)
{
    ASSERT(minuend != out && subtrahend != out,
        "integer_out_subtract operand aliased with out parameter.");
    uint32_t*int_arena_savepoint =
        int_stack_slot_new(stack_cursor, 1 + max(minuend->value_count, subtrahend->value_count));
    struct Integer negative;
    integer_copy(stack_cursor, &negative, subtrahend);
    negative.sign *= -1;
    integer_out_add(&int_arena_savepoint, out, minuend, &negative);
    rewind_stack_cursor(stack_cursor, int_arena_savepoint);
}

struct Integer integer_subtract(void**stack_cursor, struct Integer*minuend,
    struct Integer*subtrahend)
{
    struct Integer difference;
    integer_out_subtract(stack_cursor, &difference, minuend, subtrahend);
    return difference;
}

bool integer_less_than(void*stack_cursor, struct Integer*a, struct Integer*b)
{
    struct Integer difference;
    integer_out_subtract(&stack_cursor, &difference, b, a);
    return difference.sign > 0;
}

void integer_out_multiply(void**output_stack_cursor, void*scratch_stack_cursor, struct Integer*out,
    struct Integer*a, struct Integer*b)
{
    ASSERT(a != out && b != out, "integer_out_multiply operand aliased with out parameter.");
    out->value_count = 0;
    out->values = int_stack_slot_new(output_stack_cursor, a->value_count + b->value_count);
    out->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)a->values[i] * b->values[j];
            size_t shift = i + j;
            struct Integer integer_component =
                { shift, int_stack_slot_new(&scratch_stack_cursor, shift + 2), 0 };
            memset(integer_component.values, 0, shift * sizeof(uint32_t));
            if (product_component > 0)
            {
                integer_component.values[shift] = product_component;
                integer_component.value_count += 1;
                integer_component.sign = 1;
                uint32_t high_bytes = (product_component & 0xffffffff00000000) >> 32;
                if (high_bytes > 0)
                {
                    integer_component.values[shift + 1] = high_bytes;
                    integer_component.value_count += 1;
                }
            }
            integer_add_to_a_in_place(out, &integer_component);
        }
    }
    out->sign = a->sign * b->sign;
    rewind_stack_cursor(output_stack_cursor, out->values + out->value_count);
}

struct Integer integer_multiply(void**output_stack_cursor, void*scratch_stack_cursor,
    struct Integer*a, struct Integer*b)
{
    struct Integer product;
    integer_out_multiply(output_stack_cursor, scratch_stack_cursor, &product, a, b);
    return product;
}

struct Division
{
    struct Integer quotient;
    struct Integer remainder;
};

int leading_digit_place(struct Integer*a)
{
    if (!a->value_count)
    {
        return 0;
    }
    uint32_t divisor_leading_digit = 0x80000000;
    uint32_t last_value = a->values[a->value_count - 1];
    for (int i = 31; i >= 0; --i)
    {
        if ((last_value & divisor_leading_digit) != 0)
        {
            return i;
        }
        divisor_leading_digit = divisor_leading_digit >> 1;
    }
}

void calculate_division_values(void*stack_cursor, struct Integer*positive_divisor,
    struct Division*division, size_t quotient_value_index, uint32_t quotient_digit)
{
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer difference;
            integer_out_subtract(&stack_cursor, &difference, &division->remainder,
                positive_divisor);
            if (difference.sign >= 0)
            {
                division->quotient.values[quotient_value_index] |= quotient_digit;
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
            for (int j = 0; j < positive_divisor->value_count - 1; ++j)
            {
                positive_divisor->values[j] = positive_divisor->values[j] >> 1 |
                    positive_divisor->values[j + 1] << 31;
            }
            positive_divisor->values[positive_divisor->value_count - 1] =
                positive_divisor->values[positive_divisor->value_count - 1] >> 1;
        }
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
        struct Integer positive_divisor = { dividend->value_count,
            int_stack_slot_new(&scratch_stack_cursor, dividend->value_count), 1 };
        size_t quotient_value_index = dividend->value_count - divisor->value_count;
        memset(positive_divisor.values, 0, quotient_value_index * sizeof(uint32_t));
        memcpy(positive_divisor.values + quotient_value_index, divisor->values,
            divisor->value_count * sizeof(uint32_t));
        int shift = dividend_leading_digit_place - divisor_leading_digit_place;
        uint32_t quotient_digit;
        if (shift >= 0)
        {
            positive_divisor.values[quotient_value_index] =
                positive_divisor.values[quotient_value_index] << shift;
            for (int i = quotient_value_index + 1; i < positive_divisor.value_count; ++i)
            {
                positive_divisor.values[i] = positive_divisor.values[i] << shift |
                    positive_divisor.values[i - 1] >> 32 - shift;
            }
            quotient_digit = 1 << shift;
        }
        else
        {
            shift *= -1;
            for (int i = quotient_value_index; i < positive_divisor.value_count; ++i)
            {
                positive_divisor.values[i - 1] = positive_divisor.values[i] << 32 - shift |
                    positive_divisor.values[i - 1] >> shift;
            }
            positive_divisor.values[positive_divisor.value_count - 1] =
                positive_divisor.values[positive_divisor.value_count - 1] >> shift;
            quotient_digit = 1 << 32 - shift;
            quotient_value_index -= 1;
        }
        int8_t dividend_sign = dividend->sign;
        dividend->sign = 1;
        out->quotient.value_count = dividend->value_count;
        out->quotient.values = *output_stack_cursor;
        out->quotient.sign = quotient_sign;
        memset(out->quotient.values, 0, out->quotient.value_count * sizeof(uint32_t));
        out->remainder = *dividend;
        calculate_division_values(scratch_stack_cursor, &positive_divisor, out,
            quotient_value_index, quotient_digit);
        trim_leading_zeroes(&out->quotient);
        *(uint32_t**)output_stack_cursor += out->quotient.value_count;
        integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->remainder);
        if (out->remainder.sign != 0)
        {
            out->remainder.sign = dividend->sign;
        }
        dividend->sign = dividend_sign;
        return;
    }
    out->quotient.value_count = 0;
    out->quotient.values = 0;
    out->quotient.sign = 0;
    integer_copy(output_stack_cursor, &out->remainder, dividend);
}

void integer_exponentiate(uint32_t*(allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, uint32_t*scratch_stack_cursor, struct Integer*out,
    struct Integer*base, struct Integer*exponent)
{
    integer_new(int_stack_slot_new, &scratch_stack_cursor, out, 1);
    struct Integer base_to_a_power_of_two = *base;
    struct Integer remaining_exponent = *exponent;
    struct Integer two;
    integer_new(int_stack_slot_new, &scratch_stack_cursor, &two, 2);
    while (remaining_exponent.sign > 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, output_stack_cursor, &division,
            &remaining_exponent, &two);
        if (division.remainder.sign > 0)
        {
            *out = integer_multiply(&scratch_stack_cursor, output_stack_cursor, out,
                &base_to_a_power_of_two);
        }
        base_to_a_power_of_two = integer_multiply(&scratch_stack_cursor, output_stack_cursor,
            &base_to_a_power_of_two, &base_to_a_power_of_two);
        remaining_exponent = division.quotient;
    }
    integer_copy_values_to_allocator(allocator, memory, out);
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
    struct Integer quotient = *a;
    struct Integer power;
    integer_new(int_stack_slot_new, &scratch_stack_cursor, &power, 10);
    while (quotient.sign != 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, *output_stack_cursor, &division, &quotient,
            &power);
        quotient = division.quotient;
        next_char -= 1;
        if (division.remainder.sign != 0)
        {
            *next_char = division.remainder.values[0] + '0';
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

void get_gcd(void*stack_a_cursor, void*stack_b_cursor, struct Integer*out, struct Integer*a,
    struct Integer*b)
{
    while (b->value_count > 0)
    {
        struct Integer*c = b;
        struct Division*division = stack_slot_new(&stack_a_cursor, sizeof(struct Division));
        integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, division, a, b);
        b = &division->remainder;
        a = c;
    }
    *out = *a;
}

struct ExtendedGCDInfo
{
    struct Integer gcd;
    struct Integer a_coefficient;
    struct Integer b_coefficient;
    struct Integer a_over_gcd;
    struct Integer b_over_gcd;
};

void extended_gcd(void**output_stack_cursor, void*scratch_stack_cursor, struct ExtendedGCDInfo*out,
    struct Integer*a, struct Integer*b)
{
    struct Integer zero = { 0, 0, 0 };
    struct Integer*a_coefficient = &zero;
    struct Integer*b_coefficient = stack_slot_new(&scratch_stack_cursor, sizeof(struct Integer));
    integer_new(int_stack_slot_new, &scratch_stack_cursor, b_coefficient, 1);
    struct Integer*a_over_gcd = &zero;
    struct Integer*b_over_gcd = b_coefficient;
    while (!integer_equals(a, &zero))
    {
        struct Division*division = stack_slot_new(&scratch_stack_cursor, sizeof(struct Division));
        integer_euclidean_divide(output_stack_cursor, scratch_stack_cursor, division, b, a);
        struct Integer product;
        integer_out_multiply(&scratch_stack_cursor, *output_stack_cursor, &product, b_over_gcd,
            &division->quotient);
        struct Integer*m = stack_slot_new(&scratch_stack_cursor, sizeof(struct Integer));
        integer_out_subtract(&scratch_stack_cursor, m, a_coefficient, &product);
        integer_out_multiply(&scratch_stack_cursor, *output_stack_cursor, &product,
            a_over_gcd, &division->quotient);
        struct Integer*n = stack_slot_new(&scratch_stack_cursor, sizeof(struct Integer));
        integer_out_subtract(&scratch_stack_cursor, n, b_coefficient, &product);
        b = a;
        a = &division->remainder;
        a_coefficient = b_over_gcd;
        b_coefficient = a_over_gcd;
        b_over_gcd = m;
        a_over_gcd = n;
    }
    out->gcd = *b;
    integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->gcd);
    out->a_coefficient = *a_coefficient;
    integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->a_coefficient);
    out->b_coefficient = *b_coefficient;
    integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->b_coefficient);
    out->a_over_gcd = *a_over_gcd;
    integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->a_over_gcd);
    out->b_over_gcd = *b_over_gcd;
    out->b_over_gcd.sign *= -1;
    integer_copy_values_to_allocator(int_stack_slot_new, output_stack_cursor, &out->b_over_gcd);
}

struct Rational
{
    struct Integer numerator;
    struct Integer denominator;
};

bool rational_new(uint32_t*(allocator)(void*, size_t), void*memory, void*output_stack_cursor,
    void*scratch_stack_cursor, struct Rational*out, struct Integer*numerator,
    struct Integer*denominator)
{
    if (denominator->sign == 0)
    {
        printf("Tried to divide by 0.");
        return false;
    }
    if (numerator->sign == 0)
    {
        out->numerator.value_count = 0;
        out->numerator.sign = 0;
        integer_new(allocator, memory, &out->denominator, 1);
        return true;
    }
    struct ExtendedGCDInfo gcd_info;
    extended_gcd(&scratch_stack_cursor, output_stack_cursor, &gcd_info, numerator, denominator);
    out->numerator = gcd_info.a_over_gcd;
    out->denominator = gcd_info.b_over_gcd;
    if (out->denominator.sign < 0)
    {
        out->numerator.sign *= -1;
        out->denominator.sign *= -1;
    }
    integer_copy_values_to_allocator(allocator, memory, &out->numerator);
    integer_copy_values_to_allocator(allocator, memory, &out->denominator);
    return true;
}

bool rational_out_add(uint32_t*(allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, uint32_t*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Rational*b)
{
    struct Integer numerator_term_a;
    integer_out_multiply(&scratch_stack_cursor, output_stack_cursor, &numerator_term_a,
        &a->numerator, &b->denominator);
    struct Integer numerator_term_b;
    integer_out_multiply(&scratch_stack_cursor, output_stack_cursor, &numerator_term_b,
        &b->numerator, &a->denominator);
    struct Integer numerator;
    integer_out_add(&scratch_stack_cursor, &numerator, &numerator_term_a, &numerator_term_b);
    struct Integer denominator;
    integer_out_multiply(&scratch_stack_cursor, output_stack_cursor, &denominator, &a->denominator,
        &b->denominator);
    return rational_new(allocator, memory, output_stack_cursor, scratch_stack_cursor, out,
        &numerator, &denominator);
}

bool rational_out_multiply(uint32_t*(allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Rational*out, struct Rational*a,
    struct Rational*b)
{
    struct Integer numerator;
    integer_out_multiply(&scratch_stack_cursor, output_stack_cursor, &numerator, &a->numerator,
        &b->numerator);
    struct Integer denominator;
    integer_out_multiply(&scratch_stack_cursor, output_stack_cursor, &denominator, &a->denominator,
        &b->denominator);
    return rational_new(allocator, memory, output_stack_cursor, scratch_stack_cursor, out,
        &numerator, &denominator);
}

struct Rational rational_multiply(uint32_t*(allocator)(void*, size_t), void*memory,
    void*stack_a_cursor, void*stack_b_cursor, struct Rational*a, struct Rational*b)
{
    struct Rational product;
    rational_out_multiply(allocator, memory, stack_a_cursor, stack_b_cursor, &product, a, b);
    return product;
}

struct Number
{
    union
    {
        struct Rational value;
        struct Node
        {
            struct Number*left;
            struct Number*right;
        };
    };
    char operation;
};

struct NumberSlot
{
    struct Number number;
    struct NumberSlot*previous;
    struct NumberSlot*next;
};

struct NumberSlot*number_slot_new(struct NumberSlot**pool_cursor, struct NumberSlot**free_list)
{
    struct NumberSlot*allocation = *free_list;
    if (allocation)
    {
        *free_list = allocation->next;
    }
    else
    {
        allocation = stack_slot_new(pool_cursor, sizeof(struct NumberSlot));
        allocation->next = 0;
    }
    return allocation;
}

void number_slot_free(struct NumberSlot**number_slot_free_list, struct NumberSlot*slot,
    struct IntPool*int_pool)
{
    if (slot->number.operation == 'r')
    {
        if (slot->number.value.numerator.value_count)
        {
            int_pool_slot_free(int_pool, &slot->number.value.numerator);
        }
        int_pool_slot_free(int_pool, &slot->number.value.denominator);
    }
    slot->next = *number_slot_free_list;
    *number_slot_free_list = slot;
}

bool get_input(struct NumberSlot**number_slot_pool_cursor, struct NumberSlot**number_free_list,
    struct IntPool*int_pool, void*stack_a_cursor, void*stack_b_cursor)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    struct NumberSlot*previous_slot = 0;
    while (true)
    {
        struct NumberSlot*slot = number_slot_new(number_slot_pool_cursor, number_free_list);
        slot->previous = previous_slot;
        slot->next = *number_slot_pool_cursor;
        if (isdigit(next_char))
        {
            slot->number.operation = 'r';
            integer_new(int_pool_slot_new, int_pool, &slot->number.value.denominator, 1);
            integer_new(int_stack_slot_new, &stack_a_cursor, &slot->number.value.numerator,
                next_char - '0');
            struct Integer ten;
            integer_new(int_stack_slot_new, &stack_a_cursor, &ten, 10);
            next_char = getchar();
            while (isdigit(next_char))
            {
                slot->number.value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
                    &slot->number.value.numerator, &ten);
                struct Integer digit;
                integer_new(int_stack_slot_new, &stack_a_cursor, &digit, next_char - '0');
                slot->number.value.numerator =
                    integer_add(&stack_a_cursor, &slot->number.value.numerator, &digit);
                next_char = getchar();
            }
            integer_copy_values_to_allocator(int_pool_slot_new, int_pool,
                &slot->number.value.numerator);
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
            slot->number.operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            slot->next = 0;
            return true;
        }
        previous_slot = slot;
    }
}

bool is_numeric(struct NumberSlot*slot)
{
    if (!slot)
    {
        return false;
    }
    return slot->number.operation == 'r' || slot->number.left;
}

bool is_unparsed_operation(struct NumberSlot*slot, char operation)
{
    return slot->number.operation == operation && !slot->number.left;
}

bool parse_binary_operation(struct NumberSlot*slot, char operation)
{
    if (!is_unparsed_operation(slot, operation))
    {
        return true;
    }
    if (!slot->previous || !is_numeric(slot->previous))
    {
        printf("%c missing left operand.", operation);
        return false;
    }
    if (slot->next && is_numeric(slot->next))
    {
        slot->number.left = &slot->previous->number;
        slot->number.right = &slot->next->number;
        slot->previous = slot->previous->previous;
        if (slot->previous)
        {
            slot->previous->next = slot;
        }
        slot->next = slot->next->next;
        if (slot->next)
        {
            slot->next->previous = slot;
        }
    }
    else
    {
        printf("%c missing right operand.", operation);
        return false;
    }
    return true;
}

void rewind_to_first_slot(struct NumberSlot**slot)
{
    while ((*slot)->previous)
    {
        *slot = (*slot)->previous;
    }
}

bool parse_binary_operation_pair(struct NumberSlot**slot, char operation_a, char operation_b)
{
    rewind_to_first_slot(slot);
    while (true)
    {
        if (is_unparsed_operation(*slot, operation_a) &&
            !parse_binary_operation(*slot, operation_a))
        {
            return false;
        }
        else if (is_unparsed_operation(*slot, operation_b) &&
            !parse_binary_operation(*slot, operation_b))
        {
            return false;
        }
        if (!(*slot)->next)
        {
            return true;
        }
        *slot = (*slot)->next;
    }
}

struct NumberSlot*parse_input(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntPool*int_pool,
    struct NumberSlot*number_slot)
{
    if (!number_slot)
    {
        printf("Empty expression.");
        return 0;
    }
    while (true)
    {
        if (number_slot->number.operation == ')')
        {
            printf("Unmatched ).");
            return 0;
        }
        if (number_slot->number.operation == '(')
        {
            int unmatched_paren_count = 1;
            struct NumberSlot*nested_number_slot = number_slot;
            while (unmatched_paren_count > 0)
            {
                nested_number_slot = nested_number_slot->next;
                if (!nested_number_slot)
                {
                    printf("Unmatched (.");
                    return 0;
                }
                if (nested_number_slot->number.operation == '(')
                {
                    unmatched_paren_count += 1;
                }
                else if (nested_number_slot->number.operation == ')')
                {
                    unmatched_paren_count -= 1;
                }
            }
            struct NumberSlot*previous = number_slot->previous;
            struct NumberSlot*next = nested_number_slot->next;
            number_slot->next->previous = 0;
            nested_number_slot->previous->next = 0;
            struct NumberSlot*nested_expression = parse_input(number_slot_pool_cursor,
                number_slot_free_list, int_pool, number_slot->next);
            if (nested_expression)
            {
                nested_expression->previous = previous;
                if (previous)
                {
                    previous->next = nested_expression;
                }
                nested_expression->next = next;
                if (next)
                {
                    next->previous = nested_expression;
                }
                number_slot_free(number_slot_free_list, number_slot, int_pool);
                number_slot_free(number_slot_free_list, nested_number_slot, int_pool);
                number_slot = nested_expression;
            }
            else
            {
                return 0;
            }
        }
        if (!number_slot->next)
        {
            break;
        }
        number_slot = number_slot->next;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        struct NumberSlot*next_slot = number_slot->next;
        if (is_unparsed_operation(number_slot, '+') && !is_numeric(number_slot->previous))
        {
            if (!next_slot || (!is_numeric(next_slot) && !is_unparsed_operation(next_slot, '+') &&
                !is_unparsed_operation(next_slot, '-')))
            {
                printf("+ missing right operand.");
                return 0;
            }
            next_slot->previous = number_slot->previous;
            if (number_slot->previous)
            {
                number_slot->previous->next = next_slot;
            }
            number_slot_free(number_slot_free_list, number_slot, int_pool);
        }
        if (!next_slot)
        {
            break;
        }
        number_slot = next_slot;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        if (!parse_binary_operation(number_slot, '^'))
        {
            return 0;
        }
        if (!number_slot->next)
        {
            break;
        }
        number_slot = number_slot->next;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        if (is_unparsed_operation(number_slot, '-') && !is_numeric(number_slot->previous))
        {
            if (!number_slot->next)
            {
                printf("- missing right operand.");
                return 0;
            }
            if (is_numeric(number_slot->next))
            {
                number_slot->number.operation = 'r';
                integer_new(int_pool_slot_new, int_pool, &number_slot->number.value.numerator, -1);
                integer_new(int_pool_slot_new, int_pool, &number_slot->number.value.denominator, 1);
                struct NumberSlot*times =
                    number_slot_new(number_slot_pool_cursor, number_slot_free_list);
                times->number.operation = '*';
                times->number.left = &number_slot->number;
                times->number.right = &number_slot->next->number;
                times->next = number_slot->next->next;
                if (number_slot->previous)
                {
                    number_slot->previous->next = times;
                }
                number_slot = times;
                if (number_slot->next)
                {
                    number_slot->next->previous = number_slot;
                }
                else
                {
                    break;
                }
                number_slot = number_slot->next;
            }
            else if (is_unparsed_operation(number_slot->next, '-'))
            {
                struct NumberSlot*new_next = number_slot->next->next;
                if (!new_next)
                {
                    printf("- missing right operand.");
                    return 0;
                }
                new_next->previous = number_slot->previous;
                if (number_slot->previous)
                {
                    number_slot->previous->next = new_next;
                }
                number_slot_free(number_slot_free_list, number_slot->next, int_pool);
                number_slot_free(number_slot_free_list, number_slot, int_pool);
                number_slot = new_next;
            }
            else
            {
                printf("- missing left operand.");
                return 0;
            }
        }
        else
        {
            if (!number_slot->next)
            {
                break;
            }
            number_slot = number_slot->next;
        }
    }
    if (!parse_binary_operation_pair(&number_slot, '*', '/'))
    {
        return 0;
    }
    if (!parse_binary_operation_pair(&number_slot, '+', '-'))
    {
        return 0;
    }
    return number_slot;
}

bool number_add(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntPool*int_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*a, struct Number*b)
{
    if (a->left->operation == 'r')
    {
        if (a->right->operation == 'r')
        {
            out->operation = 'r';
            if (!rational_out_add(int_pool_slot_new, int_pool, stack_a_cursor, stack_b_cursor,
                &out->value, &a->value, &b->value))
            {
                return false;
            }
            number_slot_free(number_slot_free_list, a, int_pool);
            number_slot_free(number_slot_free_list, b, int_pool);
        }
    }
    else
    {
        return number_add(number_slot_pool_cursor, number_slot_free_list, int_pool, stack_a_cursor,
            stack_b_cursor, out, b, a);
    }
    return true;
}

bool number_multiply(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntPool*int_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        if (b->operation == 'r')
        {
            out->operation = 'r';
            if (!rational_out_multiply(int_pool_slot_new, int_pool, stack_a_cursor,
                stack_b_cursor, &out->value, &a->value, &b->value))
            {
                return false;
            }
            number_slot_free(number_slot_free_list, a, int_pool);
            number_slot_free(number_slot_free_list, b, int_pool);
            return true;
        }
    }
    else
    {
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, out, b, a);
    }
    out->operation = '*';
    out->left = a;
    out->right = b;
    return true;
}

void number_reciprocal(struct Number*a)
{
    if (a->operation == 'r')
    {
        struct Integer old_numerator = a->value.numerator;
        a->value.numerator = a->value.denominator;
        a->value.denominator = old_numerator;
    }
    else
    {
        ASSERT(0, "number_reciprocal called on Number type for which it is not yet implemented.");
    }
}

struct Integer*primes;

struct Integer*next_prime_slot(struct Integer*prime)
{
    return (size_t)(prime + 1) + prime->value_count * sizeof(uint32_t);
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
    integer_out_add(&stack_b_cursor, next_slot, prime, primes);
    while (true)
    {
        struct Integer*potential_factor = primes;
        while (true)
        {
            struct Division division;
            integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division, next_slot,
                potential_factor);
            if (division.remainder.value_count == 0)
            {
                break;
            }
            potential_factor = next_prime_slot(potential_factor);
            if (potential_factor == next_slot)
            {
                memcpy(next_slot + 1, next_slot->values, next_slot->value_count * sizeof(uint32_t));
                next_slot->values = next_slot + 1;
                return next_slot;
            }
        }
        integer_add_to_a_in_place(next_slot, primes);
    }
}

struct Factor
{
    struct Integer*value;
    struct Integer multiplicity;
};

bool number_exponentiate(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntPool*int_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*base, struct Number*exponent)
{
    if (exponent->operation == 'r')
    {
        if (base->operation == 'r')
        {
            if (exponent->value.numerator.sign < 0)
            {
                number_reciprocal(base);
                exponent->value.numerator.sign = 1;
                return number_exponentiate(number_slot_pool_cursor, number_slot_free_list, int_pool,
                    stack_a_cursor, stack_b_cursor, out, base, exponent);
            }
            if (integer_equals_one(&base->value.denominator))
            {
                if (integer_equals_one(&exponent->value.denominator))
                {
                    out->operation = 'r';
                    integer_new(int_pool_slot_new, int_pool, &out->value.denominator, 1);
                    integer_exponentiate(int_pool_slot_new, int_pool, stack_a_cursor,
                        stack_b_cursor, &out->value, &base->value.numerator,
                        &exponent->value.numerator);
                    return true;
                }
                struct Integer radicand;
                integer_exponentiate(int_stack_slot_new, &stack_a_cursor, stack_a_cursor,
                    stack_b_cursor, &radicand, &base->value.numerator, &exponent->value.numerator);
                size_t factor_count = 0;
                struct Factor*factors = stack_b_cursor;
                struct Integer*prime = primes;
                struct Integer one;
                integer_new(int_stack_slot_new, &stack_a_cursor, &one, 1);
                while (integer_less_than(stack_a_cursor, prime, &radicand))
                {
                    struct Division division;
                    integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division, &radicand,
                        prime);
                    if (division.remainder.value_count == 0)
                    {
                        factor_count += 1;
                        struct Factor*factor =
                            stack_slot_new(&stack_b_cursor, sizeof(struct Factor));
                        factor->value = prime;
                        factor->multiplicity.value_count = 0;
                        factor->multiplicity.sign = 0;
                        while (division.remainder.value_count == 0)
                        {
                            factor->multiplicity =
                                integer_add(&stack_a_cursor, &factor->multiplicity, &one);
                            radicand = division.quotient;
                            integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                                &radicand, factor->value);
                        }
                    }
                    prime = next_prime(stack_a_cursor, stack_b_cursor, prime);
                }
                if (!integer_equals_one(&radicand))
                {
                    factor_count += 1;
                    struct Factor*factor = stack_slot_new(&stack_b_cursor, sizeof(struct Factor));
                    factor->value = &radicand;
                    factor->multiplicity = one;
                }
                struct Integer coefficient;
                integer_new(int_stack_slot_new, &stack_a_cursor, &coefficient, 1);
                struct Integer multiplicity_gcd = exponent->value.denominator;
                for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
                {
                    struct Division division;
                    integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                        &factors[factor_index].multiplicity, &exponent->value.denominator);
                    factors[factor_index].multiplicity = division.remainder;
                    struct Integer exponentiation;
                    integer_exponentiate(int_stack_slot_new, &stack_a_cursor, stack_a_cursor,
                        stack_b_cursor, &exponentiation, factors[factor_index].value,
                        &division.quotient);
                    coefficient = integer_multiply(&stack_a_cursor, stack_b_cursor, &coefficient,
                        &exponentiation);
                    get_gcd(stack_a_cursor, stack_b_cursor, &multiplicity_gcd, &multiplicity_gcd,
                        &factors[factor_index].multiplicity);
                }
                struct Division division;
                integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                    &exponent->value.denominator, &multiplicity_gcd);
                number_slot_free(number_slot_free_list, base, int_pool);
                if (integer_equals_one(&division.quotient))
                {
                    number_slot_free(number_slot_free_list, exponent, int_pool);
                    out->operation = 'r';
                    integer_new(int_pool_slot_new, int_pool, &out->value.denominator, 1);
                    out->value.numerator = coefficient;
                    integer_copy_values_to_allocator(int_pool_slot_new, int_pool,
                        &out->value.numerator);
                    return true;
                }
                int_pool_slot_free(int_pool, &exponent->value.denominator);
                exponent->value.denominator = division.quotient;
                integer_copy_values_to_allocator(int_pool_slot_new, int_pool,
                    &exponent->value.denominator);
                struct Number*number_coefficient =
                    number_slot_new(number_slot_pool_cursor, number_slot_free_list);
                number_coefficient->operation = 'r';
                number_coefficient->value.numerator = coefficient;
                integer_new(int_pool_slot_new, int_pool, &number_coefficient->value.denominator, 1);
                struct Number*number_radicand =
                    number_slot_new(number_slot_pool_cursor, number_slot_free_list);
                number_radicand->operation = 'r';
                integer_new(int_pool_slot_new, int_pool, &number_radicand->value.denominator, 1);
                integer_new(int_stack_slot_new, &stack_a_cursor, &number_radicand->value.numerator,
                    1);
                for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
                {
                    struct Division division;
                    integer_euclidean_divide(&stack_a_cursor, stack_b_cursor, &division,
                        &factors[factor_index].multiplicity, &multiplicity_gcd);
                    struct Integer exponentiation;
                    integer_exponentiate(int_stack_slot_new, &stack_a_cursor, stack_a_cursor,
                        stack_b_cursor, &exponentiation, factors[factor_index].value,
                        &division.quotient);
                    number_radicand->value.numerator = integer_multiply(&stack_a_cursor,
                        stack_b_cursor, &number_radicand->value.numerator, &exponentiation);
                }
                integer_copy_values_to_allocator(int_pool_slot_new, int_pool,
                    &number_radicand->value.numerator);
                struct Number*surd =
                    number_slot_new(number_slot_pool_cursor, number_slot_free_list);
                surd->operation = '^';
                surd->left = number_radicand;
                surd->right = exponent;
                return number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
                    stack_a_cursor, stack_b_cursor, out, number_coefficient, surd);
            }
            struct Number*new_denominator =
                number_slot_new(number_slot_pool_cursor, number_slot_free_list);
            new_denominator->operation = 'r';
            integer_new(int_pool_slot_new, int_pool, &new_denominator->value.denominator, 1);
            integer_exponentiate(int_pool_slot_new, int_pool, stack_a_cursor, stack_b_cursor,
                &new_denominator->value.numerator, &base->value.denominator,
                &exponent->value.numerator);
            struct Integer one;
            integer_new(int_stack_slot_new, &stack_a_cursor, &one, 1);
            struct Integer old_denominator_exponent;
            integer_out_subtract(&stack_a_cursor, &old_denominator_exponent,
                &exponent->value.denominator, &one);
            struct Integer old_denominator_power;
            integer_exponentiate(int_stack_slot_new, &stack_a_cursor, stack_a_cursor,
                stack_b_cursor, &old_denominator_power, &base->value.denominator,
                &old_denominator_exponent);
            struct Number*new_numerator_base =
                number_slot_new(number_slot_pool_cursor, number_slot_free_list);
            new_numerator_base->operation = 'r';
            integer_new(int_pool_slot_new, int_pool, &new_numerator_base->value.denominator, 1);
            integer_out_multiply(&stack_a_cursor, stack_b_cursor,
                &new_numerator_base->value.numerator, &base->value.numerator,
                &old_denominator_power);
            integer_copy_values_to_allocator(int_pool_slot_new, int_pool,
                &new_numerator_base->value.numerator);
            struct Number*new_numerator =
                number_slot_new(number_slot_pool_cursor, number_slot_free_list);
            number_exponentiate(number_slot_pool_cursor, number_slot_free_list, int_pool,
                stack_a_cursor, stack_b_cursor, new_numerator, new_numerator_base, exponent);
            number_reciprocal(new_denominator);
            return number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
                stack_a_cursor, stack_b_cursor, out, new_numerator, new_denominator);
        }
        return true;
    }
    printf("The input expression contains an exponentiation whose exponent is not both real "
        "and rational; this program doesn't handle transcendental numbers.");
    return false;
}

bool evaluate_root(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntPool*int_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*a)
{
    switch (a->operation)
    {
    case '+':
        return number_add(number_slot_pool_cursor, number_slot_free_list, int_pool, stack_a_cursor,
            stack_b_cursor, a, a->left, a->right);
    case '-':
    {
        struct NumberSlot*negative =
            number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        negative->number.operation = 'r';
        integer_new(int_pool_slot_new, int_pool, &negative->number.value.numerator, -1);
        integer_new(int_pool_slot_new, int_pool, &negative->number.value.denominator, 1);
        struct NumberSlot*product = number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        if (!number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, product, negative, a->right))
        {
            return false;
        }
        return number_add(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, product);
    }
    case '*':
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    case '/':
        number_reciprocal(a->right);
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    case '^':
        return number_exponentiate(number_slot_pool_cursor, number_slot_free_list, int_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    }
}

bool evaluate(struct NumberSlot**number_slot_pool_cursor, struct NumberSlot**number_slot_free_list,
    struct IntPool*int_pool, uint32_t*stack_a_cursor, uint32_t*stack_b_cursor,
    struct Number*a)
{
    if (a->operation == 'r')
    {
        return true;
    }
    if (!evaluate(number_slot_pool_cursor, number_slot_free_list, int_pool, stack_a_cursor,
        stack_b_cursor, a->left))
    {
        return false;
    }
    if (!evaluate(number_slot_pool_cursor, number_slot_free_list, int_pool, stack_a_cursor,
        stack_b_cursor, a->right))
    {
        return false;
    }
    if (!evaluate_root(number_slot_pool_cursor, number_slot_free_list, int_pool, stack_a_cursor,
        stack_b_cursor, a))
    {
        return false;
    }
}

void print_number(void*stack_a_cursor, void*stack_b_cursor, struct Number*number)
{
    if (number->operation == 'r')
    {
        char*string = stack_a_cursor;
        integer_string(&stack_a_cursor, stack_b_cursor, &number->value.numerator);
        if (!integer_equals_one(&number->value.denominator))
        {
            *(char*)stack_slot_new(&stack_a_cursor, 1) = '/';
            integer_string(&stack_a_cursor, stack_b_cursor, &number->value.denominator);
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

int main()
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwPageSize;
    primes = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    VirtualAlloc(primes, page_size, MEM_COMMIT, PAGE_READWRITE);
    primes->value_count = 1;
    primes->values = primes + 1;
    primes->sign = 1;
    primes->values[0] = 2;
    struct Integer*three = next_prime_slot(primes);
    three->value_count = 1;
    three->values = three + 1;
    three->sign = 1;
    three->values[0] = 3;
    struct NumberSlot*number_slot_pool = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    struct NumberSlot*number_slot_pool_cursor = number_slot_pool;
    struct NumberSlot*number_slot_free_list = 0;
    void*stack_a = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_a_cursor = stack_a;
    void*stack_b = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_b_cursor = stack_b;
    while (true)
    {
        struct IntPool*int_pool = sized_int_pool_new();
        if (get_input(&number_slot_pool_cursor, &number_slot_free_list, int_pool, stack_a_cursor,
            stack_b_cursor))
        {
            struct NumberSlot*number_slot = number_slot_pool;
            while (number_slot->next)
            {
                if (number_slot->next->number.operation == '(' &&
                    (number_slot->number.operation == 'r' || number_slot->number.operation == ')'))
                {
                    struct NumberSlot*times =
                        number_slot_new(&number_slot_pool_cursor, &number_slot_free_list);
                    times->number.operation = '*';
                    times->previous = number_slot;
                    times->next = number_slot->next;
                    number_slot->next->previous = times;
                    number_slot->next = times;
                }
                number_slot = number_slot->next;
            }
            struct NumberSlot*input = parse_input(&number_slot_pool_cursor, &number_slot_free_list,
                int_pool, number_slot_pool);
            if (input && evaluate(&number_slot_pool_cursor, &number_slot_free_list, int_pool,
                stack_a_cursor, stack_b_cursor, &input->number))
            {
                printf("=\n");
                print_number(stack_a_cursor, stack_b_cursor, &input->number);
            }
        }
        printf("\n\n");
        rewind_stack_cursor(&number_slot_pool_cursor, number_slot_pool);
        number_slot_free_list = 0;
        rewind_stack_cursor(&stack_a_cursor, stack_a);
        rewind_stack_cursor(&stack_b_cursor, stack_b);
        while (int_pool)
        {
            struct IntPool*next_value_count_pool = int_pool->next_value_count_pool;
            struct SizedIntPage*page = &int_pool->page.next_page;
            VirtualFree(int_pool, page_size, MEM_RELEASE);
            while (page)
            {
                struct SizedIntPage*next_page = page->next_page;
                VirtualFree(page, page_size, MEM_RELEASE);
                page = next_page;
            }
            int_pool = next_value_count_pool;
        }
    }
    return 0;
}