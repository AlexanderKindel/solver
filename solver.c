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
    struct IntegerPool*value_count_pool = get_value_count_pool(pool, a->value_count);
    struct IntegerSlot*slot = (size_t)a - sizeof(struct IntegerSlot*);
    slot->next_slot = value_count_pool->free_list;
    value_count_pool->free_list = slot;
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
    return integer_subtract(&stack_a_cursor, stack_b_cursor, b, a)->sign;
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
            for (int j = 0; j < positive_divisor->value_count - 1; ++j)
            {
                (&positive_divisor->value)[j] = (&positive_divisor->value)[j] >> 1 |
                    (&positive_divisor->value)[j + 1] << 31;
            }
            (&positive_divisor->value)[positive_divisor->value_count - 1] =
                (&positive_divisor->value)[positive_divisor->value_count - 1] >> 1;
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
        if (shift >= 0)
        {
            (&positive_divisor->value)[quotient_value_index] =
                (&positive_divisor->value)[quotient_value_index] << shift;
            for (size_t i = quotient_value_index + 1; i < positive_divisor->value_count; ++i)
            {
                (&positive_divisor->value)[i] = (&positive_divisor->value)[i] << shift |
                    (&positive_divisor->value)[i - 1] >> 32 - shift;
            }
            quotient_digit = 1 << shift;
        }
        else
        {
            shift *= -1;
            for (size_t i = quotient_value_index; i < positive_divisor->value_count; ++i)
            {
                (&positive_divisor->value)[i - 1] = (&positive_divisor->value)[i] << 32 - shift |
                    (&positive_divisor->value)[i - 1] >> shift;
            }
            (&positive_divisor->value)[positive_divisor->value_count - 1] =
                (&positive_divisor->value)[positive_divisor->value_count - 1] >> shift;
            quotient_digit = 1 << 32 - shift;
            quotient_value_index -= 1;
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

struct Integer*integer_exponentiate(uint32_t*(integer_slot_allocator)(void*, size_t), void*memory,
    void*output_stack_cursor, void*scratch_stack_cursor, struct Integer*base,
    struct Integer*exponent)
{
    struct Integer*exponentiation =
        integer_new(stack_integer_slot_new, &scratch_stack_cursor, 1, 1);
    struct Integer*base_to_a_power_of_two = base;
    struct Integer*remaining_exponent = exponent;
    struct Integer two = { 1, 1, 2 };
    while (remaining_exponent->sign > 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack_cursor, output_stack_cursor, &division,
            remaining_exponent, &two);
        if (division.remainder->sign > 0)
        {
            exponentiation = integer_multiply(&scratch_stack_cursor, output_stack_cursor,
                exponentiation, base_to_a_power_of_two);
        }
        base_to_a_power_of_two = integer_multiply(&scratch_stack_cursor, output_stack_cursor,
            base_to_a_power_of_two, base_to_a_power_of_two);
        remaining_exponent = division.quotient;
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
    while (!integer_equals(a, &zero))
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

struct Rational
{
    struct Integer*numerator;
    struct Integer*denominator;
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
        out->numerator->value_count = 0;
        out->numerator->sign = 0;
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
    struct NumberSlot*slot = *free_list;
    if (slot)
    {
        *free_list = slot->next;
    }
    else
    {
        slot = stack_slot_new(pool_cursor, sizeof(struct NumberSlot));
        slot->next = 0;
    }
    return slot;
}

void number_slot_free(struct NumberSlot**number_slot_free_list, struct NumberSlot*slot,
    struct IntegerPool*integer_pool)
{
    if (slot->number.operation == 'r')
    {
        if (slot->number.value.numerator->value_count)
        {
            pool_integer_free(integer_pool, slot->number.value.numerator);
        }
        pool_integer_free(integer_pool, slot->number.value.denominator);
    }
    slot->next = *number_slot_free_list;
    *number_slot_free_list = slot;
}

bool get_input(struct NumberSlot**number_slot_pool_cursor, struct NumberSlot**number_free_list,
    struct IntegerPool*integer_pool, void*stack_a_cursor, void*stack_b_cursor)
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
            slot->number.value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            slot->number.value.numerator =
                integer_new(stack_integer_slot_new, &stack_a_cursor, next_char - '0', 1);
            struct Integer ten = { 1, 1, 10 };
            next_char = getchar();
            while (isdigit(next_char))
            {
                slot->number.value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
                    slot->number.value.numerator, &ten);
                struct Integer digit = { 1, 1, next_char - '0' };
                slot->number.value.numerator =
                    integer_add(&stack_a_cursor, slot->number.value.numerator, &digit);
                next_char = getchar();
            }
            integer_move(pool_integer_slot_new, integer_pool, &slot->number.value.numerator);
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
                {
                }
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
    struct NumberSlot**number_slot_free_list, struct IntegerPool*integer_pool,
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
                number_slot_free_list, integer_pool, number_slot->next);
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
                number_slot_free(number_slot_free_list, number_slot, integer_pool);
                number_slot_free(number_slot_free_list, nested_number_slot, integer_pool);
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
            number_slot_free(number_slot_free_list, number_slot, integer_pool);
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
                number_slot->number.value.numerator =
                    integer_new(pool_integer_slot_new, integer_pool, 1, -1);
                number_slot->number.value.denominator =
                    integer_new(pool_integer_slot_new, integer_pool, 1, 1);
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
                number_slot_free(number_slot_free_list, number_slot->next, integer_pool);
                number_slot_free(number_slot_free_list, number_slot, integer_pool);
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
    struct NumberSlot**number_slot_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        if (b->operation == 'r')
        {
            out->operation = 'r';
            if (!rational_out_add(pool_integer_slot_new, integer_pool, stack_a_cursor,
                stack_b_cursor, &out->value, &a->value, &b->value))
            {
                return false;
            }
            number_slot_free(number_slot_free_list, a, integer_pool);
            number_slot_free(number_slot_free_list, b, integer_pool);
        }
    }
    else
    {
        return number_add(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, out, b, a);
    }
    return true;
}

bool number_multiply(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        if (b->operation == 'r')
        {
            out->operation = 'r';
            if (!rational_out_multiply(pool_integer_slot_new, integer_pool, stack_a_cursor,
                stack_b_cursor, &out->value, &a->value, &b->value))
            {
                return false;
            }
            number_slot_free(number_slot_free_list, a, integer_pool);
            number_slot_free(number_slot_free_list, b, integer_pool);
            return true;
        }
    }
    else
    {
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
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
        struct Integer*old_numerator = a->value.numerator;
        a->value.numerator = a->value.denominator;
        a->value.denominator = old_numerator;
    }
    else
    {
        ABORT("number_reciprocal called on Number type for which it is not yet implemented.");
    }
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

struct Factor
{
    struct Integer*value;
    struct Integer*multiplicity;
};

bool number_exponentiate(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*out, struct Number*base, struct Number*exponent)
{
    if (exponent->operation != 'r')
    {
        printf("The input expression contains an exponentiation whose exponent is not both real "
            "and rational; this program doesn't handle transcendental numbers.");
        return false;
    }
    if (base->operation == 'r')
    {
        if (base->value.numerator->value_count == 0)
        {
            number_slot_free(number_slot_free_list, base, integer_pool);
            number_slot_free(number_slot_free_list, exponent, integer_pool);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            return true;
        }
        if (exponent->value.numerator->sign < 0)
        {
            number_reciprocal(base);
            exponent->value.numerator->sign = 1;
            return number_exponentiate(number_slot_pool_cursor, number_slot_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, out, base, exponent);
        }
        if (integer_equals_one(base->value.denominator))
        {
            if (integer_equals_one(exponent->value.denominator))
            {
                out->operation = 'r';
                out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
                out->value.numerator = integer_exponentiate(pool_integer_slot_new, integer_pool,
                    stack_a_cursor, stack_b_cursor, base->value.numerator,
                    exponent->value.numerator);
                return true;
            }
            struct Integer*radicand = integer_exponentiate(stack_integer_slot_new, &stack_a_cursor,
                stack_a_cursor, stack_b_cursor, base->value.numerator, exponent->value.numerator);
            int8_t radicand_sign = radicand->sign;
            radicand->sign = 1;
            size_t factor_count = 0;
            struct Factor*factors = stack_b_cursor;
            struct Integer*prime = primes;
            struct Integer one = { 1, 1, 1 };
            while (integer_compare(stack_a_cursor, stack_b_cursor, prime, radicand) >= 0)
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
                            integer_add(&stack_a_cursor, factor->multiplicity, &one);
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
                number_slot_free(number_slot_free_list, base, integer_pool);
                if (integer_equals_one(reduced_degree.quotient))
                {
                    number_slot_free(number_slot_free_list, exponent, integer_pool);
                    out->operation = 'r';
                    out->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move(pool_integer_slot_new, integer_pool, &out->value.numerator);
                    return true;
                }
                pool_integer_free(integer_pool, exponent->value.denominator);
                exponent->value.denominator = reduced_degree.quotient;
                integer_move(pool_integer_slot_new, integer_pool, &exponent->value.denominator);
            }
            struct Number*number_coefficient =
                number_slot_new(number_slot_pool_cursor, number_slot_free_list);
            number_coefficient->operation = 'r';
            number_coefficient->value.numerator = coefficient;
            number_coefficient->value.denominator =
                integer_new(pool_integer_slot_new, integer_pool, 1, 1);
            struct Number*number_radicand =
                number_slot_new(number_slot_pool_cursor, number_slot_free_list);
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
            struct Number*surd = number_slot_new(number_slot_pool_cursor, number_slot_free_list);
            surd->operation = '^';
            surd->left = number_radicand;
            surd->right = exponent;
            return number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
                stack_a_cursor, stack_b_cursor, out, number_coefficient, surd);
        }
        struct Number*new_denominator =
            number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        new_denominator->operation = 'r';
        new_denominator->value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        new_denominator->value.numerator = integer_exponentiate(pool_integer_slot_new, integer_pool,
            stack_a_cursor, stack_b_cursor, base->value.denominator, exponent->value.numerator);
        struct Integer one = { 1, 1, 1 };
        struct Number*new_numerator_base =
            number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        new_numerator_base->operation = 'r';
        new_numerator_base->value.denominator =
            integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        new_numerator_base->value.numerator = integer_multiply(&stack_a_cursor, stack_b_cursor,
            base->value.numerator, integer_exponentiate(stack_integer_slot_new, &stack_a_cursor,
                stack_a_cursor, stack_b_cursor, base->value.denominator,
                integer_subtract(&stack_a_cursor, stack_b_cursor, exponent->value.denominator,
                    &one)));
        integer_move(pool_integer_slot_new, integer_pool, &new_numerator_base->value.numerator);
        struct Number*new_numerator =
            number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        number_exponentiate(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, new_numerator, new_numerator_base, exponent);
        number_reciprocal(new_denominator);
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, out, new_numerator, new_denominator);
    }
    return true;
}

bool evaluate_root(struct NumberSlot**number_slot_pool_cursor,
    struct NumberSlot**number_slot_free_list, struct IntegerPool*integer_pool, void*stack_a_cursor,
    void*stack_b_cursor, struct Number*a)
{
    switch (a->operation)
    {
    case '+':
        return number_add(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    case '-':
    {
        struct NumberSlot*negative =
            number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        negative->number.operation = 'r';
        negative->number.value.numerator = integer_new(pool_integer_slot_new, integer_pool, 1, -1);
        negative->number.value.denominator = integer_new(pool_integer_slot_new, integer_pool, 1, 1);
        struct NumberSlot*product = number_slot_new(number_slot_pool_cursor, number_slot_free_list);
        if (!number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, product, negative, a->right))
        {
            return false;
        }
        return number_add(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, product);
    }
    case '*':
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    case '/':
        number_reciprocal(a->right);
        return number_multiply(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    case '^':
        return number_exponentiate(number_slot_pool_cursor, number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor, a, a->left, a->right);
    }
}

bool evaluate(struct NumberSlot**number_slot_pool_cursor, struct NumberSlot**number_slot_free_list,
    struct IntegerPool*integer_pool, uint32_t*stack_a_cursor, uint32_t*stack_b_cursor,
    struct Number*a)
{
    if (a->operation == 'r')
    {
        return true;
    }
    if (!evaluate(number_slot_pool_cursor, number_slot_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, a->left))
    {
        return false;
    }
    if (!evaluate(number_slot_pool_cursor, number_slot_free_list, integer_pool, stack_a_cursor,
        stack_b_cursor, a->right))
    {
        return false;
    }
    if (!evaluate_root(number_slot_pool_cursor, number_slot_free_list, integer_pool, stack_a_cursor,
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
}

int main()
{
    init();
    struct NumberSlot*number_slot_pool = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    struct NumberSlot*number_slot_pool_cursor = number_slot_pool;
    struct NumberSlot*number_slot_free_list = 0;
    void*stack_a = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_a_cursor = stack_a;
    void*stack_b = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    void*stack_b_cursor = stack_b;
    while (true)
    {
        struct IntegerPool*integer_pool =
            VirtualAlloc(0, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
        integer_pool_new(integer_pool);
        if (get_input(&number_slot_pool_cursor, &number_slot_free_list, integer_pool,
            stack_a_cursor, stack_b_cursor))
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
                integer_pool, number_slot_pool);
            if (input && evaluate(&number_slot_pool_cursor, &number_slot_free_list, integer_pool,
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
        while (integer_pool->first_page)
        {
            struct IntegerPool*next_value_count_pool = integer_pool + 1;
            struct IntegerSlotPage*page = integer_pool->first_page;
            while (page)
            {
                struct IntegerSlotPage*next_page = page->next_page;
                VirtualFree(page, page_size, MEM_RELEASE);
                page = next_page;
            }
            integer_pool = next_value_count_pool;
        }
    }
    return 0;
}