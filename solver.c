#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#ifdef _DEBUG
#define ABORT(message) printf(message); abort();
#define ASSERT(condition, message) if (!(condition)) { printf(message); abort(); }
#else
#define ABORT(message)
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

enum MemoryType
{
    Stack,
    Pool
};

struct Stack
{
    enum MemoryType memory_type;
    void*cursor;
    void*end;
};

void stack_new(struct Stack*out)
{
    out->memory_type = Stack;
    out->cursor = VirtualAlloc(0, arena_size, MEM_RESERVE, PAGE_READWRITE);
    out->end = (void*)((size_t)out->cursor + arena_size);
}

void*stack_slot_new(struct Stack*stack, size_t slot_size)
{
    void*slot = stack->cursor;
    size_t boundary = next_page_boundary((size_t)slot);
    stack->cursor = (void*)((size_t)slot + slot_size);
    if (stack->cursor > stack->end)
    {
        printf("Stack ran out of memory.");
        abort();
    }
    if ((size_t)stack->cursor > boundary)
    {
        VirtualAlloc((LPVOID)boundary, (size_t)stack->cursor - boundary, MEM_COMMIT,
            PAGE_READWRITE);
    }
    return slot;
}

void rewind_stack_cursor(struct Stack*stack, void*cursor_target)
{
    size_t boundary = next_page_boundary((size_t)cursor_target);
    if ((size_t)stack->cursor > boundary)
    {
        VirtualFree((LPVOID)boundary, (size_t)stack->cursor - boundary, MEM_DECOMMIT);
    }
    stack->cursor = cursor_target;
}

struct Integer
{
    size_t value_count;
    int8_t sign;
    uint32_t value;//The first element of an array of length value_count.
};

struct Integer zero = { 0, 0, 0 };
#define INT(value, sign) (struct Integer){ 1, sign 1, value }

size_t integer_size(size_t value_count)
{
    return sizeof(struct Integer) + (value_count - 1) * sizeof(uint32_t);
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
    enum MemoryType memory_type;
    struct IntegerSlotPage*first_page;
    struct IntegerSlot*cursor;
    struct IntegerSlot*free_list;
};

struct IntegerPool*permanent_integer_pool;

void integer_pool_new(struct IntegerPool*out)
{
    out->memory_type = Pool;
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

void pool_integer_free(struct IntegerPool*pool, struct Integer*a)
{
    if (a->value_count)
    {
        struct IntegerPool*value_count_pool = get_value_count_pool(pool, a->value_count);
        struct IntegerSlot*slot = (struct IntegerSlot*)((size_t)a - sizeof(struct IntegerSlot*));
        slot->next_slot = value_count_pool->free_list;
        value_count_pool->free_list = slot;
    }
}

struct Integer*integer_slot_new(void*memory, size_t value_count)
{
    if (*(enum MemoryType*)memory == Stack)
    {
        return stack_slot_new(memory, integer_size(value_count));
    }
    if (!value_count)
    {
        return &zero;
    }
    struct IntegerPool*pool = memory;
    struct IntegerPool*value_count_pool = get_value_count_pool(pool, value_count);
    struct IntegerSlot*out = pool->free_list;
    if (out)
    {
        pool->free_list = out->next_slot;
    }
    else
    {
        out = pool->cursor;
        pool->cursor = (struct IntegerSlot*)((size_t)(pool->cursor + 1) +
            (value_count - 1) * sizeof(uint32_t));
        size_t boundary = next_page_boundary((size_t)out);
        if ((size_t)pool->cursor > boundary)
        {
            struct IntegerSlotPage*new_page =
                VirtualAlloc(0, page_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
            pool->cursor = &new_page->memory;
            ((struct IntegerSlotPage*)(boundary - page_size))->next_page = new_page;
        }
    }
    return &out->integer;
}

struct Integer*integer_new(void*memory, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_slot_new(memory, 1);
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

struct Integer*integer_copy(void*memory, struct Integer*a)
{
    struct Integer*copy = integer_slot_new(memory, a->value_count);
    memcpy(copy, a, integer_size(a->value_count));
    return copy;
}

void integer_move(void*memory, struct Integer**a)
{
    *a = integer_copy(memory, *a);
}

void integer_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack, struct Integer**a)
{
    struct Integer*pool_copy = *a;
    integer_move(output_stack, a);
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

struct Integer*integer_magnitude(void*memory, struct Integer*a)
{
    struct Integer*out = integer_copy(memory, a);
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
    for (size_t i = a->value_count; i > 0; --i)
    {
        if ((&a->value)[i - 1] != 0)
        {
            a->value_count = i;
            return;
        }
    }
    a->value_count = 0;
    a->sign = 0;
}

struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b)
{
    if (a->sign == 0)
    {
        return integer_copy(output_stack, b);
    }
    if (b->sign == 0)
    {
        return integer_copy(output_stack, a);
    }
    struct Integer*long_integer = a;
    struct Integer*short_integer = b;
    if (a->value_count < b->value_count)
    {
        long_integer = b;
        short_integer = a;
    }
    struct Integer*sum = stack_slot_new(output_stack, integer_size(long_integer->value_count + 1));
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
    rewind_stack_cursor(output_stack, &sum->value + sum->value_count);
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

struct Integer*integer_negative(void*memory, struct Integer*a)
{
    struct Integer*negative = integer_copy(memory, a);
    negative->sign *= -1;
    return negative;
}

struct Integer*integer_subtract(struct Stack*output_stack, struct Stack scratch_stack,
    struct Integer*minuend, struct Integer*subtrahend)
{
    return integer_add(output_stack, minuend, integer_negative(&scratch_stack, subtrahend));
}

int8_t integer_compare(struct Stack*stack_a, struct Stack*stack_b, struct Integer*a,
    struct Integer*b)
{
    return integer_subtract(stack_a, *stack_b, a, b)->sign;
}

struct Integer*integer_multiply(struct Stack*output_stack, struct Stack scratch_stack,
    struct Integer*a, struct Integer*b)
{
    struct Integer*product =
        stack_slot_new(output_stack, integer_size(a->value_count + b->value_count));
    product->value_count = 0;
    product->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)(&a->value)[i] * (&b->value)[j];
            size_t shift = i + j;
            struct Integer*integer_component =
                stack_slot_new(&scratch_stack, integer_size(shift + 2));
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
    rewind_stack_cursor(output_stack, &product->value + product->value_count);
    return product;
}

struct Division
{
    struct Integer*quotient;
    struct Integer*remainder;
};

void integer_upshift(struct Integer*a, uint8_t shift)
{
    for (int64_t i = a->value_count - 1; i > 0; --i)
    {
        (&a->value)[i] = (&a->value)[i] << shift | (&a->value)[i - 1] >> (32 - shift);
    }
    a->value = a->value << shift;
}

void integer_downshift(struct Integer*a, uint8_t shift)
{
    for (size_t i = 0; i < a->value_count - 1; ++i)
    {
        (&a->value)[i] = (&a->value)[i] >> shift | (&a->value)[i + 1] << (32 - shift);
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

void calculate_division_values(struct Stack stack_a, struct Stack*stack_b,
    struct Integer*divisor_magnitude, struct Division*division, size_t quotient_value_index,
    uint32_t quotient_digit)
{
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer*difference =
                integer_subtract(&stack_a, *stack_b, division->remainder, divisor_magnitude);
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
            integer_halve(divisor_magnitude);
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
    ABORT("Most significant Integer value was 0.");
}

void integer_euclidean_divide(struct Stack*output_stack, struct Stack scratch_stack,
    struct Division*out, struct Integer*dividend, struct Integer*divisor)
{
    int8_t quotient_sign = dividend->sign * divisor->sign;
    int dividend_leading_digit_place = leading_digit_place(dividend);
    int divisor_leading_digit_place = leading_digit_place(divisor);
    if (dividend->value_count > divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
            dividend_leading_digit_place >= divisor_leading_digit_place))
    {
        struct Integer*divisor_magnitude =
            stack_slot_new(&scratch_stack, integer_size(dividend->value_count));
        divisor_magnitude->value_count = dividend->value_count;
        divisor_magnitude->sign = 1;
        size_t quotient_value_index = dividend->value_count - divisor->value_count;
        memset(&divisor_magnitude->value, 0, quotient_value_index * sizeof(uint32_t));
        memcpy(&divisor_magnitude->value + quotient_value_index, &divisor->value,
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
        int8_t dividend_sign = dividend->sign;
        dividend->sign = 1;
        out->quotient = stack_slot_new(output_stack, integer_size(dividend->value_count));
        out->quotient->value_count = dividend->value_count;
        out->quotient->sign = quotient_sign;
        memset(&out->quotient->value, 0, out->quotient->value_count * sizeof(uint32_t));
        out->remainder = dividend;
        calculate_division_values(scratch_stack, output_stack, divisor_magnitude, out,
            quotient_value_index, quotient_digit);
        trim_leading_zeroes(out->quotient);
        rewind_stack_cursor(output_stack, &out->quotient->value + out->quotient->value_count);
        integer_move(output_stack, &out->remainder);
        if (out->remainder->sign != 0)
        {
            out->remainder->sign = dividend->sign;
        }
        dividend->sign = dividend_sign;
    }
    else
    {
        out->quotient = &zero;
        out->remainder = integer_copy(output_stack, dividend);
    }
}

struct Integer*integer_euclidean_quotient(struct Stack*output_stack, struct Stack*scratch_stack,
    struct Integer*dividend, struct Integer*divisor)
{
    struct Division division;
    integer_euclidean_divide(output_stack, *scratch_stack, &division, dividend, divisor);
    return division.quotient;
}

struct Integer*integer_euclidean_remainder(struct Stack*output_stack, struct Stack*scratch_stack,
    struct Integer*dividend, struct Integer*divisor)
{
    struct Division division;
    integer_euclidean_divide(output_stack, *scratch_stack, &division, dividend, divisor);
    return division.remainder;
}

struct Integer*integer_doubled(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out;
    if ((&a->value)[a->value_count - 1] & 0x80000000)
    {
        out = stack_slot_new(output_stack, integer_size(a->value_count + 1));
        out->value_count = a->value_count + 1;
        (&out->value)[a->value_count] = 0;
    }
    else
    {
        out = stack_slot_new(output_stack, integer_size(a->value_count));
        out->value_count = a->value_count;
    }
    out->sign = a->sign;
    memcpy(&out->value, &a->value, a->value_count * sizeof(uint32_t));
    integer_upshift(out, 1);
    return out;
}

struct Integer*integer_half(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy(output_stack, a);
    integer_halve(out);
    return out;
}

//Requires that out be initialized to the multiplicative identity. Leaves garbage allocations on
//output_stack along with the final value of out.
void exponentiate(void(multiply)(struct Stack*, struct Stack*, void*, void*, void*),
    struct Stack*output_stack, struct Stack*scratch_stack, void*out, void*base,
    struct Integer*exponent)
{
    void*base_to_a_power_of_two = base;
    while (exponent->sign > 0)
    {
        if (exponent->value & 1)
        {
            multiply(output_stack, scratch_stack, out, out, base_to_a_power_of_two);
        }
        multiply(output_stack, scratch_stack, base_to_a_power_of_two, base_to_a_power_of_two,
            base_to_a_power_of_two);
        exponent = integer_half(output_stack, exponent);
    }
}

void integer_multiply_for_exponentiate(struct Stack*output_stack, struct Stack*scratch_stack,
    struct Integer**out, struct Integer**a, struct Integer**b)
{
    *out = integer_multiply(output_stack, *scratch_stack, *a, *b);
}

struct Integer*integer_exponentiate(void*memory, struct Stack*output_stack,
    struct Stack scratch_stack, struct Integer*base, struct Integer*exponent)
{
    struct Integer*out = integer_new(&scratch_stack, 1, 1);
    exponentiate(integer_multiply_for_exponentiate, &scratch_stack, output_stack, &out, &base,
        exponent);
    integer_move(memory, &out);
    return out;
}

void integer_string(struct Stack*output_stack, struct Stack scratch_stack, struct Integer*a)
{
    if (a->sign == 0)
    {
        *(char*)stack_slot_new(output_stack, 1) = '0';
        return;
    }
    if (a->sign < 0)
    {
        *(char*)stack_slot_new(output_stack, 1) = '-';
    }
    char*buffer_start = stack_slot_new(output_stack, 10 * a->value_count + 1);
    char*next_char = output_stack->cursor;
    struct Integer*quotient = a;
    struct Integer power = { 1, 1, 10 };
    while (quotient->sign != 0)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack, *output_stack, &division, quotient, &power);
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
    size_t char_count = (char*)output_stack->cursor - next_char;
    memcpy(buffer_start, next_char, char_count);
    rewind_stack_cursor(output_stack, buffer_start + char_count);
}

struct Integer*get_gcd(struct Stack*output_stack, struct Stack scratch_stack, struct Integer*a,
    struct Integer*b)
{
    while (b->value_count > 0)
    {
        struct Integer*c = b;
        b = integer_euclidean_remainder(&scratch_stack, output_stack, a, b);
        a = c;
    }
    return integer_copy(output_stack, a);
}

struct ExtendedGCDInfo
{
    struct Integer*gcd;
    struct Integer*a_coefficient;
    struct Integer*b_coefficient;
    struct Integer*a_over_gcd;
    struct Integer*b_over_gcd;
};

void extended_gcd(struct Stack*output_stack, struct Stack scratch_stack, struct ExtendedGCDInfo*out,
    struct Integer*a, struct Integer*b)
{
    out->a_coefficient = &zero;
    out->b_coefficient = integer_new(&scratch_stack, 1, 1);
    out->a_over_gcd = &zero;
    out->b_over_gcd = out->b_coefficient;
    while (a->value_count)
    {
        struct Division division;
        integer_euclidean_divide(&scratch_stack, *output_stack, &division, b, a);
        struct Integer*m = integer_subtract(&scratch_stack, *output_stack, out->a_coefficient,
            integer_multiply(&scratch_stack, *output_stack, out->b_over_gcd, division.quotient));
        struct Integer*n = integer_subtract(&scratch_stack, *output_stack, out->b_coefficient,
            integer_multiply(&scratch_stack, *output_stack, out->a_over_gcd, division.quotient));
        b = a;
        a = division.remainder;
        out->a_coefficient = out->b_over_gcd;
        out->b_coefficient = out->a_over_gcd;
        out->b_over_gcd = m;
        out->a_over_gcd = n;
    }
    out->gcd = integer_copy(output_stack, b);
    integer_move(output_stack, &out->a_coefficient);
    integer_move(output_stack, &out->b_coefficient);
    integer_move(output_stack, &out->a_over_gcd);
    integer_move(output_stack, &out->b_over_gcd);
    out->b_over_gcd->sign *= -1;
}

struct Stack prime_stack;
struct Integer*primes;

struct Integer*next_prime_slot(struct Integer*prime)
{
    return (struct Integer*)((size_t)prime + integer_size(prime->value_count));
}

struct Integer*next_prime(struct Stack stack_a, struct Stack*stack_b, struct Integer*prime)
{
    struct Integer*next_slot = next_prime_slot(prime);
    size_t boundary = next_page_boundary((size_t)prime);
    if ((size_t)next_slot == boundary)
    {
        VirtualAlloc((LPVOID)boundary, page_size, MEM_COMMIT, PAGE_READWRITE);
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
            if (integer_euclidean_remainder(&stack_a, stack_b, next_slot,
                potential_factor)->value_count == 0)
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

#define RATIONAL(numerator, denominator, sign)\
(struct Rational){ &INT(numerator, sign), &INT(denominator, +) }

struct RationalInterval
{
    struct Rational min;
    struct Rational max;
};

bool rational_new(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Integer*numerator, struct Integer*denominator)
{
    if (denominator->sign == 0)
    {
        printf("Tried to divide by 0.");
        return false;
    }
    if (numerator->sign == 0)
    {
        out->numerator = &zero;
        out->denominator = integer_new(memory, 1, 1);
        return true;
    }
    struct ExtendedGCDInfo gcd_info;
    extended_gcd(&scratch_stack, *output_stack, &gcd_info, numerator, denominator);
    out->numerator = gcd_info.a_over_gcd;
    out->denominator = gcd_info.b_over_gcd;
    if (out->denominator->sign < 0)
    {
        out->numerator->sign *= -1;
        out->denominator->sign *= -1;
    }
    integer_move(memory, &out->numerator);
    integer_move(memory, &out->denominator);
    return true;
}

void rational_free(struct IntegerPool*integer_pool, struct Rational*a)
{
    pool_integer_free(integer_pool, a->numerator);
    pool_integer_free(integer_pool, a->denominator);
}

void rational_move(void*memory, struct Rational*a)
{
    integer_move(memory, &a->numerator);
    integer_move(memory, &a->denominator);
}

void rational_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack, struct Rational*a)
{
    integer_move_from_pool(pool, output_stack, &a->numerator);
    integer_move_from_pool(pool, output_stack, &a->denominator);
}

void rational_magnitude(void*memory, struct Rational*out, struct Rational*a)
{
    out->numerator = integer_magnitude(memory, a->numerator);
    out->denominator = a->denominator;
}

bool rational_add(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*a, struct Rational*b)
{
    return rational_new(memory, output_stack, scratch_stack, out,
        integer_add(&scratch_stack,
            integer_multiply(&scratch_stack, *output_stack, a->numerator, b->denominator),
            integer_multiply(&scratch_stack, *output_stack, b->numerator, a->denominator)),
        integer_multiply(&scratch_stack, *output_stack, a->denominator, b->denominator));
}

void rational_integer_add(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*a, struct Integer*b)
{
    out->numerator = integer_add(memory, a->numerator,
        integer_multiply(&scratch_stack, *output_stack, b, a->denominator));
    out->denominator = integer_copy(memory, a->denominator);
}

bool rational_subtract(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*minuend, struct Rational*subtrahend)
{
    return rational_new(memory, output_stack, scratch_stack, out,
        integer_subtract(&scratch_stack, *output_stack,
            integer_multiply(&scratch_stack, *output_stack, minuend->numerator,
                subtrahend->denominator),
            integer_multiply(&scratch_stack, *output_stack, subtrahend->numerator,
                minuend->denominator)),
        integer_multiply(&scratch_stack, *output_stack, minuend->denominator,
            subtrahend->denominator));
}

bool rational_multiply(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*a, struct Rational*b)
{
    return rational_new(memory, output_stack, scratch_stack, out,
        integer_multiply(&scratch_stack, *output_stack, a->numerator, b->numerator),
        integer_multiply(&scratch_stack, *output_stack, a->denominator, b->denominator));
}

void rational_doubled(struct Stack*output_stack, struct Stack scratch_stack, struct Rational*out,
    struct Rational*a)
{
    rational_new(output_stack, output_stack, scratch_stack, out,
        integer_doubled(&scratch_stack, a->numerator), a->denominator);
}

bool rational_integer_multiply(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*a, struct Integer*b)
{
    return rational_new(memory, output_stack, scratch_stack, out,
        integer_multiply(&scratch_stack, *output_stack, a->numerator, b), a->denominator);
}

void rational_unreduced_multiply(struct Stack*output_stack, struct Stack*scratch_stack,
    struct Rational*out, struct Rational*a, struct Rational*b)
{
    out->numerator = integer_multiply(output_stack, *scratch_stack, a->numerator, b->numerator);
    out->denominator =
        integer_multiply(output_stack, *scratch_stack, a->denominator, b->denominator);
}

bool rational_divide(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*dividend, struct Rational*divisor)
{
    return rational_new(memory, output_stack, scratch_stack, out,
        integer_multiply(&scratch_stack, *output_stack, dividend->numerator, divisor->denominator),
        integer_multiply(&scratch_stack, *output_stack, dividend->denominator, divisor->numerator));
}

bool rational_integer_divide(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*dividend, struct Integer*divisor)
{
    return rational_new(memory, output_stack, scratch_stack, out, dividend->numerator,
        integer_multiply(&scratch_stack, *output_stack, dividend->denominator, divisor));
}

void rational_exponentiate(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Rational*out, struct Rational*base, struct Integer*exponent)
{
    out->numerator = integer_new(&scratch_stack, 1, 1);
    out->denominator = integer_new(&scratch_stack, 1, 1);
    struct Rational base_copy = *base;
    exponentiate(rational_unreduced_multiply, &scratch_stack, output_stack, out, &base_copy,
        exponent);
    rational_move(memory, out);
}

int8_t rational_compare(struct Stack stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    return integer_compare(&stack_a, stack_b,
        integer_multiply(&stack_a, *stack_b, a->numerator, b->denominator),
        integer_multiply(&stack_a, *stack_b, a->denominator, b->numerator));
}

struct Rational*rational_min(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    if (rational_compare(*stack_a, stack_b, a, b) < 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

struct Rational*rational_max(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    if (rational_compare(*stack_a, stack_b, a, b) > 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

void rational_estimate_trig_function(struct Stack*output_stack, struct Stack scratch_stack,
    struct RationalInterval*out, struct Rational*a_squared, struct Integer*factorial_component,
    struct Rational*delta, struct Rational*interval_size)
{
    while (true)
    {
        struct Rational delta_magnitude;
        rational_magnitude(&scratch_stack, &delta_magnitude, delta);
        if (rational_compare(*output_stack, &scratch_stack, &delta_magnitude, interval_size) <= 0)
        {
            break;
        }
        rational_add(&scratch_stack, &scratch_stack, *output_stack, &out->min, &out->min, delta);
        factorial_component = integer_add(&scratch_stack, factorial_component, &INT(1, +));
        struct Rational delta_factor;
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack, &delta_factor,
            a_squared, factorial_component);
        rational_multiply(&scratch_stack, &scratch_stack, *output_stack, delta, delta,
            &delta_factor);
        factorial_component = integer_add(&scratch_stack, factorial_component, &INT(1, +));
        struct Integer*negative_factorial_component =
            integer_negative(&scratch_stack, factorial_component);
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack, delta, delta,
            negative_factorial_component);
    }
    rational_move(output_stack, &out->min);
    if (delta->numerator->value_count > 0)
    {
        rational_add(output_stack, output_stack, scratch_stack, &out->max, &out->min, delta);
    }
    else
    {
        out->max = out->min;
        rational_add(output_stack, output_stack, scratch_stack, &out->min, &out->max, delta);
    }
}

void rational_estimate_cosine(struct Stack*output_stack, struct Stack scratch_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    out->min.numerator = &INT(1, +);
    out->min.denominator = &INT(1, +);
    struct Rational a_squared;
    rational_multiply(&scratch_stack, &scratch_stack, *output_stack, &a_squared, a, a);
    struct Integer*factorial_component = &INT(2, +);
    struct Rational delta;
    rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack, &delta, &a_squared,
        &INT(2, -));
    rational_estimate_trig_function(output_stack, scratch_stack, out, &a_squared,
        factorial_component, &delta, interval_size);
}

void rational_estimate_sine(struct Stack*output_stack, struct Stack scratch_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    out->min = *a;
    struct Rational a_squared;
    rational_multiply(&scratch_stack, &scratch_stack, *output_stack, &a_squared, a, a);
    struct Integer*factorial_component = &INT(3, +);
    struct Rational a_cubed;
    rational_multiply(&scratch_stack, &scratch_stack, *output_stack, &a_cubed, &a_squared, a);
    struct Rational delta;
    rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack, &delta, &a_cubed,
        &INT(6, -));
    rational_estimate_trig_function(output_stack, scratch_stack, out, &a_squared,
        factorial_component, &delta, interval_size);
}

struct Rational pi_estimate_min;
struct Rational pi_interval_size;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

void pi_refine_interval(struct Stack stack_a, struct Stack*stack_b)
{
    struct Rational term_a = { &INT(4, +), integer_add(&stack_a, pi_eight_k, &INT(1, +)) };
    struct Rational term_b;
    rational_new(&stack_a, &stack_a, *stack_b, &term_b, &INT(2, +),
        integer_add(&stack_a, pi_eight_k, &INT(4, +)));
    struct Rational sum;
    rational_subtract(&stack_a, &stack_a, *stack_b, &sum, &term_a, &term_b);
    struct Rational term_c = { &INT(1, +), integer_add(&stack_a, pi_eight_k, &INT(5, +)) };
    rational_subtract(&stack_a, &stack_a, *stack_b, &sum, &sum, &term_c);
    struct Rational term_d = { &INT(1, +), integer_add(&stack_a, pi_eight_k, &INT(6, +)) };
    rational_subtract(&stack_a, &stack_a, *stack_b, &sum, &sum, &term_d);
    struct Rational term;
    rational_integer_divide(&stack_a, &stack_a, *stack_b, &term, &sum, pi_sixteen_to_the_k);
    rational_add(&stack_a, &stack_a, *stack_b, &pi_estimate_min, &pi_estimate_min, &term);
    rational_multiply(&stack_a, &stack_a, *stack_b, &pi_interval_size, &pi_interval_size,
        &RATIONAL(1, 16, +));
    pi_eight_k = integer_add(&stack_a, pi_eight_k, &INT(8, +));
    pi_sixteen_to_the_k = integer_multiply(&stack_a, *stack_b, pi_sixteen_to_the_k, &INT(16, +));
}

void pi_estimate(struct Stack stack_a, struct Stack*stack_b, struct Rational*interval_size)
{
    int8_t interval_size_comparison =
        rational_compare(stack_a, stack_b, &pi_interval_size, interval_size);
    if (interval_size_comparison > 0)
    {
        rational_move_from_pool(permanent_integer_pool, &stack_a, &pi_estimate_min);
        rational_move_from_pool(permanent_integer_pool, &stack_a, &pi_interval_size);
        integer_move_from_pool(permanent_integer_pool, &stack_a, &pi_sixteen_to_the_k);
        integer_move_from_pool(permanent_integer_pool, &stack_a, &pi_eight_k);
        while (interval_size_comparison > 0)
        {
            pi_refine_interval(stack_a, stack_b);
            interval_size_comparison =
                rational_compare(stack_a, stack_b, &pi_interval_size, interval_size);
        }
        rational_move(permanent_integer_pool, &pi_estimate_min);
        rational_move(permanent_integer_pool, &pi_interval_size);
        integer_move(permanent_integer_pool, &pi_sixteen_to_the_k);
        integer_move(permanent_integer_pool, &pi_eight_k);
    }
}

void pi_shrink_interval_to_one_side_of_value(struct Stack stack_a, struct Stack*stack_b,
    struct Rational*value)
{
    struct Rational pi_estimate_max;
    rational_add(&stack_a, &stack_a, *stack_b, &pi_estimate_max, &pi_estimate_min,
        &pi_interval_size);
    if (rational_compare(stack_a, stack_b, &pi_estimate_min, value) <= 0 &&
        rational_compare(stack_a, stack_b, value, &pi_estimate_max) <= 0)
    {
        rational_move_from_pool(permanent_integer_pool, &stack_a, &pi_estimate_min);
        rational_move_from_pool(permanent_integer_pool, &stack_a, &pi_interval_size);
        integer_move_from_pool(permanent_integer_pool, &stack_a, &pi_sixteen_to_the_k);
        integer_move_from_pool(permanent_integer_pool, &stack_a, &pi_eight_k);
        while (true)
        {
            pi_refine_interval(stack_a, stack_b);
            rational_add(&stack_a, &stack_a, *stack_b, &pi_estimate_max, &pi_estimate_min,
                &pi_interval_size);
            if (rational_compare(stack_a, stack_b, &pi_estimate_min, value) > 0 ||
                rational_compare(stack_a, stack_b, value, &pi_estimate_max) > 0)
            {
                break;
            }
        }
        rational_move(permanent_integer_pool, &pi_estimate_min);
        rational_move(permanent_integer_pool, &pi_interval_size);
        integer_move(permanent_integer_pool, &pi_sixteen_to_the_k);
        integer_move(permanent_integer_pool, &pi_eight_k);
    }
}

struct Float
{
    struct Integer*significand;
    struct Integer*exponent;
};

void float_new(void*memory, struct Stack*output_stack, struct Stack scratch_stack, struct Float*out,
    struct Integer*significand, struct Integer*exponent)
{
    while (!(significand->value & 1))
    {
        significand = integer_half(&scratch_stack, significand);
        exponent = integer_add(&scratch_stack, exponent, &INT(1, +));
    }
    integer_move(memory, &significand);
    integer_move(memory, &exponent);
    out->significand = significand;
    out->exponent = exponent;
}

void float_move(void*memory, struct Float*a)
{
    integer_move(memory, &a->significand);
    integer_move(memory, &a->exponent);
}

void float_move_from_pool(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Float*a)
{
    integer_move_from_pool(integer_pool, output_stack, &a->significand);
    integer_move_from_pool(integer_pool, output_stack, &a->exponent);
}

void float_free(struct IntegerPool*integer_pool, struct Float*a)
{
    pool_integer_free(integer_pool, a->significand);
    pool_integer_free(integer_pool, a->exponent);
}

void float_add(void*memory, struct Stack*output_stack, struct Stack scratch_stack, struct Float*out,
    struct Float*a, struct Float*b)
{
    struct Integer*exponent_difference =
        integer_subtract(&scratch_stack, *output_stack, a->exponent, b->exponent);
    if (exponent_difference->sign > 0)
    {
        out->significand = integer_add(output_stack, b->significand,
            integer_multiply(&scratch_stack, *output_stack, a->significand,
                integer_exponentiate(&scratch_stack, &scratch_stack, *output_stack, &INT(2, +),
                    exponent_difference)));
        out->exponent = b->exponent;
    }
    else if (exponent_difference->sign < 0)
    {
        float_add(memory, output_stack, scratch_stack, out, b, a);
    }
    else
    {
        float_new(memory, output_stack, scratch_stack, out,
            integer_add(&scratch_stack, a->significand, b->significand), a->exponent);
    }
}

void float_negative(struct Stack*output_stack, struct Stack scratch_stack, struct Float*out,
    struct Float*a)
{
    out->significand = integer_negative(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
}

void float_subtract(struct Stack*output_stack, struct Stack scratch_stack, struct Float*out,
    struct Float*minuend, struct Float*subtrahend)
{
    struct Float negative;
    float_negative(&scratch_stack, *output_stack, &negative, subtrahend);
    float_add(output_stack, output_stack, scratch_stack, out, minuend, &negative);
}

void float_multiply(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Float*out, struct Float*a, struct Float*b)
{
    float_new(memory, output_stack, scratch_stack, out,
        integer_multiply(&scratch_stack, *output_stack, a->significand, b->significand),
        integer_add(&scratch_stack, a->exponent, b->exponent));
}

void float_multiply_for_exponentiate(struct Stack*output_stack, struct Stack*scratch_stack,
    struct Float*out, struct Float*a, struct Float*b)
{
    float_multiply(output_stack, output_stack, *scratch_stack, out, a, b);
}

void float_exponentiate(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct Float*out, struct Float*base, struct Integer*exponent)
{
    out->significand = integer_new(&scratch_stack, 1, 1);
    out->exponent = &zero;
    struct Float base_copy = *base;
    exponentiate(float_multiply_for_exponentiate, &scratch_stack, output_stack, out, &base_copy,
        exponent);
    float_move(memory, out);
}

void float_to_rational(struct Stack*output_stack, struct Stack scratch_stack, struct Rational*out,
    struct Float*a)
{
    if (a->exponent->sign < 0)
    {
        out->numerator = integer_copy(output_stack, a->significand);
        struct Integer*exponent_magnitude = integer_magnitude(&scratch_stack, a->exponent);
        out->denominator = integer_exponentiate(output_stack, output_stack, scratch_stack,
            &INT(2, +), exponent_magnitude);
    }
    else
    {
        out->numerator = integer_multiply(output_stack, scratch_stack,
            integer_exponentiate(&scratch_stack, &scratch_stack, *output_stack, &INT(2, +),
                a->exponent), a->significand);
        out->denominator = integer_new(output_stack, 1, 1);
    }
}

struct FloatInterval
{
    struct Float min;
    struct Float max;
};

void float_interval_move(void*memory, struct FloatInterval*a)
{
    if (a->min.significand)
    {
        float_move(memory, &a->min);
        float_move(memory, &a->max);
    }
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*scratch_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    float_to_rational(output_stack, *scratch_stack, &out->min, &a->min);
    float_to_rational(output_stack, *scratch_stack, &out->max, &a->max);
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
            struct FloatInterval real_part_estimate;
            struct FloatInterval imaginary_part_estimate;
        };
    };
    struct FloatInterval argument_estimate;
    struct FloatInterval magnitude_estimate;
    struct Number*previous;
    struct Number*next;
    char operation;
};

struct Number*number_slot_new(struct Stack*number_stack, struct Number**free_list)
{
    struct Number*out = *free_list;
    if (out)
    {
        *free_list = out->next;
    }
    else
    {
        out = stack_slot_new(number_stack, sizeof(struct Number));
    }
    memset(out, 0, sizeof(struct Number));
    return out;
}

struct Number*number_copy(struct Stack*number_stack, struct Number**free_list,
    struct IntegerPool*integer_pool, struct Number*a)
{
    struct Number*out = number_slot_new(number_stack, free_list);
    memcpy(out, a, sizeof(struct Number));
    if (a->operation == 'r')
    {
        rational_move(integer_pool, &out->value);
    }
    else
    {
        out->left = number_copy(number_stack, free_list, integer_pool, a->left);
        out->right = number_copy(number_stack, free_list, integer_pool, a->right);
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

bool get_input(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack stack_a, struct Stack*stack_b)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    struct Number*previous = 0;
    while (true)
    {
        struct Number*number = number_slot_new(number_stack, number_free_list);
        number->previous = previous;
        number->next = number_stack->cursor;
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.denominator = integer_new(integer_pool, 1, 1);
            number->value.numerator = integer_new(&stack_a, next_char - '0', 1);
            struct Integer ten = { 1, 1, 10 };
            next_char = getchar();
            while (isdigit(next_char))
            {
                number->value.numerator =
                    integer_multiply(&stack_a, *stack_b, number->value.numerator, &ten);
                struct Integer digit = { 1, 1, next_char - '0' };
                number->value.numerator = integer_add(&stack_a, number->value.numerator, &digit);
                next_char = getchar();
            }
            integer_move(integer_pool, &number->value.numerator);
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

struct Number*parse_input(struct Stack*number_stack, struct Number**number_free_list,
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
                parse_input(number_stack, number_free_list, integer_pool, input->next);
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
                input->value.numerator = integer_new(integer_pool, 1, -1);
                input->value.denominator = integer_new(integer_pool, 1, 1);
                struct Number*times = number_slot_new(number_stack, number_free_list);
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

bool estimate_needs_refinement(struct Stack stack_a, struct Stack*stack_b,
    struct FloatInterval*estimate, struct Rational*interval_size)
{
    if (estimate->min.significand)
    {
        struct Float current_interval_size;
        float_subtract(&stack_a, *stack_b, &current_interval_size, &estimate->max, &estimate->min);
        struct Rational rational_current_interval_size;
        float_to_rational(&stack_a, *stack_b, &rational_current_interval_size,
            &current_interval_size);
        if (rational_compare(stack_a, stack_b, interval_size, &rational_current_interval_size) >= 0)
        {
            return false;
        }
    }
    return true;
}

//Leaves garbage allocations on scratch_stack along with the final values of estimate_denominator
//and estimate_remainder. The value of out goes on memory.
void rational_continue_float_estimate(void*memory, struct Stack*output_stack,
    struct Stack*scratch_stack, struct Float*out_min, struct Float*out_max,
    struct Integer**estimate_denominator, struct Integer**estimate_remainder, struct Rational*a,
    struct Rational*interval_size)
{
    while (integer_compare(output_stack, scratch_stack, interval_size->denominator,
        integer_multiply(scratch_stack, *output_stack, *estimate_denominator,
            interval_size->numerator)) > 0)
    {
        if ((*estimate_remainder)->value_count == 0)
        {
            float_move(memory, out_min);
            *out_max = *out_min;
            return;
        }
        struct Division division;
        integer_euclidean_divide(scratch_stack, *output_stack, &division,
            integer_doubled(scratch_stack, *estimate_remainder), a->denominator);
        *estimate_denominator = integer_doubled(scratch_stack, *estimate_denominator);
        *estimate_remainder = division.remainder;
        out_min->significand = integer_add(scratch_stack, division.quotient,
            integer_doubled(scratch_stack, out_min->significand));
        out_min->exponent = integer_add(scratch_stack, out_min->exponent, &INT(1, -));
    }
    float_new(memory, output_stack, *scratch_stack, out_min, out_min->significand,
        out_min->exponent);
    float_new(memory, output_stack, *scratch_stack, out_max,
        integer_add(output_stack, out_min->significand, &INT(1, +)), out_min->exponent);
}

void rational_float_estimate(void*memory, struct Stack*output_stack, struct Stack scratch_stack,
    struct FloatInterval*out, struct Rational*a, struct Rational*interval_size)
{
    struct Integer*numerator_magnitude = integer_magnitude(&scratch_stack, a->numerator);
    struct Division division;
    integer_euclidean_divide(&scratch_stack, *output_stack, &division, numerator_magnitude,
        a->denominator);
    struct Integer*estimate_remainder = division.remainder;
    struct Integer*estimate_denominator = integer_new(&scratch_stack, 1, 1);
    if (a->numerator->sign < 0)
    {
        out->max.significand = division.quotient;
        out->max.exponent = &zero;
        rational_continue_float_estimate(memory, output_stack, &scratch_stack, &out->max, &out->min,
            &estimate_denominator, &estimate_remainder, a, interval_size);
        out->min.significand->sign = -1;
        out->max.significand->sign = -1;
    }
    else
    {
        out->min.significand = division.quotient;
        out->min.exponent = &zero;
        rational_continue_float_estimate(memory, output_stack, &scratch_stack, &out->min, &out->max,
            &estimate_denominator, &estimate_remainder, a, interval_size);
    }
}

void float_estimate_root(struct Stack*output_stack, struct Stack scratch_stack,
    struct FloatInterval*out, struct Float*a, struct Rational*interval_size, struct Integer*index)
{
    ASSERT(a->significand->sign >= 0, "float_estimate_root was called on a negative a value.");
    if (a->exponent < 0)
    {
        out->max.significand = integer_new(&scratch_stack, 1, 1);
        out->max.exponent = &zero;
    }
    else
    {
        out->max = *a;
        float_move(&scratch_stack, &out->max);
    }
    struct Rational rational_radicand;
    float_to_rational(&scratch_stack, *output_stack, &rational_radicand, a);
    struct Integer*index_minus_one = integer_add(&scratch_stack, index, &INT(1, -));
    while (true)
    {
        struct Float float_exponentiation;
        float_exponentiate(&scratch_stack, &scratch_stack, *output_stack, &float_exponentiation,
            &out->max, index_minus_one);
        struct Rational exponentiation;
        float_to_rational(&scratch_stack, *output_stack, &exponentiation, &float_exponentiation);
        struct Rational delta_numerator_minuend;
        rational_divide(&scratch_stack, &scratch_stack, *output_stack, &delta_numerator_minuend,
            &rational_radicand, &exponentiation);
        struct Rational rational_out;
        float_to_rational(&scratch_stack, *output_stack, &rational_out, &out->max);
        struct Rational delta_numerator;
        rational_subtract(&scratch_stack, &scratch_stack, *output_stack, &delta_numerator,
            &delta_numerator_minuend, &rational_out);
        struct Rational delta;
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack, &delta,
            &delta_numerator, index);
        struct Rational delta_float_estimate_interval_size;
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack,
            &delta_float_estimate_interval_size, &delta, &INT(2, -));
        struct FloatInterval delta_float_estimate;
        rational_float_estimate(&scratch_stack, &scratch_stack, *output_stack,
            &delta_float_estimate, &delta, &delta_float_estimate_interval_size);
        float_add(&scratch_stack, &scratch_stack, *output_stack, &out->max, &out->max,
            &delta_float_estimate.max);
        struct Rational delta_estimate_max;
        float_to_rational(&scratch_stack, *output_stack, &delta_estimate_max,
            &delta_float_estimate.max);
        struct Rational twice_delta;
        rational_new(&scratch_stack, &scratch_stack, *output_stack, &twice_delta,
            integer_doubled(&scratch_stack, delta.numerator), delta.denominator);
        struct Rational interval_size_comparison_value;
        rational_subtract(&scratch_stack, &scratch_stack, *output_stack,
            &interval_size_comparison_value, &delta_estimate_max, &twice_delta);
        if (rational_compare(*output_stack, &scratch_stack, &interval_size_comparison_value,
            interval_size) <= 0)
        {
            float_add(&scratch_stack, &scratch_stack, *output_stack, &out->min, &out->max,
                &delta_float_estimate.max);
            return;
        }
    }
}

struct FloatInterval*number_float_magnitude_estimate(struct Stack stack_a, struct Stack*stack_b,
    struct IntegerPool*integer_pool, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(stack_a, stack_b, &a->magnitude_estimate, interval_size))
    {
        return &a->magnitude_estimate;
    }
    if (a->operation == 'r')
    {
        if (!a->magnitude_estimate.min.significand)
        {
            struct Division division;
            integer_euclidean_divide(&stack_a, *stack_b, &division,
                integer_magnitude(&stack_a, a->value.numerator), a->value.denominator);
            a->magnitude_estimate_denominator = integer_new(&stack_a, 1, 1);
            float_new(&stack_a, &stack_a, *stack_b, &a->magnitude_estimate.min, division.quotient,
                &zero);
            a->magnitude_estimate_remainder = division.remainder;
        }
        else
        {
            float_move_from_pool(integer_pool, &stack_a, &a->magnitude_estimate.min);
            float_move_from_pool(integer_pool, &stack_a, &a->magnitude_estimate.max);
            integer_move_from_pool(integer_pool, &stack_a, &a->magnitude_estimate_denominator);
            integer_move_from_pool(integer_pool, &stack_a, &a->magnitude_estimate_remainder);
        }
        rational_continue_float_estimate(integer_pool, stack_b, &stack_a,
            &a->magnitude_estimate.min, &a->magnitude_estimate.max,
            &a->magnitude_estimate_denominator, &a->magnitude_estimate_remainder, &a->value,
            interval_size);
        integer_move(integer_pool, &a->magnitude_estimate_denominator);
        integer_move(integer_pool, &a->magnitude_estimate_remainder);
        return &a->magnitude_estimate;
    }
    if (a->magnitude_estimate.min.significand)
    {
        float_free(integer_pool, &a->magnitude_estimate.min);
        float_free(integer_pool, &a->magnitude_estimate.max);
    }
    if (a->operation == '^')
    {
        struct Rational radicand_magnitude_interval_size;
        rational_exponentiate(&stack_a, &stack_a, *stack_b, &radicand_magnitude_interval_size,
            interval_size, a->right->value.denominator);
        rational_integer_divide(&stack_a, &stack_a, *stack_b, &radicand_magnitude_interval_size,
            &radicand_magnitude_interval_size, &INT(3, +));
        struct FloatInterval*radicand_magnitude_estimate = number_float_magnitude_estimate(stack_a,
            stack_b, integer_pool, a->left, &radicand_magnitude_interval_size);
        struct Rational bound_interval_size;
        rational_integer_divide(&stack_a, &stack_a, *stack_b, &bound_interval_size, interval_size,
            &INT(3, +));
        float_estimate_root(&stack_a, *stack_b, &a->magnitude_estimate,
            &radicand_magnitude_estimate->min, interval_size, a->right->value.denominator);
        struct FloatInterval max_estimate;
        float_estimate_root(&stack_a, *stack_b, &max_estimate, &radicand_magnitude_estimate->max,
            interval_size, a->right->value.denominator);
        a->magnitude_estimate.max = max_estimate.max;
        float_interval_move(integer_pool, &a->magnitude_estimate);
        return &a->magnitude_estimate;
    }
    ABORT("number_float_magnitude_estimate case not yet implemented.");
}

void number_rational_magnitude_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    float_interval_to_rational_interval(output_stack, &scratch_stack, out,
        number_float_magnitude_estimate(scratch_stack, output_stack, integer_pool, a,
            interval_size));
}

void number_float_estimate_from_rational(void (get_rational_estimate)(struct Stack*, struct Stack,
    struct IntegerPool*, struct RationalInterval*, struct Number*, struct Rational*),
    struct Stack stack_a, struct Stack*stack_b, struct IntegerPool*integer_pool,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    struct Rational rational_estimate_interval_size = RATIONAL(1, 1, +);
    while (rational_compare(stack_a, stack_b, &rational_estimate_interval_size, interval_size) > 0)
    {
        rational_estimate_interval_size.denominator =
            integer_doubled(&stack_a, rational_estimate_interval_size.denominator);
    }
    rational_estimate_interval_size.denominator =
        integer_doubled(&stack_a, rational_estimate_interval_size.denominator);
    struct RationalInterval rational_estimate;
    get_rational_estimate(&stack_a, *stack_b, integer_pool, &rational_estimate, a,
        &rational_estimate_interval_size);
    rational_float_estimate(&stack_a, &stack_a, *stack_b, out, &rational_estimate.min,
        &rational_estimate_interval_size);
    struct FloatInterval max_estimate;
    rational_float_estimate(&stack_a, &stack_a, *stack_b, &max_estimate, &rational_estimate.max,
        &rational_estimate_interval_size);
    out->max = max_estimate.max;
    float_interval_move(integer_pool, out);
}

struct FloatInterval*number_float_argument_estimate(struct Stack*stack_a, struct Stack*stack_b,
    struct IntegerPool*integer_pool, struct Number*a, struct Rational*interval_size);

void number_rational_argument_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->sign < 0)
        {
            pi_estimate(*output_stack, &scratch_stack, interval_size);
            out->min = pi_estimate_min;
            rational_add(output_stack, output_stack, scratch_stack, &out->max, &pi_estimate_min,
                &pi_interval_size);
        }
        else
        {
            out->min.numerator = &zero;
            out->min.denominator = integer_new(output_stack, 1, 1);
            out->max.numerator = &zero;
            out->max.denominator = out->min.denominator;
        }
        return;
    case '^':
    {
        struct Rational radicand_argument_interval_size;
        rational_integer_multiply(&scratch_stack, &scratch_stack, *output_stack,
            &radicand_argument_interval_size, interval_size, a->right->value.denominator);
        struct RationalInterval radicand_rational_argument_estimate;
        number_rational_argument_estimate(&scratch_stack, *output_stack, integer_pool,
            &radicand_rational_argument_estimate, a->left, &radicand_argument_interval_size);
        rational_integer_divide(output_stack, output_stack, scratch_stack, &out->min,
            &radicand_rational_argument_estimate.min, a->right->value.denominator);
        rational_integer_divide(output_stack, output_stack, scratch_stack, &out->max,
            &radicand_rational_argument_estimate.max, a->right->value.denominator);
        return;
    }
    case '*':
    {
        float_interval_to_rational_interval(output_stack, &scratch_stack, out,
            number_float_argument_estimate(&scratch_stack, output_stack, integer_pool, a,
                interval_size));
    }
    case '+':
        ABORT("number_rational_argument_estimate case not yet implemented.");
    }
}

void number_estimate_part_sum(struct FloatInterval*(estimate_term)(void*, void*,
    struct IntegerPool*, struct Number*, struct Rational*), struct Stack stack_a,
    struct Stack*stack_b, struct IntegerPool*integer_pool, struct FloatInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    struct Rational term_estimate_interval_size;
    rational_integer_divide(&stack_a, &stack_a, *stack_b, &term_estimate_interval_size,
        interval_size, &INT(2, +));
    struct FloatInterval*left_term_estimate = estimate_term(&stack_a, &stack_a, integer_pool,
        a->left, &term_estimate_interval_size);
    struct FloatInterval*right_term_estimate = estimate_term(&stack_a, stack_b, integer_pool,
        a->right, &term_estimate_interval_size);
    float_add(integer_pool, &stack_a, *stack_b, &out->min, &left_term_estimate->min,
        &right_term_estimate->min);
    float_add(integer_pool, &stack_a, *stack_b, &out->max, &left_term_estimate->max,
        &right_term_estimate->max);
}

struct FloatInterval*number_float_argument_estimate(struct Stack*stack_a, struct Stack*stack_b,
    struct IntegerPool*integer_pool, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(*stack_a, stack_b, &a->argument_estimate, interval_size))
    {
        return &a->magnitude_estimate;
    }
    if (a->operation == '*')
    {
        number_estimate_part_sum(number_float_argument_estimate, *stack_a, stack_b, integer_pool,
            &a->argument_estimate, a, interval_size);
        return &a->argument_estimate;
    }
    else
    {
        number_float_estimate_from_rational(number_rational_argument_estimate, *stack_a, stack_b,
            integer_pool, &a->argument_estimate, a, interval_size);
        return &a->argument_estimate;
    }
}

void number_argument_cosine_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    struct Rational argument_estimate_interval_size;
    rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack,
        &argument_estimate_interval_size, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(&scratch_stack, *output_stack, integer_pool,
        &argument_estimate, a, &argument_estimate_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(&scratch_stack, *output_stack, &argument_min_cosine,
        &argument_estimate.min, &argument_estimate_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(&scratch_stack, *output_stack, &argument_max_cosine,
        &argument_estimate.max, &argument_estimate_interval_size);
    out->max = *rational_max(output_stack, &scratch_stack, &argument_min_cosine.max,
        &argument_max_cosine.max);
    rational_move(output_stack, &out->max);
    pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack, &argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack, &argument_estimate.max);
    struct Rational pi_estimate_max;
    rational_add(&scratch_stack, &scratch_stack, *output_stack, &pi_estimate_max, &pi_estimate_min,
        &pi_interval_size);
    if (rational_compare(*output_stack, &scratch_stack, &argument_estimate.min,
        &pi_estimate_min) <= 0 && rational_compare(*output_stack, &scratch_stack, &pi_estimate_max,
            &argument_estimate.max) <= 0)
    {
        out->min.numerator = integer_new(output_stack, 1, -1);
        out->min.denominator = integer_new(output_stack, 1, 1);
    }
    else
    {
        out->min = *rational_max(output_stack, &scratch_stack, &argument_min_cosine.min,
            &argument_max_cosine.min);
        rational_move(output_stack, &out->min);
    }
}

void number_argument_sine_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    struct Rational argument_estimate_interval_size;
    rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack,
        &argument_estimate_interval_size, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(&scratch_stack, *output_stack, integer_pool,
        &argument_estimate, a, &argument_estimate_interval_size);
    struct RationalInterval argument_min_sine;
    rational_estimate_sine(&scratch_stack, *output_stack, &argument_min_sine,
        &argument_estimate.min, &argument_estimate_interval_size);
    struct RationalInterval argument_max_sine;
    rational_estimate_sine(&scratch_stack, *output_stack, &argument_max_sine,
        &argument_estimate.max, &argument_estimate_interval_size);
    struct RationalInterval argument_estimate_multiple;
    rational_doubled(&scratch_stack, *output_stack, &argument_estimate_multiple.min,
        &argument_estimate.min);
    rational_doubled(&scratch_stack, *output_stack, &argument_estimate_multiple.max,
        &argument_estimate.max);
    pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack,
        &argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack,
        &argument_estimate_multiple.max);
    struct Rational pi_estimate_max;
    rational_add(&scratch_stack, &scratch_stack, *output_stack, &pi_estimate_max, &pi_estimate_min,
        &pi_interval_size);
    if (rational_compare(*output_stack, &scratch_stack, &argument_estimate_multiple.min,
        &pi_estimate_min) <= 0 && rational_compare(*output_stack, &scratch_stack, &pi_estimate_max,
            &argument_estimate_multiple.max) <= 0)
    {
        out->max.numerator = integer_new(output_stack, 1, 1);
        out->max.denominator = integer_new(output_stack, 1, 1);
        out->min = *rational_max(output_stack, &scratch_stack, &argument_min_sine.min,
            &argument_max_sine.min);
        rational_move(output_stack, &out->min);
    }
    else
    {
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack,
            &argument_estimate_multiple.min, &argument_estimate_multiple.min, &INT(3, +));
        rational_integer_divide(&scratch_stack, &scratch_stack, *output_stack,
            &argument_estimate_multiple.max, &argument_estimate_multiple.max, &INT(3, +));
        pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack,
            &argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(*output_stack, &scratch_stack,
            &argument_estimate_multiple.max);
        rational_add(&scratch_stack, &scratch_stack, *output_stack, &pi_estimate_max,
            &pi_estimate_min, &pi_interval_size);
        if (rational_compare(*output_stack, &scratch_stack, &argument_estimate_multiple.min,
            &pi_estimate_min) <= 0 && rational_compare(*output_stack, &scratch_stack,
                &pi_estimate_max, &argument_estimate_multiple.max) <= 0)
        {
            out->min.numerator = integer_new(output_stack, 1, -1);
            out->min.denominator = integer_new(output_stack, 1, 1);
            out->max = *rational_max(output_stack, &scratch_stack, &argument_min_sine.max,
                &argument_max_sine.max);
            rational_move(output_stack, &out->max);
        }
        else
        {
            out->min = *rational_max(output_stack, &scratch_stack, &argument_min_sine.min,
                &argument_max_sine.min);
            rational_move(output_stack, &out->min);
            out->max = *rational_max(output_stack, &scratch_stack, &argument_min_sine.max,
                &argument_max_sine.max);
            rational_move(output_stack, &out->max);
        }
    }
}

void number_rectangular_part_from_polar_form(void (trig_function)(struct Stack*, struct Stack,
    struct IntegerPool*, struct RationalInterval*, struct Number*, struct Rational*interval_size),
    struct Stack*output_stack, struct Stack scratch_stack, struct IntegerPool*integer_pool,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    struct RationalInterval magnitude_estimate;
    number_rational_magnitude_estimate(&scratch_stack, *output_stack, integer_pool,
        &magnitude_estimate, a, &RATIONAL(1, 1, +));
    struct Rational interval_scale;
    rational_integer_add(&scratch_stack, &scratch_stack, *output_stack, &interval_scale,
        &magnitude_estimate.max, &INT(2, +));
    struct Rational factor_interval_size;
    rational_divide(&scratch_stack, &scratch_stack, *output_stack, &factor_interval_size,
        interval_size, &interval_scale);
    number_rational_magnitude_estimate(&scratch_stack, *output_stack, integer_pool,
        &magnitude_estimate, a, &factor_interval_size);
    struct RationalInterval trig_value;
    trig_function(&scratch_stack, *output_stack, integer_pool, &trig_value, a,
        &factor_interval_size);
    if (trig_value.min.numerator->sign >= 0)
    {
        rational_multiply(output_stack, output_stack, scratch_stack, &out->min, &trig_value.min,
            &magnitude_estimate.min);
        rational_multiply(output_stack, output_stack, scratch_stack, &out->max, &trig_value.max,
            &magnitude_estimate.max);
    }
    else if (trig_value.max.numerator->sign <= 0)
    {
        rational_multiply(output_stack, output_stack, scratch_stack, &out->min, &trig_value.min,
            &magnitude_estimate.max);
        rational_multiply(output_stack, output_stack, scratch_stack, &out->max, &trig_value.max,
            &magnitude_estimate.min);
    }
    else
    {
        rational_multiply(output_stack, output_stack, scratch_stack, &out->min, &trig_value.min,
            &magnitude_estimate.max);
        rational_multiply(output_stack, output_stack, scratch_stack, &out->max, &trig_value.max,
            &magnitude_estimate.max);
    }
}

struct FloatInterval*number_float_real_part_estimate(struct Stack*output_stack,
    struct Stack*scratch_stack, struct IntegerPool*integer_pool, struct Number*a,
    struct Rational*interval_size);

void number_rational_real_part_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_cosine_estimate, output_stack,
            scratch_stack, integer_pool, out, a, interval_size);
    }
    else
    {
        struct FloatInterval*float_estimate = number_float_real_part_estimate(output_stack,
            &scratch_stack, integer_pool, a, interval_size);
        float_to_rational(output_stack, scratch_stack, &out->min, &float_estimate->min);
        float_to_rational(output_stack, scratch_stack, &out->max, &float_estimate->max);
    }
}

//When a takes the union variant with the real_part_estimate field, the return value is a pointer to
//that field, whose Integer components are allocated in integer_pool. Otherwise, the return value
//and its Integer components are allocated on output_stack.
struct FloatInterval*number_float_real_part_estimate(struct Stack*output_stack,
    struct Stack*scratch_stack, struct IntegerPool*integer_pool, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*magnitude_estimate = number_float_magnitude_estimate(*output_stack,
            scratch_stack, integer_pool, a, interval_size);
        if (a->value.numerator->sign < 0)
        {
            struct FloatInterval*real_part_estimate =
                stack_slot_new(output_stack, sizeof(struct FloatInterval));
            float_negative(output_stack, *scratch_stack, &real_part_estimate->min,
                &magnitude_estimate->max);
            float_negative(output_stack, *scratch_stack, &real_part_estimate->max,
                &magnitude_estimate->min);
            return real_part_estimate;
        }
        else
        {
            return magnitude_estimate;
        }
    }
    if (!estimate_needs_refinement(*output_stack, scratch_stack, &a->real_part_estimate,
        interval_size))
    {
        return &a->real_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
        number_float_estimate_from_rational(number_rational_real_part_estimate, *output_stack,
            scratch_stack, integer_pool, &a->real_part_estimate, a, interval_size);
        return &a->real_part_estimate;
    case '*':
        ABORT("number_float_real_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_real_part_estimate, *output_stack, scratch_stack,
            integer_pool, &a->real_part_estimate, a, interval_size);
        return &a->real_part_estimate;
    default:
        ABORT("Number operation not recognized.");
    }
}

struct FloatInterval*number_float_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*scratch_stack, struct IntegerPool*integer_pool, struct Number*a,
    struct Rational*interval_size);

void number_rational_imaginary_part_estimate(struct Stack*output_stack, struct Stack scratch_stack,
    struct IntegerPool*integer_pool, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_sine_estimate, output_stack,
            scratch_stack, integer_pool, out, a, interval_size);
    }
    else
    {
        struct FloatInterval*float_estimate = number_float_imaginary_part_estimate(output_stack,
            &scratch_stack, integer_pool, a, interval_size);
        float_to_rational(output_stack, scratch_stack, &out->min, &float_estimate->min);
        float_to_rational(output_stack, scratch_stack, &out->max, &float_estimate->max);
    }
}

//When a takes the union variant with the imaginary_part_estimate field, the return value is a
//pointer to that field, whose Integer components are allocated in integer_pool. Otherwise, the
//return value and its Integer components are allocated on output_stack.
struct FloatInterval*number_float_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*scratch_stack, struct IntegerPool*integer_pool, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*imaginary_part_estimate =
            stack_slot_new(output_stack, sizeof(struct FloatInterval));
        imaginary_part_estimate->min.significand = &zero;
        imaginary_part_estimate->min.exponent = &zero;
        imaginary_part_estimate->max.significand = &zero;
        imaginary_part_estimate->max.exponent = &zero;
        return imaginary_part_estimate;
    }
    if (!estimate_needs_refinement(*output_stack, scratch_stack, &a->imaginary_part_estimate,
        interval_size))
    {
        return &a->imaginary_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
        number_float_estimate_from_rational(number_rational_imaginary_part_estimate, *output_stack,
            scratch_stack, integer_pool, &a->imaginary_part_estimate, a, interval_size);
        return &a->imaginary_part_estimate;
    case '*':
        ABORT("number_float_imaginary_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_imaginary_part_estimate, *output_stack, scratch_stack,
            integer_pool, &a->imaginary_part_estimate, a, interval_size);
        return &a->imaginary_part_estimate;
    default:
        ABORT("Number operation not recognized.");
    }
}

struct Number*number_add(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            number_free(number_free_list, integer_pool, a);
            return b;
        }
        if (b->operation == 'r')
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = 'r';
            if (!rational_add(integer_pool, &stack_a, *stack_b, &out->value, &a->value, &b->value))
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
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
        break;
    case '*':
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
        break;
    case '+':
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '+':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    return number_add(number_stack, number_free_list, integer_pool, stack_a, stack_b, b, a);
}

struct Number*number_multiply(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = integer_new(integer_pool, 1, 1);
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
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = 'r';
            if (!rational_multiply(integer_pool, &stack_a, *stack_b, &out->value, &a->value,
                &b->value))
            {
                return 0;
            }
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            return out;
        }
        break;
    case '^':
    {
        switch (b->operation)
        {
        case 'r':
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '*';
            out->left = b;
            out->right = a;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return out;
        }
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case 'r':
        {
            if (a->left->operation == 'r')
            {
                struct Number*coefficient = number_multiply(number_stack, number_free_list,
                    integer_pool, stack_a, stack_b, b, a->left);
                if (!coefficient)
                {
                    return 0;
                }
                number_free(number_free_list, integer_pool, b);
                struct Number*out = number_multiply(number_stack, number_free_list, integer_pool,
                    stack_a, stack_b, coefficient, a->right);
                number_root_free(number_free_list, integer_pool, a);
                return out;
            }
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '*';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            struct Number*out = number_multiply(number_stack, number_free_list, integer_pool,
                stack_a, stack_b, a->right, b);
            if (!out)
            {
                return 0;
            }
            out = number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                a->left, out);
            number_root_free(number_free_list, integer_pool, a);
            return out;
        }
        }
        break;
    case '+':
    {
        struct Number*b_copy = number_copy(number_stack, number_free_list, integer_pool, b);
        if (b->operation == 'r')
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = '+';
            out->left = number_slot_new(number_stack, number_free_list);
            out->left = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
                stack_b, a->left, b);
            if (!out->left)
            {
                return 0;
            }
            out->right = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
                stack_b, a->right, b_copy);
            if (!out->right)
            {
                return 0;
            }
            number_root_free(number_free_list, integer_pool, a);
            return out;
        }
        struct Number*left = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
            stack_b, a->left, b);
        if (!left)
        {
            return 0;
        }
        struct Number*right = number_multiply(number_stack, number_free_list, integer_pool,
            stack_a, stack_b, a->right, b_copy);
        if (!right)
        {
            return 0;
        }
        number_root_free(number_free_list, integer_pool, a);
        return number_add(number_stack, number_free_list, integer_pool, stack_a, stack_b, left,
            right);
    }
    }
    return number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b, b, a);
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

struct Number*number_exponentiate(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack stack_a, struct Stack stack_b, struct Number*base,
    struct Number*exponent)
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
        return number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            base, exponent);
    }
    if (base->operation == 'r')
    {
        if (base->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            number_free(number_free_list, integer_pool, base);
            number_free(number_free_list, integer_pool, exponent);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = integer_new(integer_pool, 1, 1);
            return out;
        }
        if (integer_equals_one(base->value.denominator))
        {
            if (integer_equals_one(exponent->value.denominator))
            {
                struct Number*out = number_slot_new(number_stack, number_free_list);
                out->operation = 'r';
                out->value.denominator = integer_new(integer_pool, 1, 1);
                out->value.numerator = integer_exponentiate(integer_pool, &stack_a, stack_b,
                    base->value.numerator, exponent->value.numerator);
                number_free(number_free_list, integer_pool, base);
                number_free(number_free_list, integer_pool, exponent);
                return out;
            }
            struct Integer*radicand = integer_exponentiate(&stack_a, &stack_a, stack_b,
                base->value.numerator, exponent->value.numerator);
            number_free(number_free_list, integer_pool, base);
            int8_t radicand_sign = radicand->sign;
            radicand->sign = 1;
            size_t factor_count = 0;
            struct Factor*factors = stack_b.cursor;
            struct Integer*prime = primes;
            while (integer_compare(&stack_a, &stack_b, prime, radicand) <= 0)
            {
                struct Division division;
                integer_euclidean_divide(&stack_a, stack_b, &division, radicand, prime);
                if (division.remainder->value_count == 0)
                {
                    factor_count += 1;
                    struct Factor*factor = stack_slot_new(&stack_b, sizeof(struct Factor));
                    factor->value = prime;
                    factor->multiplicity = &zero;
                    while (division.remainder->value_count == 0)
                    {
                        factor->multiplicity =
                            integer_add(&stack_a, factor->multiplicity, &INT(1, +));
                        radicand = division.quotient;
                        integer_euclidean_divide(&stack_a, stack_b, &division,
                            radicand, factor->value);
                    }
                }
                prime = next_prime(stack_a, &stack_b, prime);
            }
            struct Integer*coefficient = integer_new(&stack_a, 1, 1);
            struct Integer*multiplicity_gcd = exponent->value.denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Division multiplicity_reduction;
                integer_euclidean_divide(&stack_a, stack_b, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->value.denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(&stack_a, stack_b, coefficient,
                    integer_exponentiate(&stack_a, &stack_a, stack_b, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = get_gcd(&stack_a, stack_b, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            if (radicand_sign > 0)
            {
                struct Integer*reduced_degree = integer_euclidean_quotient(&stack_a, &stack_b,
                    exponent->value.denominator, multiplicity_gcd);
                if (integer_equals_one(reduced_degree))
                {
                    number_free(number_free_list, integer_pool, exponent);
                    struct Number*out = number_slot_new(number_stack, number_free_list);
                    out->operation = 'r';
                    out->value.denominator = integer_new(integer_pool, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move(integer_pool, &out->value.numerator);
                    return out;
                }
                pool_integer_free(integer_pool, exponent->value.denominator);
                exponent->value.denominator = reduced_degree;
                integer_move(integer_pool, &exponent->value.denominator);
            }
            struct Number*number_coefficient = number_slot_new(number_stack, number_free_list);
            number_coefficient->operation = 'r';
            number_coefficient->value.numerator = coefficient;
            integer_move(integer_pool, &number_coefficient->value.numerator);
            number_coefficient->value.denominator = integer_new(integer_pool, 1, 1);
            struct Number*number_radicand = number_slot_new(number_stack, number_free_list);
            number_radicand->operation = 'r';
            number_radicand->value.denominator = integer_new(integer_pool, 1, 1);
            number_radicand->value.numerator = integer_new(&stack_a, 1, radicand_sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(&stack_a, &stack_b,
                    factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(&stack_a, &stack_a, stack_b,
                    factors[factor_index].value, reduced_multiplicity);
                number_radicand->value.numerator = integer_multiply(&stack_a, stack_b,
                    number_radicand->value.numerator, exponentiation);
            }
            integer_move(integer_pool, &number_radicand->value.numerator);
            struct Number*surd = number_slot_new(number_stack, number_free_list);
            surd->operation = '^';
            surd->left = number_radicand;
            surd->right = exponent;
            return number_multiply(number_stack, number_free_list, integer_pool, stack_a, &stack_b,
                number_coefficient, surd);
        }
        struct Number*new_denominator = number_slot_new(number_stack, number_free_list);
        new_denominator->operation = 'r';
        new_denominator->value.denominator = integer_new(integer_pool, 1, 1);
        new_denominator->value.numerator = integer_exponentiate(integer_pool, &stack_a, stack_b,
            base->value.denominator, exponent->value.numerator);
        struct Number*new_numerator_base = number_slot_new(number_stack, number_free_list);
        new_numerator_base->operation = 'r';
        new_numerator_base->value.denominator = integer_new(integer_pool, 1, 1);
        new_numerator_base->value.numerator = integer_multiply(&stack_a, stack_b,
            base->value.numerator,
            integer_exponentiate(&stack_a, &stack_a, stack_b, base->value.denominator,
                integer_subtract(&stack_a, stack_b, exponent->value.denominator, &INT(1, +))));
        integer_move(integer_pool, &new_numerator_base->value.numerator);
        struct Number*new_numerator = number_exponentiate(number_stack, number_free_list,
            integer_pool, stack_a, stack_b, new_numerator_base, exponent);
        number_free(number_free_list, integer_pool, new_numerator_base);
        number_reciprocal(new_denominator);
        struct Number*out = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
            &stack_b, new_numerator, new_denominator);
        number_free(number_free_list, integer_pool, new_numerator);
        number_free(number_free_list, integer_pool, new_denominator);
        return out;
    }
    ABORT("number_exponentiate case not yet implemented.");
}

struct Number*number_evaluate_root(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    switch (a->operation)
    {
    case '+':
        return number_add(number_stack, number_free_list, integer_pool, *stack_a, stack_b, a->left,
            a->right);
    case '-':
    {
        struct Number*negative = number_slot_new(number_stack, number_free_list);
        negative->operation = 'r';
        negative->value.numerator = integer_new(integer_pool, 1, -1);
        negative->value.denominator = integer_new(integer_pool, 1, 1);
        struct Number*product = number_multiply(number_stack, number_free_list, integer_pool,
            *stack_a, stack_b, negative, a->right);
        if (!product)
        {
            return 0;
        }
        return number_add(number_stack, number_free_list, integer_pool, *stack_a, stack_b, a->left,
            product);
    }
    case '*':
        return number_multiply(number_stack, number_free_list, integer_pool, *stack_a, stack_b,
            a->left, a->right);
    case '/':
        number_reciprocal(a->right);
        return number_multiply(number_stack, number_free_list, integer_pool, *stack_a, stack_b,
            a->left, a->right);
    case '^':
        return number_exponentiate(number_stack, number_free_list, integer_pool, *stack_a, *stack_b,
            a->left, a->right);
    default:
        ABORT("Number operation not recognized.");
    }
}

struct Number*number_evaluate(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    if (a->operation == 'r')
    {
        return a;
    }
    a->left =
        number_evaluate(number_stack, number_free_list, integer_pool, stack_a, stack_b, a->left);
    if (!a->left)
    {
        return 0;
    }
    a->right =
        number_evaluate(number_stack, number_free_list, integer_pool, stack_a, stack_b, a->right);
    if (!a->right)
    {
        return 0;
    }
    return number_evaluate_root(number_stack, number_free_list, integer_pool, stack_a, stack_b, a);
}

void print_number(struct Stack stack_a, struct Stack*stack_b, struct Number*number)
{
    if (number->operation == 'r')
    {
        char*string = stack_a.cursor;
        integer_string(&stack_a, *stack_b, number->value.numerator);
        if (!integer_equals_one(number->value.denominator))
        {
            *(char*)stack_slot_new(&stack_a, 1) = '/';
            integer_string(&stack_a, *stack_b, number->value.denominator);
        }
        *(char*)stack_slot_new(&stack_a, 1) = 0;
        printf("%s", string);
    }
    else
    {
        printf("(");
        print_number(stack_a, stack_b, number->left);
        printf(")");
        printf("%c", number->operation);
        printf("(");
        print_number(stack_a, stack_b, number->right);
        printf(")");
    }
}

void init()
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwPageSize;
    stack_new(&prime_stack);
    primes = prime_stack.cursor;
    integer_new(&prime_stack, 2, 1);
    integer_new(&prime_stack, 3, 1);
    permanent_integer_pool = VirtualAlloc(0, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    integer_pool_new(permanent_integer_pool);
    pi_estimate_min.numerator = integer_new(permanent_integer_pool, 47, 1);
    pi_estimate_min.denominator = integer_new(permanent_integer_pool, 15, 1);
    pi_interval_size.numerator = integer_new(permanent_integer_pool, 1696, 1);
    pi_interval_size.denominator = integer_new(permanent_integer_pool, 12285, 1);
    pi_sixteen_to_the_k = integer_new(permanent_integer_pool, 16, 1);
    pi_eight_k = integer_new(permanent_integer_pool, 8, 1);
}

int main()
{
    init();
    struct Stack number_stack;
    stack_new(&number_stack);
    void*number_stack_start = number_stack.cursor;
    struct Number*number_free_list = 0;
    struct Stack stack_a;
    stack_new(&stack_a);
    struct Stack stack_b;
    stack_new(&stack_b);
    while (true)
    {
        struct IntegerPool*integer_pool =
            VirtualAlloc(0, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
        integer_pool_new(integer_pool);
        if (get_input(&number_stack, &number_free_list, integer_pool, stack_a, &stack_b))
        {
            struct Number*number = number_stack_start;
            while (number->next)
            {
                if ((number->operation == 'r' && number->next->operation == '(') ||
                    (number->operation == ')' &&
                    (number->next->operation == 'r' || number->next->operation == '(')))
                {
                    struct Number*times = number_slot_new(&number_stack, &number_free_list);
                    times->operation = '*';
                    times->previous = number;
                    times->next = number->next;
                    number->next->previous = times;
                    number->next = times;
                }
                number = number->next;
            }
            struct Number*input =
                parse_input(&number_stack, &number_free_list, integer_pool, number_stack_start);
            if (input)
            {
                struct Number*evaluation = number_evaluate(&number_stack, &number_free_list,
                    integer_pool, &stack_a, &stack_b, input);
                if (evaluation)
                {
                    printf("=\n");
                    print_number(stack_a, &stack_b, evaluation);
                }
            }
        }
        printf("\n\n");
        rewind_stack_cursor(&number_stack, number_stack_start);
        number_free_list = 0;
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