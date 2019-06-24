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

struct Stack
{
    size_t start;
    size_t end;
    void*cursor;
    size_t allocation_cursor;
};

struct Integer
{
    size_t value_count;
    int8_t sign;
    uint32_t value;//The first element of an array of length value_count.
};

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

struct SizedIntegerPool
{
    struct IntegerSlotPage*first_page;
    struct IntegerSlot*cursor;
    struct IntegerSlot*free_list;
};

struct IntegerPool
{
    struct SizedIntegerPool*sized_pools;
    struct Stack slot_page_stack;
};

struct RingOperations
{
    void*(*copy)(struct Stack*, void*);
    bool(*equals_zero)(void*);
    void*additive_identity;
    void*multiplicative_identity;
    void*(*add)(struct Stack*, struct Stack*, void*, void*, void*);
    void*(*negative)(struct Stack*, struct Stack*, void*, void*);
    void*(*multiply)(struct Stack*, struct Stack*, void*, void*, void*);
};

struct Division
{
    void*quotient;
    void*remainder;
};

struct ExtendedGCDInfo
{
    void*gcd;
    void*a_coefficient;
    void*b_coefficient;
    void*a_over_gcd;
    void*b_over_gcd;
};

struct IntegerDivision
{
    struct Integer*quotient;
    struct Integer*remainder;
};

struct Rational
{
    struct Integer*numerator;
    struct Integer*denominator;
};

struct RationalInterval
{
    struct Rational*min;
    struct Rational*max;
};

struct Polynomial
{
    size_t coefficient_count;
    void*coefficients;//The first element of an array of void*'s of length coefficient_count.
};

struct PolynomialDivision
{
    struct Polynomial*quotient;
    struct Polynomial*remainder;
};

struct IntegerPolynomial
{
    size_t coefficient_count;
    struct Integer*coefficients;
};

struct RationalPolynomial
{
    size_t coefficient_count;
    struct Rational*coefficients;
};

struct Float
{
    struct Integer*significand;
    struct Integer*exponent;
};

struct FloatInterval
{
    struct Float min;
    struct Float max;
};

struct Factor
{
    struct Integer*value;
    struct Integer*multiplicity;
};

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
    struct RationalPolynomial*minimal_polynomial;
    struct Number*previous;
    struct Number*next;
    char operation;
};

size_t page_size;

struct Integer zero = { 0, 0, 0 };
struct Integer one = { 1, 1, 1 };

struct Rational rational_zero = { &zero, &one };
struct Rational rational_one = { &one, &one };

size_t pool_memory_start;
size_t pool_memory_end;

struct IntegerPool permanent_integer_pool;
struct Stack prime_stack;

struct Integer*primes;

struct Rational pi_estimate_min;
struct Rational pi_interval_size;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

struct Polynomial polynomial_zero = { 0, 0 };
struct IntegerPolynomial integer_polynomial_one = { 1, &one };
struct RationalPolynomial rational_polynomial_one = { 1, &rational_one };

struct Float float_zero = { &zero, &zero };
struct Float float_one = { &one, &zero };

size_t next_page_boundary(size_t address)
{
    size_t distance_from_page_start = address % page_size;
    if (distance_from_page_start)
    {
        return address + page_size - distance_from_page_start;
    }
    return address;
}

size_t next_aligned_address(size_t address, size_t alignment)
{
    return ((address + alignment - 1) / alignment) * alignment;
}

void stack_new(struct Stack*out, size_t start, size_t end)
{
    out->start = start;
    out->end = end;
    out->cursor = (void*)start;
    out->allocation_cursor = out->start;
}

void stack_free(struct Stack*out)
{
    while (out->allocation_cursor > out->start)
    {
        out->allocation_cursor -= page_size;
        VirtualFree((void*)out->allocation_cursor, 0, MEM_RELEASE);
    }
    out->cursor = (void*)out->start;
}

void*stack_slot_new(struct Stack*output_stack, size_t slot_size, size_t alignment)
{
    size_t slot = next_aligned_address((size_t)output_stack->cursor, alignment);
    output_stack->cursor = (void*)(slot + slot_size);
    if ((size_t)output_stack->cursor > output_stack->end)
    {
        printf("Stack ran out of virtual address space.");
        abort();
    }
    while ((size_t)output_stack->cursor > output_stack->allocation_cursor)
    {
        if (!VirtualAlloc((void*)output_stack->allocation_cursor, page_size,
            MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE))
        {
            printf("Ran out of physical memory.");
            abort();
        }
        output_stack->allocation_cursor = output_stack->allocation_cursor + page_size;
    }
    return (void*)slot;
}

size_t integer_size(size_t value_count)
{
    return sizeof(struct Integer) + (value_count - 1) * sizeof(uint32_t);
}

struct Integer*integer_stack_slot_new(struct Stack*output_stack, size_t value_count)
{
    return stack_slot_new(output_stack, integer_size(value_count), _Alignof(struct Integer));
}

struct IntegerSlotPage*integer_pool_allocate_slot_page(struct IntegerPool*pool)
{
    pool->slot_page_stack.allocation_cursor = pool->slot_page_stack.allocation_cursor + page_size;
    if (pool->slot_page_stack.allocation_cursor > pool->slot_page_stack.end)
    {
        printf("Integer pool ran out of virtual address space.");
        abort();
    }
    void*new_page_address = pool->slot_page_stack.cursor;
    pool->slot_page_stack.cursor = (void*)pool->slot_page_stack.allocation_cursor;
    return VirtualAlloc(new_page_address, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
}

struct SizedIntegerPool*get_sized_integer_pool(struct IntegerPool*pool, size_t value_count)
{
    return &pool->sized_pools[(value_count - 1) / 2];
}

struct Integer*integer_pool_slot_new(struct IntegerPool*pool, size_t value_count)
{
    if (!value_count)
    {
        return &zero;
    }
    struct SizedIntegerPool*sized_pool = get_sized_integer_pool(pool, value_count);
    if (!sized_pool->first_page)
    {
        sized_pool->first_page = integer_pool_allocate_slot_page(pool);
        sized_pool->cursor = &sized_pool->first_page->memory;
    }
    struct IntegerSlot*out = sized_pool->free_list;
    if (out)
    {
        sized_pool->free_list = out->next_slot;
    }
    else
    {
        out = (struct IntegerSlot*)next_aligned_address((size_t)sized_pool->cursor,
            _Alignof(struct IntegerSlot));
        sized_pool->cursor = (struct IntegerSlot*)((size_t)out +
            offsetof(struct IntegerSlot, integer) + integer_size(value_count));
        size_t boundary = next_page_boundary((size_t)out);
        if ((size_t)sized_pool->cursor > boundary)
        {
            struct IntegerSlotPage*new_page = integer_pool_allocate_slot_page(pool);
            ((struct IntegerSlotPage*)(boundary - page_size))->next_page = new_page;
            out = &new_page->memory;
        }
    }
    return &out->integer;
}

void integer_pool_new(struct IntegerPool*out, size_t start, size_t end)
{
    out->sized_pools =
        VirtualAlloc((void*)start, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    stack_new(&out->slot_page_stack, start + page_size, end);
}

void*generic_gcd(struct RingOperations*operations,
    void(*euclidean_divide)(struct Stack*, struct Stack*, struct Division*, void*, void*, void*),
    struct Stack*output_stack, struct Stack*local_stack, void*a, void*b, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!operations->equals_zero(b))
    {
        void*c = b;
        struct Division division;
        euclidean_divide(local_stack, output_stack, &division, a, b, misc);
        b = division.remainder;
        a = c;
    }
    a = operations->copy(output_stack, a);
    local_stack->cursor = local_stack_savepoint;
    return a;
}

void generic_extended_gcd(struct RingOperations*operations,
    void(*euclidean_divide)(struct Stack*, struct Stack*, struct Division*, void*, void*, void*),
    struct Stack*output_stack, struct Stack*local_stack, struct ExtendedGCDInfo*out, void*a, void*b,
    void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->a_coefficient = operations->additive_identity;
    out->b_coefficient = operations->multiplicative_identity;
    out->a_over_gcd = operations->additive_identity;
    out->b_over_gcd = operations->multiplicative_identity;
    while (!operations->equals_zero(a))
    {
        struct Division division;
        euclidean_divide(local_stack, output_stack, &division, b, a, misc);
        struct Integer*m = operations->add(local_stack, output_stack, out->a_coefficient,
            operations->negative(local_stack, output_stack,
                operations->multiply(local_stack, output_stack, out->b_over_gcd, division.quotient,
                    misc), misc), misc);
        struct Integer*n = operations->add(local_stack, output_stack, out->b_coefficient,
            operations->negative(local_stack, output_stack,
                operations->multiply(local_stack, output_stack, out->a_over_gcd, division.quotient,
                    misc), misc), misc);
        b = a;
        a = division.remainder;
        out->a_coefficient = out->b_over_gcd;
        out->b_coefficient = out->a_over_gcd;
        out->b_over_gcd = m;
        out->a_over_gcd = n;
    }
    out->a_coefficient = operations->copy(output_stack, out->a_coefficient);
    out->a_over_gcd = operations->copy(output_stack, out->a_over_gcd);
    out->b_coefficient = operations->copy(output_stack, out->b_coefficient);
    out->b_over_gcd = operations->negative(output_stack, local_stack, out->b_over_gcd, misc);
    out->gcd = operations->copy(output_stack, b);
    local_stack->cursor = local_stack_savepoint;
}

#define INT(value, sign) (struct Integer){ 1, sign 1, value }

struct Integer*stack_integer_new(struct Stack*stack, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_stack_slot_new(stack, 1);
    out->value_count = 1;
    out->sign = sign;
    out->value = value;
    return out;
}

struct Integer*pool_integer_new(struct IntegerPool*pool, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_pool_slot_new(pool, 1);
    out->value_count = 1;
    out->sign = sign;
    out->value = value;
    return out;
}

void pool_integer_free(struct IntegerPool*pool, struct Integer*a)
{
    ASSERT(pool_memory_start <= (size_t)a && (size_t)a <= pool_memory_end,
        "pool_integer_free argument wasn't in a pool.");
    if (a->value_count)
    {
        struct SizedIntegerPool*sized_pool = get_sized_integer_pool(pool, a->value_count);
        struct IntegerSlot*slot =
            (struct IntegerSlot*)((size_t)a - offsetof(struct IntegerSlot, integer));
        slot->next_slot = sized_pool->free_list;
        sized_pool->free_list = slot;
    }
}

struct Integer*integer_from_char(struct Stack*output_stack, char value)
{
    uint32_t uint_value = value - '0';
    if (value)
    {
        struct Integer*out = integer_stack_slot_new(output_stack, 1);
        out->value_count = 1;
        out->sign = 1;
        out->value = uint_value;
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
        if ((&a->value)[i - 1] != 0)
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
    struct Integer*out = integer_stack_slot_new(output_stack, value_count);
    out->value_count = value_count;
    *(size_t*)&out->value = value;
    integer_trim_leading_zeroes(out);
    output_stack->cursor = &out->value + out->value_count;
    return out;
}

struct Integer*integer_copy_to_stack(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*copy = integer_stack_slot_new(output_stack, a->value_count);
    memcpy(copy, a, integer_size(a->value_count));
    return copy;
}

struct Integer*integer_copy_to_pool(struct IntegerPool*pool, struct Integer*a)
{
    struct Integer*out = integer_pool_slot_new(pool, a->value_count);
    memcpy(out, a, integer_size(a->value_count));
    return out;
}

void integer_move_to_pool(struct IntegerPool*pool, struct Integer**a)
{
    struct Integer*pool_copy = integer_pool_slot_new(pool, (*a)->value_count);
    memcpy(pool_copy, *a, integer_size((*a)->value_count));
    *a = pool_copy;
}

void integer_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack, struct Integer**a)
{
    ASSERT(pool_memory_start <= (size_t)*a && (size_t)*a <= pool_memory_end,
        "integer_move_from_pool argument wasn't in a pool.");
    struct Integer*pool_copy = *a;
    *a = integer_copy_to_stack(output_stack, *a);
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

bool integer_equals_zero(struct Integer*a)
{
    return a->value_count == 0;
}

struct Integer*integer_magnitude(void*output_arena, struct Integer*a)
{
    struct Integer*out = integer_copy_to_stack(output_arena, a);
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

//Assumes that, at a->values, there is enough allocated but unused memory for a->value_count to
//increase to max(a->value_count, b->value_count) + 1.
void integer_add_to_a_in_place(struct Integer*a, struct Integer*b)
{
    ASSERT(pool_memory_start > (size_t)a || (size_t)a > pool_memory_end,
        "The a argument of integer_add_to_a_in_place was a pool Integer.");
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
        (&a->value)[a->value_count] = 0;
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
        if ((&a->value)[a->value_count - 1] != 0)
        {
            twos_complement(a);
            a->sign *= -1;
        }
    }
    integer_trim_leading_zeroes(a);
}

struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b)
{
    struct Integer*out = integer_copy_to_stack(output_stack, a);
    stack_slot_new(output_stack, (b->value_count + 1) * sizeof(uint32_t), _Alignof(struct Integer));
    integer_add_to_a_in_place(out, b);
    output_stack->cursor = &out->value + out->value_count;
    return out;
}

struct Integer*integer_generic_add(struct Stack*output_stack, struct Stack*unused_stack,
    struct Integer*a, struct Integer*b, void*unused)
{
    return integer_add(output_stack, a, b);
}

struct Integer*integer_negative(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy_to_stack(output_stack, a);
    out->sign = -out->sign;
    return out;
}

struct Integer*integer_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Integer*a, void*unused)
{
    return integer_negative(output_stack, a);
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

int8_t integer_compare(struct Stack*stack_a, struct Stack*stack_b, struct Integer*a,
    struct Integer*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = integer_subtract(stack_a, stack_b, a, b)->sign;
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Integer*integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*product = integer_stack_slot_new(output_stack, a->value_count + b->value_count);
    product->value_count = 0;
    product->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)(&a->value)[i] * (&b->value)[j];
            size_t shift = i + j;
            struct Integer*integer_component = integer_stack_slot_new(local_stack, shift + 2);
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
    output_stack->cursor = &product->value + product->value_count;
    local_stack->cursor = local_stack_savepoint;
    return product;
}

struct Integer*integer_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b, void*unused)
{
    return integer_multiply(output_stack, local_stack, a, b);
}

void integer_upshift(struct Integer*a, uint8_t shift)
{
    for (size_t i = a->value_count; i-- > 1;)
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
    ASSERT(pool_memory_start > (size_t)a || (size_t)a > pool_memory_end,
        "integer_halve was called on a pool Integer.");
    integer_downshift(a, 1);
    integer_trim_leading_zeroes(a);
}

void calculate_division_values(struct Stack*stack_a, struct Stack*stack_b,
    struct Integer*divisor_magnitude, struct IntegerDivision*division, size_t quotient_value_index,
    uint32_t quotient_digit)
{
    void*stack_a_savepoint = stack_a->cursor;
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer*difference =
                integer_subtract(stack_a, stack_b, division->remainder, divisor_magnitude);
            if (difference->sign >= 0)
            {
                (&division->quotient->value)[quotient_value_index] |= quotient_digit;
                division->remainder = difference;
            }
            if (quotient_digit == 1)
            {
                if (quotient_value_index == 0)
                {
                    stack_a->cursor = stack_a_savepoint;
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

void integer_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerDivision*out, struct Integer*dividend, struct Integer*divisor)
{
    int8_t quotient_sign = dividend->sign * divisor->sign;
    int dividend_leading_digit_place = leading_digit_place(dividend);
    int divisor_leading_digit_place = leading_digit_place(divisor);
    if (dividend->value_count > divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
            dividend_leading_digit_place >= divisor_leading_digit_place))
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Integer*divisor_magnitude =
            integer_stack_slot_new(local_stack, dividend->value_count);
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
        out->quotient = integer_stack_slot_new(output_stack, dividend->value_count);
        out->quotient->value_count = dividend->value_count;
        out->quotient->sign = quotient_sign;
        memset(&out->quotient->value, 0, out->quotient->value_count * sizeof(uint32_t));
        out->remainder = dividend;
        calculate_division_values(output_stack, local_stack, divisor_magnitude, out,
            quotient_value_index, quotient_digit);
        integer_trim_leading_zeroes(out->quotient);
        output_stack->cursor = &out->quotient->value + out->quotient->value_count;
        out->remainder = integer_copy_to_stack(output_stack, out->remainder);
        if (out->remainder->sign != 0)
        {
            out->remainder->sign = dividend->sign;
        }
        dividend->sign = dividend_sign;
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        out->quotient = &zero;
        out->remainder = integer_copy_to_stack(output_stack, dividend);
    }
}

void integer_generic_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerDivision*out, struct Integer*dividend, struct Integer*divisor, void*unused)
{
    integer_euclidean_divide(output_stack, local_stack, out, dividend, divisor);
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
    if ((&a->value)[a->value_count - 1] & 0x80000000)
    {
        out = integer_stack_slot_new(output_stack, a->value_count + 1);
        out->value_count = a->value_count + 1;
        (&out->value)[a->value_count] = 0;
    }
    else
    {
        out = integer_stack_slot_new(output_stack, a->value_count);
        out->value_count = a->value_count;
    }
    out->sign = a->sign;
    memcpy(&out->value, &a->value, a->value_count * sizeof(uint32_t));
    integer_upshift(out, 1);
    return out;
}

struct Integer*integer_half(struct Stack*output_stack, struct Integer*a)
{
    struct Integer*out = integer_copy_to_stack(output_stack, a);
    integer_halve(out);
    return out;
}

void*generic_exponentiate(struct RingOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, void*base, struct Integer*exponent, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = operations->multiplicative_identity;
    void*base_to_a_power_of_two = base;
    while (exponent->sign > 0)
    {
        if (exponent->value & 1)
        {
            out = operations->multiply(local_stack, output_stack, out, base_to_a_power_of_two, 0);
        }
        base_to_a_power_of_two = operations->multiply(local_stack, output_stack,
            base_to_a_power_of_two, base_to_a_power_of_two, misc);
        exponent = integer_half(local_stack, exponent);
    }
    out = operations->copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RingOperations integer_operations = { integer_copy_to_stack, integer_equals_zero, &zero,
    &one, integer_generic_add, integer_generic_negative, integer_generic_multiply };

struct Integer*integer_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*base, struct Integer*exponent)
{
    return generic_exponentiate(&integer_operations, output_stack, local_stack, base, exponent, 0);
}

void integer_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct ExtendedGCDInfo*out, struct Integer*a, struct Integer*b)
{
    generic_extended_gcd(&integer_operations, integer_generic_euclidean_divide, output_stack,
        local_stack, out, a, b, 0);
}

struct Integer*integer_gcd(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a,
    struct Integer*b)
{
    return generic_gcd(&integer_operations, integer_generic_euclidean_divide, output_stack,
        local_stack, a, b, 0);
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

void integer_string(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a)
{
    if (a->sign == 0)
    {
        *(char*)stack_slot_new(output_stack, 1, 1) = '0';
        return;
    }
    if (a->sign < 0)
    {
        *(char*)stack_slot_new(output_stack, 1, 1) = '-';
    }
    void*local_stack_savepoint = local_stack->cursor;
    char*buffer_start = stack_slot_new(output_stack, 10 * a->value_count + 1, 1);
    char*next_char = output_stack->cursor;
    struct Integer*quotient = a;
    struct Integer power = { 1, 1, 10 };
    while (quotient->sign != 0)
    {
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division, quotient, &power);
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
    output_stack->cursor = buffer_start + char_count;
    local_stack->cursor = local_stack_savepoint;
}

struct Integer*next_prime_slot(struct Integer*prime)
{
    return (struct Integer*)((size_t)prime + integer_size(prime->value_count));
}

struct Integer*next_prime(struct Stack*stack_a, struct Stack*stack_b, struct Integer*prime)
{
    struct Integer*next_slot = next_prime_slot(prime);
    size_t boundary = next_page_boundary((size_t)prime);
    if ((size_t)next_slot == boundary)
    {
        VirtualAlloc((void*)boundary, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    }
    if (next_slot->value_count)
    {
        return next_slot;
    }
    void*stack_a_savepoint = stack_a->cursor;
    memcpy(next_slot, prime, integer_size(prime->value_count));
    integer_add_to_a_in_place(next_slot, primes);
    while (true)
    {
        struct Integer*potential_factor = primes;
        while (true)
        {
            if ((integer_euclidean_remainder(stack_a, stack_b, next_slot,
                potential_factor))->value_count == 0)
            {
                break;
            }
            potential_factor = next_prime_slot(potential_factor);
            if (potential_factor == next_slot)
            {
                stack_a->cursor = stack_a_savepoint;
                return next_slot;
            }
        }
        integer_add_to_a_in_place(next_slot, primes);
    }
}

struct Rational*rational_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*numerator, struct Integer*denominator)
{
    if (denominator->sign == 0)
    {
        printf("Tried to divide by 0.");
        return 0;
    }
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    if (numerator->sign == 0)
    {
        out->numerator = &zero;
        out->denominator = &one;
        return out;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo gcd_info;
    generic_extended_gcd(&integer_operations, integer_generic_euclidean_divide, local_stack,
        output_stack, &gcd_info, numerator, denominator, 0);
    out->numerator = gcd_info.a_over_gcd;
    out->denominator = gcd_info.b_over_gcd;
    if (out->denominator->sign < 0)
    {
        out->numerator->sign = -out->numerator->sign;
        out->denominator->sign = -out->denominator->sign;
    }
    out->numerator = integer_copy_to_stack(output_stack, out->numerator);
    out->denominator = integer_copy_to_stack(output_stack, out->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void rational_free(struct IntegerPool*integer_pool, struct Rational*a)
{
    pool_integer_free(integer_pool, a->numerator);
    pool_integer_free(integer_pool, a->denominator);
}

struct Rational*rational_copy_to_stack(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_copy_to_stack(output_stack, a->numerator);
    out->denominator = integer_copy_to_stack(output_stack, a->denominator);
    return out;
}

void rational_move_to_pool(struct IntegerPool*pool, struct Rational*a)
{
    integer_move_to_pool(pool, &a->numerator);
    integer_move_to_pool(pool, &a->denominator);
}

void rational_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack, struct Rational*a)
{
    integer_move_from_pool(pool, output_stack, &a->numerator);
    integer_move_from_pool(pool, output_stack, &a->denominator);
}

bool rational_equals_zero(struct Rational*a)
{
    return integer_equals_zero(a->numerator) == 0;
}

struct Rational*rational_magnitude(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_magnitude(output_stack, a->numerator);
    out->denominator = a->denominator;
    return out;
}

struct Rational*rational_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack, integer_add(local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->denominator),
        integer_multiply(local_stack, output_stack, b->numerator, a->denominator)),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_generic_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b, void*unused)
{
    return rational_add(output_stack, local_stack, a, b);
}

struct Rational*rational_integer_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_add(output_stack, a->numerator,
        integer_multiply(local_stack, output_stack, b, a->denominator));
    out->denominator = integer_copy_to_stack(output_stack, a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_negative(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_negative(output_stack, a->numerator);
    out->denominator = integer_copy_to_stack(output_stack, a->denominator);
    return out;
}

struct Rational*rational_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Rational*a, void*unused)
{
    return rational_negative(output_stack, a);
}

struct Rational*rational_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*minuend, struct Rational*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_subtract(local_stack, output_stack,
            integer_multiply(local_stack, output_stack, minuend->numerator,
                subtrahend->denominator),
            integer_multiply(local_stack, output_stack, subtrahend->numerator,
                minuend->denominator)),
        integer_multiply(local_stack, output_stack, minuend->denominator,
            subtrahend->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->numerator),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b, void*unused)
{
    return rational_multiply(output_stack, local_stack, a, b);
}

struct Rational*rational_doubled(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_doubled(local_stack, a->numerator), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_unreduced_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b, void*unused)
{
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_multiply(output_stack, local_stack, a->numerator, b->numerator);
    out->denominator = integer_multiply(output_stack, local_stack, a->denominator, b->denominator);
    return out;
}

struct RingOperations rational_ring_operations = { rational_copy_to_stack, rational_equals_zero,
    &rational_zero, &rational_one, rational_generic_add, rational_generic_negative,
    rational_generic_multiply };

struct Rational*rational_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, void*unused)
{
    struct Rational*out =
        stack_slot_new(output_stack, sizeof(struct Rational), _Alignof(struct Rational));
    out->numerator = integer_copy_to_stack(output_stack, a->denominator);
    out->denominator = integer_copy_to_stack(output_stack, a->numerator);
    return out;
}

struct Rational*rational_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Rational*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, dividend->numerator, divisor->denominator),
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor->numerator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_generic_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Rational*divisor, void*misc)
{
    return rational_divide(output_stack, local_stack, dividend, divisor);
}

struct Rational*rational_integer_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Integer*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack, dividend->numerator,
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*base, struct Integer*exponent)
{
    return generic_exponentiate(&(struct RingOperations){ rational_copy_to_stack,
        rational_equals_zero, &rational_zero, &rational_one, rational_generic_add,
        rational_generic_negative, rational_unreduced_multiply },
        output_stack, local_stack, base, exponent, 0);
}

int8_t rational_compare(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = integer_compare(stack_a, stack_b,
        integer_multiply(stack_a, stack_b, a->numerator, b->denominator),
        integer_multiply(stack_a, stack_b, a->denominator, b->numerator));
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Rational*rational_min(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    if (rational_compare(stack_a, stack_b, a, b) < 0)
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
    if (rational_compare(stack_a, stack_b, a, b) > 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

void rational_estimate_trig_function(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a_squared, struct Integer*factorial_component,
    struct Rational*delta, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (true)
    {
        if (rational_compare(output_stack, local_stack, rational_magnitude(local_stack, delta),
            interval_size) <= 0)
        {
            break;
        }
        out->min = rational_add(local_stack, output_stack, out->min, delta);
        factorial_component = integer_add(local_stack, factorial_component, &one);
        delta = rational_multiply(local_stack, output_stack, delta,
            rational_integer_divide(local_stack, output_stack, a_squared, factorial_component));
        factorial_component = integer_add(local_stack, factorial_component, &one);
        struct Integer*negative_factorial_component =
            integer_negative(local_stack, factorial_component);
        delta =
            rational_integer_divide(local_stack, output_stack, delta, negative_factorial_component);
    }
    out->min = rational_copy_to_stack(output_stack, out->min);
    if (delta->numerator->value_count > 0)
    {
        out->max = rational_add(output_stack, local_stack, out->min, delta);
    }
    else
    {
        out->max = out->min;
        out->min = rational_add(output_stack, local_stack, out->max, delta);
    }
    local_stack->cursor = local_stack_savepoint;
}

void rational_estimate_cosine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min->numerator = &one;
    out->min->denominator = &one;
    struct Rational*a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = &INT(2, +);
    struct Rational*delta =
        rational_integer_divide(local_stack, output_stack, a_squared, &INT(2, -));
    rational_estimate_trig_function(output_stack, local_stack, out, a_squared, factorial_component,
        delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
}

void rational_estimate_sine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min = a;
    struct Rational*a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = &INT(3, +);
    struct Rational*delta = rational_integer_divide(local_stack, output_stack,
        rational_multiply(local_stack, output_stack, a_squared, a), &INT(6, -));
    rational_estimate_trig_function(output_stack, local_stack, out, a_squared, factorial_component,
        delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
}

void pi_refine_interval(struct Stack*stack_a, struct Stack*stack_b)
{
    void*stack_a_savepoint = stack_a->cursor;
    pi_estimate_min = *rational_add(stack_a, stack_b, &pi_estimate_min,
        rational_integer_divide(stack_a, stack_b, rational_subtract(stack_a, stack_b,
            &(struct Rational){ &INT(4, +), integer_add(stack_a, pi_eight_k, &one) },
            rational_add(stack_a, stack_b,
                rational_reduced(stack_a, stack_b, &INT(2, +),
                    integer_add(stack_a, pi_eight_k, &INT(4, +))),
                rational_add(stack_a, stack_b,
                    &(struct Rational){ &INT(1, +), integer_add(stack_a, pi_eight_k, &INT(5, +)) },
                    &(struct Rational){ &INT(1, +), integer_add(stack_a, pi_eight_k, &INT(6, +)) })
            )), pi_sixteen_to_the_k));
    pi_interval_size = *rational_multiply(stack_a, stack_b, &pi_interval_size,
        &(struct Rational){ &one, &INT(16, +) });
    pi_eight_k = integer_add(stack_a, pi_eight_k, &INT(8, +));
    pi_sixteen_to_the_k = integer_multiply(stack_a, stack_b, pi_sixteen_to_the_k, &INT(16, +));
    stack_a->cursor = stack_a_savepoint;
}

void pi_estimate(struct Stack*stack_a, struct Stack*stack_b, struct Rational*interval_size)
{
    int8_t interval_size_comparison =
        rational_compare(stack_a, stack_b, &pi_interval_size, interval_size);
    if (interval_size_comparison > 0)
    {
        void*stack_a_savepoint = stack_a->cursor;
        integer_move_from_pool(&permanent_integer_pool, stack_a, &pi_sixteen_to_the_k);
        integer_move_from_pool(&permanent_integer_pool, stack_a, &pi_eight_k);
        while (interval_size_comparison > 0)
        {
            pi_refine_interval(stack_a, stack_b);
            interval_size_comparison =
                rational_compare(stack_a, stack_b, &pi_interval_size, interval_size);
        }
        pi_sixteen_to_the_k = integer_copy_to_pool(&permanent_integer_pool, pi_sixteen_to_the_k);
        pi_eight_k = integer_copy_to_pool(&permanent_integer_pool, pi_eight_k);
        stack_a->cursor = stack_a_savepoint;
    }
}

void pi_shrink_interval_to_one_side_of_value(struct Stack*stack_a, struct Stack*stack_b,
    struct Rational*value)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*pi_estimate_max =
        rational_add(stack_a, stack_b, &pi_estimate_min, &pi_interval_size);
    if (rational_compare(stack_a, stack_b, &pi_estimate_min, value) <= 0 &&
        rational_compare(stack_a, stack_b, value, pi_estimate_max) <= 0)
    {
        integer_move_from_pool(&permanent_integer_pool, stack_a, &pi_sixteen_to_the_k);
        integer_move_from_pool(&permanent_integer_pool, stack_a, &pi_eight_k);
        while (true)
        {
            pi_refine_interval(stack_a, stack_b);
            pi_estimate_max =
                rational_add(stack_a, stack_b, &pi_estimate_min, &pi_interval_size);
            if (rational_compare(stack_a, stack_b, &pi_estimate_min, value) > 0 ||
                rational_compare(stack_a, stack_b, value, pi_estimate_max) > 0)
            {
                break;
            }
        }
        pi_sixteen_to_the_k = integer_copy_to_pool(&permanent_integer_pool, pi_sixteen_to_the_k);
        pi_eight_k = integer_copy_to_pool(&permanent_integer_pool, pi_eight_k);
    }
    stack_a->cursor = stack_a_savepoint;
}

void*polynomial_slot_new(struct Stack*output_stack, size_t coefficient_count)
{
    struct Polynomial*out = stack_slot_new(output_stack, (coefficient_count + 1) * sizeof(size_t),
        _Alignof(struct Polynomial));
    out->coefficient_count = coefficient_count;
    return out;
}

void polynomial_copy_coefficients(void*(coefficient_copy)(struct Stack*, void*),
    struct Stack*output_stack, struct Polynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        (&a->coefficients)[i] = coefficient_copy(output_stack, (&a->coefficients)[i]);
    }
}

void*polynomial_copy(void*(coefficient_copy)(struct Stack*, void*), struct Stack*output_stack,
    struct Polynomial*a)
{
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    memcpy(&out->coefficients, &a->coefficients, a->coefficient_count * sizeof(void*));
    polynomial_copy_coefficients(coefficient_copy, output_stack, a);
    return out;
}

void polynomial_trim_leading_zeroes(bool(coefficient_equals_zero)(void*), struct Polynomial*a)
{
    for (size_t i = *(size_t*)a; i-- > 0;)
    {
        if (coefficient_equals_zero((&a->coefficients)[i]))
        {
            *(size_t*)a -= 1;
        }
        else
        {
            return;
        }
    }
}

bool polynomial_equals_zero(struct Polynomial*a)
{
    return a->coefficient_count == 0;
}

void*polynomial_add(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        void*temp = a;
        a = b;
        b = temp;
    }
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        (&out->coefficients)[i] = coefficient_operations->add(output_stack, local_stack,
            (&a->coefficients)[i], (&b->coefficients)[i], 0);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        (&a->coefficients)[i] = coefficient_operations->copy(output_stack, (&a->coefficients)[i]);
    }
    polynomial_trim_leading_zeroes(coefficient_operations->equals_zero, out);
    return out;
}

void*polynomial_negative(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Polynomial*a)
{
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        (&out->coefficients)[i] =
            coefficient_operations->negative(output_stack, 0, (&a->coefficients)[i], 0);
    }
    return out;
}

void*polynomial_subtract(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = polynomial_add(coefficient_operations, output_stack, local_stack, a,
        polynomial_negative(coefficient_operations, local_stack, b));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*polynomial_multiply(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_slot_new(output_stack, 0);
    }
    struct Polynomial*out =
        polynomial_slot_new(output_stack, a->coefficient_count + b->coefficient_count - 1);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        (&out->coefficients)[i] = coefficient_operations->additive_identity;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            (&out->coefficients)[i + j] = coefficient_operations->add(output_stack, local_stack,
                (&out->coefficients)[i + j], coefficient_operations->multiply(local_stack,
                    output_stack, (&a->coefficients)[i], (&b->coefficients)[j], 0), 0);
        }
    }
    polynomial_trim_leading_zeroes(coefficient_operations->equals_zero, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*polynomial_multiply_by_coefficient(struct RingOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a, void*b)
{
    if (coefficient_operations->equals_zero(b))
    {
        return &polynomial_zero;
    }
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        (&out->coefficients)[i] = coefficient_operations->multiply(output_stack, local_stack,
            (&a->coefficients)[i], b, 0);
    }
    return out;
}

void polynomial_euclidean_divide(struct RingOperations*coefficient_operations,
    void*(coefficient_reciprocal)(struct Stack*, struct Stack*, void*, void*),
    struct Stack*output_stack, struct Stack*local_stack, struct PolynomialDivision*out,
    struct Polynomial*dividend, struct Polynomial*divisor, void*misc)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        out->quotient = polynomial_slot_new(output_stack, 0);
        out->quotient->coefficient_count = 0;
        out->remainder = polynomial_copy(coefficient_operations->copy, output_stack, dividend);
        return;
    }
    void*output_stack_savepoint = output_stack->cursor;
    void*local_stack_savepoint = local_stack->cursor;
    size_t quotient_coefficient_count =
        1 + dividend->coefficient_count - divisor->coefficient_count;
    out->quotient = polynomial_slot_new(output_stack, quotient_coefficient_count);
    out->remainder = polynomial_slot_new(output_stack, dividend->coefficient_count);
    memcpy(&out->remainder->coefficients, &dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    void*leading_coefficient_reciprocal = coefficient_reciprocal(local_stack, output_stack,
        (&divisor->coefficients)[divisor->coefficient_count - 1], misc);
    for (size_t i = 1; i <= quotient_coefficient_count; ++i)
    {
        void*quotient = coefficient_operations->multiply(local_stack, output_stack,
            (&out->remainder->coefficients)[out->remainder->coefficient_count - i],
            leading_coefficient_reciprocal, misc);
        memcpy(&out->quotient->coefficients + quotient_coefficient_count - i, quotient,
            sizeof(void*));
        for (size_t j = 1; j < divisor->coefficient_count; ++j)
        {
            (&out->remainder)[out->remainder->coefficient_count - i - j] =
                coefficient_operations->add(local_stack, output_stack,
                    (&out->remainder)[out->remainder->coefficient_count - i - j],
                    coefficient_operations->negative(local_stack, output_stack,
                        coefficient_operations->multiply(local_stack, output_stack, quotient,
                            (&divisor->coefficients)[divisor->coefficient_count - 1 - j], misc),
                        misc), misc);
        }
    }
    polynomial_copy_coefficients(coefficient_operations->copy, output_stack, out->quotient);
    polynomial_trim_leading_zeroes(coefficient_operations->equals_zero, out->remainder);
    polynomial_copy_coefficients(coefficient_operations->copy, output_stack, out->remainder);
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_copy(integer_copy_to_stack, output_stack, (struct Polynomial*)a);
}

struct IntegerPolynomial*integer_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    return polynomial_add(&integer_operations, output_stack, local_stack, (struct Polynomial*)a,
        (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_generic_add(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b, void*unused)
{
    return polynomial_add(&integer_operations, output_stack, local_stack, (struct Polynomial*)a,
        (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_negative(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_negative(&integer_operations, output_stack, (struct Polynomial*)a);
}

struct IntegerPolynomial*integer_polynomial_generic_negative(struct Stack*output_stack,
    struct Stack*unused_stack, struct IntegerPolynomial*a, void*unused)
{
    return polynomial_negative(&integer_operations, output_stack, (struct Polynomial*)a);
}

struct IntegerPolynomial*integer_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*minuend, struct IntegerPolynomial*subtrahend)
{
    return polynomial_subtract(&integer_operations, output_stack, local_stack,
        (struct Polynomial*)minuend, (struct Polynomial*)subtrahend);
}

struct IntegerPolynomial*integer_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    return polynomial_multiply(&integer_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_generic_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b, void*unused)
{
    return polynomial_multiply(&integer_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b)
{
    return polynomial_multiply_by_coefficient(&integer_operations, output_stack, local_stack,
        (struct Polynomial*)a, b);
}

struct RingOperations integer_polynomial_ring_operations = { integer_polynomial_copy,
    polynomial_equals_zero, &polynomial_zero, &integer_polynomial_one,
    integer_polynomial_generic_add, integer_polynomial_generic_negative,
    integer_polynomial_generic_multiply };

//Sets the fields of out to 0 when the quotient wouldn't have integer coefficients.
void integer_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        out->quotient = polynomial_slot_new(output_stack, 0);
        out->quotient->coefficient_count = 0;
        out->remainder = (struct Polynomial*)integer_polynomial_copy(output_stack, dividend);
        return;
    }
    void*output_stack_savepoint = output_stack->cursor;
    void*local_stack_savepoint = local_stack->cursor;
    size_t quotient_coefficient_count =
        1 + dividend->coefficient_count - divisor->coefficient_count;
    out->quotient = polynomial_slot_new(output_stack, quotient_coefficient_count);
    out->remainder = polynomial_slot_new(output_stack, dividend->coefficient_count);
    memcpy(&out->remainder->coefficients, &dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    for (size_t i = 1; i <= quotient_coefficient_count; ++i)
    {
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division,
            (&out->remainder->coefficients)[out->remainder->coefficient_count - i],
            (&divisor->coefficients)[divisor->coefficient_count - 1]);
        if (!division.remainder->value_count)
        {
            memcpy(&out->quotient->coefficients + quotient_coefficient_count - i, division.quotient,
                sizeof(struct Integer*));
            for (size_t j = 1; j < divisor->coefficient_count; ++j)
            {
                (&out->remainder->coefficients)[out->remainder->coefficient_count - i - j] =
                    integer_subtract(local_stack, output_stack,
                        (&out->remainder->coefficients)[out->remainder->coefficient_count - i - j],
                        integer_multiply(local_stack, output_stack, division.quotient,
                            (&divisor->coefficients)[divisor->coefficient_count - 1 - j]));
            }
        }
        else
        {
            out->quotient = 0;
            out->remainder = 0;
            return;
        }
    }
    polynomial_copy_coefficients(integer_copy_to_stack, output_stack, out->quotient);
    polynomial_trim_leading_zeroes(integer_equals_zero, out->remainder);
    polynomial_copy_coefficients(integer_copy_to_stack, output_stack, out->remainder);
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*integer_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor)
{
    struct PolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct IntegerPolynomial*)division.remainder;
}

struct IntegerPolynomial*integer_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = polynomial_slot_new(output_stack, a->coefficient_count - 1);
    struct Integer*multiplier = &zero;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        multiplier = integer_add(local_stack, multiplier, &one);
        (&out->coefficients)[i - 1] =
            integer_multiply(output_stack, local_stack, multiplier, (&a->coefficients)[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*integer_polynomial_content(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return &one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = a->coefficients;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        out = integer_gcd(local_stack, output_stack, out, (&a->coefficients)[i]);
    }
    out = integer_copy_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_integer_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct Integer*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = polynomial_slot_new(output_stack, dividend->coefficient_count);
    for (size_t i = 0; i < dividend->coefficient_count; ++i)
    {
        (&out->coefficients)[i] = integer_euclidean_quotient(output_stack, local_stack,
            (&dividend->coefficients)[i], divisor);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_integer_divide(output_stack, local_stack, a,
        integer_polynomial_content(local_stack, output_stack, a));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (a->coefficient_count == 0)
    {
        return integer_polynomial_copy(output_stack, b);
    }
    if (b->coefficient_count == 0)
    {
        return integer_polynomial_copy(output_stack, a);
    }
    if (b->coefficient_count > a->coefficient_count)
    {
        struct IntegerPolynomial*t = a;
        a = b;
        b = t;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*a_content = integer_polynomial_content(local_stack, output_stack, a);
    struct Integer*b_content = integer_polynomial_content(local_stack, output_stack, b);
    struct IntegerPolynomial*a_primitive_part =
        integer_polynomial_integer_divide(local_stack, output_stack, a, a_content);
    struct IntegerPolynomial*b_primitive_part =
        integer_polynomial_integer_divide(local_stack, output_stack, b, b_content);
    struct Integer*d = integer_gcd(local_stack, output_stack, a_content, b_content);
    struct Integer*g = &one;
    struct Integer*h = &one;
    struct IntegerPolynomial*remainder;
    while (true)
    {
        struct Integer*degree =
            integer_from_size_t(local_stack, a->coefficient_count - b->coefficient_count);
        remainder = integer_polynomial_euclidean_remainder(local_stack, output_stack,
            integer_polynomial_integer_multiply(local_stack, output_stack, a_primitive_part,
                integer_exponentiate(local_stack, output_stack,
                    (&b_primitive_part->coefficients)[b_primitive_part->coefficient_count - 1],
                    integer_add(local_stack, degree, &one))), b_primitive_part);
        if (remainder->coefficient_count <= 1)
        {
            break;
        }
        a = b;
        b = integer_polynomial_integer_divide(local_stack, output_stack, remainder,
            integer_multiply(local_stack, output_stack, g,
                integer_exponentiate(local_stack, output_stack, h, degree)));
        g = (&a->coefficients)[a->coefficient_count - 1];
        struct Integer*h_exponent = integer_subtract(local_stack, output_stack, &one, degree);
        if (h_exponent->sign >= 0)
        {
            h = integer_multiply(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, h, h_exponent),
                integer_exponentiate(local_stack, output_stack, g, degree));
        }
        else
        {
            h_exponent->sign = -h_exponent->sign;
            h = integer_euclidean_quotient(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, g, degree),
                integer_exponentiate(local_stack, output_stack, h, h_exponent));
        }
    }
    if (remainder->coefficient_count == 1)
    {
        struct IntegerPolynomial*out = polynomial_slot_new(output_stack, 1);
        out->coefficients = integer_copy_to_stack(output_stack, d);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct IntegerPolynomial*out = integer_polynomial_integer_multiply(output_stack, local_stack,
        integer_polynomial_primitive_part(local_stack, output_stack, b), d);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo info;
    integer_extended_gcd(local_stack, output_stack, &info, a, characteristic);
    struct Integer*out;
    if (((struct Integer*)info.a_coefficient)->sign < 0)
    {
        out = integer_add(output_stack, info.a_coefficient, characteristic);
    }
    else
    {
        out = integer_copy_to_stack(output_stack, info.a_coefficient);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_reduced(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    struct IntegerPolynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        (&out->coefficients)[i] = integer_euclidean_remainder(local_stack, output_stack,
            (&a->coefficients)[i], characteristic);
        if ((&out->coefficients)[i]->sign < 0)
        {
            (&out->coefficients)[i] =
                integer_add(local_stack, (&out->coefficients)[i], characteristic);
        }
    }
    polynomial_trim_leading_zeroes(integer_equals_zero, (struct Polynomial*)out);
    polynomial_copy_coefficients(integer_copy_to_stack, output_stack, (struct Polynomial*)out);
    return out;
}

struct IntegerPolynomial*modded_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_add(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_negative(local_stack, a), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_subtract(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_multiply(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_multiply_by_coefficient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        polynomial_multiply_by_coefficient(&integer_operations, local_stack, output_stack,
        (struct Polynomial*)a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RingOperations modded_polynomial_operations = { integer_polynomial_copy,
    polynomial_equals_zero, &polynomial_zero, &integer_polynomial_one, modded_polynomial_add,
    modded_polynomial_negative, modded_polynomial_multiply };

void modded_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    polynomial_euclidean_divide(&integer_operations, modded_integer_reciprocal, local_stack,
        output_stack, out, (struct Polynomial*)dividend, (struct Polynomial*)divisor,
        characteristic);
    out->quotient = (struct Polynomial*)modded_polynomial_reduced(output_stack, local_stack,
        (struct IntegerPolynomial*)out->quotient, characteristic);
    out->remainder = (struct Polynomial*)modded_polynomial_reduced(output_stack, local_stack,
        (struct IntegerPolynomial*)out->remainder, characteristic);
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*modded_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic)
{
    struct PolynomialDivision division;
    modded_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor,
        characteristic);
    return (struct IntegerPolynomial*)division.quotient;
}

struct IntegerPolynomial*modded_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic)
{
    struct PolynomialDivision division;
    modded_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor,
        characteristic);
    return (struct IntegerPolynomial*)division.remainder;
}

struct IntegerPolynomial*modded_polynomial_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*base, struct Integer*exponent,
    struct Integer*characteristic)
{
    return generic_exponentiate(&modded_polynomial_operations, output_stack, local_stack, base,
        exponent, characteristic);
}

struct IntegerPolynomial*modded_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic)
{
    return generic_gcd(&modded_polynomial_operations, modded_polynomial_euclidean_divide,
        output_stack, local_stack, (struct Polynomial*)a, (struct Polynomial*)b, characteristic);
}

void cantor_zassenhaus_split(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, size_t degree,
    struct IntegerPolynomial**factor_cursor)
{
    void*local_stack_savepoint = local_stack->cursor;
    if ((a->coefficient_count - 1) / degree == 1)
    {
        *factor_cursor = modded_polynomial_multiply_by_coefficient(output_stack, local_stack, a,
            modded_integer_reciprocal(local_stack, output_stack,
            (&a->coefficients)[a->coefficient_count - 1], characteristic), characteristic);
        *factor_cursor += 1;
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    struct IntegerPolynomial*b;
    if (integer_equals(characteristic, &INT(2, +)))
    {
        struct IntegerPolynomial*x_squared = polynomial_slot_new(local_stack, 3);
        x_squared->coefficients = &zero;
        (&x_squared->coefficients)[1] = &zero;
        (&x_squared->coefficients)[2] = &one;
        struct IntegerPolynomial*t = &integer_polynomial_one;
        struct IntegerPolynomial*c = &integer_polynomial_one;
        while (true)
        {
            for (size_t i = 1; i < degree; ++i)
            {
                c = modded_polynomial_add(local_stack, output_stack, t,
                    modded_polynomial_euclidean_remainder(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack, c, c, characteristic),
                        a, characteristic), characteristic);
            }
            b = modded_polynomial_gcd(local_stack, output_stack, a, c, characteristic);
            if (b->coefficient_count > 1 && b->coefficient_count != a->coefficient_count)
            {
                break;
            }
            t = modded_polynomial_multiply(local_stack, output_stack, t, x_squared, characteristic);
            c = t;
        }
    }
    else
    {
        struct Integer*p_minus_one = integer_add(local_stack, characteristic, &INT(1, -));
        struct IntegerPolynomial*t = polynomial_slot_new(local_stack, 2 * degree);
        for (size_t i = 0; i < t->coefficient_count - 1; ++i)
        {
            (&t->coefficients)[i] = p_minus_one;
        }
        (&t->coefficients)[t->coefficient_count - 1] = &one;
        struct Integer*power = integer_add(local_stack,
            integer_exponentiate(local_stack, output_stack, characteristic,
                integer_from_size_t(local_stack, degree)),
            &INT(1, -));
        integer_halve(power);
        while (true)
        {
            b = modded_polynomial_gcd(local_stack, output_stack, a,
                modded_polynomial_subtract(local_stack, output_stack,
                    modded_polynomial_exponentiate(local_stack, output_stack, t, power,
                        characteristic), &integer_polynomial_one, characteristic), characteristic);
            if (b->coefficient_count > 1 && b->coefficient_count != a->coefficient_count)
            {
                break;
            }
        }
        size_t i = t->coefficient_count - 2;
        while (true)
        {
            if ((&t->coefficients)[i]->value_count)
            {
                (&t->coefficients)[i] = integer_add(local_stack, (&t->coefficients)[i], &INT(1, -));
                break;
            }
            if (i == 0)
            {
                t->coefficient_count -= 1;
                for (size_t i = 0; i < t->coefficient_count - 1; ++i)
                {
                    (&t->coefficients)[i] = p_minus_one;
                }
                (&t->coefficients)[t->coefficient_count - 1] = &one;
                break;
            }
            (&t->coefficients)[i] = p_minus_one;
            --i;
        }
    }
    cantor_zassenhaus_split(output_stack, local_stack, b, characteristic, degree, factor_cursor);
    cantor_zassenhaus_split(output_stack, local_stack,
        modded_polynomial_euclidean_quotient(local_stack, output_stack, a, b, characteristic),
        characteristic, degree, factor_cursor);
    local_stack->cursor = local_stack_savepoint;
}

void squarefree_modded_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic,
    struct IntegerPolynomial**factor_cursor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*x = polynomial_slot_new(local_stack, 2);
    (&x->coefficients)[0] = &zero;
    (&x->coefficients)[1] = &one;
    struct IntegerPolynomial*v = a;
    struct IntegerPolynomial*w = x;
    size_t d = 0;
    while (2 * d + 2 < a->coefficient_count)
    {
        ++d;
        w = modded_polynomial_euclidean_remainder(local_stack, output_stack,
            modded_polynomial_exponentiate(local_stack, output_stack, w, characteristic,
                characteristic),
            a, characteristic);
        struct IntegerPolynomial*degree_d_factor_product =
            modded_polynomial_gcd(local_stack, output_stack,
                modded_polynomial_subtract(local_stack, output_stack, w, x, characteristic),
                v, characteristic);
        if (degree_d_factor_product->coefficient_count > 1)
        {
            cantor_zassenhaus_split(output_stack, local_stack, degree_d_factor_product,
                characteristic, d, factor_cursor);
            v = modded_polynomial_euclidean_quotient(local_stack, output_stack, v,
                degree_d_factor_product, characteristic);
            w = modded_polynomial_euclidean_remainder(local_stack, output_stack, w, v,
                characteristic);
        }
    }
    if (v->coefficient_count > 1)
    {
        cantor_zassenhaus_split(output_stack, local_stack, v, characteristic,
            v->coefficient_count - 1, factor_cursor);
    }
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*bound_coefficients(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic_power)
{
    struct IntegerPolynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        struct Integer*remainder = integer_euclidean_remainder(local_stack, output_stack,
            (&a->coefficients)[i], characteristic_power);
        if (remainder->sign < 0)
        {
            remainder = integer_add(local_stack, remainder, characteristic_power);
        }
        if (integer_compare(output_stack, local_stack, integer_doubled(local_stack, remainder),
            characteristic_power) > 0)
        {
            (&a->coefficients)[i] =
                integer_subtract(output_stack, local_stack, remainder, characteristic_power);
        }
        else
        {
            (&a->coefficients)[i] = integer_copy_to_stack(output_stack, remainder);
        }
    }
    return out;
}

struct IntegerPolynomial*find_valid_combination(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*combination, size_t combination_size,
    struct IntegerPolynomial*factor_candidates, struct IntegerPolynomial**factor_candidates_end,
    struct IntegerPolynomial**a, struct Integer*characteristic_power)
{
    if (combination_size == 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct IntegerPolynomial*bounded_combination =
            bound_coefficients(local_stack, output_stack, combination, characteristic_power);
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                (&(*a)->coefficients)[(*a)->coefficient_count - 1]), bounded_combination);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            do
            {
                *a = (struct IntegerPolynomial*)division.quotient;
                integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
                    integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                        (&(*a)->coefficients)[(*a)->coefficient_count - 1]),
                    bounded_combination);
            } while (division.remainder && division.remainder->coefficient_count == 0);
            struct IntegerPolynomial*factor =
                integer_polynomial_primitive_part(output_stack, local_stack, bounded_combination);
            local_stack->cursor = local_stack_savepoint;
            return factor;
        }
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        for (struct IntegerPolynomial**candidate = &factor_candidates;
            *candidate <= *factor_candidates_end - combination_size; ++*candidate)
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
                integer_polynomial_multiply(local_stack, output_stack, combination, *candidate),
                combination_size - 1, factor_candidates + 1, factor_candidates_end, a,
                characteristic_power);
            local_stack->cursor = local_stack_savepoint;
            if (factor)
            {
                *factor_candidates_end -= 1;
                *candidate = *factor_candidates_end;
                return factor;
            }
        }
    }
    return 0;
}

struct RationalPolynomial*rational_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_add(&rational_ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_multiply(&rational_ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

void rational_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor)
{
    polynomial_euclidean_divide(&rational_ring_operations, rational_reciprocal, output_stack,
        local_stack, out, (struct Polynomial*)dividend, (struct Polynomial*)divisor, 0);
}

struct IntegerPolynomial*rational_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*denominator_lcm = &one;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        denominator_lcm = integer_multiply(local_stack, output_stack, denominator_lcm,
            (&a->coefficients[i])->denominator);
    }
    struct IntegerPolynomial*integer_polynomial =
        polynomial_slot_new(local_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        (&integer_polynomial->coefficients)[i] = integer_multiply(local_stack, output_stack,
            (&a->coefficients)[i]->numerator, integer_euclidean_quotient(local_stack, output_stack,
                denominator_lcm, (&a->coefficients)[i]->denominator));
    }
    struct IntegerPolynomial*out =
        integer_polynomial_primitive_part(output_stack, local_stack, integer_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*significand, struct Integer*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!(significand->value & 1))
    {
        significand = integer_half(local_stack, significand);
        exponent = integer_add(local_stack, exponent, &one);
    }
    struct Float*out = stack_slot_new(output_stack, sizeof(struct Float), _Alignof(struct Float));
    out->significand = integer_copy_to_stack(output_stack, significand);
    out->exponent = integer_copy_to_stack(output_stack, exponent);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_copy_to_stack(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = stack_slot_new(output_stack, sizeof(struct Float), _Alignof(struct Float));
    out->significand = integer_copy_to_stack(output_stack, a->significand);
    out->exponent = integer_copy_to_stack(output_stack, a->exponent);
    return out;
}

void float_copy_to_pool(struct IntegerPool*pool, struct Float*a)
{
    a->significand = integer_copy_to_pool(pool, a->significand);
    a->exponent = integer_copy_to_pool(pool, a->exponent);
}

void float_move_to_pool(struct IntegerPool*pool, struct Float*a)
{
    integer_move_to_pool(pool, &a->significand);
    integer_move_to_pool(pool, &a->exponent);
}

void float_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack, struct Float*a)
{
    integer_move_from_pool(pool, output_stack, &a->significand);
    integer_move_from_pool(pool, output_stack, &a->exponent);
}

void float_free(struct IntegerPool*pool, struct Float*a)
{
    pool_integer_free(pool, a->significand);
    pool_integer_free(pool, a->exponent);
}

bool float_equals_zero(struct Float*a)
{
    return integer_equals_zero(a->significand) == 0;
}

struct Float*float_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*exponent_difference =
        integer_subtract(local_stack, output_stack, a->exponent, b->exponent);
    struct Float*out;
    if (exponent_difference->sign > 0)
    {
        out = stack_slot_new(output_stack, sizeof(struct Float), _Alignof(struct Float));
        out->significand = integer_add(output_stack, b->significand,
            integer_multiply(local_stack, output_stack, a->significand,
                integer_exponentiate(local_stack, output_stack, &INT(2, +),
                    exponent_difference)));
        out->exponent = integer_copy_to_stack(output_stack, b->exponent);
    }
    else if (exponent_difference->sign < 0)
    {
        out = float_add(output_stack, local_stack, b, a);
    }
    else
    {
        out = float_reduced(output_stack, local_stack,
            integer_add(local_stack, a->significand, b->significand), a->exponent);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_generic_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b, void*unused)
{
    return float_add(output_stack, local_stack, a, b);
}

struct Float*float_negative(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = stack_slot_new(output_stack, sizeof(struct Float), _Alignof(struct Float));
    out->significand = integer_negative(output_stack, a->significand);
    out->exponent = integer_copy_to_stack(output_stack, a->exponent);
    return out;
}

struct Float*float_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Float*a, void*unused)
{
    return float_negative(output_stack, a);
}

struct Float*float_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*minuend, struct Float*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = stack_slot_new(output_stack, sizeof(struct Float), _Alignof(struct Float));
    out = float_add(output_stack, local_stack, minuend, float_negative(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = float_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->significand, b->significand),
        integer_add(local_stack, a->exponent, b->exponent));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a, struct Float*b, void*unused)
{
    return float_multiply(output_stack, local_stack, a, b);
}

struct RingOperations float_ring_operations = { float_copy_to_stack, float_equals_zero,
    &float_zero, &float_one, float_generic_add, float_generic_negative, float_generic_multiply };

struct Float*float_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*base, struct Integer*exponent)
{
    return generic_exponentiate(&float_ring_operations, output_stack, local_stack, base, exponent,
        0);
}

struct Rational*float_to_rational(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = stack_slot_new(output_stack, sizeof(struct Rational*),
        _Alignof(struct Rational*));
    if (a->exponent->sign < 0)
    {
        out->numerator = integer_copy_to_stack(output_stack, a->significand);
        out->denominator = integer_exponentiate(output_stack, local_stack, &INT(2, +),
            integer_magnitude(local_stack, a->exponent));
    }
    else
    {
        out->numerator = integer_multiply(output_stack, local_stack,
            integer_exponentiate(local_stack, output_stack, &INT(2, +), a->exponent),
            a->significand);
        out->denominator = integer_stack_slot_new(output_stack, 1);
        out->denominator->value_count = 1;
        out->denominator->value = 1;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void float_interval_move_to_pool(struct IntegerPool*pool, struct FloatInterval*a)
{
    float_move_to_pool(pool, &a->min);
    float_move_to_pool(pool, &a->max);
}

void float_interval_move_from_pool(struct IntegerPool*pool, struct Stack*output_stack,
    struct FloatInterval*a)
{
    float_move_from_pool(pool, output_stack, &a->min);
    float_move_from_pool(pool, output_stack, &a->max);
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    out->min = float_to_rational(output_stack, local_stack, &a->min);
    out->max = float_to_rational(output_stack, local_stack, &a->max);
}

struct Number*number_slot_new(struct Stack*number_stack, struct Number**free_list)
{
    struct Number*out = *free_list;
    if (out)
    {
        *free_list = out->next;
    }
    else
    {
        out = stack_slot_new(number_stack, sizeof(struct Number), _Alignof(struct Number));
    }
    memset(out, 0, sizeof(struct Number));
    return out;
}

void number_node_free(struct Number**free_list, struct IntegerPool*integer_pool, struct Number*a)
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
    if (a->operation != 'r')
    {
        number_free(free_list, integer_pool, a->left);
        number_free(free_list, integer_pool, a->right);
    }
    number_node_free(free_list, integer_pool, a);
}

struct Number*number_copy(struct Stack*number_stack, struct Number**free_list,
    struct IntegerPool*integer_pool, struct Number*a)
{
    struct Number*out = number_slot_new(number_stack, free_list);
    memcpy(out, a, sizeof(struct Number));
    if (a->operation == 'r')
    {
        rational_move_to_pool(integer_pool, &out->value);
    }
    else
    {
        out->left = number_copy(number_stack, free_list, integer_pool, a->left);
        out->right = number_copy(number_stack, free_list, integer_pool, a->right);
        float_interval_move_to_pool(integer_pool, &a->imaginary_part_estimate);
    }
    float_interval_move_to_pool(integer_pool, &a->real_part_estimate);
    float_interval_move_to_pool(integer_pool, &a->argument_estimate);
    float_interval_move_to_pool(integer_pool, &a->magnitude_estimate);
    return out;
}

bool get_input(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    void*stack_a_savepoint = stack_a->cursor;
    struct Number*previous = 0;
    while (true)
    {
        struct Number*number = number_slot_new(number_stack, number_free_list);
        number->previous = previous;
        number->next = number_stack->cursor;
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.denominator = pool_integer_new(integer_pool, 1, 1);
            number->value.numerator = integer_from_char(stack_a, next_char);
            next_char = getchar();
            while (isdigit(next_char))
            {
                number->value.numerator =
                    integer_multiply(stack_a, stack_b, number->value.numerator, &INT(10, +));
                number->value.numerator = integer_add(stack_a, number->value.numerator,
                    integer_from_char(stack_a, next_char));
                next_char = getchar();
            }
            integer_move_to_pool(integer_pool, &number->value.numerator);
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
                stack_a->cursor = stack_a_savepoint;
                return false;
            }
            number->operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            number->next = 0;
            stack_a->cursor = stack_a_savepoint;
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
                number_node_free(number_free_list, integer_pool, input);
                number_node_free(number_free_list, integer_pool, nested_number);
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
            number_node_free(number_free_list, integer_pool, input);
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
                input->value.numerator = pool_integer_new(integer_pool, 1, -1);
                input->value.denominator = pool_integer_new(integer_pool, 1, 1);
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
                number_node_free(number_free_list, integer_pool, input->next);
                number_node_free(number_free_list, integer_pool, input);
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

bool estimate_needs_refinement(struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*estimate, struct Rational*interval_size)
{
    if (estimate->min.significand)
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Float*current_interval_size =
            float_subtract(stack_a, stack_b, &estimate->max, &estimate->min);
        struct Rational*rational_current_interval_size =
            float_to_rational(stack_a, stack_b, current_interval_size);
        if (rational_compare(stack_a, stack_b, interval_size, rational_current_interval_size) >= 0)
        {
            stack_a->cursor = stack_a_savepoint;
            return false;
        }
        stack_a->cursor = stack_a_savepoint;
    }
    return true;
}

//Leaves garbage allocations on local_stack along with the final values of estimate_denominator and
//estimate_remainder. The values of out_min and out_max go on output_stack.
void rational_continue_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*out_min, struct Float*out_max, struct Integer**estimate_denominator,
    struct Integer**estimate_remainder, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (integer_compare(output_stack, local_stack, interval_size->denominator,
        integer_multiply(local_stack, output_stack, *estimate_denominator,
            interval_size->numerator)) > 0)
    {
        if ((*estimate_remainder)->value_count == 0)
        {
            out_min = float_copy_to_stack(output_stack, out_min);
            *out_max = *out_min;
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division,
            integer_doubled(local_stack, *estimate_remainder), a->denominator);
        *estimate_denominator = integer_doubled(local_stack, *estimate_denominator);
        *estimate_remainder = division.remainder;
        out_min->significand = integer_add(local_stack, division.quotient,
            integer_doubled(local_stack, out_min->significand));
        out_min->exponent = integer_add(local_stack, out_min->exponent, &INT(1, -));
    }
    out_min = float_reduced(output_stack, local_stack, out_min->significand, out_min->exponent);
    out_max = float_reduced(output_stack, local_stack,
        integer_add(local_stack, out_min->significand, &one), out_min->exponent);
    local_stack->cursor = local_stack_savepoint;
}

void rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*out_min, struct Float*out_max, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerDivision division;
    integer_euclidean_divide(local_stack, output_stack, &division,
        integer_magnitude(local_stack, a->numerator), a->denominator);
    struct Integer*estimate_remainder = division.remainder;
    struct Integer*estimate_denominator = &one;
    if (a->numerator->sign < 0)
    {
        out_max->significand = division.quotient;
        out_max->exponent = &zero;
        rational_continue_float_estimate(output_stack, local_stack, out_max, out_min,
            &estimate_denominator, &estimate_remainder, a, interval_size);
        out_min->significand->sign = -1;
        out_max->significand->sign = -1;
    }
    else
    {
        out_min->significand = division.quotient;
        out_min->exponent = &zero;
        rational_continue_float_estimate(output_stack, local_stack, out_min, out_max,
            &estimate_denominator, &estimate_remainder, a, interval_size);
    }
    local_stack->cursor = local_stack_savepoint;
}

void float_estimate_root(struct Stack*output_stack, struct Stack*local_stack, struct Float*out_min,
    struct Float*out_max, struct Float*a, struct Rational*interval_size, struct Integer*index)
{
    ASSERT(a->significand->sign >= 0, "float_estimate_root was called on a negative a value.");
    void*local_stack_savepoint = local_stack->cursor;
    if (a->exponent < 0)
    {
        out_max->significand = &one;
        out_max->exponent = &zero;
    }
    else
    {
        out_max = float_copy_to_stack(local_stack, a);
    }
    struct Rational*rational_radicand = float_to_rational(local_stack, output_stack, a);
    struct Integer*index_minus_one = integer_add(local_stack, index, &INT(1, -));
    while (true)
    {
        struct Rational*delta = rational_integer_divide(local_stack, output_stack,
            rational_subtract(local_stack, output_stack,
                rational_divide(local_stack, output_stack, rational_radicand,
                    float_to_rational(local_stack, output_stack,
                        float_exponentiate(local_stack, output_stack, out_max, index_minus_one))),
                float_to_rational(local_stack, output_stack, out_max)), index);
        struct FloatInterval delta_float_estimate;
        rational_float_estimate(local_stack, output_stack, &delta_float_estimate.min,
            &delta_float_estimate.max, delta,
            rational_integer_divide(local_stack, output_stack, delta, &INT(2, -)));
        out_max = float_add(local_stack, output_stack, out_max, &delta_float_estimate.max);
        struct Rational*interval_size_comparison_value =
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, &delta_float_estimate.max),
                rational_reduced(local_stack, output_stack,
                    integer_doubled(local_stack, delta->numerator), delta->denominator));
        if (rational_compare(output_stack, local_stack,
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, &delta_float_estimate.max),
                rational_reduced(local_stack, output_stack,
                    integer_doubled(local_stack, delta->numerator), delta->denominator)),
            interval_size) <= 0)
        {
            out_min = float_add(output_stack, local_stack, out_max, &delta_float_estimate.max);
            out_max = float_copy_to_stack(output_stack, out_max);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
    }
}

struct FloatInterval*number_float_magnitude_estimate(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(stack_a, stack_b, &a->magnitude_estimate, interval_size))
    {
        return &a->magnitude_estimate;
    }
    void*stack_a_savepoint = stack_a->cursor;
    if (a->operation == 'r')
    {
        if (!a->magnitude_estimate.min.significand)
        {
            struct IntegerDivision division;
            integer_euclidean_divide(stack_a, stack_b, &division,
                integer_magnitude(stack_a, a->value.numerator), a->value.denominator);
            a->magnitude_estimate_denominator = &one;
            a->magnitude_estimate.min.significand = division.quotient;
            a->magnitude_estimate.min.exponent = &zero;
            a->magnitude_estimate_remainder = division.remainder;
        }
        else
        {
            float_interval_move_from_pool(integer_pool, stack_a, &a->magnitude_estimate);
            integer_move_from_pool(integer_pool, stack_a, &a->magnitude_estimate_denominator);
            integer_move_from_pool(integer_pool, stack_a, &a->magnitude_estimate_remainder);
        }
        rational_continue_float_estimate(stack_a, stack_b, &a->magnitude_estimate.min,
            &a->magnitude_estimate.max, &a->magnitude_estimate_denominator,
            &a->magnitude_estimate_remainder, &a->value, interval_size);
        float_interval_move_to_pool(integer_pool, &a->magnitude_estimate);
        integer_move_to_pool(integer_pool, &a->magnitude_estimate_denominator);
        integer_move_to_pool(integer_pool, &a->magnitude_estimate_remainder);
        stack_a->cursor = stack_a_savepoint;
        return &a->magnitude_estimate;
    }
    if (a->magnitude_estimate.min.significand)
    {
        float_free(integer_pool, &a->magnitude_estimate.min);
        float_free(integer_pool, &a->magnitude_estimate.max);
    }
    if (a->operation == '^')
    {
        struct FloatInterval*radicand_magnitude_estimate =
            number_float_magnitude_estimate(integer_pool, stack_a, stack_b, a->left,
                rational_integer_divide(stack_a, stack_b,
                    rational_exponentiate(stack_a, stack_b, interval_size,
                        a->right->value.denominator),
                    &INT(3, +)));
        struct Rational*bound_interval_size =
            rational_integer_divide(stack_a, stack_b, interval_size, &INT(3, +));
        struct Float unused_estimate_bound;
        float_estimate_root(stack_a, stack_b, &a->magnitude_estimate.min, &unused_estimate_bound,
            &radicand_magnitude_estimate->min, bound_interval_size, a->right->value.denominator);
        float_estimate_root(stack_a, stack_b, &unused_estimate_bound, &a->magnitude_estimate.max,
            &radicand_magnitude_estimate->max, bound_interval_size, a->right->value.denominator);
        float_interval_move_to_pool(integer_pool, &a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return &a->magnitude_estimate;
    }
    ABORT("number_float_magnitude_estimate case not yet implemented.");
}

void number_rational_magnitude_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    float_interval_to_rational_interval(output_stack, local_stack, out,
        number_float_magnitude_estimate(integer_pool, output_stack, local_stack, a, interval_size));
}

void number_float_estimate_from_rational(void(get_rational_estimate)(struct IntegerPool*,
    struct Stack*, struct Stack*, struct RationalInterval*, struct Number*, struct Rational*),
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b,
    struct Float*out_min, struct Float*out_max, struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*rational_estimate_interval_size = &rational_one;
    while (rational_compare(stack_a, stack_b, rational_estimate_interval_size, interval_size) > 0)
    {
        rational_estimate_interval_size->denominator =
            integer_doubled(stack_a, rational_estimate_interval_size->denominator);
    }
    rational_estimate_interval_size->denominator =
        integer_doubled(stack_a, rational_estimate_interval_size->denominator);
    struct RationalInterval rational_estimate;
    get_rational_estimate(integer_pool, stack_a, stack_b, &rational_estimate, a,
        rational_estimate_interval_size);
    struct Float unused_estimate_bound;
    rational_float_estimate(stack_a, stack_b, out_min, &unused_estimate_bound,
        rational_estimate.min, rational_estimate_interval_size);
    rational_float_estimate(stack_a, stack_b, &unused_estimate_bound, out_max,
        rational_estimate.max, rational_estimate_interval_size);
    float_move_to_pool(integer_pool, out_min);
    float_move_to_pool(integer_pool, out_max);
    stack_a->cursor = stack_a_savepoint;
}

struct FloatInterval*number_float_argument_estimate(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Rational*interval_size);

void number_rational_argument_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->sign < 0)
        {
            pi_estimate(output_stack, local_stack, interval_size);
            out->min = &pi_estimate_min;
            out->max = rational_add(output_stack, local_stack, &pi_estimate_min, &pi_interval_size);
        }
        else
        {
            out->min->numerator = &zero;
            out->min->denominator = &one;
            out->max->numerator = &zero;
            out->max->denominator = &one;
        }
        return;
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval radicand_rational_argument_estimate;
        number_rational_argument_estimate(integer_pool, local_stack, output_stack,
            &radicand_rational_argument_estimate, a->left,
            rational_integer_multiply(local_stack, output_stack, interval_size,
                a->right->value.denominator));
        out->min = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate.min, a->right->value.denominator);
        out->max = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate.max, a->right->value.denominator);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        float_interval_to_rational_interval(output_stack, local_stack, out,
            number_float_argument_estimate(integer_pool, output_stack, local_stack, a,
                interval_size));
        local_stack->cursor = local_stack_savepoint;
    }
    case '+':
        ABORT("number_rational_argument_estimate case not yet implemented.");
    }
}

void number_estimate_part_sum(struct FloatInterval*(estimate_term)(struct IntegerPool*,
    struct Stack*, struct Stack*, struct Number*, struct Rational*),
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b,
    struct Float*out_min, struct Float*out_max, struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*term_estimate_interval_size =
        rational_integer_divide(stack_a, stack_b, interval_size, &INT(2, +));
    struct FloatInterval*left_term_estimate =
        estimate_term(integer_pool, stack_a, stack_b, a->left, term_estimate_interval_size);
    struct FloatInterval*right_term_estimate =
        estimate_term(integer_pool, stack_a, stack_b, a->right, term_estimate_interval_size);
    out_min = float_add(stack_a, stack_b, &left_term_estimate->min, &right_term_estimate->min);
    out_max = float_add(stack_a, stack_b, &left_term_estimate->max, &right_term_estimate->max);
    float_move_to_pool(integer_pool, out_min);
    float_move_to_pool(integer_pool, out_max);
}

struct FloatInterval*float_argument_estimate_for_part_sum(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    return number_float_argument_estimate(integer_pool, stack_a, stack_b, a, interval_size);
}

struct FloatInterval*number_float_argument_estimate(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(stack_a, stack_b, &a->argument_estimate, interval_size))
    {
        return &a->magnitude_estimate;
    }
    if (a->operation == '*')
    {
        number_estimate_part_sum(float_argument_estimate_for_part_sum, integer_pool, stack_a,
            stack_b, &a->argument_estimate.min, &a->argument_estimate.max, a, interval_size);
        return &a->argument_estimate;
    }
    else
    {
        number_float_estimate_from_rational(number_rational_argument_estimate, integer_pool,
            stack_a, stack_b, &a->argument_estimate.min, &a->argument_estimate.max, a,
            interval_size);
        return &a->argument_estimate;
    }
}

void number_argument_cosine_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(integer_pool, local_stack, output_stack, &argument_estimate,
        a, argument_estimate_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_min_cosine, argument_estimate.min,
        argument_estimate_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_max_cosine, argument_estimate.max,
        argument_estimate_interval_size);
    out->max = rational_copy_to_stack(output_stack,
        rational_max(output_stack, local_stack, argument_min_cosine.max, argument_max_cosine.max));
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack, argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack, argument_estimate.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate.min, &pi_estimate_min) <= 0 &&
        rational_compare(output_stack, local_stack, pi_estimate_max, argument_estimate.max) <= 0)
    {
        out->min->numerator = stack_integer_new(output_stack, 1, -1);
        out->min->denominator = &one;
    }
    else
    {
        out->min = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
            argument_min_cosine.min, argument_max_cosine.min));
        out->min = rational_copy_to_stack(output_stack, out->min);
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_argument_sine_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(integer_pool, local_stack, output_stack, &argument_estimate,
        a, argument_estimate_interval_size);
    struct RationalInterval argument_min_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_min_sine, argument_estimate.min,
        argument_estimate_interval_size);
    struct RationalInterval argument_max_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_max_sine, argument_estimate.max,
        argument_estimate_interval_size);
    struct RationalInterval argument_estimate_multiple =
        { rational_doubled(local_stack, output_stack, argument_estimate.min),
        rational_doubled(local_stack, output_stack, argument_estimate.max) };
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
        argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
        argument_estimate_multiple.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
        &pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
            argument_estimate_multiple.max) <= 0)
    {
        out->max->numerator = stack_integer_new(output_stack, 1, 1);
        out->max->denominator = stack_integer_new(output_stack, 1, 1);
        out->min = rational_copy_to_stack(output_stack,
            rational_max(output_stack, local_stack, argument_min_sine.min, argument_max_sine.min));
    }
    else
    {
        argument_estimate_multiple.min = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.min, &INT(3, +));
        argument_estimate_multiple.max = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.max, &INT(3, +));
        pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
            argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
            argument_estimate_multiple.max);
        pi_estimate_max =
            rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
        if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
            &pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
                argument_estimate_multiple.max) <= 0)
        {
            out->min->numerator = stack_integer_new(output_stack, 1, -1);
            out->min->denominator = stack_integer_new(output_stack, 1, 1);
            out->max = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
        else
        {
            out->min = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.min, argument_max_sine.min));
            out->max = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_rectangular_part_from_polar_form(void(trig_function)(struct IntegerPool*, struct Stack*,
    struct Stack*, struct RationalInterval*, struct Number*, struct Rational*),
    struct IntegerPool*integer_pool, struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval magnitude_estimate;
    number_rational_magnitude_estimate(integer_pool, local_stack, output_stack, &magnitude_estimate,
        a, &rational_one);
    struct Rational*factor_interval_size = rational_divide(local_stack, output_stack, interval_size,
        rational_integer_add(local_stack, output_stack, magnitude_estimate.max, &INT(2, +)));
    number_rational_magnitude_estimate(integer_pool, local_stack, output_stack, &magnitude_estimate,
        a, factor_interval_size);
    struct RationalInterval trig_value;
    trig_function(integer_pool, local_stack, output_stack, &trig_value, a, factor_interval_size);
    if (trig_value.min->numerator->sign >= 0)
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.min);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.max);
    }
    else if (trig_value.max->numerator->sign <= 0)
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.max);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.min);
    }
    else
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.max);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.max);
    }
    local_stack->cursor = local_stack_savepoint;
}

struct FloatInterval*number_float_real_part_estimate(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size);

void number_rational_real_part_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_cosine_estimate, integer_pool,
            output_stack, local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval*float_estimate = number_float_real_part_estimate(integer_pool,
            local_stack, output_stack, a, interval_size);
        out->min = float_to_rational(output_stack, local_stack, &float_estimate->min);
        out->max = float_to_rational(output_stack, local_stack, &float_estimate->max);
        local_stack->cursor = local_stack_savepoint;
    }
}

//When a takes the union variant with the real_part_estimate field, the return value is a pointer to
//that field, whose Integer components are allocated in integer_pool. Otherwise, the return value
//and its Integer components are allocated on output_stack.
struct FloatInterval*number_float_real_part_estimate(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*magnitude_estimate = number_float_magnitude_estimate(integer_pool,
            output_stack, local_stack, a, interval_size);
        if (a->value.numerator->sign < 0)
        {
            struct FloatInterval*real_part_estimate = stack_slot_new(output_stack,
                sizeof(struct FloatInterval), _Alignof(struct FloatInterval));
            real_part_estimate->min = *float_negative(output_stack, &magnitude_estimate->max);
            real_part_estimate->max = *float_negative(output_stack, &magnitude_estimate->min);
            return real_part_estimate;
        }
        else
        {
            return magnitude_estimate;
        }
    }
    if (!estimate_needs_refinement(output_stack, local_stack, &a->real_part_estimate,
        interval_size))
    {
        return &a->real_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
        number_float_estimate_from_rational(number_rational_real_part_estimate, integer_pool,
            output_stack, local_stack, &a->real_part_estimate.min, &a->real_part_estimate.max, a,
            interval_size);
        return &a->real_part_estimate;
    case '*':
        ABORT("number_float_real_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_real_part_estimate, integer_pool, output_stack,
            local_stack, &a->real_part_estimate.min, &a->real_part_estimate.max, a, interval_size);
        return &a->real_part_estimate;
    default:
        ABORT("Number operation not recognized.");
    }
}

struct FloatInterval*number_float_imaginary_part_estimate(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size);

void number_rational_imaginary_part_estimate(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct RationalInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_sine_estimate, integer_pool,
            output_stack, local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval*float_estimate = number_float_imaginary_part_estimate(integer_pool,
            local_stack, output_stack, a, interval_size);
        out->min = float_to_rational(output_stack, local_stack, &float_estimate->min);
        out->max = float_to_rational(output_stack, local_stack, &float_estimate->max);
        local_stack->cursor = local_stack_savepoint;
    }
}

//When a takes the union variant with the imaginary_part_estimate field, the return value is a
//pointer to that field, whose Integer components are allocated in integer_pool. Otherwise, the
//return value and its Integer components are allocated on output_stack.
struct FloatInterval*number_float_imaginary_part_estimate(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*imaginary_part_estimate = stack_slot_new(output_stack,
            sizeof(struct FloatInterval), _Alignof(struct FloatInterval));
        imaginary_part_estimate->min.significand = &zero;
        imaginary_part_estimate->min.exponent = &zero;
        imaginary_part_estimate->max.significand = &zero;
        imaginary_part_estimate->max.exponent = &zero;
        return imaginary_part_estimate;
    }
    if (!estimate_needs_refinement(output_stack, local_stack, &a->imaginary_part_estimate,
        interval_size))
    {
        return &a->imaginary_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
        number_float_estimate_from_rational(number_rational_imaginary_part_estimate, integer_pool,
            output_stack, local_stack, &a->imaginary_part_estimate.min,
            &a->imaginary_part_estimate.max, a, interval_size);
        return &a->imaginary_part_estimate;
    case '*':
        ABORT("number_float_imaginary_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_imaginary_part_estimate, integer_pool, output_stack,
            local_stack, &a->imaginary_part_estimate.min, &a->imaginary_part_estimate.max, a,
            interval_size);
        return &a->imaginary_part_estimate;
    default:
        ABORT("Number operation not recognized.");
    }
}

struct RationalPoynomial*number_minimal_polynomial(struct Stack*polynomial_stack,
    struct IntegerPool*integer_pool, struct Number*a)
{
    if (a->minimal_polynomial)
    {
        return a->minimal_polynomial;
    }
    switch (a->operation)
    {
    case 'r':
        a->minimal_polynomial = polynomial_slot_new(polynomial_stack, 2);
        a->minimal_polynomial->coefficients->numerator =
            integer_copy_to_pool(integer_pool, a->value.numerator);
        a->minimal_polynomial->coefficients->numerator->sign =
            -a->minimal_polynomial->coefficients->numerator->sign;
        a->minimal_polynomial->coefficients->denominator =
            integer_copy_to_pool(integer_pool, a->value.denominator);
        (&a->minimal_polynomial->coefficients)[1]->numerator = pool_integer_new(integer_pool, 1, 1);
        (&a->minimal_polynomial->coefficients)[1]->denominator =
            pool_integer_new(integer_pool, 1, 1);
        return a->minimal_polynomial;
    default:
        ABORT("number_minimal_polynomial case not yet implemented.");
    }
}

struct Number*number_add(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
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
            void*stack_a_savepoint = stack_a->cursor;
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = 'r';
            struct Rational*sum = rational_add(stack_a, stack_b, &a->value, &b->value);
            if (!sum)
            {
                return 0;
            }
            out->value = *sum;
            rational_move_to_pool(integer_pool, &out->value);
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            stack_a->cursor = stack_a_savepoint;
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
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
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
            out->value.denominator = pool_integer_new(integer_pool, 1, 1);
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            return out;
        }
        if (integer_equals(a->value.numerator, &one) && integer_equals(a->value.denominator, &one))
        {
            number_free(number_free_list, integer_pool, a);
            return b;
        }
        if (b->operation == 'r')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct Number*out = number_slot_new(number_stack, number_free_list);
            out->operation = 'r';
            struct Rational*product = rational_multiply(stack_a, stack_b, &a->value, &b->value);
            if (!product)
            {
                return 0;
            }
            out->value = *product;
            rational_move_to_pool(integer_pool, &out->value);
            number_free(number_free_list, integer_pool, a);
            number_free(number_free_list, integer_pool, b);
            stack_a->cursor = stack_a_savepoint;
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
                number_node_free(number_free_list, integer_pool, a);
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
                stack_a, stack_b, a->left, b);
            if (!out)
            {
                return 0;
            }
            out = number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                out, a->right);
            number_node_free(number_free_list, integer_pool, a);
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
            number_node_free(number_free_list, integer_pool, a);
            return out;
        }
        struct Number*left = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
            stack_b, a->left, b);
        if (!left)
        {
            return 0;
        }
        struct Number*right = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
            stack_b, a->right, b_copy);
        if (!right)
        {
            return 0;
        }
        number_node_free(number_free_list, integer_pool, a);
        return number_add(number_stack, number_free_list, integer_pool, stack_a, stack_b, left,
            right);
    }
    }
    return number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b, b, a);
}

struct Number*number_reciprocal(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a);

struct Number*number_divide(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*dividend, struct Number*divisor)
{
    return number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b, dividend,
        number_reciprocal(number_stack, number_free_list, integer_pool, stack_a, stack_b, divisor));
}

struct Number*number_exponentiate(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*base,
    struct Number*exponent);

struct Number*number_reciprocal(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        if (!a->value.numerator->value_count)
        {
            printf("Tried to divide by 0.");
            return 0;
        }
        struct Integer*old_numerator = a->value.numerator;
        a->value.numerator = a->value.denominator;
        a->value.denominator = old_numerator;
        return a;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        pool_integer_free(integer_pool, a->right->value.numerator);
        a->right->value.numerator = integer_add(stack_a, a->right->value.denominator, &INT(1, -));
        integer_move_to_pool(integer_pool, &a->right->value.numerator);
        stack_a->cursor = stack_a_savepoint;
        struct Number*radicand_copy =
            number_copy(number_stack, number_free_list, integer_pool, a->left);
        return number_divide(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                a->left, a->right),
            radicand_copy);
    }
    case '*':
    {
        struct Number*out =
            number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                number_reciprocal(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                    a->left),
                number_reciprocal(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                    a->right));
        number_node_free(number_free_list, integer_pool, a);
        return out;
    }
    case '+':
        ABORT("number_reciprocal case not yet implemented.");
    default:
        ABORT("Number operation not recognized.");
    }
}

struct Number*number_exponentiate(struct Stack*number_stack, struct Number**number_free_list,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*base,
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
        base =
            number_reciprocal(number_stack, number_free_list, integer_pool, stack_a, stack_b, base);
        if (!base)
        {
            return 0;
        }
        exponent->value.numerator->sign = 1;
        return number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            base, exponent);
    }
    switch (base->operation)
    {
    case 'r':
        if (base->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_stack, number_free_list);
            number_free(number_free_list, integer_pool, base);
            number_free(number_free_list, integer_pool, exponent);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = pool_integer_new(integer_pool, 1, 1);
            return out;
        }
        void*stack_a_savepoint = stack_a->cursor;
        if (integer_equals(base->value.denominator, &one))
        {
            if (integer_equals(exponent->value.denominator, &one))
            {
                struct Number*out = number_slot_new(number_stack, number_free_list);
                out->operation = 'r';
                out->value.denominator = pool_integer_new(integer_pool, 1, 1);
                out->value.numerator = integer_exponentiate(stack_a, stack_b, base->value.numerator,
                    exponent->value.numerator);
                integer_move_to_pool(integer_pool, &out->value.numerator);
                number_free(number_free_list, integer_pool, base);
                number_free(number_free_list, integer_pool, exponent);
                stack_a->cursor = stack_a_savepoint;
                return out;
            }
            void*stack_b_savepoint = stack_b->cursor;
            struct Integer*radicand = integer_exponentiate(stack_a, stack_b, base->value.numerator,
                exponent->value.numerator);
            pool_integer_free(integer_pool, exponent->value.numerator);
            exponent->value.numerator = pool_integer_new(integer_pool, 1, 1);
            number_free(number_free_list, integer_pool, base);
            int8_t radicand_sign = radicand->sign;
            radicand->sign = 1;
            size_t factor_count = 0;
            struct Factor*factors = stack_b->cursor;
            struct Integer*prime = primes;
            while (integer_compare(stack_a, stack_b, prime, radicand) <= 0)
            {
                struct IntegerDivision division;
                integer_euclidean_divide(stack_a, stack_b, &division, radicand, prime);
                if (division.remainder->value_count == 0)
                {
                    factor_count += 1;
                    struct Factor*factor =
                        stack_slot_new(stack_b, sizeof(struct Factor), _Alignof(struct Factor));
                    factor->value = prime;
                    factor->multiplicity = &zero;
                    while (division.remainder->value_count == 0)
                    {
                        factor->multiplicity = integer_add(stack_a, factor->multiplicity, &one);
                        radicand = division.quotient;
                        integer_euclidean_divide(stack_a, stack_b, &division, radicand,
                            factor->value);
                    }
                }
                prime = next_prime(stack_a, stack_b, prime);
            }
            struct Integer*coefficient = &one;
            struct Integer*multiplicity_gcd = exponent->value.denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct IntegerDivision multiplicity_reduction;
                integer_euclidean_divide(stack_a, stack_b, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->value.denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(stack_a, stack_b, coefficient,
                    integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = integer_gcd(stack_a, stack_b, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            if (radicand_sign > 0)
            {
                struct Integer*reduced_degree = integer_euclidean_quotient(stack_a, stack_b,
                    exponent->value.denominator, multiplicity_gcd);
                if (integer_equals(reduced_degree, &one))
                {
                    number_free(number_free_list, integer_pool, exponent);
                    struct Number*out = number_slot_new(number_stack, number_free_list);
                    out->operation = 'r';
                    out->value.denominator = pool_integer_new(integer_pool, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move_to_pool(integer_pool, &out->value.numerator);
                    stack_a->cursor = stack_a_savepoint;
                    stack_b->cursor = stack_b_savepoint;
                    return out;
                }
                pool_integer_free(integer_pool, exponent->value.denominator);
                exponent->value.denominator = reduced_degree;
                integer_move_to_pool(integer_pool, &exponent->value.denominator);
            }
            struct Number*number_coefficient = number_slot_new(number_stack, number_free_list);
            number_coefficient->operation = 'r';
            number_coefficient->value.numerator = coefficient;
            integer_move_to_pool(integer_pool, &number_coefficient->value.numerator);
            number_coefficient->value.denominator = pool_integer_new(integer_pool, 1, 1);
            struct Number*number_radicand = number_slot_new(number_stack, number_free_list);
            number_radicand->operation = 'r';
            number_radicand->value.denominator = pool_integer_new(integer_pool, 1, 1);
            number_radicand->value.numerator = stack_integer_new(stack_a, 1, radicand_sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(stack_a, stack_b,
                    factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(stack_a, stack_b,
                    factors[factor_index].value, reduced_multiplicity);
                number_radicand->value.numerator = integer_multiply(stack_a, stack_b,
                    number_radicand->value.numerator, exponentiation);
            }
            integer_move_to_pool(integer_pool, &number_radicand->value.numerator);
            struct Number*surd = number_slot_new(number_stack, number_free_list);
            surd->operation = '^';
            surd->left = number_radicand;
            surd->right = exponent;
            stack_a->cursor = stack_a_savepoint;
            stack_b->cursor = stack_b_savepoint;
            return number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                number_coefficient, surd);
        }
        struct Number*new_denominator = number_slot_new(number_stack, number_free_list);
        new_denominator->operation = 'r';
        new_denominator->value.denominator = pool_integer_new(integer_pool, 1, 1);
        new_denominator->value.numerator = integer_exponentiate(stack_a, stack_b,
            base->value.denominator, exponent->value.numerator);
        integer_move_to_pool(integer_pool, &new_denominator->value.numerator);
        struct Number*new_numerator_base = number_slot_new(number_stack, number_free_list);
        new_numerator_base->operation = 'r';
        new_numerator_base->value.denominator = pool_integer_new(integer_pool, 1, 1);
        new_numerator_base->value.numerator =
            integer_multiply(stack_a, stack_b, base->value.numerator,
                integer_exponentiate(stack_a, stack_b, base->value.denominator,
                    integer_add(stack_a, exponent->value.denominator, &INT(1, -))));
        integer_move_to_pool(integer_pool, &new_numerator_base->value.numerator);
        struct Number*new_numerator = number_exponentiate(number_stack, number_free_list,
            integer_pool, stack_a, stack_b, new_numerator_base, exponent);
        number_free(number_free_list, integer_pool, new_numerator_base);
        new_denominator = number_reciprocal(number_stack, number_free_list, integer_pool, stack_a,
            stack_b, new_denominator);
        struct Number*out = number_multiply(number_stack, number_free_list, integer_pool, stack_a,
            stack_b, new_numerator, new_denominator);
        number_free(number_free_list, integer_pool, new_numerator);
        number_free(number_free_list, integer_pool, new_denominator);
        stack_a->cursor = stack_a_savepoint;
        return out;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        rational_move_from_pool(integer_pool, stack_a, &exponent->value);
        exponent->value =
            *rational_multiply(stack_a, stack_b, &exponent->value, &base->right->value);
        rational_move_to_pool(integer_pool, &exponent->value);
        stack_a->cursor = stack_a_savepoint;
        number_free(number_free_list, integer_pool, base->right);
        struct Number*out = number_exponentiate(number_stack, number_free_list, integer_pool,
            stack_a, stack_b, base->left, exponent);
        number_node_free(number_free_list, integer_pool, base);
        return out;
    }
    case '*':
    {
        struct Number*exponent_copy =
            number_copy(number_stack, number_free_list, integer_pool, exponent);
        struct Number*out =
            number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                    base->left, exponent),
                number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
                    base->right, exponent_copy));
        number_node_free(number_free_list, integer_pool, base);
        return out;
    }
    case '+':
        ABORT("number_exponentiate case not yet implemented.");
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
    switch (a->operation)
    {
    case '+':
        return number_add(number_stack, number_free_list, integer_pool, stack_a, stack_b, a->left,
            a->right);
    case '-':
    {
        struct Number*negative = number_slot_new(number_stack, number_free_list);
        negative->operation = 'r';
        negative->value.numerator = pool_integer_new(integer_pool, 1, -1);
        negative->value.denominator = pool_integer_new(integer_pool, 1, 1);
        struct Number*product = number_multiply(number_stack, number_free_list, integer_pool,
            stack_a, stack_b, negative, a->right);
        if (!product)
        {
            return 0;
        }
        return number_add(number_stack, number_free_list, integer_pool, stack_a, stack_b, a->left,
            product);
    }
    case '*':
        return number_multiply(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            a->left, a->right);
    case '/':
        return number_divide(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            a->left, a->right);
    case '^':
        return number_exponentiate(number_stack, number_free_list, integer_pool, stack_a, stack_b,
            a->left, a->right);
    default:
        ABORT("Number operation not recognized.");
    }
}

void print_number(struct Stack*stack_a, struct Stack*stack_b, struct Number*number)
{
    if (number->operation == 'r')
    {
        void*stack_a_savepoint = stack_a->cursor;
        char*string = stack_a->cursor;
        integer_string(stack_a, stack_b, number->value.numerator);
        if (!integer_equals(number->value.denominator, &one))
        {
            *(char*)stack_slot_new(stack_a, 1, 1) = '/';
            integer_string(stack_a, stack_b, number->value.denominator);
        }
        *(char*)stack_slot_new(stack_a, 1, 1) = 0;
        printf("%s", string);
        stack_a->cursor = stack_a_savepoint;
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

void init(struct IntegerPool*integer_pool, struct Stack*polynomial_stack, struct Stack*number_stack,
    struct Stack*stack_a, struct Stack*stack_b)
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwAllocationGranularity;
    size_t arena_size = page_size * (((size_t)system_info.lpMaximumApplicationAddress -
        (size_t)system_info.lpMinimumApplicationAddress) / (7 * page_size));
    pool_memory_start = (size_t)system_info.lpMinimumApplicationAddress;
    integer_pool_new(&permanent_integer_pool, pool_memory_start, pool_memory_start + arena_size);
    pool_memory_end = permanent_integer_pool.slot_page_stack.end + arena_size;
    integer_pool_new(integer_pool, permanent_integer_pool.slot_page_stack.end, pool_memory_end);
    stack_new(&prime_stack, integer_pool->slot_page_stack.end,
        integer_pool->slot_page_stack.end + arena_size);
    stack_new(polynomial_stack, prime_stack.end, prime_stack.end + arena_size);
    stack_new(number_stack, polynomial_stack->end, polynomial_stack->end + arena_size);
    stack_new(stack_a, number_stack->end, number_stack->end + arena_size);
    stack_new(stack_b, stack_a->end, stack_a->end + arena_size);
    primes = prime_stack.cursor;
    stack_integer_new(&prime_stack, 2, 1);
    stack_integer_new(&prime_stack, 3, 1);
    pi_estimate_min.numerator = pool_integer_new(&permanent_integer_pool, 47, 1);
    pi_estimate_min.denominator = pool_integer_new(&permanent_integer_pool, 15, 1);
    pi_interval_size.numerator = pool_integer_new(&permanent_integer_pool, 1696, 1);
    pi_interval_size.denominator = pool_integer_new(&permanent_integer_pool, 12285, 1);
    pi_sixteen_to_the_k = pool_integer_new(&permanent_integer_pool, 16, 1);
    pi_eight_k = pool_integer_new(&permanent_integer_pool, 8, 1);
}

int main()
{
    struct IntegerPool integer_pool;
    struct Stack number_stack;
    struct Number*number_free_list = 0;
    struct Stack stack_a;
    struct Stack stack_b;
    struct Stack polynomial_stack;
    init(&integer_pool, &polynomial_stack, &number_stack, &stack_a, &stack_b);
    while (true)
    {
        if (get_input(&number_stack, &number_free_list, &integer_pool, &stack_a, &stack_b))
        {
            struct Number*number = (struct Number*)number_stack.start;
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
            struct Number*input = parse_input(&number_stack, &number_free_list, &integer_pool,
                (struct Number*)number_stack.start);
            if (input)
            {
                struct Number*evaluation = number_evaluate(&number_stack, &number_free_list,
                    &integer_pool, &stack_a, &stack_b, input);
                if (evaluation)
                {
                    printf("=\n");
                    print_number(&stack_a, &stack_b, evaluation);
                }
            }
        }
        printf("\n\n");
        stack_free(&number_stack);
        number_free_list = 0;
        stack_free(&stack_a);
        stack_free(&stack_b);
        stack_free(&polynomial_stack);
        memset(integer_pool.sized_pools, 0, page_size);
        stack_free(&integer_pool.slot_page_stack);
    }
    return 0;
}