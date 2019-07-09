#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#ifdef _DEBUG
#define ASSERT(condition, message) if (!(condition)) { printf(message); abort(); }
#else
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
    uint32_t value[];
};

struct Factor
{
    struct Integer*value;
    struct Integer*multiplicity;
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
    void*coefficients[];
};

struct PolynomialDivision
{
    struct Polynomial*quotient;
    struct Polynomial*remainder;
};

struct IntegerPolynomial
{
    size_t coefficient_count;
    struct Integer*coefficients[];
};

struct RationalPolynomial
{
    size_t coefficient_count;
    struct Rational*coefficients[];
};

struct Matrix
{
    struct Rational***rows;
    size_t width;
    size_t height;
};

struct AlgebraicNumber
{
    struct AlgebraicNumber*next_term;
    struct Rational*term_coefficient;
    size_t generator_degrees[2];
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
    union
    {
        struct RationalPolynomial*minimal_polynomial;//During evaluation.
        struct Number*previous;//During parsing, when the input is stored as a doubly linked list.
    };

    //Points to the next element of the doubly linked list representing the input during parsing.
    //When on a free list during evaluation, points to the next element of the free list. While not
    //on a free list during evaluation, points to the next element in a singly linked list
    //consisting of the conjugates of the first element. If a Number is not part of a list of
    //conjugates and doesn't have operation == 'r', next == 0 means the conjugates have not yet been
    //calculated.
    struct Number*next;
    char operation;
};

struct NumberPool
{
    struct Stack stack;
    struct Number*free_list;
};

size_t page_size;

struct Integer zero = { 0, 0, 0 };
struct Integer one = { 1, 1, 1 };

struct Rational rational_zero = { &zero, &one };
struct Rational rational_one = { &one, &one };

size_t integer_pool_memory_start;
size_t integer_pool_memory_end;

struct NumberPool permanent_number_pool;
struct Stack permanent_stack;
struct IntegerPool permanent_integer_pool;
struct Stack polynomial_stack;

struct Number**roots_of_unity;
struct Integer**primes;

struct Rational pi_estimate_min;
struct Rational pi_interval_size;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

struct Polynomial polynomial_zero = { 0, 0 };
struct IntegerPolynomial integer_polynomial_one = { 1, &one };
struct RationalPolynomial rational_polynomial_one = { 1, &rational_one };

struct Float float_zero = { &zero, &zero };
struct Float float_one = { &one, &zero };

struct Number*number_divide_by_zero_error = 0;
struct Number*number_product_consolidation_failed = (struct Number*)1;

#define POINTER_SWAP(a, b) { void*temp = a; a = b; b = temp; }

__declspec(noreturn) void crash(char*message)
{
    printf(message);
    abort();
}

void array_reverse(void**a, size_t element_count)
{
    for (size_t i = 0; i < element_count / 2; ++i)
    {
        a[i] = a[element_count - 1 - i];
    }
}

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
        crash("Stack ran out of virtual address space.");
    }
    while ((size_t)output_stack->cursor > output_stack->allocation_cursor)
    {
        if (!VirtualAlloc((void*)output_stack->allocation_cursor, page_size,
            MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE))
        {
            crash("Ran out of physical memory.");
        }
        output_stack->allocation_cursor = output_stack->allocation_cursor + page_size;
    }
    return (void*)slot;
}

#define STACK_SLOT_NEW(stack, type) stack_slot_new(stack, sizeof(type), _Alignof(type))

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
        crash("Integer pool ran out of virtual address space.");
    }
    void*new_page_address = pool->slot_page_stack.cursor;
    pool->slot_page_stack.cursor = (void*)pool->slot_page_stack.allocation_cursor;
    return VirtualAlloc(new_page_address, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
}

struct SizedIntegerPool*get_sized_integer_pool(struct IntegerPool*pool, size_t value_count)
{
    return &pool->sized_pools[value_count / 2];
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
        void*m = operations->add(local_stack, output_stack, out->a_coefficient,
            operations->negative(local_stack, output_stack,
                operations->multiply(local_stack, output_stack, out->b_over_gcd, division.quotient,
                    misc), misc), misc);
        void*n = operations->add(local_stack, output_stack, out->b_coefficient,
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
    out->value[0] = value;
    return out;
}

struct Integer*pool_integer_new(struct IntegerPool*pool, uint32_t value, int8_t sign)
{
    struct Integer*out = integer_pool_slot_new(pool, 1);
    out->value_count = 1;
    out->sign = sign;
    out->value[0] = value;
    return out;
}

void pool_integer_free(struct IntegerPool*pool, struct Integer*a)
{
    ASSERT(integer_pool_memory_start <= (size_t)a && (size_t)a <= integer_pool_memory_end,
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
    if (uint_value)
    {
        struct Integer*out = integer_stack_slot_new(output_stack, 1);
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
    struct Integer*out = integer_stack_slot_new(output_stack, value_count);
    out->value_count = value_count;
    *(size_t*)&out->value = value;
    integer_trim_leading_zeroes(out);
    output_stack->cursor = &out->value[out->value_count];
    return out;
}

size_t integer_to_size_t(struct Integer*a)
{
    switch (a->value_count)
    {
    case 0:
        return 0;
    case 1:
    case 2:
        return *(size_t*)&a->value;
    default:
        crash("integer_to_size_t called on an Integer too large to fit in a size_t.");
    }
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
    ASSERT(integer_pool_memory_start <= (size_t)*a && (size_t)*a <= integer_pool_memory_end,
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
    for (int i = 0; i < a->value_count; ++i)
    {
        a->value[i] = ~a->value[i];
    }
    for (int i = 0; i < a->value_count; ++i)
    {
        uint32_t power = 1;
        for (int j = 0; j < 32; ++j)
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
    ASSERT(integer_pool_memory_start > (size_t)a || (size_t)a > integer_pool_memory_end,
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
    struct Integer*out = integer_copy_to_stack(output_stack, a);
    stack_slot_new(output_stack, (b->value_count + 1) * sizeof(uint32_t), _Alignof(struct Integer));
    integer_add_to_a_in_place(out, b);
    output_stack->cursor = &out->value[out->value_count];
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
    struct Integer*out = integer_stack_slot_new(output_stack, a->value_count + b->value_count);
    out->value_count = 0;
    out->sign = 0;
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)a->value[i] * b->value[j];
            size_t shift = i + j;
            struct Integer*integer_component = integer_stack_slot_new(local_stack, shift + 2);
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

struct Integer*integer_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b, void*unused)
{
    return integer_multiply(output_stack, local_stack, a, b);
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

void integer_halve(struct Integer*a)
{
    ASSERT(integer_pool_memory_start > (size_t)a || (size_t)a > integer_pool_memory_end,
        "integer_halve was called on a pool Integer.");
    integer_downshift(a, 1);
    integer_trim_leading_zeroes(a);
}

int leading_digit_place(struct Integer*a)
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
        out->quotient = integer_stack_slot_new(output_stack, dividend->value_count);
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
        out->remainder = integer_copy_to_stack(output_stack, out->remainder);
        if (out->remainder->sign != 0)
        {
            out->remainder->sign = dividend->sign;
        }
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
    if (a->value[a->value_count - 1] & 0x80000000)
    {
        out = integer_stack_slot_new(output_stack, a->value_count + 1);
        out->value_count = a->value_count + 1;
        out->value[a->value_count] = 0;
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

struct Integer*n_choose_k(struct Stack*output_stack, struct Stack*local_stack, struct Integer*n,
    struct Integer*k)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*numerator = &one;
    for (struct Integer*i =
        integer_add(local_stack, &one, integer_subtract(local_stack, output_stack, n, k));
        integer_compare(output_stack, local_stack, i, n) <= 0;
        i = integer_add(local_stack, i, &one))
    {
        numerator = integer_multiply(local_stack, output_stack, numerator, i);
    }
    struct Integer*denominator = &one;
    for (struct Integer*i = &INT(2, +); integer_compare(output_stack, local_stack, i, k) <= 0;
        i = integer_add(local_stack, i, &one))
    {
        denominator = integer_multiply(local_stack, output_stack, denominator, i);
    }
    struct Integer*out =
        integer_euclidean_quotient(output_stack, local_stack, numerator, denominator);
    local_stack->cursor = local_stack_savepoint;
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
        if (exponent->value[0] & 1)
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
            *next_char = division.remainder->value[0] + '0';
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
                primes[index] = integer_copy_to_stack(&permanent_stack, primes[index]);
                return primes[index];
            }
        }
        primes[index] = integer_add(stack_a, primes[index], primes[0]);
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
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
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
    out->numerator->sign *= out->denominator->sign;
    out->denominator->sign *= out->denominator->sign;
    out->numerator = integer_copy_to_stack(output_stack, out->numerator);
    out->denominator = integer_copy_to_stack(output_stack, out->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void rational_free(struct IntegerPool*pool, struct Rational*a)
{
    pool_integer_free(pool, a->numerator);
    pool_integer_free(pool, a->denominator);
}

struct Rational*rational_copy_to_stack(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
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
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
    out->numerator = integer_magnitude(output_stack, a->numerator);
    out->denominator = integer_copy_to_stack(output_stack, a->denominator);
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
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
    out->numerator = integer_add(output_stack, a->numerator,
        integer_multiply(local_stack, output_stack, b, a->denominator));
    out->denominator = integer_copy_to_stack(output_stack, a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_negative(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
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
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
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
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
    out->numerator = integer_copy_to_stack(output_stack, a->denominator);
    out->denominator = integer_copy_to_stack(output_stack, a->numerator);
    out->numerator->sign *= out->denominator->sign;
    out->denominator->sign = 1;
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

void rational_interval_expand_bounds(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalInterval*a, struct Rational*bound_candidate)
{
    if (rational_compare(stack_a, stack_b, bound_candidate, a->min) < 0)
    {
        a->min = bound_candidate;
    }
    else if (rational_compare(stack_a, stack_b, bound_candidate, a->max) > 0)
    {
        a->max = bound_candidate;
    }
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
                    &(struct Rational){ &one, integer_add(stack_a, pi_eight_k, &INT(5, +)) },
                    &(struct Rational){ &one, integer_add(stack_a, pi_eight_k, &INT(6, +)) }))),
            pi_sixteen_to_the_k));
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

struct Float*float_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*significand, struct Integer*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!(significand->value[0] & 1))
    {
        significand = integer_half(local_stack, significand);
        exponent = integer_add(local_stack, exponent, &one);
    }
    struct Float*out = STACK_SLOT_NEW(output_stack, struct Float);
    out->significand = integer_copy_to_stack(output_stack, significand);
    out->exponent = integer_copy_to_stack(output_stack, exponent);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_copy_to_stack(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = STACK_SLOT_NEW(output_stack, struct Float);
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

struct Float*float_magnitude(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = STACK_SLOT_NEW(output_stack, struct Float);
    out->significand = integer_magnitude(output_stack, a->significand);
    out->exponent = integer_copy_to_stack(output_stack, a->exponent);
    return out;
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
        out = STACK_SLOT_NEW(output_stack, struct Float);
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
    struct Float*out = STACK_SLOT_NEW(output_stack, struct Float);
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
    struct Float*out = STACK_SLOT_NEW(output_stack, struct Float);
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

int8_t float_compare(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = float_subtract(stack_a, stack_b, a, b)->significand->sign;
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Float*float_max(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b)
{
    if (float_compare(stack_a, stack_b, a, b) > 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

struct Rational*float_to_rational(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = STACK_SLOT_NEW(output_stack, struct Rational);
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
        out->denominator = &one;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
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
    struct Float**out_min, struct Float**out_max, struct Integer**estimate_denominator,
    struct Integer**estimate_remainder, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (true)
    {
        if ((*estimate_remainder)->value_count == 0)
        {
            *out_min = float_copy_to_stack(output_stack, *out_min);
            *out_max = *out_min;
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        if (integer_compare(output_stack, local_stack, interval_size->denominator,
            integer_multiply(local_stack, output_stack, *estimate_denominator,
                interval_size->numerator)) <= 0)
        {
            break;
        }
        struct IntegerDivision division;
        integer_euclidean_divide(local_stack, output_stack, &division,
            integer_doubled(local_stack, *estimate_remainder), a->denominator);
        *estimate_denominator = integer_doubled(local_stack, *estimate_denominator);
        *estimate_remainder = division.remainder;
        (*out_min)->significand = integer_add(local_stack, division.quotient,
            integer_doubled(local_stack, (*out_min)->significand));
        (*out_min)->exponent = integer_add(local_stack, (*out_min)->exponent, &INT(1, -));
    }
    *out_min =
        float_reduced(output_stack, local_stack, (*out_min)->significand, (*out_min)->exponent);
    *out_max = float_reduced(output_stack, local_stack,
        integer_add(local_stack, (*out_min)->significand, &one), (*out_min)->exponent);
    local_stack->cursor = local_stack_savepoint;
}

void rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerDivision division;
    integer_euclidean_divide(local_stack, output_stack, &division,
        integer_magnitude(local_stack, a->numerator), a->denominator);
    struct Integer*estimate_remainder = division.remainder;
    struct Integer*estimate_denominator = &one;
    if (a->numerator->sign < 0)
    {
        *out_max = STACK_SLOT_NEW(local_stack, struct Float);
        (*out_max)->significand = division.quotient;
        (*out_max)->exponent = &zero;
        rational_continue_float_estimate(output_stack, local_stack, out_max, out_min,
            &estimate_denominator, &estimate_remainder, rational_negative(local_stack, a),
            interval_size);
        (*out_min)->significand->sign = -1;
        (*out_max)->significand->sign = -1;
    }
    else
    {
        *out_min = STACK_SLOT_NEW(local_stack, struct Float);
        (*out_min)->significand = division.quotient;
        (*out_min)->exponent = &zero;
        rational_continue_float_estimate(output_stack, local_stack, out_min, out_max,
            &estimate_denominator, &estimate_remainder, a, interval_size);
    }
    local_stack->cursor = local_stack_savepoint;
}

void rational_interval_to_float_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct RationalInterval*a,
    struct Rational*bound_interval_size)
{
    struct Float*unused_bound;
    rational_float_estimate(output_stack, local_stack, out_min, &unused_bound, a->min,
        bound_interval_size);
    rational_float_estimate(output_stack, local_stack, &unused_bound, out_max, a->max,
        bound_interval_size);
}

void float_estimate_root(struct Stack*output_stack, struct Stack*local_stack, struct Float**out_min,
    struct Float**out_max, struct Float*a, struct Rational*interval_size, struct Integer*index)
{
    ASSERT(a->significand->sign >= 0, "float_estimate_root was called on a negative a value.");
    void*local_stack_savepoint = local_stack->cursor;
    if (a->exponent < 0)
    {
        (*out_max)->significand = &one;
        (*out_max)->exponent = &zero;
    }
    else
    {
        *out_max = float_copy_to_stack(local_stack, a);
    }
    struct Rational*rational_radicand = float_to_rational(local_stack, output_stack, a);
    struct Integer*index_minus_one = integer_add(local_stack, index, &INT(1, -));
    while (true)
    {
        struct Rational*delta = rational_integer_divide(local_stack, output_stack,
            rational_subtract(local_stack, output_stack,
                rational_divide(local_stack, output_stack, rational_radicand,
                    float_to_rational(local_stack, output_stack,
                        float_exponentiate(local_stack, output_stack, *out_max, index_minus_one))),
                float_to_rational(local_stack, output_stack, *out_max)), index);
        struct Float*delta_float_estimate_min;
        struct Float*delta_float_estimate_max;
        rational_float_estimate(local_stack, output_stack, &delta_float_estimate_min,
            &delta_float_estimate_max, delta,
            rational_integer_divide(local_stack, output_stack, delta, &INT(2, -)));
        *out_max = float_add(local_stack, output_stack, *out_max, delta_float_estimate_max);
        struct Rational*interval_size_comparison_value =
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, delta_float_estimate_max),
                rational_reduced(local_stack, output_stack,
                    integer_doubled(local_stack, delta->numerator), delta->denominator));
        if (rational_compare(output_stack, local_stack,
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, delta_float_estimate_max),
                rational_reduced(local_stack, output_stack,
                    integer_doubled(local_stack, delta->numerator), delta->denominator)),
            interval_size) <= 0)
        {
            *out_min = float_add(output_stack, local_stack, *out_max, delta_float_estimate_max);
            *out_max = float_copy_to_stack(output_stack, *out_max);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
    }
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

void float_interval_free(struct IntegerPool*integer_pool, struct FloatInterval*a)
{
    if (a->min.significand)
    {
        float_free(integer_pool, &a->min);
        float_free(integer_pool, &a->max);
    }
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    out->min = float_to_rational(output_stack, local_stack, &a->min);
    out->max = float_to_rational(output_stack, local_stack, &a->max);
}

//Rounded up to the nearest integer.
struct Integer*integer_square_root(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*float_square_root_min;
    struct Float*float_square_root_max;
    float_estimate_root(local_stack, output_stack, &float_square_root_min, &float_square_root_max,
        &(struct Float){ a, &zero }, &rational_one, &INT(2, +));
    struct Rational*rational_square_root =
        float_to_rational(local_stack, output_stack, float_square_root_max);
    struct Integer*out;
    if (integer_equals(rational_square_root->denominator, &one))
    {
        out = integer_copy_to_stack(output_stack, rational_square_root->numerator);
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

size_t size_t_factor(struct Stack*output_stack, struct Stack*local_stack, size_t**out, size_t a)
{
    void*local_stack_savepoint = local_stack->cursor;
    *out = STACK_SLOT_NEW(output_stack, size_t);
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
            STACK_SLOT_NEW(output_stack, size_t);
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
    *out = STACK_SLOT_NEW(output_stack, struct Factor);
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
            struct Factor*factor = *out + factor_count;
            factor->value = prime;
            factor->multiplicity = &zero;
            factor_count += 1;
            STACK_SLOT_NEW(output_stack, struct Factor);
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
        (*out)[factor_count].value = a;
        (*out)[factor_count].multiplicity = &one;
        factor_count += 1;
    }
    else
    {
        output_stack->cursor = *out + factor_count;
    }
    for (size_t i = 0; i < factor_count; ++i)
    {
        (*out)[i].value = integer_copy_to_stack(output_stack, (*out)[i].value);
        (*out)[i].multiplicity = integer_copy_to_stack(output_stack, (*out)[i].multiplicity);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
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
        a->coefficients[i] = coefficient_copy(output_stack, a->coefficients[i]);
    }
}

void*polynomial_copy(void*(coefficient_copy)(struct Stack*, void*), struct Stack*output_stack,
    struct Polynomial*a)
{
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    memcpy(out->coefficients, a->coefficients, a->coefficient_count * sizeof(void*));
    polynomial_copy_coefficients(coefficient_copy, output_stack, a);
    return out;
}

void polynomial_trim_leading_zeroes(bool(coefficient_equals_zero)(void*), struct Polynomial*a)
{
    for (size_t i = *(size_t*)a; i-- > 0;)
    {
        if (coefficient_equals_zero(a->coefficients[i]))
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
        POINTER_SWAP(a, b);
    }
    struct Polynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = coefficient_operations->add(output_stack, local_stack,
            a->coefficients[i], b->coefficients[i], 0);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        a->coefficients[i] = coefficient_operations->copy(output_stack, a->coefficients[i]);
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
        out->coefficients[i] =
            coefficient_operations->negative(output_stack, 0, a->coefficients[i], 0);
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
        out->coefficients[i] = coefficient_operations->additive_identity;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            out->coefficients[i + j] = coefficient_operations->add(output_stack, local_stack,
                out->coefficients[i + j], coefficient_operations->multiply(local_stack,
                    output_stack, a->coefficients[i], b->coefficients[j], 0), 0);
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
        out->coefficients[i] = coefficient_operations->multiply(output_stack, local_stack,
            a->coefficients[i], b, 0);
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
        divisor->coefficients[divisor->coefficient_count - 1], misc);
    for (size_t i = 1; i <= quotient_coefficient_count; ++i)
    {
        void*quotient = coefficient_operations->multiply(local_stack, output_stack,
            out->remainder->coefficients[out->remainder->coefficient_count - i],
            leading_coefficient_reciprocal, misc);
        memcpy(&out->quotient->coefficients[quotient_coefficient_count - i], quotient,
            sizeof(void*));
        for (size_t j = 1; j < divisor->coefficient_count; ++j)
        {
            out->remainder->coefficients[out->remainder->coefficient_count - i - j] =
                coefficient_operations->add(local_stack, output_stack,
                    out->remainder->coefficients[out->remainder->coefficient_count - i - j],
                    coefficient_operations->negative(local_stack, output_stack,
                        coefficient_operations->multiply(local_stack, output_stack, quotient,
                            divisor->coefficients[divisor->coefficient_count - 1 - j], misc),
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
            out->remainder->coefficients[out->remainder->coefficient_count - i],
            divisor->coefficients[divisor->coefficient_count - 1]);
        if (!division.remainder->value_count)
        {
            memcpy(&out->quotient->coefficients[quotient_coefficient_count - i], division.quotient,
                sizeof(struct Integer*));
            for (size_t j = 1; j < divisor->coefficient_count; ++j)
            {
                out->remainder->coefficients[out->remainder->coefficient_count - i - j] =
                    integer_subtract(local_stack, output_stack,
                        out->remainder->coefficients[out->remainder->coefficient_count - i - j],
                        integer_multiply(local_stack, output_stack, division.quotient,
                            divisor->coefficients[divisor->coefficient_count - 1 - j]));
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

struct IntegerPolynomial*integer_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor)
{
    struct PolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct IntegerPolynomial*)division.quotient;
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
        out->coefficients[i - 1] =
            integer_multiply(output_stack, local_stack, multiplier, a->coefficients[i]);
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
    struct Integer*out = a->coefficients[0];
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        out = integer_gcd(local_stack, output_stack, out, a->coefficients[i]);
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
        out->coefficients[i] = integer_euclidean_quotient(output_stack, local_stack,
            dividend->coefficients[i], divisor);
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
                    b_primitive_part->coefficients[b_primitive_part->coefficient_count - 1],
                    integer_add(local_stack, degree, &one))), b_primitive_part);
        if (remainder->coefficient_count <= 1)
        {
            break;
        }
        a = b;
        b = integer_polynomial_integer_divide(local_stack, output_stack, remainder,
            integer_multiply(local_stack, output_stack, g,
                integer_exponentiate(local_stack, output_stack, h, degree)));
        g = a->coefficients[a->coefficient_count - 1];
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
        out->coefficients[0] = integer_copy_to_stack(output_stack, d);
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
        out->coefficients[i] = integer_euclidean_remainder(local_stack, output_stack,
            a->coefficients[i], characteristic);
        if (out->coefficients[i]->sign < 0)
        {
            out->coefficients[i] = integer_add(local_stack, out->coefficients[i], characteristic);
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

struct IntegerPolynomial*modded_polynomial_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_multiply_by_coefficient(output_stack,
        local_stack, a, modded_integer_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
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

void modded_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct ExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    generic_extended_gcd(&modded_polynomial_operations, modded_polynomial_euclidean_divide,
        output_stack, local_stack, out, (struct Polynomial*)a, (struct Polynomial*)b,
        characteristic);
}

size_t cantor_zassenhaus_split(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, size_t degree,
    struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    if ((a->coefficient_count - 1) / degree == 1)
    {
        out[0] = modded_polynomial_multiply_by_coefficient(output_stack, local_stack, a,
            modded_integer_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
        local_stack->cursor = local_stack_savepoint;
        return 1;
    }
    struct IntegerPolynomial*b;
    if (integer_equals(characteristic, &INT(2, +)))
    {
        struct IntegerPolynomial*x_squared = polynomial_slot_new(local_stack, 3);
        x_squared->coefficients[0] = &zero;
        x_squared->coefficients[1] = &zero;
        x_squared->coefficients[2] = &one;
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
            t->coefficients[i] = p_minus_one;
        }
        t->coefficients[t->coefficient_count - 1] = &one;
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
            if (t->coefficients[i]->value_count)
            {
                t->coefficients[i] = integer_add(local_stack, t->coefficients[i], &INT(1, -));
                break;
            }
            if (i == 0)
            {
                t->coefficient_count -= 1;
                for (size_t i = 0; i < t->coefficient_count - 1; ++i)
                {
                    t->coefficients[i] = p_minus_one;
                }
                t->coefficients[t->coefficient_count - 1] = &one;
                break;
            }
            t->coefficients[i] = p_minus_one;
            --i;
        }
    }
    size_t factor_count =
        cantor_zassenhaus_split(output_stack, local_stack, b, characteristic, degree, out);
    factor_count += cantor_zassenhaus_split(output_stack, local_stack,
        modded_polynomial_euclidean_quotient(local_stack, output_stack, a, b, characteristic),
        characteristic, degree, out + factor_count);
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t squarefree_modded_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    struct IntegerPolynomial*x = polynomial_slot_new(local_stack, 2);
    x->coefficients[0] = &zero;
    x->coefficients[1] = &one;
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
            factor_count += cantor_zassenhaus_split(output_stack, local_stack,
                degree_d_factor_product, characteristic, d, out + factor_count);
            v = modded_polynomial_euclidean_quotient(local_stack, output_stack, v,
                degree_d_factor_product, characteristic);
            w = modded_polynomial_euclidean_remainder(local_stack, output_stack, w, v,
                characteristic);
        }
    }
    if (v->coefficient_count > 1)
    {
        factor_count += cantor_zassenhaus_split(output_stack, local_stack, v, characteristic,
            v->coefficient_count - 1, out + factor_count);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

struct IntegerPolynomial*bound_coefficients(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic_power)
{
    struct IntegerPolynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        struct Integer*remainder = integer_euclidean_remainder(local_stack, output_stack,
            a->coefficients[i], characteristic_power);
        if (remainder->sign < 0)
        {
            remainder = integer_add(local_stack, remainder, characteristic_power);
        }
        if (integer_compare(output_stack, local_stack, integer_doubled(local_stack, remainder),
            characteristic_power) > 0)
        {
            a->coefficients[i] =
                integer_subtract(output_stack, local_stack, remainder, characteristic_power);
        }
        else
        {
            a->coefficients[i] = integer_copy_to_stack(output_stack, remainder);
        }
    }
    return out;
}

struct IntegerPolynomial*find_valid_combination(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*combination, size_t combination_size,
    struct IntegerPolynomial**factor_candidates_start,
    struct IntegerPolynomial**factor_candidates_end, struct IntegerPolynomial**a,
    struct Integer*characteristic_power)
{
    if (combination_size == 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct IntegerPolynomial*bounded_combination =
            bound_coefficients(local_stack, output_stack, combination, characteristic_power);
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                (*a)->coefficients[(*a)->coefficient_count - 1]), bounded_combination);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            do
            {
                *a = (struct IntegerPolynomial*)division.quotient;
                integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
                    integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                        (*a)->coefficients[(*a)->coefficient_count - 1]),
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
        for (struct IntegerPolynomial**candidate = factor_candidates_start;
            *candidate <= *factor_candidates_end - combination_size; ++*candidate)
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
                integer_polynomial_multiply(local_stack, output_stack, combination, *candidate),
                combination_size - 1, factor_candidates_start + 1, factor_candidates_end, a,
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

size_t squarefree_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*modded_a;
    struct IntegerPolynomial*gcd;
    size_t prime_index = 0;
    struct Integer*prime = primes[0];
    while (true)
    {
        modded_a = modded_polynomial_reduced(local_stack, output_stack, a, prime);
        if (modded_a->coefficient_count == a->coefficient_count)
        {
            gcd = modded_polynomial_gcd(local_stack, output_stack, modded_a,
                modded_polynomial_reduced(local_stack, output_stack,
                    integer_polynomial_derivative(local_stack, output_stack, modded_a), prime),
                prime);
            if (gcd->coefficient_count == 1)
            {
                break;
            }
        }
        ++prime_index;
        prime = get_prime(output_stack, local_stack, prime_index);
    }
    modded_a = modded_polynomial_monic(local_stack, output_stack, modded_a, prime);
    struct IntegerPolynomial**modded_a_factors = stack_slot_new(local_stack,
        (modded_a->coefficient_count - 1) * sizeof(struct IntegerPolynomial*),
        _Alignof(struct IntegerPolynomial*));
    size_t modded_a_factor_count = squarefree_modded_polynomial_factor(local_stack, output_stack,
        modded_a, prime, modded_a_factors);
    struct Integer*coefficient_bound = &zero;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficient_bound = integer_add(local_stack, coefficient_bound,
            integer_multiply(local_stack, output_stack, a->coefficients[i], a->coefficients[i]));
    }
    struct Integer*a_degree_minus_one = integer_from_size_t(local_stack, a->coefficient_count - 2);
    struct Integer*k =
        integer_euclidean_quotient(local_stack, output_stack, a_degree_minus_one, &INT(2, +));
    struct Integer*leading_coefficient_magnitude =
        integer_magnitude(local_stack, a->coefficients[a->coefficient_count - 1]);
    coefficient_bound = integer_doubled(local_stack, integer_multiply(local_stack, output_stack,
        leading_coefficient_magnitude, integer_add(local_stack,
            integer_multiply(local_stack, output_stack,
                integer_square_root(local_stack, output_stack, coefficient_bound),
                n_choose_k(local_stack, output_stack, a_degree_minus_one, k)),
            integer_multiply(local_stack, output_stack, leading_coefficient_magnitude,
                n_choose_k(local_stack, output_stack, a_degree_minus_one,
                    integer_add(local_stack, k, &INT(1, -)))))));
    struct Integer*e = &one;
    struct Integer*prime_power = prime;
    while (integer_compare(output_stack, local_stack, prime_power, coefficient_bound) < 0)
    {
        e = integer_add(local_stack, e, &one);
        prime_power = integer_multiply(local_stack, output_stack, prime_power, prime);
    }
    struct IntegerPolynomial*product_of_unlifted_factors = modded_a;
    struct IntegerPolynomial*unmodded_b_times_c = modded_polynomial_monic(local_stack, output_stack,
        modded_polynomial_reduced(local_stack, output_stack, a, prime_power), prime_power);
    struct IntegerPolynomial**lifted_factors = stack_slot_new(local_stack,
        modded_a_factor_count * sizeof(struct IntegerPolynomial*),
        _Alignof(struct IntegerPolynomial*));
    size_t lifted_factor_count = 0;
    modded_a_factor_count -= 1;
    for (size_t i = 0; i < modded_a_factor_count; ++i)
    {
        product_of_unlifted_factors = modded_polynomial_euclidean_quotient(local_stack,
            output_stack, product_of_unlifted_factors, modded_a_factors[i], prime);
        struct Integer*lift_characteristic = prime;
        struct IntegerPolynomial*b = modded_a_factors[i];
        struct IntegerPolynomial*c = product_of_unlifted_factors;
        while (integer_compare(output_stack, local_stack, lift_characteristic, prime_power) < 0)
        {
            struct IntegerPolynomial*b_mod_prime =
                modded_polynomial_reduced(local_stack, output_stack, b, prime);
            struct IntegerPolynomial*c_mod_prime =
                modded_polynomial_reduced(local_stack, output_stack, c, prime);
            struct ExtendedGCDInfo gcd_info;
            modded_polynomial_extended_gcd(local_stack, output_stack, &gcd_info, b_mod_prime,
                c_mod_prime, prime);
            struct Integer*reciprocal = modded_integer_reciprocal(local_stack, output_stack,
                ((struct IntegerPolynomial*)gcd_info.gcd)->coefficients[0], prime);
            gcd_info.a_coefficient = modded_polynomial_multiply_by_coefficient(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.a_coefficient, reciprocal, prime);
            gcd_info.b_coefficient = modded_polynomial_multiply_by_coefficient(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.b_coefficient, reciprocal, prime);
            struct IntegerPolynomial*f = modded_polynomial_reduced(local_stack, output_stack,
                integer_polynomial_integer_divide(local_stack, output_stack,
                    integer_polynomial_subtract(local_stack, output_stack, unmodded_b_times_c,
                        integer_polynomial_multiply(local_stack, output_stack, b, c)),
                    lift_characteristic), prime);
            struct PolynomialDivision division;
            modded_polynomial_euclidean_divide(local_stack, output_stack, &division,
                modded_polynomial_multiply(local_stack, output_stack, f,
                    (struct IntegerPolynomial*)gcd_info.b_coefficient, prime), b_mod_prime, prime);
            b = integer_polynomial_add(local_stack, output_stack, b,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    (struct IntegerPolynomial*)division.remainder, lift_characteristic));
            c = integer_polynomial_add(local_stack, output_stack, c,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    modded_polynomial_add(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack,
                            (struct IntegerPolynomial*)gcd_info.a_coefficient, f, prime),
                        modded_polynomial_multiply(local_stack, output_stack, c_mod_prime,
                            (struct IntegerPolynomial*)division.quotient, prime), prime),
                    lift_characteristic));
            lift_characteristic =
                integer_multiply(local_stack, output_stack, lift_characteristic, prime);
        }
        lifted_factors[lifted_factor_count] = b;
        ++lifted_factor_count;
        unmodded_b_times_c = c;
    }
    lifted_factors[lifted_factor_count] = unmodded_b_times_c;
    ++lifted_factor_count;
    struct IntegerPolynomial*a_times_leading_coefficient =
        integer_polynomial_integer_multiply(local_stack, output_stack, a,
            a->coefficients[a->coefficient_count - 1]);
    size_t factor_count = 0;
    size_t combination_size = 1;
    struct IntegerPolynomial*unfactored_component_of_a = a;
    struct IntegerPolynomial**lifted_factors_end = lifted_factors + lifted_factor_count;
    while (2 * combination_size < lifted_factor_count)
    {
        while (true)
        {
            struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
                &integer_polynomial_one, combination_size, lifted_factors, lifted_factors_end,
                &unfactored_component_of_a, prime_power);
            if (!factor)
            {
                break;
            }
            out[factor_count] = factor;
            ++factor_count;
        }
        ++combination_size;
    }
    if (2 * combination_size == lifted_factor_count)
    {
        lifted_factors_end -= 1;
        struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
            *lifted_factors_end, combination_size - 1, lifted_factors, lifted_factors_end,
            &unfactored_component_of_a, prime_power);
        if (factor)
        {
            out[factor_count] = factor;
            ++factor_count;
        }
    }
    struct IntegerPolynomial*final_factor = polynomial_slot_new(local_stack, 1);
    final_factor->coefficients[0] =
        unfactored_component_of_a->coefficients[unfactored_component_of_a->coefficient_count - 1];
    for (struct IntegerPolynomial**factor = lifted_factors; factor < lifted_factors_end;
        factor += 1)
    {
        final_factor =
            integer_polynomial_multiply(local_stack, output_stack, final_factor, *factor);
    }
    if (final_factor->coefficient_count > 1)
    {
        out[factor_count] = integer_polynomial_primitive_part(output_stack, local_stack,
            bound_coefficients(local_stack, output_stack, final_factor, prime_power));
        ++factor_count;
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

bool list_contains_a(struct IntegerPolynomial**list, size_t list_length,
    struct IntegerPolynomial*a)
{
    for (size_t i = 0; i < list_length; ++i)
    {
        if (a->coefficient_count == list[i]->coefficient_count)
        {
            for (size_t j = 0; j < a->coefficient_count; ++j)
            {
                if (!integer_equals(a->coefficients[j], list[i]->coefficients[j]))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

size_t primitive_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    if (!a->coefficients[0]->value_count)
    {
        out[0] = polynomial_slot_new(output_stack, 2);
        out[0]->coefficients[0] = &zero;
        out[0]->coefficients[1] = &one;
        factor_count = 1;
        do
        {
            a = integer_polynomial_euclidean_quotient(local_stack, output_stack, a, *out);
        } while (!a->coefficients[0]->value_count);
    }
    if (a->coefficients[0]->value_count < 2)
    {
        local_stack->cursor = local_stack_savepoint;
        return factor_count;
    }
    struct IntegerPolynomial**squarefree_factors =
        stack_slot_new(local_stack, (a->coefficient_count - 1) * sizeof(struct IntegerPolynomial*),
            _Alignof(struct IntegerPolynomial*));
    size_t squarefree_factor_count = 0;
    struct IntegerPolynomial*d = integer_polynomial_derivative(local_stack, output_stack, a);
    struct IntegerPolynomial*b =
        integer_polynomial_gcd(local_stack, output_stack, a, d);
    struct IntegerPolynomial*c =
        integer_polynomial_euclidean_quotient(local_stack, output_stack, a, b);
    while (c->coefficient_count != 1 ||
        !integer_equals(integer_magnitude(local_stack, c->coefficients[0]), &one))
    {
        d = integer_polynomial_subtract(local_stack, output_stack,
            integer_polynomial_euclidean_quotient(local_stack, output_stack, d, b),
            integer_polynomial_derivative(local_stack, output_stack, c));
        b = integer_polynomial_gcd(local_stack, output_stack, c, d);
        if (!list_contains_a(squarefree_factors, squarefree_factor_count, b))
        {
            squarefree_factors[squarefree_factor_count] = b;
            ++squarefree_factor_count;
        }
        c = integer_polynomial_euclidean_quotient(local_stack, output_stack, c, b);
    }
    for (size_t i = 0; i < squarefree_factor_count; ++i)
    {
        if (integer_compare(output_stack, local_stack,
            integer_magnitude(local_stack, squarefree_factors[i]->coefficients[0]),
            integer_magnitude(local_stack, squarefree_factors[i]->coefficients
                [squarefree_factors[i]->coefficient_count - 1])) < 0)
        {
            array_reverse(squarefree_factors[i]->coefficients,
                squarefree_factors[i]->coefficient_count);
            struct IntegerPolynomial**reversed_irreducible_factors = stack_slot_new(local_stack,
                (squarefree_factor_count - 1) * sizeof(struct IntegerPolynomial*),
                _Alignof(struct IntegerPolynomial*));
            size_t reversed_factor_count = squarefree_integer_polynomial_factor(local_stack,
                output_stack, squarefree_factors[i], reversed_irreducible_factors);
            for (size_t j = 0; j < reversed_factor_count; ++j)
            {
                array_reverse(reversed_irreducible_factors[j]->coefficients,
                    reversed_irreducible_factors[j]->coefficient_count);
                out[factor_count] =
                    integer_polynomial_copy(output_stack, reversed_irreducible_factors[j]);
                ++factor_count;
            }
        }
        else
        {
            factor_count += squarefree_integer_polynomial_factor(output_stack, local_stack,
                squarefree_factors[i], out);
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

bool rational_polynomial_equals(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (rational_compare(stack_a, stack_b, a->coefficients[i], b->coefficients[i]) != 0)
        {
            return false;
        }
    }
    return true;
}

struct RationalPolynomial*rational_polynomial_copy_to_number_memory(
    struct IntegerPool*integer_pool, struct RationalPolynomial*a)
{
    struct RationalPolynomial*out = polynomial_slot_new(&polynomial_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = a->coefficients[i];
        rational_move_to_pool(integer_pool, out->coefficients[i]);
    }
    return out;
}

struct RationalPolynomial*rational_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_add(&rational_ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_subtract(&rational_ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_multiply(&rational_ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Rational*b)
{
    return polynomial_multiply_by_coefficient(&integer_operations, local_stack, output_stack,
        (struct Polynomial*)a, b);
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
            a->coefficients[i]->denominator);
    }
    struct IntegerPolynomial*integer_polynomial =
        polynomial_slot_new(local_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        integer_polynomial->coefficients[i] = integer_multiply(local_stack, output_stack,
            a->coefficients[i]->numerator, integer_euclidean_quotient(local_stack, output_stack,
                denominator_lcm, a->coefficients[i]->denominator));
    }
    struct IntegerPolynomial*out =
        integer_polynomial_primitive_part(output_stack, local_stack, integer_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*integer_polynomial_to_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    struct RationalPolynomial*out = polynomial_slot_new(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_reduced(output_stack, local_stack, a->coefficients[i],
            a->coefficients[a->coefficient_count - 1]);
    }
    return out;
}

size_t rational_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = primitive_integer_polynomial_factor(local_stack, output_stack,
        rational_polynomial_primitive_part(local_stack, output_stack, a),
        (struct IntegerPolynomial**)out);
    for (size_t i = 0; i < factor_count; ++i)
    {
        out[i] = integer_polynomial_to_monic(local_stack, output_stack,
            (struct IntegerPolynomial*)out[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

void matrix_row_echelon_form(struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a,
    struct RationalPolynomial**augmentation)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t small_dimension = min(a->width, a->height);
    for (size_t i = 0; i < small_dimension; ++i)
    {
        for (size_t j = i + 1; j < a->height; ++j)
        {
            if (a->rows[i][i]->numerator->value_count == 0)
            {
                POINTER_SWAP(a->rows[i], a->rows[j]);
                POINTER_SWAP(augmentation[i], augmentation[j]);
            }
            else
            {
                for (size_t k = i + 1; k < a->height; ++k)
                {
                    struct Rational*scalar =
                        rational_divide(local_stack, output_stack, a->rows[k][i], a->rows[i][i]);
                    for (size_t l = i; l < a->width; ++l)
                    {
                        a->rows[k][l] = rational_subtract(local_stack, output_stack, a->rows[k][l],
                            rational_multiply(local_stack, output_stack, a->rows[i][l], scalar));
                    }
                    augmentation[k] =
                        rational_polynomial_subtract(local_stack, output_stack, augmentation[k],
                            rational_polynomial_rational_multiply(local_stack, output_stack,
                                augmentation[i], scalar));
                }
                break;
            }
        }
    }
    for (size_t i = 0; i < a->height; ++i)
    {
        for (size_t j = 0; j < a->width; ++j)
        {
            a->rows[i][j] = rational_copy_to_stack(output_stack, a->rows[i][j]);
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

int8_t degree_compare(size_t a[2], size_t b[2])
{
    for (size_t i = 0; i < 2; ++i)
    {
        if (a[i] < b[i])
        {
            return -1;
        }
        else if (a[i] > b[i])
        {
            return 1;
        }
    }
    return 0;
}

void algebraic_number_move_coefficients_to_stack(struct Stack*output_stack,
    struct AlgebraicNumber*a)
{
    while (a)
    {
        rational_copy_to_stack(output_stack, a->term_coefficient);
        a = a->next_term;
    }
}

struct AlgebraicNumber*pool_algebraic_number_slot_new(struct Stack*output_stack,
    struct AlgebraicNumber*free_list)
{
    if (free_list)
    {
        struct AlgebraicNumber*out = free_list;
        free_list = free_list->next_term;
        return out;
    }
    return STACK_SLOT_NEW(output_stack, struct AlgebraicNumber);
}

//There must not be anything on output_stack after a. Leaves the coefficients on local_stack.
void algebraic_number_add_to_a_in_place(struct Stack*output_stack, struct Stack*local_stack,
    struct AlgebraicNumber**a, struct AlgebraicNumber*b)
{
    struct AlgebraicNumber*free_list = 0;
    struct AlgebraicNumber*a_term = *a;
    struct AlgebraicNumber*previous_term = 0;
    while (b)
    {
        int8_t degree_comparison = degree_compare(a_term->generator_degrees, b->generator_degrees);
        while (degree_comparison < 0)
        {
            if (a_term->next_term)
            {
                previous_term = a_term;
                a_term = a_term->next_term;
                degree_comparison = degree_compare(a_term->generator_degrees, b->generator_degrees);
            }
            else
            {
                a_term->next_term = pool_algebraic_number_slot_new(output_stack, free_list);
                while (b)
                {
                    a_term = a_term->next_term;
                    memcpy(a_term, b, sizeof(struct AlgebraicNumber));
                    b = b->next_term;
                    a_term->next_term = pool_algebraic_number_slot_new(output_stack, free_list);
                }
                output_stack->cursor = a_term->next_term;
                a_term->next_term = 0;
                return;
            }
        }
        if (degree_comparison == 0)
        {
            a_term->term_coefficient = rational_add(local_stack, output_stack,
                a_term->term_coefficient, b->term_coefficient);
            if (!a_term->term_coefficient->numerator->value_count)
            {
                if (previous_term)
                {
                    previous_term->next_term = a_term->next_term;
                }
                else
                {
                    *a = a_term->next_term;
                }
                a_term->next_term = free_list;
                free_list = a_term;
            }
        }
        else
        {
            struct AlgebraicNumber*term_to_insert =
                pool_algebraic_number_slot_new(output_stack, free_list);
            memcpy(term_to_insert, b, sizeof(struct AlgebraicNumber));
            term_to_insert->next_term = a_term;
            if (previous_term)
            {
                previous_term->next_term = term_to_insert;
            }
            else
            {
                *a = term_to_insert;
            }
        }
        b = b->next_term;
    }
}

struct AlgebraicNumber*algebraic_number_add(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct AlgebraicNumber*out = STACK_SLOT_NEW(output_stack, struct AlgebraicNumber);
    struct AlgebraicNumber*out_term = out;
    out_term->next_term = out_term;
    while (a)
    {
        out_term = out_term->next_term;
        memcpy(out_term, a, sizeof(struct AlgebraicNumber));
        a = a->next_term;
        out_term->next_term = STACK_SLOT_NEW(output_stack, struct AlgebraicNumber);
    }
    output_stack->cursor = out_term->next_term;
    out_term->next_term = 0;    
    algebraic_number_add_to_a_in_place(output_stack, local_stack, &out, b);
    algebraic_number_move_coefficients_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct AlgebraicNumber*algebraic_number_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b,
    struct RationalPolynomial*generator_annulling_polynomials[2])
{
    void*local_stack_savepoint = local_stack->cursor;
    struct AlgebraicNumber*out = 0;
    while (a)
    {
        while (b)
        {
            struct Rational*coefficient = rational_multiply(local_stack, output_stack,
                a->term_coefficient, b->term_coefficient);
            if (coefficient->numerator->value_count)
            {
                struct AlgebraicNumber*product_component =
                    STACK_SLOT_NEW(local_stack, struct AlgebraicNumber);
                memset(product_component, 0, sizeof(struct AlgebraicNumber));
                product_component->term_coefficient = coefficient;
                for (size_t i = 0; i < 2; ++i)
                {
                    size_t component_generator_degree =
                        a->generator_degrees[i] + b->generator_degrees[i];
                    if (component_generator_degree <
                        generator_annulling_polynomials[i]->coefficient_count - 1)
                    {
                        struct AlgebraicNumber*product_component_term = product_component;
                        while (product_component_term)
                        {
                            product_component_term->generator_degrees[i] =
                                component_generator_degree;
                        }
                    }
                    else
                    {
                        size_t reduced_degree = component_generator_degree -
                            generator_annulling_polynomials[i]->coefficient_count + 1;
                        struct AlgebraicNumber*component_factor =
                            STACK_SLOT_NEW(local_stack, struct AlgebraicNumber);
                        struct AlgebraicNumber*component_factor_term = component_factor;
                        component_factor_term->next_term = component_factor_term;
                        for (size_t j = 0;
                            j < generator_annulling_polynomials[i]->coefficient_count - 1; ++j)
                        {
                            coefficient = generator_annulling_polynomials[i]->coefficients[j];
                            if (coefficient->numerator->value_count)
                            {
                                component_factor_term = component_factor_term->next_term;
                                memset(component_factor_term, 0, sizeof(struct AlgebraicNumber));
                                component_factor_term->generator_degrees[i] = j + reduced_degree;
                                component_factor_term->term_coefficient =
                                    rational_negative(local_stack, coefficient);
                                component_factor_term->next_term =
                                    STACK_SLOT_NEW(local_stack, struct AlgebraicNumber);
                            }
                        }
                        component_factor_term->next_term = 0;
                        product_component = algebraic_number_multiply(local_stack, output_stack,
                            product_component, component_factor, generator_annulling_polynomials);
                    }
                }
                algebraic_number_add_to_a_in_place(local_stack, output_stack, &out,
                    product_component);
            }
            b = b->next_term;
        }
        a = a->next_term;
    }
    algebraic_number_move_coefficients_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_pool_new(struct NumberPool*out, size_t start, size_t end)
{
    stack_new(&out->stack, start, end);
    out->free_list = 0;
}

struct Number*number_slot_new(struct NumberPool*pool)
{
    struct Number*out = pool->free_list;
    if (out)
    {
        pool->free_list = out->next;
    }
    else
    {
        out = STACK_SLOT_NEW(&pool->stack, struct Number);
    }
    memset(out, 0, sizeof(struct Number));
    return out;
}

void number_free(struct NumberPool*number_pool, struct IntegerPool*integer_pool, struct Number*a);

void number_node_free_during_parse(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Number*a)
{
    if (a->operation == 'r')
    {
        rational_free(integer_pool, &a->value);
    }
    a->next = number_pool->free_list;
    number_pool->free_list = a;
}

void number_node_free(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Number*a)
{
    if (a->operation != 'r')
    {
        float_interval_free(integer_pool, &a->real_part_estimate);
        float_interval_free(integer_pool, &a->imaginary_part_estimate);
    }
    float_interval_free(integer_pool, &a->argument_estimate);
    float_interval_free(integer_pool, &a->magnitude_estimate);
    if (a->next)
    {
        number_free(number_pool, integer_pool, a->next);
    }
    number_node_free_during_parse(number_pool, integer_pool, a);
}

void number_free(struct NumberPool*number_pool, struct IntegerPool*integer_pool, struct Number*a)
{
    if (a->operation != 'r')
    {
        number_free(number_pool, integer_pool, a->left);
        number_free(number_pool, integer_pool, a->right);
    }
    number_node_free(number_pool, integer_pool, a);
}

struct Number*number_copy(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Number*a)
{
    struct Number*out = number_slot_new(number_pool);
    memcpy(out, a, sizeof(struct Number));
    if (a->operation == 'r')
    {
        rational_move_to_pool(integer_pool, &out->value);
    }
    else
    {
        out->left = number_copy(number_pool, integer_pool, a->left);
        out->right = number_copy(number_pool, integer_pool, a->right);
        if (a->real_part_estimate.min.significand)
        {
            float_interval_move_to_pool(integer_pool, &out->real_part_estimate);
        }
        if (a->imaginary_part_estimate.min.significand)
        {
            float_interval_move_to_pool(integer_pool, &out->imaginary_part_estimate);
        }
    }
    if (a->argument_estimate.min.significand)
    {
        float_interval_move_to_pool(integer_pool, &out->argument_estimate);
    }
    if (a->magnitude_estimate.min.significand)
    {
        float_interval_move_to_pool(integer_pool, &out->magnitude_estimate);
    }
    return out;
}

bool get_input(struct NumberPool*transient_number_pool, struct IntegerPool*transient_integer_pool,
    struct Stack*stack_a, struct Stack*stack_b)
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
        struct Number*number = number_slot_new(transient_number_pool);
        number->previous = previous;
        number->next = transient_number_pool->stack.cursor;
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.denominator = pool_integer_new(transient_integer_pool, 1, 1);
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
            integer_move_to_pool(transient_integer_pool, &number->value.numerator);
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

struct Number*parse_input(struct NumberPool*transient_number_pool,
    struct IntegerPool*transient_integer_pool, struct Number*input)
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
                parse_input(transient_number_pool, transient_integer_pool, input->next);
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
                number_node_free_during_parse(transient_number_pool, transient_integer_pool, input);
                number_node_free_during_parse(transient_number_pool, transient_integer_pool,
                    nested_number);
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
            number_node_free_during_parse(transient_number_pool, transient_integer_pool, input);
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
                input->value.numerator = pool_integer_new(transient_integer_pool, 1, -1);
                input->value.denominator = pool_integer_new(transient_integer_pool, 1, 1);
                struct Number*times = number_slot_new(transient_number_pool);
                times->operation = '*';
                times->left = input;
                times->right = input->next;
                times->next = input->next->next;
                if (input->previous)
                {
                    input->previous->next = times;
                    times->previous = input;
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
                number_node_free_during_parse(transient_number_pool, transient_integer_pool,
                    input->next);
                number_node_free_during_parse(transient_number_pool, transient_integer_pool, input);
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

void number_rational_magnitude_estimate(struct IntegerPool*integer_pool, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size);

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
        struct Float*magnitude_estimate_min;
        struct Float*magnitude_estimate_max;
        rational_continue_float_estimate(stack_a, stack_b, &magnitude_estimate_min,
            &magnitude_estimate_max, &a->magnitude_estimate_denominator,
            &a->magnitude_estimate_remainder, &a->value, interval_size);
        a->argument_estimate.min = *magnitude_estimate_min;
        a->argument_estimate.max = *magnitude_estimate_max;
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
    switch (a->operation)
    {
    case '^':
    {
        struct FloatInterval*radicand_magnitude_estimate =
            number_float_magnitude_estimate(integer_pool, stack_a, stack_b, a->left,
                rational_integer_divide(stack_a, stack_b,
                    rational_exponentiate(stack_a, stack_b, interval_size,
                        a->right->value.denominator),
                    &INT(3, +)));
        struct Rational*bound_interval_size =
            rational_integer_divide(stack_a, stack_b, interval_size, &INT(3, +));
        struct Float*magnitude_estimate_min;
        struct Float*magnitude_estimate_max;
        struct Float*unused_estimate_bound;
        float_estimate_root(stack_a, stack_b, &magnitude_estimate_min, &unused_estimate_bound,
            &radicand_magnitude_estimate->min, bound_interval_size, a->right->value.denominator);
        float_estimate_root(stack_a, stack_b, &unused_estimate_bound, &magnitude_estimate_max,
            &radicand_magnitude_estimate->max, bound_interval_size, a->right->value.denominator);
        a->argument_estimate.min = *magnitude_estimate_min;
        a->argument_estimate.max = *magnitude_estimate_max;
        float_interval_move_to_pool(integer_pool, &a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return &a->magnitude_estimate;
    }
    case '*':
    {
        struct RationalInterval left_factor_magnitude_estimate;
        number_rational_magnitude_estimate(integer_pool, stack_a, stack_b,
            &left_factor_magnitude_estimate, a->left, &rational_one);
        struct RationalInterval right_factor_magnitude_estimate;
        number_rational_magnitude_estimate(integer_pool, stack_a, stack_b,
            &right_factor_magnitude_estimate, a->right, &rational_one);
        interval_size = rational_divide(stack_a, stack_b, interval_size,
            rational_multiply(stack_a, stack_b,
                rational_add(stack_a, stack_b, left_factor_magnitude_estimate.max, &rational_one),
                rational_add(stack_a, stack_b, right_factor_magnitude_estimate.max, 
                    &rational_one)));
        number_float_magnitude_estimate(integer_pool, stack_a, stack_b, a->left, interval_size);
        number_float_magnitude_estimate(integer_pool, stack_a, stack_b, a->right, interval_size);
        a->argument_estimate.min = *float_multiply(stack_a, stack_b,
            &a->left->magnitude_estimate.min, &a->right->magnitude_estimate.min);
        a->argument_estimate.max = *float_multiply(stack_a, stack_b,
            &a->left->magnitude_estimate.max, &a->right->magnitude_estimate.max);
        float_interval_move_to_pool(integer_pool, &a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return &a->magnitude_estimate;
    }
    }
    crash("number_float_magnitude_estimate case not yet implemented.");
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
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*rational_estimate_interval_size =
        rational_place_value(stack_a, stack_b, interval_size);
    struct RationalInterval rational_estimate;
    get_rational_estimate(integer_pool, stack_a, stack_b, &rational_estimate, a,
        rational_estimate_interval_size);
    struct Float*out_min;
    struct Float*out_max;
    rational_interval_to_float_interval(stack_a, stack_b, &out_min, &out_max, &rational_estimate,
        rational_estimate_interval_size);
    float_move_to_pool(integer_pool, out_min);
    float_move_to_pool(integer_pool, out_max);
    out->min = *out_min;
    out->max = *out_max;
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
            out->min = STACK_SLOT_NEW(output_stack, struct Rational);
            out->min->numerator = &zero;
            out->min->denominator = &one;
            out->max = STACK_SLOT_NEW(output_stack, struct Rational);
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
        crash("number_rational_argument_estimate case not yet implemented.");
    }
}

void number_estimate_part_sum(struct FloatInterval*(estimate_term)(struct IntegerPool*,
    struct Stack*, struct Stack*, struct Number*, struct Rational*),
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*term_estimate_interval_size =
        rational_integer_divide(stack_a, stack_b, interval_size, &INT(2, +));
    struct FloatInterval*left_term_estimate =
        estimate_term(integer_pool, stack_a, stack_b, a->left, term_estimate_interval_size);
    struct FloatInterval*right_term_estimate =
        estimate_term(integer_pool, stack_a, stack_b, a->right, term_estimate_interval_size);
    out->min = *float_add(stack_a, stack_b, &left_term_estimate->min, &right_term_estimate->min);
    out->max = *float_add(stack_a, stack_b, &left_term_estimate->max, &right_term_estimate->max);
    float_move_to_pool(integer_pool, &out->min);
    float_move_to_pool(integer_pool, &out->max);
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
            stack_b, &a->argument_estimate, a, interval_size);
        return &a->argument_estimate;
    }
    else
    {
        number_float_estimate_from_rational(number_rational_argument_estimate, integer_pool,
            stack_a, stack_b, &a->argument_estimate, a, interval_size);
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
            struct FloatInterval*real_part_estimate =
                STACK_SLOT_NEW(output_stack, struct FloatInterval);
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
            output_stack, local_stack, &a->real_part_estimate, a, interval_size);
        return &a->real_part_estimate;
    case '*':
        crash("number_float_real_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_real_part_estimate, integer_pool, output_stack,
            local_stack, &a->real_part_estimate, a, interval_size);
    }
    return &a->real_part_estimate;
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
        struct FloatInterval*imaginary_part_estimate =
            STACK_SLOT_NEW(output_stack, struct FloatInterval);
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
    {
        number_float_estimate_from_rational(number_rational_imaginary_part_estimate, integer_pool,
            output_stack, local_stack, &a->imaginary_part_estimate, a, interval_size);
    }
    case '*':
        crash("number_float_imaginary_part_estimate case not yet implemented.");
    case '+':
        number_estimate_part_sum(number_float_imaginary_part_estimate, integer_pool, output_stack,
            local_stack, &a->imaginary_part_estimate, a, interval_size);
    }
    return &a->imaginary_part_estimate;
}

void rational_polynomial_estimate_evaluation(struct IntegerPool*integer_pool,
    struct Stack*output_stack, struct Stack*local_stack, struct Float**real_estimate_min_out,
    struct Float**real_estimate_max_out, struct Float**imaginary_estimate_min_out,
    struct Float**imaginary_estimate_max_out, struct RationalPolynomial*a, struct Number*argument,
    struct Rational*interval_size)
{
    if (a->coefficient_count == 0)
    {
        *real_estimate_min_out = &float_zero;
        *real_estimate_max_out = &float_zero;
        *imaginary_estimate_min_out = &float_zero;
        *imaginary_estimate_max_out = &float_zero;
        return;
    }
    if (a->coefficient_count == 1)
    {
        rational_float_estimate(output_stack, local_stack, real_estimate_min_out,
            real_estimate_max_out, a->coefficients[0], interval_size);
        *imaginary_estimate_min_out = &float_zero;
        *imaginary_estimate_max_out = &float_zero;
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*evaluation_real_term_count;
    struct Integer*evaluation_imaginary_term_count;
    struct IntegerDivision division;
    integer_euclidean_divide(local_stack, output_stack, &division,
        integer_from_size_t(local_stack, a->coefficient_count), &INT(2, +));
    if (!division.remainder->value_count)
    {
        evaluation_imaginary_term_count =
            integer_multiply(local_stack, output_stack, division.quotient, division.quotient);
        evaluation_real_term_count = integer_add(local_stack, evaluation_imaginary_term_count,
            integer_add(local_stack, division.quotient, &INT(1, -)));
    }
    else
    {
        evaluation_imaginary_term_count = integer_add(local_stack, division.quotient,
            integer_multiply(local_stack, output_stack, division.quotient, division.quotient));
        evaluation_real_term_count =
            integer_add(local_stack, evaluation_imaginary_term_count, division.quotient);
    }
    struct Rational*real_interval_scale = rational_magnitude(local_stack, a->coefficients[0]);
    struct Rational*imaginary_interval_scale = &rational_one;
    struct FloatInterval*argument_real_part_estimate = number_float_real_part_estimate(integer_pool,
        local_stack, output_stack, argument, &rational_one);
    struct FloatInterval*argument_imaginary_part_estimate =
        number_float_imaginary_part_estimate(integer_pool, local_stack, output_stack, argument,
            &rational_one);
    struct Rational*real_magnitude_bound = float_to_rational(local_stack, output_stack,
        float_max(output_stack, local_stack,
            float_magnitude(local_stack, &argument_real_part_estimate->min),
            float_magnitude(local_stack, &argument_real_part_estimate->max)));
    struct Rational*imaginary_magnitude_bound = float_to_rational(local_stack, output_stack,
        float_max(output_stack, local_stack,
            float_magnitude(local_stack, &argument_imaginary_part_estimate->min),
            float_magnitude(local_stack, &argument_imaginary_part_estimate->max)));
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct Rational*coefficient_magnitude = rational_magnitude(local_stack, a->coefficients[i]);
        struct Rational*real_magnitude_power = rational_exponentiate(local_stack, output_stack,
            real_magnitude_bound, integer_from_size_t(local_stack, i - 1));
        real_interval_scale = rational_max(output_stack, local_stack, real_interval_scale,
            rational_doubled(local_stack, output_stack,
                rational_multiply(local_stack, output_stack, coefficient_magnitude,
                    real_magnitude_power)));
        struct Rational*imaginary_magnitude_power = &rational_one;
        struct Integer*active_term_count = evaluation_imaginary_term_count;
        struct Integer*inactive_term_count = evaluation_real_term_count;
        for (size_t j = 1; j < i; ++j)
        {
            struct Rational*next_imaginary_magnitude_power = rational_multiply(local_stack,
                output_stack, imaginary_magnitude_power, imaginary_magnitude_bound);
            struct Rational*shared_scale_component = rational_doubled(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    rational_integer_multiply(local_stack, output_stack,
                        rational_multiply(local_stack, output_stack, coefficient_magnitude,
                            rational_add(local_stack, output_stack, &rational_one,
                                rational_add(local_stack, output_stack, real_magnitude_power,
                                    next_imaginary_magnitude_power))),
                        n_choose_k(local_stack, output_stack, integer_from_size_t(local_stack, i),
                            integer_from_size_t(local_stack, j))), active_term_count));
            real_magnitude_power = rational_divide(local_stack, output_stack, real_magnitude_power,
                real_magnitude_bound);
            real_interval_scale = rational_max(output_stack, local_stack, real_interval_scale,
                rational_multiply(local_stack, output_stack, shared_scale_component,
                    real_magnitude_power));
            imaginary_interval_scale = rational_max(output_stack, local_stack,
                imaginary_interval_scale, rational_multiply(local_stack, output_stack,
                    shared_scale_component, imaginary_magnitude_power));
            imaginary_magnitude_power = next_imaginary_magnitude_power;
            POINTER_SWAP(active_term_count, inactive_term_count);
        }
        imaginary_interval_scale = rational_max(output_stack, local_stack, imaginary_interval_scale,
            rational_doubled(local_stack, output_stack,
                rational_multiply(local_stack, output_stack, imaginary_magnitude_power,
                    rational_integer_multiply(local_stack, output_stack, coefficient_magnitude,
                        active_term_count))));
    }
    interval_size = rational_place_value(local_stack, output_stack, interval_size);
    argument_real_part_estimate =
        number_float_real_part_estimate(integer_pool, local_stack, output_stack, argument,
            rational_divide(local_stack, output_stack, interval_size, real_interval_scale));
    argument_imaginary_part_estimate =
        number_float_imaginary_part_estimate(integer_pool, local_stack, output_stack, argument,
            rational_divide(local_stack, output_stack, interval_size, imaginary_interval_scale));
    struct RationalInterval evaluation_real_part = { a->coefficients[0], a->coefficients[0] };
    struct RationalInterval evaluation_imaginary_part = { &rational_zero, &rational_zero };
    struct Rational*real_part_min =
        float_to_rational(local_stack, output_stack, &argument_real_part_estimate->min);
    struct Rational*real_part_max =
        float_to_rational(local_stack, output_stack, &argument_real_part_estimate->max);
    struct Rational*imaginary_part_min =
        float_to_rational(local_stack, output_stack, &argument_imaginary_part_estimate->min);
    struct Rational*imaginary_part_max =
        float_to_rational(local_stack, output_stack, &argument_imaginary_part_estimate->max);
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct Integer*term_degree = integer_from_size_t(local_stack, i);
        struct Rational*real_part_min_power =
            rational_exponentiate(local_stack, output_stack, real_part_min, term_degree);
        struct Rational*real_part_max_power =
            rational_exponentiate(local_stack, output_stack, real_part_max, term_degree);
        struct Rational*imaginary_part_min_power = &rational_one;
        struct Rational*imaginary_part_max_power = &rational_one;
        for (size_t j = 0; j <= i; ++j)
        {
            struct Rational*coefficient = rational_integer_multiply(local_stack, output_stack,
                a->coefficients[i], n_choose_k(local_stack, output_stack, term_degree,
                    integer_from_size_t(local_stack, j)));
            struct Rational*term_bound_candidate_a =
                rational_multiply(local_stack, output_stack, coefficient, real_part_min_power);
            struct Rational*term_bound_candidate_b =
                rational_multiply(local_stack, output_stack, coefficient, real_part_max_power);
            struct Rational*term_bound_candidate_c = rational_multiply(local_stack, output_stack,
                term_bound_candidate_a, imaginary_part_min_power);
            struct Rational*term_bound_candidate_d = rational_multiply(local_stack, output_stack,
                term_bound_candidate_b, imaginary_part_min_power);
            term_bound_candidate_a = rational_multiply(local_stack, output_stack,
                term_bound_candidate_a, imaginary_part_max_power);
            term_bound_candidate_b = rational_multiply(local_stack, output_stack,
                term_bound_candidate_b, imaginary_part_max_power);
            struct RationalInterval term_interval =
                { term_bound_candidate_a, term_bound_candidate_a };
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_b);
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_c);
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_d);
            size_t degree_remainder = j % 4;
            switch (degree_remainder)
            {
            case 0:
                rational_add(local_stack, output_stack, evaluation_real_part.min,
                    term_interval.min);
                rational_add(local_stack, output_stack, evaluation_real_part.max,
                    term_interval.max);
                break;
            case 1:
                rational_add(local_stack, output_stack, evaluation_imaginary_part.min,
                    term_interval.min);
                rational_add(local_stack, output_stack, evaluation_imaginary_part.max,
                    term_interval.max);
                break;
            case 2:
                rational_subtract(local_stack, output_stack, evaluation_real_part.min,
                    term_interval.max);
                rational_subtract(local_stack, output_stack, evaluation_real_part.max,
                    term_interval.min);
                break;
            case 3:
                rational_subtract(local_stack, output_stack, evaluation_imaginary_part.min,
                    term_interval.max);
                rational_subtract(local_stack, output_stack, evaluation_imaginary_part.max,
                    term_interval.min);
                break;
            }
            real_part_min_power =
                rational_divide(local_stack, output_stack, real_part_min_power, real_part_min);
            real_part_max_power =
                rational_divide(local_stack, output_stack, real_part_max_power, real_part_max);
            imaginary_part_min_power = rational_multiply(local_stack, output_stack,
                imaginary_part_min_power, imaginary_part_min);
            imaginary_part_max_power = rational_multiply(local_stack, output_stack,
                imaginary_part_max_power, imaginary_part_max);
        }
    }
    rational_interval_to_float_interval(output_stack, local_stack, real_estimate_min_out,
        real_estimate_max_out, &evaluation_real_part, interval_size);
    rational_interval_to_float_interval(output_stack, local_stack, imaginary_estimate_min_out,
        imaginary_estimate_max_out, &evaluation_imaginary_part, interval_size);
    local_stack->cursor = local_stack_savepoint;
}

void number_calculate_minimal_polynomial_from_annulling_polynomial(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct RationalPolynomial*annulling_polynomial,
    struct Number*a)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct RationalPolynomial**candidates = stack_slot_new(stack_a,
        annulling_polynomial->coefficient_count * sizeof(struct RationalPolynomial*),
        _Alignof(struct RationalPolynomial*));
    size_t candidate_count =
        rational_polynomial_factor(stack_a, stack_b, annulling_polynomial, candidates);
    if (candidate_count == 1)
    {
        a->minimal_polynomial =
            rational_polynomial_copy_to_number_memory(integer_pool, candidates[0]);
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    struct Rational*interval_size = &rational_one;
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct Float*real_estimate_min;
            struct Float*real_estimate_max;
            struct Float*imaginary_estimate_min;
            struct Float*imaginary_estimate_max;
            rational_polynomial_estimate_evaluation(integer_pool, stack_a, stack_b,
                &real_estimate_min, &real_estimate_max, &imaginary_estimate_min,
                &imaginary_estimate_max, candidates[i], a, interval_size);
            if (float_compare(stack_a, stack_b, &float_zero, real_estimate_min) < 0 ||
                float_compare(stack_a, stack_b, real_estimate_max, &float_zero) < 0 ||
                float_compare(stack_a, stack_b, &float_zero, imaginary_estimate_min) < 0 ||
                float_compare(stack_a, stack_b, imaginary_estimate_max, &float_zero) < 0)
            {
                candidate_count -= 1;
                candidates[i] = candidates[candidate_count];
                if (candidate_count == 1)
                {
                    a->minimal_polynomial =
                        rational_polynomial_copy_to_number_memory(integer_pool, candidates[0]);
                    stack_a->cursor = stack_a_savepoint;
                    return;
                }
            }
            else
            {
                ++i;
            }
        }
        interval_size->denominator = integer_doubled(stack_a, interval_size->denominator);
    }
}

void number_calculate_minimal_polynomial_from_algebraic_form(struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct AlgebraicNumber*algebraic_form,
    struct RationalPolynomial*generator_minimal_polynomials[2])
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t term_count_upper_bound = (generator_minimal_polynomials[0]->coefficient_count - 1) *
        (generator_minimal_polynomials[1]->coefficient_count - 1) + 1;
    struct Matrix matrix;
    matrix.rows = stack_slot_new(stack_a, term_count_upper_bound * sizeof(struct Rational**),
        _Alignof(struct Rational**));
    matrix.height = 0;
    size_t(**terms_present)[2] = stack_slot_new(stack_a,
        term_count_upper_bound * sizeof(size_t(*)[2]), _Alignof(size_t(*)[2]));
    size_t terms_present_count = 1;
    struct AlgebraicNumber*power = STACK_SLOT_NEW(stack_a, struct AlgebraicNumber);
    power->term_coefficient = &rational_one;
    power->generator_degrees[0] = 0;
    power->generator_degrees[1] = 0;
    power->next_term = 0;
    struct AlgebraicNumber**powers =
        stack_slot_new(stack_a, term_count_upper_bound * sizeof(struct AlgebraicNumber*),
            _Alignof(struct AlgebraicNumber*));
    bool constant_is_present = false;
    while (!constant_is_present || matrix.height < terms_present_count)
    {
        power = algebraic_number_multiply(stack_a, stack_b, power, algebraic_form,
            generator_minimal_polynomials);
        powers[matrix.height] = power;
        struct AlgebraicNumber*power_term = powers[matrix.height];
        matrix.rows[matrix.height] = stack_slot_new(stack_a,
            term_count_upper_bound * sizeof(struct Rational*), _Alignof(struct Rational*));
        for (size_t i = 0; i < term_count_upper_bound; ++i)
        {
            matrix.rows[matrix.height][i] = &rational_zero;
        }
        while (power_term)
        {
            for (size_t i = 0; i < terms_present_count; ++i)
            {
                if (degree_compare(*terms_present[i], power_term->generator_degrees) == 0)
                {
                    matrix.rows[matrix.height][i] = power_term->term_coefficient;
                    goto degree_already_present;
                }
            }
            terms_present[terms_present_count] = &power_term->generator_degrees;
            matrix.rows[matrix.height][matrix.height] = power_term->term_coefficient;
        degree_already_present:
            power_term = power_term->next_term;
        }
        if (matrix.rows[matrix.height][0]->numerator->value_count)
        {
            constant_is_present = true;
        }
        ++matrix.height;
    }
    matrix.width = matrix.height;
    struct RationalPolynomial**augmentation = stack_slot_new(stack_a,
        matrix.height * sizeof(struct RationalPolynomial*), _Alignof(struct RationalPolynomial*));
    for (size_t i = 0; i < matrix.height; ++i)
    {
        array_reverse(matrix.rows[i], matrix.width);
        augmentation[i] = polynomial_slot_new(stack_a, i + 2);
        for (size_t j = 0; j <= i; ++j)
        {
            augmentation[i]->coefficients[j] = &rational_zero;
        }
        augmentation[i]->coefficients[i + 1] = &rational_one;
    }
    matrix_row_echelon_form(stack_a, stack_b, &matrix, augmentation);
    struct RationalPolynomial*annulling_polynomial = augmentation[matrix.height - 1];
    annulling_polynomial->coefficients[0] = rational_subtract(stack_a, stack_b,
        annulling_polynomial->coefficients[0], matrix.rows[matrix.height - 1][matrix.width - 1]);
    number_calculate_minimal_polynomial_from_annulling_polynomial(integer_pool, stack_b, stack_a,
        annulling_polynomial, a);
    stack_a->cursor = stack_a_savepoint;
}

struct Number*number_add(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            number_free(number_pool, integer_pool, a);
            return b;
        }
        if (b->operation == 'r')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct Number*out = number_slot_new(number_pool);
            out->operation = 'r';
            out->value = *rational_add(stack_a, stack_b, &a->value, &b->value);
            rational_move_to_pool(integer_pool, &out->value);
            number_free(number_pool, integer_pool, a);
            number_free(number_pool, integer_pool, b);
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
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
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
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
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
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '+':
        {
            //Placeholder.
            struct Number*out = number_slot_new(number_pool);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    return number_add(number_pool, integer_pool, stack_a, stack_b, b, a);
}

struct Number*number_multiply(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Number*b);

struct Number*number_rational_multiply(struct NumberPool*number_pool,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct Rational*b)
{
    if (b->numerator->value_count == 0)
    {
        struct Number*out = number_slot_new(number_pool);
        out->operation = 'r';
        out->value.numerator = &zero;
        out->value.denominator = pool_integer_new(integer_pool, 1, 1);
        number_free(number_pool, integer_pool, a);
        return out;
    }
    if (integer_equals(b->numerator, &one) && integer_equals(b->denominator, &one))
    {
        return a;
    }
    switch (a->operation)
    {
    case 'r':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_slot_new(number_pool);
        out->operation = 'r';
        out->value = *rational_multiply(stack_a, stack_b, &a->value, b);
        rational_move_to_pool(integer_pool, &out->value);
        number_free(number_pool, integer_pool, a);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '^':
    {
        struct Number*out = number_slot_new(number_pool);
        out->operation = '*';
        out->left = number_slot_new(number_pool);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(integer_pool, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(integer_pool, b->denominator);
        out->right = a;
        return out;
    }
    case '*':
    {
        if (a->left->operation == 'r')
        {
            struct Number*out = number_multiply(number_pool, integer_pool, stack_a, stack_b,
                number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, a->left, b),
                a->right);
            number_node_free(number_pool, integer_pool, a);
            return out;
        }
        struct Number*out = number_slot_new(number_pool);
        out->operation = '*';
        out->left = number_slot_new(number_pool);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(integer_pool, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(integer_pool, b->denominator);
        out->right = a;
        return out;
    }
    case '+':
    {
        struct Number*out = number_slot_new(number_pool);
        out->operation = '+';
        out->left = number_slot_new(number_pool);
        out->left =
            number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, a->left, b);
        out->right =
            number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, a->right, b);
        number_node_free(number_pool, integer_pool, a);
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_exponentiate(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*base, struct Number*exponent);

struct Number**get_roots_of_unity(struct Stack*stack_a, struct Stack*stack_b,
    struct Integer*degree);

struct Number*number_consolidate_product(struct NumberPool*number_pool,
    struct IntegerPool*integer_pool, struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b)
{
    if (b->operation == 'r')
    {
        struct Number*out =
            number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, a, &b->value);
        number_free(number_pool, integer_pool, b);
        return out;
    }
    switch (a->operation)
    {
    case 'r':
    {
        struct Number*out =
            number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, b, &a->value);
        number_free(number_pool, integer_pool, a);
        return out;
    }
    case '^':
    {
        if (b->operation == '^')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct ExtendedGCDInfo info;
            integer_extended_gcd(stack_a, stack_b, &info, a->right->value.denominator,
                b->right->value.denominator);
            struct Integer*product_index = integer_euclidean_quotient(stack_a, stack_b,
                integer_multiply(stack_a, stack_b, a->right->value.denominator,
                    b->right->value.denominator), info.gcd);
            struct Number*exponent = number_slot_new(number_pool);
            exponent->operation = 'r';
            exponent->value.numerator =
                integer_copy_to_pool(integer_pool, integer_magnitude(stack_a, info.b_over_gcd));
            exponent->value.denominator = pool_integer_new(integer_pool, 1, 1);
            struct Number*radicand_factor_a = number_exponentiate(number_pool, integer_pool,
                stack_a, stack_b, number_copy(number_pool, integer_pool, a->left), exponent);
            exponent = number_slot_new(number_pool);
            exponent->operation = 'r';
            exponent->value.numerator =
                integer_copy_to_pool(integer_pool, integer_magnitude(stack_a, info.a_over_gcd));
            exponent->value.denominator = pool_integer_new(integer_pool, 1, 1);
            struct Number*radicand_factor_b = number_exponentiate(number_pool, integer_pool,
                stack_a, stack_b, number_copy(number_pool, integer_pool, b->left), exponent);
            struct RationalInterval radicand_factor_a_argument_estimate;
            number_rational_argument_estimate(integer_pool, stack_a, stack_b,
                &radicand_factor_a_argument_estimate, radicand_factor_a, &rational_one);
            struct RationalInterval radicand_factor_b_argument_estimate;
            number_rational_argument_estimate(integer_pool, stack_a, stack_b,
                &radicand_factor_b_argument_estimate, radicand_factor_b, &rational_one);
            if (!integer_equals(a->right->value.denominator, b->right->value.denominator))
            {
                struct Rational*radicand_factor_argument_interval_size =
                    rational_integer_divide(stack_a, stack_b, &pi_estimate_min, product_index);
                struct RationalInterval a_argument_estimate;
                number_rational_argument_estimate(integer_pool, stack_a, stack_b,
                    &a_argument_estimate, a, radicand_factor_argument_interval_size);
                struct RationalInterval b_argument_estimate;
                number_rational_argument_estimate(integer_pool, stack_a, stack_b,
                    &b_argument_estimate, b, radicand_factor_argument_interval_size);
                if (rational_compare(stack_a, stack_b,
                    rational_integer_divide(stack_a, stack_b,
                        radicand_factor_a_argument_estimate.max, product_index),
                    a_argument_estimate.min) < 0 ||
                    rational_compare(stack_a, stack_b,
                        rational_integer_divide(stack_a, stack_b,
                            radicand_factor_b_argument_estimate.max, product_index),
                        b_argument_estimate.min) < 0)
                {
                    number_free(number_pool, integer_pool, radicand_factor_a);
                    number_free(number_pool, integer_pool, radicand_factor_b);
                    stack_a->cursor = stack_a_savepoint;
                    return number_product_consolidation_failed;
                }
            }
            struct Number*product_radicand = number_multiply(number_pool, integer_pool, stack_a,
                stack_b, radicand_factor_a, radicand_factor_b);
            exponent = number_slot_new(number_pool);
            exponent->operation = 'r';
            exponent->value.numerator = pool_integer_new(integer_pool, 1, 1);
            exponent->value.denominator = integer_copy_to_pool(integer_pool, product_index);
            struct Number*out = number_exponentiate(number_pool, integer_pool, stack_a, stack_b,
                product_radicand, exponent);
            struct RationalInterval product_radicand_rational_argument_estimate;
            number_rational_argument_estimate(integer_pool, stack_a, stack_b,
                &product_radicand_rational_argument_estimate, product_radicand, &pi_estimate_min);
            if (rational_compare(stack_a, stack_b, product_radicand_rational_argument_estimate.max,
                rational_add(stack_a, stack_b, radicand_factor_a_argument_estimate.min,
                    radicand_factor_b_argument_estimate.min)) < 0)
            {
                out = number_multiply(number_pool, integer_pool, stack_a, stack_b, out,
                    number_copy(number_pool, integer_pool,
                        get_roots_of_unity(stack_a, stack_b, product_index)[1]));
            }
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case '^':
        {
            struct Number*factor =
                number_consolidate_product(number_pool, integer_pool, stack_a, stack_b, a->left, b);
            if (factor)
            {
                return number_multiply(number_pool, integer_pool, stack_a, stack_b, a->right,
                    factor);
            }
            factor = number_consolidate_product(number_pool, integer_pool, stack_a, stack_b,
                a->right, b);
            if (factor)
            {
                return number_multiply(number_pool, integer_pool, stack_a, stack_b, a->left,
                    factor);
            }
            return number_product_consolidation_failed;
        }
        case '*':
        {
            struct Number*out =
                number_multiply(number_pool, integer_pool, stack_a, stack_b, a->left, b);
            if (!out)
            {
                return 0;
            }
            out = number_multiply(number_pool, integer_pool, stack_a, stack_b, out, a->right);
            number_node_free(number_pool, integer_pool, a);
            return out;
        }
        }
        break;
    case '+':
    {
        struct Number*b_copy = number_copy(number_pool, integer_pool, b);
        struct Number*left =
            number_multiply(number_pool, integer_pool, stack_a, stack_b, a->left, b);
        if (!left)
        {
            return 0;
        }
        struct Number*right =
            number_multiply(number_pool, integer_pool, stack_a, stack_b, a->right, b_copy);
        if (!right)
        {
            return 0;
        }
        number_node_free(number_pool, integer_pool, a);
        return number_add(number_pool, integer_pool, stack_a, stack_b, left, right);
    }
    }
    return number_multiply(number_pool, integer_pool, stack_a, stack_b, b, a);
}

struct Number*number_multiply(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Number*b)
{
    struct Number*out =
        number_consolidate_product(number_pool, integer_pool, stack_a, stack_b, a, b);
    if (out == number_product_consolidation_failed)
    {
        out = number_slot_new(number_pool);
        out->operation = '*';
        out->left = a;
        out->right = b;
    }
    return out;
}

struct Number*number_divide(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*dividend, struct Number*divisor);

void number_calculate_minimal_polynomial(struct Stack*stack_a, struct Stack*stack_b,
    struct IntegerPool*integer_pool, struct Number*a);

void number_calculate_conjugates(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a);

struct Number*number_reciprocal(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        if (!a->value.numerator->value_count)
        {
            printf("Tried to divide by 0.");
            return number_divide_by_zero_error;
        }
        POINTER_SWAP(a->value.numerator, a->value.denominator);
        return a;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        pool_integer_free(integer_pool, a->right->value.numerator);
        a->right->value.numerator = integer_add(stack_a, a->right->value.denominator, &INT(1, -));
        integer_move_to_pool(integer_pool, &a->right->value.numerator);
        stack_a->cursor = stack_a_savepoint;
        struct Number*radicand_copy = number_copy(number_pool, integer_pool, a->left);
        return number_divide(number_pool, integer_pool, stack_a, stack_b,
            number_exponentiate(number_pool, integer_pool, stack_a, stack_b, a->left, a->right),
            radicand_copy);
    }
    case '*':
    {
        struct Number*out = number_multiply(number_pool, integer_pool, stack_a, stack_b,
            number_reciprocal(number_pool, integer_pool, stack_a, stack_b, a->left),
            number_reciprocal(number_pool, integer_pool, stack_a, stack_b, a->right));
        number_node_free(number_pool, integer_pool, a);
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        number_calculate_conjugates(number_pool, integer_pool, stack_a, stack_b, a);
        struct Number*out = number_slot_new(number_pool);
        out->operation = 'r';
        out->value.numerator = pool_integer_new(integer_pool, 1, 1);
        out->value.denominator = pool_integer_new(integer_pool, 1, 1);
        struct Number*conjugate = a->next;
        while (conjugate)
        {
            out = number_multiply(number_pool, integer_pool, stack_a, stack_b, out,
                number_copy(number_pool, integer_pool, conjugate));
            conjugate = conjugate->next;
        }
        number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a);
        struct Rational*coefficient =
            rational_reciprocal(stack_a, stack_b, a->minimal_polynomial->coefficients[0], 0);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out =
            number_rational_multiply(number_pool, integer_pool, stack_a, stack_b, out, coefficient);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_divide(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal =
        number_reciprocal(number_pool, integer_pool, stack_a, stack_b, divisor);
    if (reciprocal == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    return number_multiply(number_pool, integer_pool, stack_a, stack_b, dividend, reciprocal);
}

struct Number*number_exponentiate(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*base, struct Number*exponent)
{
    if (exponent->operation != 'r')
    {
        printf("The input expression contains an exponentiation whose exponent is not both real "
            "and rational; this program doesn't handle transcendental numbers.");
        return 0;
    }
    if (exponent->value.numerator->sign < 0)
    {
        base = number_reciprocal(number_pool, integer_pool, stack_a, stack_b, base);
        if (base == number_divide_by_zero_error)
        {
            return number_divide_by_zero_error;
        }
        exponent->value.numerator->sign = 1;
        return number_exponentiate(number_pool, integer_pool, stack_a, stack_b, base, exponent);
    }
    switch (base->operation)
    {
    case 'r':
        if (base->value.numerator->value_count == 0)
        {
            struct Number*out = number_slot_new(number_pool);
            number_free(number_pool, integer_pool, base);
            number_free(number_pool, integer_pool, exponent);
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
                struct Number*out = number_slot_new(number_pool);
                out->operation = 'r';
                out->value.denominator = pool_integer_new(integer_pool, 1, 1);
                out->value.numerator = integer_exponentiate(stack_a, stack_b, base->value.numerator,
                    exponent->value.numerator);
                integer_move_to_pool(integer_pool, &out->value.numerator);
                number_free(number_pool, integer_pool, base);
                number_free(number_pool, integer_pool, exponent);
                stack_a->cursor = stack_a_savepoint;
                return out;
            }
            struct Integer*radicand = integer_exponentiate(stack_a, stack_b, base->value.numerator,
                exponent->value.numerator);
            pool_integer_free(integer_pool, exponent->value.numerator);
            exponent->value.numerator = pool_integer_new(integer_pool, 1, 1);
            number_free(number_pool, integer_pool, base);
            int8_t radicand_sign = radicand->sign;
            radicand->sign = 1;
            struct Factor*factors;
            size_t factor_count = integer_factor(stack_a, stack_b, &factors, radicand);
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
                    number_free(number_pool, integer_pool, exponent);
                    struct Number*out = number_slot_new(number_pool);
                    out->operation = 'r';
                    out->value.denominator = pool_integer_new(integer_pool, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move_to_pool(integer_pool, &out->value.numerator);
                    stack_a->cursor = stack_a_savepoint;
                    return out;
                }
                pool_integer_free(integer_pool, exponent->value.denominator);
                exponent->value.denominator = reduced_degree;
                integer_move_to_pool(integer_pool, &exponent->value.denominator);
            }
            struct Number*number_coefficient = number_slot_new(number_pool);
            number_coefficient->operation = 'r';
            number_coefficient->value.numerator = coefficient;
            integer_move_to_pool(integer_pool, &number_coefficient->value.numerator);
            number_coefficient->value.denominator = pool_integer_new(integer_pool, 1, 1);
            struct Number*number_radicand = number_slot_new(number_pool);
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
            struct Number*surd = number_slot_new(number_pool);
            surd->operation = '^';
            surd->left = number_radicand;
            surd->right = exponent;
            stack_a->cursor = stack_a_savepoint;
            return number_multiply(number_pool, integer_pool, stack_a, stack_b, number_coefficient,
                surd);
        }
        struct Number*new_denominator = number_slot_new(number_pool);
        new_denominator->operation = 'r';
        new_denominator->value.numerator = pool_integer_new(integer_pool, 1, 1);
        new_denominator->value.denominator = integer_exponentiate(stack_a, stack_b,
            base->value.denominator, exponent->value.numerator);
        integer_move_to_pool(integer_pool, &new_denominator->value.denominator);
        struct Number*new_numerator_base = number_slot_new(number_pool);
        new_numerator_base->operation = 'r';
        new_numerator_base->value.denominator = pool_integer_new(integer_pool, 1, 1);
        new_numerator_base->value.numerator =
            integer_multiply(stack_a, stack_b, base->value.numerator,
                integer_exponentiate(stack_a, stack_b, base->value.denominator,
                    integer_add(stack_a, exponent->value.denominator, &INT(1, -))));
        integer_move_to_pool(integer_pool, &new_numerator_base->value.numerator);
        struct Number*new_numerator = number_exponentiate(number_pool, integer_pool, stack_a,
            stack_b, new_numerator_base, exponent);
        struct Number*out = number_multiply(number_pool, integer_pool, stack_a, stack_b,
            new_numerator, new_denominator);
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
        number_free(number_pool, integer_pool, base->right);
        struct Number*out =
            number_exponentiate(number_pool, integer_pool, stack_a, stack_b, base->left, exponent);
        number_node_free(number_pool, integer_pool, base);
        return out;
    }
    case '*':
    {
        struct Number*exponent_copy = number_copy(number_pool, integer_pool, exponent);
        struct Number*out = number_multiply(number_pool, integer_pool, stack_a, stack_b,
            number_exponentiate(number_pool, integer_pool, stack_a, stack_b, base->left, exponent),
            number_exponentiate(number_pool, integer_pool, stack_a, stack_b, base->right,
                exponent_copy));
        number_node_free(number_pool, integer_pool, base);
        return out;
    }
    case '+':
        crash("number_exponentiate case not yet implemented.");
    default:
        crash("Number operation not recognized.");
    }
}

void quick_sort(struct Stack*stack_a, struct Stack*stack_b, struct Number**roots,
    struct FloatInterval**argument_estimates, size_t start, size_t end)
{
    if (start < end)
    {
        size_t final_pivot_index = 0;
        for (size_t i = 0; i < end; ++i)
        {
            if (float_compare(stack_a, stack_b, &argument_estimates[i]->max,
                &argument_estimates[end]->min) < 0)
            {
                POINTER_SWAP(argument_estimates[i], argument_estimates[final_pivot_index]);
                POINTER_SWAP(roots[i], roots[final_pivot_index]);
                ++i;
            }
        }
        POINTER_SWAP(argument_estimates[final_pivot_index], argument_estimates[end]);
        POINTER_SWAP(roots[final_pivot_index], roots[end]);
        if (final_pivot_index > 1)
        {
            quick_sort(stack_a, stack_b, roots, argument_estimates, start, final_pivot_index - 1);
        }
        quick_sort(stack_a, stack_b, roots, argument_estimates, final_pivot_index + 1, end);
    }
}

void roots_of_unity_sort(struct IntegerPool*permanent_integer_pool, struct Stack*stack_a,
    struct Stack*stack_b, struct Number**roots, struct Integer*degree)
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t degree_size_t = integer_to_size_t(degree);
    struct Rational*interval_size =
        rational_integer_divide(stack_a, stack_b, &pi_estimate_min, degree);
    struct FloatInterval**argument_estimates = stack_slot_new(stack_a,
        (degree_size_t - 1) * sizeof(struct FloatInterval*), _Alignof(struct FloatInterval*));
    for (size_t i = 1; i < degree_size_t; ++i)
    {
        argument_estimates[i - 1] = number_float_argument_estimate(permanent_integer_pool, stack_a,
            stack_b, roots[i], interval_size);
    }
    quick_sort(stack_a, stack_b, &roots[1], argument_estimates, 0, degree_size_t - 1);
    stack_a->cursor = stack_a_savepoint;
}

struct Number**get_roots_of_unity(struct Stack*stack_a, struct Stack*stack_b, struct Integer*degree)
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t degree_size_t = integer_to_size_t(degree);
    size_t degree_minus_one_size_t = degree_size_t - 1;
    struct Integer*degree_minus_one = integer_from_size_t(stack_a, degree_minus_one_size_t);
    struct Number**out = roots_of_unity + integer_to_size_t(
        integer_half(stack_a, integer_multiply(stack_a, stack_b, degree, degree_minus_one)));
    if (out[0])
    {
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    size_t degree_square_root = integer_to_size_t(integer_square_root(stack_a, stack_b, degree));
    size_t prime_index = 0;
    while (true)
    {
        struct Integer*prime = get_prime(stack_a, stack_b, prime_index);
        size_t prime_size_t = integer_to_size_t(prime);
        if (prime_size_t > degree_square_root)
        {
            break;
        }
        if (degree_size_t % prime_size_t != 0)
        {
            size_t quotient = degree_size_t / prime_size_t;
            struct Number**prime_degree_roots = get_roots_of_unity(stack_a, stack_b, prime);
            struct Number**quotient_degree_roots =
                get_roots_of_unity(stack_a, stack_b, integer_from_size_t(stack_a, quotient));
            for (size_t i = 0; i < prime_size_t; ++i)
            {
                for (size_t j = 0; j < quotient; ++j)
                {
                    struct Number*exponent = number_slot_new(&permanent_number_pool);
                    exponent->operation = 'r';
                    exponent->value.numerator = pool_integer_new(&permanent_integer_pool, 1, 1);
                    exponent->value.denominator =
                        integer_copy_to_pool(&permanent_integer_pool, prime);
                    out[i] = number_multiply(&permanent_number_pool,
                        &permanent_integer_pool, stack_a, stack_b,
                        number_copy(&permanent_number_pool, &permanent_integer_pool,
                            prime_degree_roots[i]),
                        number_exponentiate(&permanent_number_pool, &permanent_integer_pool,
                            stack_a, stack_b,
                            number_copy(&permanent_number_pool, &permanent_integer_pool,
                                quotient_degree_roots[j]),
                            exponent));
                }
            }
            roots_of_unity_sort(&permanent_integer_pool, stack_a, stack_b, out, degree);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
    }
    size_t*degree_minus_one_factors;
    size_t factor_count =
        size_t_factor(stack_a, stack_b, &degree_minus_one_factors, degree_minus_one_size_t);
    struct Integer*generator = &INT(2, +);
    while (true)
    {
        generator_not_found:
        for (size_t i = 0; i < factor_count; ++i)
        {
            if (!integer_equals(integer_euclidean_remainder(stack_a, stack_b, 
                integer_exponentiate(stack_a, stack_b, generator, integer_from_size_t(stack_a,
                    degree_minus_one_size_t / degree_minus_one_factors[i])), degree), &one))
            {
                generator = integer_add(stack_a, generator, &one);
                goto generator_not_found;
            }
        }
        break;
    }
    struct RationalPolynomial*resolvent_generator_annulling_polynomials[2] =
        { polynomial_slot_new(stack_a, degree_size_t),
        polynomial_slot_new(stack_a, degree_size_t + 1) };
    resolvent_generator_annulling_polynomials[0]->coefficients[0] =
        &(struct Rational) { &INT(1, -), &one };
    resolvent_generator_annulling_polynomials[1]->coefficients[0] =
        resolvent_generator_annulling_polynomials[0]->coefficients[0];
    for (size_t i = 1; i < degree_minus_one_size_t; ++i)
    {
        resolvent_generator_annulling_polynomials[0]->coefficients[i] = &rational_zero;
        resolvent_generator_annulling_polynomials[1]->coefficients[i] = &rational_zero;
    }
    resolvent_generator_annulling_polynomials[0]->coefficients[degree_minus_one_size_t] =
        &rational_one;
    resolvent_generator_annulling_polynomials[1]->coefficients[degree_minus_one_size_t] =
        &rational_zero;
    resolvent_generator_annulling_polynomials[1]->coefficients[degree_size_t] = &rational_one;
    struct AlgebraicNumber**resolvents = stack_slot_new(stack_a,
        (degree_minus_one_size_t - 1) * sizeof(struct AlgebraicNumber*),
        _Alignof(struct AlgebraicNumber*));
    size_t resolvent_count_minus_one = degree_minus_one_size_t - 1;
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        resolvents[i] = STACK_SLOT_NEW(stack_a, struct AlgebraicNumber);
        resolvents[i]->next_term = 0;
        resolvents[i]->term_coefficient = &rational_one;
        resolvents[i]->generator_degrees[0] = 0;
        resolvents[i]->generator_degrees[1] = 1;
    }
    struct Integer*generator_power = &one;
    struct Rational**resolvent_multiples_in_terms_of_degree_minus_first_roots =
        stack_slot_new(stack_a, resolvent_count_minus_one * degree_minus_one_size_t *
            sizeof(struct Rational*), _Alignof(struct Rational*));
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        size_t resolvent_multiple_index = degree_minus_one_size_t * i;
        resolvent_multiples_in_terms_of_degree_minus_first_roots[resolvent_multiple_index] =
            &rational_zero;
        generator_power = integer_euclidean_remainder(stack_a, stack_b,
            integer_multiply(stack_a, stack_b, generator_power, generator), degree);
        size_t degree_minus_first_root_exponent = 1;
        for (size_t j = 1; j < degree_minus_one_size_t; ++j)
        {
            struct AlgebraicNumber*resolvent_term = STACK_SLOT_NEW(stack_a, struct AlgebraicNumber);
            resolvent_term->next_term = 0;
            resolvent_term->term_coefficient = &rational_one;
            resolvent_term->generator_degrees[0] = degree_minus_first_root_exponent;
            resolvent_term->generator_degrees[1] = integer_to_size_t(generator_power);
            resolvents[j - 1] =
                algebraic_number_add(stack_a, stack_b, resolvents[j - 1], resolvent_term);
            resolvent_multiples_in_terms_of_degree_minus_first_roots[resolvent_multiple_index + j] =
                &rational_zero;
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i + 1) % degree_minus_one_size_t;
        }
    }
    struct AlgebraicNumber*resolvent_power = STACK_SLOT_NEW(stack_a, struct AlgebraicNumber);
    resolvent_power->next_term = 0;
    resolvent_power->term_coefficient = &rational_one;
    resolvent_power->generator_degrees[0] = 0;
    resolvent_power->generator_degrees[1] = 0;
    for (size_t i = resolvent_count_minus_one; i-- > 0;)
    {
        resolvent_power = algebraic_number_multiply(stack_a, stack_b, resolvent_power,
            resolvents[0], resolvent_generator_annulling_polynomials);
        struct AlgebraicNumber*resolvent_product = algebraic_number_multiply(stack_a, stack_b,
            resolvent_power, resolvents[i], resolvent_generator_annulling_polynomials);
        struct AlgebraicNumber*resolvent_product_term = resolvent_product;
        while (resolvent_product_term)
        {
            switch (resolvent_product_term->generator_degrees[1])
            {
            case 0:
            {
                size_t index =
                    degree_minus_one_size_t * i + resolvent_product_term->generator_degrees[0];
                resolvent_multiples_in_terms_of_degree_minus_first_roots[index] =
                    rational_add(stack_a, stack_b, resolvent_product_term->term_coefficient,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[index]);
            }
            case 1:
            {
                size_t index =
                    degree_minus_one_size_t * i + resolvent_product_term->generator_degrees[0];
                resolvent_multiples_in_terms_of_degree_minus_first_roots[index] =
                    rational_subtract(stack_a, stack_b, resolvent_product_term->term_coefficient,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[index]);
            }
            }
            resolvent_product_term = resolvent_product_term->next_term;
        }
    }
    struct Number**degree_minus_first_roots =
        get_roots_of_unity(stack_a, stack_b, degree_minus_one);
    struct Number**resolvent_product_values = stack_slot_new(stack_a,
        resolvent_count_minus_one * sizeof(struct Number*), _Alignof(struct Number*));
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        size_t resolvent_multiple_index = degree_minus_one_size_t * i;
        resolvent_product_values[i] =
            number_copy(&permanent_number_pool, &permanent_integer_pool, roots_of_unity[0]);
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            resolvent_product_values[i] = number_add(&permanent_number_pool,
                &permanent_integer_pool, stack_a, stack_b, resolvent_product_values[i],
                number_rational_multiply(&permanent_number_pool, &permanent_integer_pool, stack_a,
                    stack_b, number_copy(&permanent_number_pool, &permanent_integer_pool,
                        degree_minus_first_roots[i]),
                    resolvent_multiples_in_terms_of_degree_minus_first_roots
                        [resolvent_multiple_index]));
        }
    }
    struct Number**resolvent_values = stack_slot_new(stack_a,
        degree_minus_one_size_t * sizeof(struct Number*), _Alignof(struct Number*));
    resolvent_values[0] = number_slot_new(&permanent_number_pool);
    resolvent_values[0]->operation = 'r';
    resolvent_values[0]->value.numerator = pool_integer_new(&permanent_integer_pool, 1, -1);
    resolvent_values[0]->value.denominator = pool_integer_new(&permanent_integer_pool, 1, 1);
    struct Number*exponent = number_slot_new(&permanent_number_pool);
    exponent->operation = 'r';
    exponent->value.numerator = pool_integer_new(&permanent_integer_pool, 1, 1);
    exponent->value.denominator = integer_copy_to_pool(&permanent_integer_pool, degree_minus_one);
    resolvent_values[1] = number_exponentiate(&permanent_number_pool, &permanent_integer_pool,
        stack_a, stack_b, resolvent_product_values[0], exponent);
    for (size_t i = 1; i < resolvent_count_minus_one; ++i)
    {
        resolvent_values[i + 1] =
            number_divide(&permanent_number_pool, &permanent_integer_pool, stack_a, stack_b,
                number_copy(&permanent_number_pool, &permanent_integer_pool, resolvent_values[1]),
                resolvent_product_values[i]);
    }
    struct Rational*degree_minus_one_reciprocal = STACK_SLOT_NEW(stack_a, struct Rational);
    degree_minus_one_reciprocal->numerator = &one;
    degree_minus_one_reciprocal->denominator = degree_minus_one;
    out[0] = roots_of_unity[0];
    for (size_t i = 1; i < degree_size_t; ++i)
    {
        out[i] = number_slot_new(&permanent_number_pool);
        out[i]->operation = 'r';
        out[i]->value.numerator = pool_integer_new(&permanent_integer_pool, 0, 0);
        out[i]->value.denominator = pool_integer_new(&permanent_integer_pool, 1, 1);
        size_t degree_minus_first_root_exponent = 0;
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            out[i] = number_add(&permanent_number_pool, &permanent_integer_pool, stack_a, stack_b,
                number_multiply(&permanent_number_pool, &permanent_integer_pool, stack_a, stack_b,
                    number_copy(&permanent_number_pool, &permanent_integer_pool,
                        degree_minus_first_roots[degree_minus_first_root_exponent]),
                    resolvent_values[j]), out[i]);
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i) % degree_minus_one_size_t;
        }
        out[i] = number_rational_multiply(&permanent_number_pool, &permanent_integer_pool, stack_a,
            stack_b, out[i], degree_minus_one_reciprocal);
    }
    roots_of_unity_sort(&permanent_integer_pool, stack_a, stack_b, out, degree);
    stack_a->cursor = stack_a_savepoint;
    return out;
}

void number_calculate_minimal_polynomial(struct Stack*stack_a, struct Stack*stack_b,
    struct IntegerPool*integer_pool, struct Number*a)
{
    if (a->minimal_polynomial)
    {
        return;
    }
    switch (a->operation)
    {
    case 'r':
        a->minimal_polynomial = polynomial_slot_new(&polynomial_stack, 2);
        a->minimal_polynomial->coefficients[0]->numerator =
            integer_copy_to_pool(integer_pool, a->value.numerator);
        a->minimal_polynomial->coefficients[0]->numerator->sign =
            -a->minimal_polynomial->coefficients[0]->numerator->sign;
        a->minimal_polynomial->coefficients[0]->denominator =
            integer_copy_to_pool(integer_pool, a->value.denominator);
        a->minimal_polynomial->coefficients[1]->numerator = pool_integer_new(integer_pool, 1, 1);
        a->minimal_polynomial->coefficients[1]->denominator =
            a->minimal_polynomial->coefficients[1]->numerator;
        return;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a->left);
        size_t surd_index = integer_to_size_t(a->left->value.denominator);
        struct RationalPolynomial*annulling_polynomial = polynomial_slot_new(stack_a,
            surd_index * (a->left->minimal_polynomial->coefficient_count - 1) + 1);
        annulling_polynomial->coefficients[0] = a->left->minimal_polynomial->coefficients[0];
        for (size_t i = 1; i < a->left->minimal_polynomial->coefficient_count; ++i)
        {
            size_t annulling_polynomial_coefficient_index = i * surd_index;
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
                a->left->minimal_polynomial->coefficients[i];
            for (size_t j = 1; j < surd_index; ++j)
            {
                annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                    &rational_zero;
            }
        }
        number_calculate_minimal_polynomial_from_annulling_polynomial(integer_pool, stack_a,
            stack_b, annulling_polynomial, a);
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    case '*':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct AlgebraicNumber*algebraic_form = STACK_SLOT_NEW(stack_a, struct AlgebraicNumber);
        algebraic_form->term_coefficient = &rational_one;
        algebraic_form->generator_degrees[0] = 1;
        algebraic_form->generator_degrees[1] = 1;
        algebraic_form->next_term = 0;
        number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a->left);
        number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a->right);
        number_calculate_minimal_polynomial_from_algebraic_form(integer_pool, stack_a, stack_b, a,
            algebraic_form, (struct RationalPolynomial*[2]) { a->left->minimal_polynomial,
                a->right->minimal_polynomial });
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    case '+':
        crash("number_minimal_polynomial case not yet implemented.");
    }
}

void number_calculate_sum_or_product_conjugates(struct Number*(operation)(struct NumberPool*,
    struct IntegerPool*, struct Stack*, struct Stack*, struct Number*, struct Number*),
    struct NumberPool*number_pool, struct IntegerPool*integer_pool, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a)
{
    number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a);
    number_calculate_conjugates(number_pool, integer_pool, stack_a, stack_b, a->left);
    number_calculate_conjugates(number_pool, integer_pool, stack_a, stack_b, a->right);
    struct Number*conjugate = a;
    struct Number*left_conjugate = a->left;
    struct Number*right_conjugate = a->right->next;
    size_t conjugate_count = 1;
    while (left_conjugate)
    {
        while (right_conjugate)
        {
            struct Number*candidate = operation(number_pool, integer_pool, stack_a, stack_b,
                number_copy(number_pool, integer_pool, left_conjugate),
                number_copy(number_pool, integer_pool, right_conjugate));
            number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, candidate);
            if (!rational_polynomial_equals(stack_a, stack_b, a->minimal_polynomial,
                candidate->minimal_polynomial))
            {
                conjugate->next = candidate;
                ++conjugate_count;
                if (conjugate_count == a->minimal_polynomial->coefficient_count - 1)
                {
                    return;
                }
                conjugate = conjugate->next;
            }
            else
            {
                number_free(number_pool, integer_pool, candidate);
            }
            right_conjugate = right_conjugate->next;
        }
        left_conjugate = left_conjugate->next;
        right_conjugate = a->right;
    }
    crash("Not enough conjugates found.");
}

void number_calculate_conjugates(struct NumberPool*number_pool, struct IntegerPool*integer_pool,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    if (a->next)
    {
        return;
    }
    switch (a->operation)
    {
    case 'r':
        return;
    case '^':
    {
        number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, a);
        number_calculate_conjugates(number_pool, integer_pool, stack_a, stack_b, a->left);
        size_t roots_of_unity_count = integer_to_size_t(a->right->value.denominator);
        struct Number**roots_of_unity =
            get_roots_of_unity(stack_a, stack_b, a->right->value.denominator);
        struct Number*radicand_conjugate = a->left->next;
        struct Number*conjugate = a;
        size_t conjugate_count = 1;        
        for (size_t i = 0; i < roots_of_unity_count; ++i)
        {
            while (radicand_conjugate)
            {
                struct Number*surd = number_slot_new(number_pool);
                surd->operation = '^';
                surd->left = number_copy(number_pool, integer_pool, radicand_conjugate);
                surd->right = number_copy(number_pool, integer_pool, a->left);
                struct Number*candidate = number_multiply(number_pool, integer_pool, stack_a,
                    stack_b, surd, number_copy(number_pool, integer_pool, radicand_conjugate));
                number_calculate_minimal_polynomial(stack_a, stack_b, integer_pool, candidate);
                if (!rational_polynomial_equals(stack_a, stack_b, a->minimal_polynomial,
                    candidate->minimal_polynomial))
                {
                    conjugate->next = candidate;
                    ++conjugate_count;
                    if (conjugate_count == a->minimal_polynomial->coefficient_count - 1)
                    {
                        return;
                    }
                    conjugate = conjugate->next;
                }
                else
                {
                    number_free(number_pool, integer_pool, candidate);
                }
                radicand_conjugate = radicand_conjugate->next;
            }
            radicand_conjugate = a->left;
        }
        crash("Not enough conjugates found.");
    }
    case '*':
    {
        number_calculate_sum_or_product_conjugates(number_multiply, number_pool, integer_pool,
            stack_a, stack_b, a);
        return;
    }
    case '+':
    {
        number_calculate_sum_or_product_conjugates(number_add, number_pool, integer_pool, stack_a,
            stack_b, a);
    }
    }
}

struct Number*number_evaluate(struct NumberPool*number_pool,
    struct IntegerPool*transient_integer_pool, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a)
{
    a->minimal_polynomial = 0;
    a->next = 0;
    if (a->operation == 'r')
    {
        return a;
    }
    a->left = number_evaluate(number_pool, transient_integer_pool, stack_a, stack_b, a->left);
    if (a->left == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    a->right = number_evaluate(number_pool, transient_integer_pool, stack_a, stack_b, a->right);
    if (a->right == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    switch (a->operation)
    {
    case '+':
        return number_add(number_pool, transient_integer_pool, stack_a, stack_b, a->left, a->right);
    case '-':
    {
        struct Number*negative = number_slot_new(number_pool);
        negative->operation = 'r';
        negative->value.numerator = pool_integer_new(transient_integer_pool, 1, -1);
        negative->value.denominator = pool_integer_new(transient_integer_pool, 1, 1);
        return number_add(number_pool, transient_integer_pool, stack_a, stack_b, a->left,
            number_multiply(number_pool, transient_integer_pool, stack_a, stack_b, negative,
                a->right));
    }
    case '*':
        return number_multiply(number_pool, transient_integer_pool, stack_a, stack_b, a->left,
            a->right);
    case '/':
        return number_divide(number_pool, transient_integer_pool, stack_a, stack_b, a->left,
            a->right);
    case '^':
        return number_exponentiate(number_pool, transient_integer_pool, stack_a, stack_b, a->left,
            a->right);
    default:
        crash("Number operation not recognized.");
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

void init(struct IntegerPool*transient_integer_pool, struct NumberPool*transient_number_pool,
    struct Stack*stack_a, struct Stack*stack_b)
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwAllocationGranularity;
    size_t arena_size = page_size * (((size_t)system_info.lpMaximumApplicationAddress -
        (size_t)system_info.lpMinimumApplicationAddress) / (8 * page_size));
    roots_of_unity = VirtualAlloc(system_info.lpMinimumApplicationAddress, page_size,
        MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    number_pool_new(&permanent_number_pool, (size_t)roots_of_unity + page_size,
        (size_t)roots_of_unity + arena_size);
    primes = VirtualAlloc((void*)permanent_number_pool.stack.end, page_size,
        MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    stack_new(&permanent_stack, permanent_number_pool.stack.end + page_size,
        permanent_number_pool.stack.end + arena_size);
    integer_pool_memory_start = permanent_stack.end;
    integer_pool_new(&permanent_integer_pool, integer_pool_memory_start,
        integer_pool_memory_start + arena_size);
    integer_pool_memory_end = permanent_integer_pool.slot_page_stack.end + arena_size;
    integer_pool_new(transient_integer_pool, permanent_integer_pool.slot_page_stack.end,
        integer_pool_memory_end);
    number_pool_new(transient_number_pool, integer_pool_memory_end,
        integer_pool_memory_end + arena_size);
    stack_new(&polynomial_stack, transient_number_pool->stack.end,
        transient_number_pool->stack.end + arena_size);
    stack_new(stack_a, polynomial_stack.end, polynomial_stack.end + arena_size);
    stack_new(stack_b, stack_a->end, stack_a->end + arena_size);
    roots_of_unity[0] = STACK_SLOT_NEW(&permanent_stack, struct Number);
    roots_of_unity[0]->operation = 'r';
    roots_of_unity[0]->value = rational_one;
    roots_of_unity[0]->magnitude_estimate.min = float_one;
    roots_of_unity[0]->magnitude_estimate.max = float_one;
    roots_of_unity[1] = roots_of_unity[0];
    roots_of_unity[2] = STACK_SLOT_NEW(&permanent_stack, struct Number);
    roots_of_unity[2]->operation = 'r';
    roots_of_unity[2]->value.numerator = pool_integer_new(&permanent_integer_pool, 1, -1);
    roots_of_unity[2]->value.denominator = &one;
    roots_of_unity[2]->magnitude_estimate.min = float_one;
    roots_of_unity[2]->magnitude_estimate.max = float_one;
    primes[0] = stack_integer_new(&permanent_stack, 2, 1);
    primes[1] = stack_integer_new(&permanent_stack, 3, 1);
    pi_estimate_min.numerator = pool_integer_new(&permanent_integer_pool, 47, 1);
    pi_estimate_min.denominator = pool_integer_new(&permanent_integer_pool, 15, 1);
    pi_interval_size.numerator = pool_integer_new(&permanent_integer_pool, 1696, 1);
    pi_interval_size.denominator = pool_integer_new(&permanent_integer_pool, 12285, 1);
    pi_sixteen_to_the_k = pool_integer_new(&permanent_integer_pool, 16, 1);
    pi_eight_k = pool_integer_new(&permanent_integer_pool, 8, 1);
}

int main()
{
    struct IntegerPool transient_integer_pool;
    struct NumberPool transient_number_pool;
    struct Stack stack_a;
    struct Stack stack_b;
    init(&transient_integer_pool, &transient_number_pool, &stack_a, &stack_b);
    while (true)
    {
        if (get_input(&transient_number_pool, &transient_integer_pool, &stack_a, &stack_b))
        {
            struct Number*number = (struct Number*)transient_number_pool.stack.start;
            while (number->next)
            {
                if ((number->operation == 'r' && number->next->operation == '(') ||
                    (number->operation == ')' &&
                    (number->next->operation == 'r' || number->next->operation == '(')))
                {
                    struct Number*times = number_slot_new(&transient_number_pool);
                    times->operation = '*';
                    times->previous = number;
                    times->next = number->next;
                    number->next->previous = times;
                    number->next = times;
                }
                number = number->next;
            }
            struct Number*input = parse_input(&transient_number_pool, &transient_integer_pool,
                (struct Number*)transient_number_pool.stack.start);
            if (input)
            {
                struct Number*evaluation = number_evaluate(&transient_number_pool,
                    &transient_integer_pool, &stack_a, &stack_b, input);
                if (evaluation != number_divide_by_zero_error)
                {
                    printf("=\n");
                    print_number(&stack_a, &stack_b, evaluation);
                }
            }
        }
        printf("\n\n");
        stack_free(&transient_number_pool.stack);
        transient_number_pool.free_list = 0;
        stack_free(&stack_a);
        stack_free(&stack_b);
        stack_free(&polynomial_stack);
        memset(transient_integer_pool.sized_pools, 0, page_size);
        stack_free(&transient_integer_pool.slot_page_stack);
    }
    return 0;
}