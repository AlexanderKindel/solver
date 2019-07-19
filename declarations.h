#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#ifdef _DEBUG
#define ASSERT(condition, message) if (!(condition)) { puts(message); abort(); }
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

struct Pool
{
    void*first_page;
    void*cursor;
    void*free_list;
};

struct PoolSet
{
    struct Pool*pool_page;
    struct Stack slot_page_stack;
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

struct RingOperations
{
    //Functions that all types that fill out this struct must implement.
    void*(*copy)(struct Stack*output_stack, void*a);
    bool(*equals)(void*a, void*b);
    void*additive_identity;
    void*multiplicative_identity;
    void*(*add)(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b, void*misc);
    void*(*negative)(struct Stack*output_stack, struct Stack*local_stack, void*a, void*misc);
    void*(*multiply)(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
        void*misc);
};

struct EuclideanDomainOperations
{
    //Functions that all types that fill out this struct must implement.
    struct RingOperations ring_operations;
    void(*euclidean_divide)(struct Stack*output_stack, struct Stack*local_stack,
        struct Division*out, void*dividend, void*divisor, void*misc);

    //Function that can be implemented for any type that implements the above, but which has
    //different optimized implementations for different types. May be left 0 if unneeded.
    void*(*gcd)(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b, void*misc);
};

struct FieldOperations
{
    struct RingOperations ring_operations;
    void*(*reciprocal)(struct Stack*, struct Stack*, void*, void*);
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

struct PolynomialExtendedGCDInfo
{
    struct Polynomial*gcd;
    struct Polynomial*a_coefficient;
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

struct NestedPolynomial
{
    size_t coefficient_count;
    struct RationalPolynomial*coefficients[];
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
    struct Float*min;
    struct Float*max;
};

struct SumFieldForm
{
    struct Number*generator;
    struct RationalPolynomial*left_term_in_terms_of_generator;
    struct RationalPolynomial*right_term_in_terms_of_generator;
};

struct Number
{
    //Points to the next element of the doubly linked list representing the input during parsing.
    //When on a free list during evaluation, points to the next element of the free list. While not
    //on a free list during evaluation, points to the next element in a singly linked list
    //consisting of the conjugates of the first element. If a Number is not part of a list of
    //conjugates and doesn't have operation == 'r', next == 0 means the conjugates have not yet been
    //calculated.
    struct Number*next;
    union
    {
        struct Number*previous;//During parsing, when the input is stored as a doubly linked list.
        struct RationalPolynomial*minimal_polynomial;//During evaluation.
    };
    struct SumFieldForm*field_form;
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
            struct FloatInterval*real_part_estimate;
            struct FloatInterval*imaginary_part_estimate;
        };
    };
    struct FloatInterval*argument_estimate;
    struct FloatInterval*magnitude_estimate;
    char operation;
};

typedef void(*rational_estimate_getter)(struct PoolSet*, struct Stack*, struct Stack*,
    struct RationalInterval*, struct Number*, struct Rational*);

typedef struct FloatInterval*(*float_estimate_getter)(struct PoolSet*, struct Stack*, struct Stack*,
    struct Number*, struct Rational*);

size_t next_page_boundary(size_t address);
size_t next_aligned_address(size_t address, size_t alignment);
void stack_initialize(struct Stack*out, size_t start, size_t end);
void stack_free(struct Stack*out);
void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment);
void pool_set_initialize(struct PoolSet*out, size_t start, size_t end);
void*pool_slot_page_allocate(struct PoolSet*pool_set);
struct Pool*get_pool(struct PoolSet*pool_set, size_t slot_size);
void*pool_slot_allocate(struct PoolSet*pool_set, size_t slot_size);
void pool_slot_free(struct PoolSet*pool_set, void*a, size_t slot_size);

#define STACK_SLOT_ALLOCATE(stack, type) stack_slot_allocate(stack, sizeof(type), _Alignof(type))

#define POOL_VALUE_ALLOCATE(pool_set, type)\
    (void*)((size_t)pool_slot_allocate(pool_set, sizeof(void*) + sizeof(type)) + sizeof(void*))

#define POOL_VALUE_FREE(pool_set, a, type)\
    pool_slot_free(pool_set, (void*)((size_t)a - sizeof(void*)), sizeof(type) + sizeof(void*))

__declspec(noreturn) void crash(char*message);
void array_reverse(void**a, size_t element_count);
void*generic_exponentiate(struct RingOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, void*base, struct Integer*exponent, void*misc);
void*generic_gcd(struct EuclideanDomainOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, void*a, void*b, void*misc);
void generic_extended_gcd(struct EuclideanDomainOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct ExtendedGCDInfo*out, void*a, void*b, void*misc);
void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count);
void polynomial_copy_coefficients(void*(coefficient_copy)(struct Stack*, void*),
    struct Stack*output_stack, struct Polynomial*a);
void*polynomial_copy(void*(coefficient_copy)(struct Stack*, void*), struct Stack*output_stack,
    struct Polynomial*a);
void polynomial_trim_leading_zeroes(struct RingOperations*coefficient_operations,
    struct Polynomial*a);
bool polynomial_equals(struct RingOperations*coefficient_operations, struct Polynomial*a,
    struct Polynomial*b);
void*polynomial_add(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b);
void*polynomial_negative(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Polynomial*a);
void*polynomial_subtract(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*minuend, struct Polynomial*subtrahend);
void*polynomial_multiply(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b, void*misc);
void*polynomial_multiply_by_coefficient(struct RingOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a, void*b, void*misc);
void polynomial_euclidean_divide(struct EuclideanDomainOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct PolynomialDivision*out,
    struct Polynomial*dividend, struct Polynomial*divisor, void*misc);
void field_polynomial_euclidean_divide(struct FieldOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct PolynomialDivision*out,
    struct Polynomial*dividend, struct Polynomial*divisor, void*misc);
void*polynomial_divide_by_coefficient(void*(coefficient_quotient)(struct Stack*, struct Stack*,
    void*, void*), struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*dividend,
    void*divisor);
void*polynomial_content(struct EuclideanDomainOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a, void*misc);
size_t polynomial_squarefree_factor(struct EuclideanDomainOperations*polynomial_operations,
    void*(derivative)(struct Stack*, struct Stack*, void*), struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial**out, void*misc);
void*polynomial_derivative(void*(coefficient_times_integer)(struct Stack*, struct Stack*, void*,
    struct Integer*), struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a);
void*integer_generic_add(struct Stack*output_stack, struct Stack*unused_stack, void*a, void*b,
    void*unused);
void*integer_generic_negative(struct Stack*output_stack, struct Stack*unused_stack, void*a,
    void*unused);
void*integer_generic_multiply(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused);
void integer_generic_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Division*out, void*dividend, void*divisor, void*unused);
void*integer_generic_gcd(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused);
void*integer_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused);
void*integer_polynomial_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    void*a, void*unused);
void*integer_polynomial_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    void*a, void*b, void*unused);
void integer_polynomial_generic_euclidean_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct Division*out, void*dividend, void*divisor, void*unused);
void*integer_polynomial_generic_gcd(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused);
void*rational_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused);
void*rational_generic_negative(struct Stack*output_stack, struct Stack*unused_stack, void*a,
    void*unused);
void*rational_generic_multiply(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused);
void*rational_generic_reciprocal(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*unused);
void*rational_generic_divide(struct Stack*output_stack, struct Stack*local_stack, void*dividend,
    void*divisor, void*unused);
void*rational_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused);
void*rational_polynomial_generic_negative(struct Stack*output_stack, struct Stack*local_stack,
    void*a, void*unused);
void*rational_polynomial_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    void*a, void*b, void*unused);
void rational_polynomial_generic_euclidean_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct Division*out, void*dividend, void*divisor, void*unused);
void*nested_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, struct NestedPolynomial*b, void*unused);
void*nested_polynomial_generic_negative(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, void*unused);

#define POINTER_SWAP(a, b) { void*temp = a; a = b; b = temp; }

struct Integer*stack_integer_allocate(struct Stack*output_stack, size_t value_count);
struct Integer*pool_integer_allocate(struct PoolSet*pool_set, size_t value_count);
void pool_integer_free(struct PoolSet*pool_set, struct Integer*a);
struct Integer*integer_copy_to_stack(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_copy_to_pool(struct PoolSet*pool_set, struct Integer*a);
void integer_move_to_pool(struct PoolSet*pool_set, struct Integer**a);
void integer_move_from_pool(struct PoolSet*pool_set, struct Stack*output_stack, struct Integer**a);
struct Integer*stack_integer_initialize(struct Stack*stack, uint32_t value, int8_t sign);
struct Integer*pool_integer_initialize(struct PoolSet*pool_set, uint32_t value, int8_t sign);
struct Integer*integer_from_char(struct Stack*output_stack, char value);
struct Integer*integer_from_size_t(struct Stack*output_stack, size_t value);
size_t integer_to_size_t(struct Integer*a);
bool integer_equals(struct Integer*a, struct Integer*b);
struct Integer*integer_magnitude(void*output_arena, struct Integer*a);
void integer_add_to_a_in_place(struct Integer*a, struct Integer*b);
struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_negative(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*minuend, struct Integer*subtrahend);
struct Integer*integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b);
void integer_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerDivision*out, struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_euclidean_quotient(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_euclidean_remainder(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_doubled(struct Stack*output_stack, struct Integer*a);
void integer_halve(struct Integer*a);
struct Integer*integer_half(struct Stack*output_stack, struct Integer*a);
int8_t integer_compare(struct Stack*stack_a, struct Stack*stack_b, struct Integer*a,
    struct Integer*b);
struct Integer*n_choose_k(struct Stack*output_stack, struct Stack*local_stack, struct Integer*n,
    struct Integer*k);
struct Integer*integer_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*base, struct Integer*exponent);
struct Integer*integer_gcd(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a,
    struct Integer*b);
void integer_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct ExtendedGCDInfo*out, struct Integer*a, struct Integer*b);
struct Integer*integer_lcm(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a,
    struct Integer*b);
struct Integer*integer_square_root(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a);
struct Integer*get_prime(struct Stack*stack_a, struct Stack*stack_b, size_t index);
size_t size_t_factor(struct Stack*output_stack, struct Stack*local_stack, size_t**out, size_t a);
size_t integer_factor(struct Stack*output_stack, struct Stack*local_stack, struct Factor**out,
    struct Integer*a);
void integer_string(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a);

#define INT(value, sign) (struct Integer){ 1, sign 1, value }

void rational_free(struct PoolSet*pool_set, struct Rational*a);
struct Rational*rational_copy_to_stack(struct Stack*output_stack, struct Rational*a);
void rational_move_to_pool(struct PoolSet*pool_set, struct Rational*a);
void rational_move_from_pool(struct PoolSet*pool_set, struct Stack*output_stack, struct Rational*a);
struct Rational*rational_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*numerator, struct Integer*denominator);
bool rational_equals(struct Rational*a, struct Rational*b);
struct Rational*rational_magnitude(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b);
struct Rational*rational_integer_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b);
struct Rational*rational_negative(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*minuend, struct Rational*subtrahend);
struct Rational*rational_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b);
struct Rational*rational_unreduced_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b, void*unused);
struct Rational*rational_integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b);
struct Rational*rational_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a);
struct Rational*rational_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Rational*divisor);
struct Rational*rational_integer_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Integer*divisor);
struct Rational*rational_doubled(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a);
struct Rational*rational_half(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a);
int8_t rational_compare(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b);
struct Rational*rational_min(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b);
struct Rational*rational_max(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b);
struct Rational*rational_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*base, struct Integer*exponent);
void rational_estimate_cosine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
void rational_estimate_sine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
void rational_estimate_arctangent(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
void estimate_atan2(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*y, struct Rational*x,
    struct Rational*interval_size);
void rational_continue_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Integer**estimate_denominator,
    struct Integer**estimate_remainder, struct Rational*a, struct Rational*interval_size);
void rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Rational*a, struct Rational*interval_size);
void rational_interval_expand_bounds(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalInterval*a, struct Rational*bound_candidate);
void rational_interval_to_float_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct RationalInterval*a,
    struct Rational*bound_interval_size);
void pi_estimate(struct Stack*stack_a, struct Stack*stack_b, struct Rational*interval_size);
void pi_shrink_interval_to_one_side_of_value(struct Stack*stack_a, struct Stack*stack_b,
    struct Rational*value);

struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a);
bool integer_polynomial_equals(struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_negative(struct Stack*output_stack,
    struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend);
struct IntegerPolynomial*integer_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b);
void integer_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_integer_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct Integer*divisor);
struct Integer*integer_polynomial_content(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
void integer_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
size_t primitive_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out);
struct IntegerPolynomial*integer_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a);
struct RationalPolynomial*integer_polynomial_to_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a);

struct Integer*modded_integer_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_reduced(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_multiply_by_coefficient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b,
    struct Integer*characteristic);
void modded_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*base, struct Integer*exponent,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic);
void modded_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct ExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
size_t squarefree_modded_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, struct IntegerPolynomial**out);

struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack,
    struct RationalPolynomial*a);
struct RationalPolynomial*rational_polynomial_copy_to_number_memory(
    struct PoolSet*pool_set, struct RationalPolynomial*a);
bool rational_polynomial_equals(struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Integer*b);
struct RationalPolynomial*rational_polynomial_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Rational*b);
void rational_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*base, struct Integer*exponent);
void rational_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
size_t rational_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial**out);
struct RationalPolynomial*rational_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a);
struct IntegerPolynomial*rational_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a);
struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a);

struct RationalPolynomial*number_field_element_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_reciprocal(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend, struct RationalPolynomial*divisor,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*base, struct Integer*exponent,
    struct RationalPolynomial*generator_minimal_polynomial);

struct NestedPolynomial*nested_polynomial_copy(struct Stack*output_stack,
    struct NestedPolynomial*a);
bool nested_polynomial_equals(struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a);
struct NestedPolynomial*nested_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*minuend, struct NestedPolynomial*subtrahend);
struct NestedPolynomial*nested_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct RationalPolynomial*b);
void nested_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor);
struct NestedPolynomial*nested_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend, struct NestedPolynomial*divisor);
struct NestedPolynomial*nested_polynomial_rational_polynomial_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend, struct RationalPolynomial*divisor);
struct NestedPolynomial*nested_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a);
struct RationalPolynomial*nested_polynomial_get_resultant(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);

struct NestedPolynomial*number_field_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomial*number_field_polynomial_element_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
void number_field_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomial*number_field_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
size_t number_field_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*generator_minimal_polynomial,
    struct NestedPolynomial**out);

void matrix_row_echelon_form(struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a,
    struct RationalPolynomial**augmentation);
struct AlgebraicNumber*pool_algebraic_number_allocate(struct Stack*output_stack,
    struct AlgebraicNumber*free_list);
void algebraic_number_move_coefficients_to_stack(struct Stack*output_stack,
    struct AlgebraicNumber*a);
struct AlgebraicNumber*algebraic_number_add(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b);
struct AlgebraicNumber*algebraic_number_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b,
    struct RationalPolynomial*generator_annulling_polynomials[2]);

void float_free(struct PoolSet*pool_set, struct Float*a);
struct Float*float_copy_to_stack(struct Stack*output_stack, struct Float*a);
void float_copy_to_pool(struct PoolSet*pool_set, struct Float*a);
void float_move_to_pool(struct PoolSet*pool_set, struct Float**a);
void float_move_from_pool(struct PoolSet*pool_set, struct Stack*output_stack, struct Float**a);
struct Float*float_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*significand, struct Integer*exponent);
bool float_equals(struct Float*a, struct Float*b);
struct Float*float_magnitude(struct Stack*output_stack, struct Float*a);
struct Float*float_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
struct Float*float_generic_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b, void*unused);
struct Float*float_negative(struct Stack*output_stack, struct Float*a);
struct Float*float_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Float*a, void*unused);
struct Float*float_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*minuend, struct Float*subtrahend);
struct Float*float_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
struct Float*float_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a, struct Float*b, void*unused);
int8_t float_compare(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b);
struct Float*float_max(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b);
struct Float*float_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*base, struct Integer*exponent);
void float_estimate_root(struct Stack*output_stack, struct Stack*local_stack, struct Float**out_min,
    struct Float**out_max, struct Float*a, struct Rational*interval_size, struct Integer*index);
struct Rational*float_to_rational(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a);
void float_interval_free(struct PoolSet*pool_set, struct FloatInterval*a);
void float_interval_move_to_pool(struct PoolSet*pool_set, struct FloatInterval*a);
void float_interval_move_from_pool(struct PoolSet*pool_set, struct Stack*output_stack,
    struct FloatInterval*a);
void float_interval_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct FloatInterval*a, struct FloatInterval*b);
void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a);

struct Number*number_allocate(struct PoolSet*pool_set);
void number_node_free(struct PoolSet*pool_set, struct Number*a);
void number_node_free_during_parse(struct PoolSet*pool_set, struct Number*a);
void number_free(struct PoolSet*pool_set, struct Number*a);
struct Number*number_copy(struct PoolSet*pool_set, struct Number*a);
struct Number*number_add(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b);
struct Number*number_rational_multiply(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*b);
struct Number*number_multiply(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b);
struct Number*number_reciprocal(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a);
struct Number*number_divide(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*dividend, struct Number*divisor);
struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a);
struct Number*number_exponentiate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*base, struct Rational*exponent);
void number_calculate_minimal_polynomial(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a);
void number_calculate_conjugates(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a);

struct FloatInterval*number_float_real_part_estimate(struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size);
void number_rational_real_part_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size);
struct FloatInterval*number_float_imaginary_part_estimate(struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size);
void number_rational_imaginary_part_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size);
struct FloatInterval*number_float_magnitude_estimate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*interval_size);
void number_rational_magnitude_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size);
struct FloatInterval*number_float_argument_estimate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*interval_size);
void number_rational_argument_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size);
void rational_polynomial_estimate_evaluation(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct Float**real_estimate_min_out,
    struct Float**real_estimate_max_out, struct Float**imaginary_estimate_min_out,
    struct Float**imaginary_estimate_max_out, struct RationalPolynomial*a, struct Number*argument,
    struct Rational*interval_size);

struct Number**get_roots_of_unity(struct Stack*stack_a, struct Stack*stack_b,
    struct Integer*degree);

//Not equal to the OS page size, but instead to whichever of the OS page size and OS allocation
//granularity is larger, so that neither the reserve nor the commit operation causes the requested
//memory address to be rounded when calling VirtualAlloc with MEM_RESERVE | MEM_COMMIT.
size_t page_size;
size_t pool_memory_start;
size_t pool_memory_end;
struct Stack permanent_stack;
struct Stack polynomial_stack;
struct PoolSet permanent_pool_set;

struct Integer zero = { 0, 0, 0 };
struct Integer one = { 1, 1, 1 };
struct EuclideanDomainOperations integer_operations = { { integer_copy_to_stack,
    integer_equals, &zero, &one, integer_generic_add, integer_generic_negative,
    integer_generic_multiply }, integer_generic_euclidean_divide, integer_generic_gcd };
struct Integer**primes;

struct Rational rational_zero = { &zero, &one };
struct Rational rational_one = { &one, &one };
struct FieldOperations rational_operations = { { rational_copy_to_stack, rational_equals,
    &rational_zero, &rational_one, rational_generic_add, rational_generic_negative,
    rational_generic_multiply }, rational_generic_reciprocal };
struct Rational pi_estimate_min;
struct Rational pi_interval_size;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

struct Float float_zero = { &zero, &zero };
struct Float float_one = { &one, &zero };
struct RingOperations float_operations = { float_copy_to_stack, float_equals, &float_zero,
    &float_one, float_generic_add, float_generic_negative, float_generic_multiply };

struct Polynomial polynomial_zero = { 0 };

struct IntegerPolynomial integer_polynomial_one = { 1, &one };
struct EuclideanDomainOperations integer_polynomial_operations = { { integer_polynomial_copy,
    integer_polynomial_equals, &polynomial_zero, &integer_polynomial_one,
    integer_polynomial_generic_add, integer_polynomial_generic_negative,
    integer_polynomial_generic_multiply }, integer_polynomial_generic_euclidean_divide,
    integer_polynomial_generic_gcd };

struct FieldOperations modded_integer_operations = { { integer_copy_to_stack, integer_equals,
    &zero, &one, integer_generic_add, integer_generic_negative, integer_generic_multiply },
    modded_integer_reciprocal };
struct EuclideanDomainOperations modded_polynomial_operations = { { integer_polynomial_copy,
    integer_polynomial_equals, &polynomial_zero, &integer_polynomial_one, modded_polynomial_add,
    modded_polynomial_negative, modded_polynomial_multiply },
    modded_polynomial_euclidean_divide, 0 };

struct RationalPolynomial rational_polynomial_one = { 1, &rational_one };
struct EuclideanDomainOperations rational_polynomial_operations = { { rational_polynomial_copy,
    rational_polynomial_equals, &polynomial_zero, &rational_polynomial_one,
    rational_polynomial_generic_add, rational_polynomial_generic_negative,
    rational_polynomial_generic_multiply }, rational_polynomial_generic_euclidean_divide, 0 };

struct FieldOperations number_field_element_operations = { { rational_polynomial_copy,
    rational_polynomial_equals, &polynomial_zero, &rational_polynomial_one,
    rational_polynomial_generic_add, rational_polynomial_generic_negative,
    number_field_element_multiply }, number_field_element_reciprocal };

struct NestedPolynomial nested_polynomial_one = { 1, &rational_polynomial_one };

struct EuclideanDomainOperations number_field_polynomial_operations =
    { { nested_polynomial_copy, nested_polynomial_equals, &polynomial_zero,
    &nested_polynomial_one, nested_polynomial_generic_add, nested_polynomial_generic_negative,
    number_field_polynomial_multiply }, number_field_polynomial_euclidean_divide,
    number_field_polynomial_gcd };

struct Number*number_divide_by_zero_error = 0;
struct Number*number_product_consolidation_failed = (struct Number*)1;

struct Number**roots_of_unity;

#endif