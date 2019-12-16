#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <ctype.h>
#include <setjmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "os.h"

#ifdef _DEBUG
#define ASSERT(condition, message) if (!(condition)) { puts(message); abort(); }
#else
#define ASSERT(condition, message)
#endif

/*Let a function's result refer to its return value together with the values to which it assigns its
out parameters, if it has any. Let the reference graph of a value of a primitive data type refer to
the value, the reference graph of a pointer refer to the reference graph of the value it points to,
and the reference graph of an aggregate data type, to the union of the reference graphs of its
elements.

If the reference graph of a function's result contains dynamically allocated elements, then the
function takes an output_stack parameter, and allocates all such elements on that Stack. This means
that functions never assign dynamically allocated values passed in by their callers directly into
the reference graphs of their results, instead deep-copying such values to output_stack first.

Functions that need to do dynamic allocations in service of calculating their results but that are
not themselves part of the result occasionally do such allocations on output_stack for efficiency
reasons, but most of the time, they take a local_stack parameter and do the allocations on that. Any
allocations a function does on a local_stack are popped before the function returns. The few
functions that do a significant amount of allocation on output_stack that isn't part of their
results are identified with comments.*/
struct Stack
{
    void*start;
    void*end;
    void*cursor;
    void*cursor_max;
};

struct Integer
{
    size_t value_count;
    int8_t sign;
    uint32_t value[];
};

struct IntegerDivision
{
    struct Integer*quotient;
    struct Integer*remainder;
};

struct Factor
{
    struct Integer*value;
    struct Integer*multiplicity;
};

struct ExtendedGCDInfo
{
    void*gcd;
    void*a_coefficient;
    void*b_coefficient;
    void*a_over_gcd;
    void*b_over_gcd;
};

struct Rational
{
    struct Integer*numerator;
    struct Integer*denominator;
};

struct Float
{
    struct Integer*significand;
    struct Integer*exponent;
};

struct RationalInterval
{
    struct Rational*min;
    struct Rational*max;
};

struct FloatInterval
{
    struct Float*min;
    struct Float*max;
};

struct PolynomialExtendedGCDInfo
{
    void*gcd;
    struct RationalPolynomial*a_coefficient;
};

struct IntegerPolynomial
{
    size_t coefficient_count;
    struct Integer*coefficients[];
};

struct IntegerPolynomialDivision
{
    struct IntegerPolynomial*quotient;
    struct IntegerPolynomial*remainder;
};

struct RationalPolynomial
{
    size_t coefficient_count;
    struct Rational*coefficients[];
};

struct RationalPolynomialDivision
{
    struct RationalPolynomial*quotient;
    struct RationalPolynomial*remainder;
};

struct GaussianRational
{
    struct Rational*real;
    struct Rational*imaginary;
};

struct GaussianRationalPolynomial
{
    size_t coefficient_count;
    struct GaussianRational*coefficients[];
};

struct NestedPolynomial
{
    size_t coefficient_count;
    struct RationalPolynomial*coefficients[];
};

struct NestedPolynomialDivision
{
    struct NestedPolynomial*quotient;
    struct NestedPolynomial*remainder;
};

struct Matrix
{
    struct Rational***rows;
    size_t width;
    size_t height;
};

struct Region
{
    struct FloatInterval*real_interval;
    struct FloatInterval*imaginary_interval;
};

struct UnevaluatedNumber
{
    union
    {
        struct Integer*value;
        struct
        {
            struct UnevaluatedNumber*left;
            struct UnevaluatedNumber*right;
        };
    };
    char operation;
};

struct Token
{
    struct Token*previous;
    struct Token*next;
    struct UnevaluatedNumber*value;
    bool is_parsed;
};

struct Number
{
    char operation;
    struct RationalPolynomial*minimal_polynomial;
    union
    {
        struct Rational*value;//When operation == 'r'.
        struct
        {//When operation == '^'.
            struct Number*radicand;
            struct Integer*index;
        };
        struct
        {//When operation == '*' || operation == '+'.
            size_t element_count;
            struct Number**elements;
        };
    };

    //When operation == '+'.
    struct Number*generator;
    struct RationalPolynomial**terms_in_terms_of_generator;
};

struct TermSplit
{
    struct Rational*rational_part;
    struct Number*irrational_part;
};

typedef struct RationalInterval*(*rational_estimate_getter)(struct Stack*restrict,
    struct Stack*restrict, struct Number*, struct Rational*);

typedef struct FloatInterval*(*float_estimate_getter)(struct Stack*restrict, struct Stack*restrict,
    struct Number*, struct Rational*);

struct EstimateGetters
{
    rational_estimate_getter rational;
    float_estimate_getter fl;
};

__declspec(noreturn) void crash(char*message);
void stack_initialize(struct Stack*out, void*start, size_t size);
void stack_reset(struct Stack*stack);
void*array_start(struct Stack*output_stack, size_t alignment);
void array_extend(struct Stack*output_stack, size_t element_size);
void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment);

#define ALLOCATE(output_stack, type) stack_slot_allocate(output_stack, sizeof(type), _Alignof(type))

#define ARRAY_ALLOCATE(output_stack, element_count, type)\
    stack_slot_allocate(output_stack, (element_count) * sizeof(type), _Alignof(type))

#define POINTER_SWAP(a, b) { void*temp = a; a = b; b = temp; }

struct SmallInteger
{
    size_t value_count;
    int8_t sign;
    uint32_t value;
};

#define INT(value, sign) (struct Integer*)&(struct SmallInteger){ 1, sign, value }

struct Poly1
{
    size_t coefficient_count;
    void*coefficients[1];
};

struct Poly2
{
    size_t coefficient_count;
    void*coefficients[2];
};

#define POLY(coefficient_type, coefficient_count, ...)\
(struct coefficient_type ## Polynomial*)&(struct Poly ## coefficient_count)\
{ coefficient_count, { __VA_ARGS__ }}

struct Integer*integer_allocate(struct Stack*output_stack, size_t value_count);
struct Integer*integer_copy(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_initialize(struct Stack*output_stack, uint32_t value, int8_t sign);
struct Integer*char_to_integer(struct Stack*output_stack, char a);
struct Integer*size_t_to_integer(struct Stack*output_stack, size_t a);
size_t integer_to_size_t(struct Integer*a);
bool integer_equals(struct Integer*a, struct Integer*b);
struct Integer*integer_get_magnitude(struct Stack*output_stack, struct Integer*a);
void integer_add_to_a_in_place(struct Integer*restrict a, struct Integer*restrict b);
struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_negate(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend);
struct Integer*integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
void integer_euclidean_divide(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerDivision*out, struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_get_quotient(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_double(struct Stack*output_stack, struct Integer*a);
void integer_halve_in_place(struct Integer*a);
struct Integer*integer_halve(struct Stack*output_stack, struct Integer*a);
int8_t integer_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Integer*a, struct Integer*b);
struct Integer*integer_get_lcm(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_take_square_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a);
struct Integer*get_next_prime(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b);
size_t size_t_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    size_t**out, size_t a);
size_t integer_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Factor**out, struct Integer*a);
size_t integer_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*a);

struct Rational*rational_copy(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*numerator, struct Integer*denominator);
bool rational_equals(struct Rational*a, struct Rational*b);
struct Rational*rational_get_magnitude(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Rational*a, struct Rational*b);
struct Rational*rational_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b);
struct Rational*rational_negate(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*minuend, struct Rational*subtrahend);
struct Rational*rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational*rational_unreduced_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational*rational_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b);
struct Rational*rational_get_reciprocal(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Rational*divisor);
struct Rational*rational_integer_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Integer*divisor);
struct Rational*rational_double(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a);
struct Rational*rational_halve(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a);
int8_t rational_get_sign(struct Rational*a);
int8_t rational_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Rational*a, struct Rational*b);
struct Rational*rational_get_min(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational*rational_get_max(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b);
struct Rational*rational_get_argument(struct Stack*output_stack, struct Rational*a);
void rational_estimate_cosine(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
void rational_estimate_sine(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
struct RationalInterval*rational_estimate_arctangent(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);
struct RationalInterval*rational_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*y, struct Rational*x,
    struct Rational*interval_size);
struct FloatInterval*rational_get_float_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size);

struct Float*float_copy(struct Stack*output_stack, struct Float*a);
struct Float*float_reduce(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*significand, struct Integer*exponent);
bool float_equals(struct Float*a, struct Float*b);
struct Float*float_get_magnitude(struct Stack*output_stack, struct Float*a);
struct Float*float_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct Float*float_negate(struct Stack*output_stack, struct Float*a);
struct Float*float_subtract(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*minuend, struct Float*subtrahend);
struct Float*float_multiply(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
int8_t float_get_sign(struct Float*a);
int8_t float_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Float*a, struct Float*b);
struct Float*float_get_min(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
struct Float*float_get_max(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b);
void float_estimate_root(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float**out_min, struct Float**out_max, struct Float*a, struct Rational*interval_size,
    struct Integer*index);
struct Rational*float_to_rational(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float*a);

struct RationalInterval*rational_interval_copy(struct Stack*output_stack,
    struct RationalInterval*a);
void rational_interval_expand_bounds(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalInterval*a,
    struct Rational*bound_candidate);
struct Rational*rational_interval_get_factor_interval_size(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*factor_a,
    struct RationalInterval*factor_b, struct Rational*product_interval_size);
struct RationalInterval*rational_interval_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size);
struct FloatInterval*rational_interval_to_float_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*a,
    struct Rational*bound_interval_size);
struct FloatInterval*float_interval_copy(struct Stack*output_stack, struct FloatInterval*a);
bool float_intervals_are_disjoint(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct FloatInterval*a, struct FloatInterval*b);
struct FloatInterval*float_interval_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*minuend,
    struct FloatInterval*subtrahend);
bool region_a_contains_b(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b);
struct RationalInterval*float_interval_to_rational_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*a);
bool regions_are_disjoint(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b);

void pi_estimate(struct Rational*interval_size);
void pi_shrink_interval_to_one_side_of_value(struct Rational*value);

void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count);
struct IntegerPolynomial*integer_polynomial_get_quotient(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
void integer_polynomial_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct PolynomialExtendedGCDInfo*out,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
size_t primitive_integer_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial**out);
struct RationalPolynomial*integer_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
struct RationalPolynomial*integer_polynomial_to_rational_polynomial(struct Stack*output_stack,
    struct IntegerPolynomial*a);

struct Integer*integer_get_residue(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic);
struct Integer*modded_integer_negate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend,
    struct Integer*characteristic);
struct Integer*modded_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic);
struct Integer*modded_integer_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_get_quotient(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
size_t squarefree_modded_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic,
    struct IntegerPolynomial**out);

struct RationalPolynomial*rational_polynomial_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Integer*b);
struct RationalPolynomial*rational_polynomial_get_quotient(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
void rational_polynomial_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct PolynomialExtendedGCDInfo*out,
    struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
size_t rational_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial**out);
size_t rational_polynomial_count_roots_in_rectangle(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a, struct RationalInterval*real,
    struct RationalInterval*imaginary);
void rational_polynomial_evaluate_at_region(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Region*out, struct RationalPolynomial*a,
    struct Region*argument, struct Rational*interval_size);
struct IntegerPolynomial*rational_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a);
struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a);

struct GaussianRational*gaussian_rational_copy(struct Stack*output_stack,
    struct GaussianRational*a);
bool gaussian_rational_equals(struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_negate(struct Stack*output_stack,
    struct GaussianRational*a);
struct GaussianRational*gaussian_rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*minuend,
    struct GaussianRational*subtrahend);
struct GaussianRational*gaussian_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct Rational*b);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct Rational*b);

struct NestedPolynomial*nested_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor);
struct RationalPolynomial*nested_polynomial_get_resultant(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct RationalPolynomial*number_field_element_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial);
size_t number_field_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial**out, struct NestedPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);

void matrix_make_row_echelon_form(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Matrix*a, struct Rational**augmentation);
void matrix_diagonalize(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Matrix*a, struct Rational**augmentation);

struct Number*number_copy(struct Stack*output_stack, struct Number*a);
struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value);
struct Number*number_surd_initialize(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*radicand, struct Integer*index);
int8_t number_formal_compare(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct Number*a, struct Number*b);
struct RationalPolynomial*number_annulling_polynomial_to_minimal_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*annulling_polynomial, struct Number*a);
struct Number*number_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);
struct Number*number_divide(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*dividend, struct Number*divisor);
struct Rational*number_get_rational_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);
struct Number*number_integer_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*exponent);
struct Number*number_rational_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Rational*exponent);
struct Number*number_take_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*index);
struct Number**number_get_conjugates(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);
size_t number_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a);

struct RationalPolynomial*sum_get_minimal_polynomial(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a,
    struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a, struct Number*b);
struct Number*number_eliminate_linear_dependencies(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);

struct RationalPolynomial*product_get_minimal_polynomial(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a,
    struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*b);
struct Number*number_multiply(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a, struct Number*b);

struct FloatInterval*number_get_float_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_get_rational_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct FloatInterval*number_get_float_imaginary_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_get_rational_imaginary_part_estimate(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct Rational*interval_size);
void number_get_region_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Region*out, struct Number*a,
    struct Rational*interval_size);
void number_halve_region_estimate_dimensions(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Region*out, struct Number*a);
struct FloatInterval*number_get_float_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_get_rational_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_estimate_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_estimate_argument_in_radians(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);

struct Number**get_roots_of_unity(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct Integer*degree);

struct Stack permanent_stack;
struct Stack pi_stack_a;
struct Stack pi_stack_b;
struct Stack prime_stack;
struct Integer**next_prime;
jmp_buf memory_error_buffer;

struct RationalInterval pi;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;

struct Number**roots_of_unity;

struct Integer zero = { 0, 0, 0 };
struct Integer one = { 1, 1, 1 };
struct Rational rational_zero = { &zero, &one };
struct Rational rational_one = { &one, &one };
struct Float float_zero = { &zero, &zero };
struct Float float_one = { &one, &zero };
void*polynomial_zero = &(size_t) { 0 };
struct IntegerPolynomial integer_polynomial_one = { 1, &one };
struct RationalPolynomial rational_polynomial_one = { 1, &rational_one };
struct GaussianRational gaussian_rational_zero = { &rational_zero, &rational_zero };
struct GaussianRational gaussian_rational_one = { &rational_one, &rational_zero };
struct GaussianRationalPolynomial gaussian_rational_polynomial_one = { 1, &gaussian_rational_one };
struct NestedPolynomial nested_polynomial_one = { 1, &rational_polynomial_one };
struct Number number_one;

struct EstimateGetters real_estimate_getters =
{ number_get_rational_real_part_estimate, number_get_float_real_part_estimate };
struct EstimateGetters imaginary_estimate_getters =
{ number_get_rational_imaginary_part_estimate, number_get_float_imaginary_part_estimate };
struct EstimateGetters magnitude_estimate_getters =
{ number_get_rational_magnitude_estimate, number_get_float_magnitude_estimate };

#endif
