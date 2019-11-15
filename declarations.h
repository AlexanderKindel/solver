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

//With the exception of functions whose names begin with "leaking_", every function that dynamically
//allocates memory that the caller is interested in has an output_stack parameter, and leaves the
//memory on that Stack. All other Stack parameters - local_stack for functions that also have an
//output_stack, or stack_a and stack_b for those that don't - are for internal use only, and have
//any memory the function allocates on them freed again before the function returns. It's common for
//a function implementation to not know whether a value it allocates is an output value until after
//the fact, for example, in a while loop where the end condition depends on the value. In such
//cases, the value is initially allocated on local_stack, then copied to output_stack. leaking_
//functions exist to handle cases where, among the values allocated by the function body, different
//callers are interested in different ones, and getting any of them onto output_stack would require
//copying them. These functions do not reset local_stack, deferring the question of which values to
//copy to the caller instead of doing it themselves.
struct Stack
{
    size_t start;
    size_t end;
    void*cursor;
    size_t cursor_max;
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

struct RectangularEstimate
{
    struct FloatInterval*real_part_estimate;
    struct FloatInterval*imaginary_part_estimate;
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

typedef struct RationalInterval*(*rational_estimate_getter)(struct Stack*, struct Stack*,
    struct Number*, struct Rational*);

typedef struct FloatInterval*(*float_estimate_getter)(struct Stack*, struct Stack*, struct Number*,
    struct Rational*);

struct EstimateGetters
{
    rational_estimate_getter rational;
    float_estimate_getter fl;
};

__declspec(noreturn) void crash(char*message);
void stack_initialize(struct Stack*out, size_t start, size_t size);
void stack_reset(struct Stack*out);
void*array_start(struct Stack*output_stack, size_t alignment);
void extend_array(struct Stack*output_stack, size_t element_size);
void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment);

#define ALLOCATE(stack, type) stack_slot_allocate(stack, sizeof(type), _Alignof(type))

#define ARRAY_ALLOCATE(stack, element_count, type)\
    stack_slot_allocate(stack, (element_count) * sizeof(type), _Alignof(type))

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

void set_page_size();

struct Integer*integer_allocate(struct Stack*output_stack, size_t value_count);
struct Integer*integer_copy(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_initialize(struct Stack*stack, uint32_t value, int8_t sign);
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
struct Integer*integer_lcm(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a,
    struct Integer*b);
struct Integer*integer_square_root(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a);
struct Integer*get_prime(struct Stack*stack_a, struct Stack*stack_b, size_t index);
size_t size_t_factor(struct Stack*output_stack, struct Stack*local_stack, size_t**out, size_t a);
size_t integer_factor(struct Stack*output_stack, struct Stack*local_stack, struct Factor**out,
    struct Integer*a);
size_t integer_string(struct Stack*output_stack, struct Stack*local_stack, struct Integer*a);

struct Rational*rational_copy(struct Stack*output_stack, struct Rational*a);
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
    struct Rational*a, struct Rational*b);
struct Rational*rational_integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b);
struct Rational*rational_reciprocal(struct Stack*output_stack, struct Rational*a);
struct Rational*rational_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Rational*divisor);
struct Rational*rational_integer_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Integer*divisor);
struct Rational*rational_doubled(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a);
struct Rational*rational_half(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a);
int8_t rational_sign(struct Rational*a);
int8_t rational_compare(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b);
struct Rational*rational_min(struct Stack*output_stack, struct Stack*local_stack, struct Rational*a,
    struct Rational*b);
struct Rational*rational_max(struct Stack*output_stack, struct Stack*local_stack, struct Rational*a,
    struct Rational*b);
struct Rational*rational_argument(struct Stack*output_stack, struct Rational*a);
void rational_estimate_cosine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
void rational_estimate_sine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size);
struct RationalInterval*rational_estimate_arctangent(struct Stack*output_stack,
    struct Stack*local_stack, struct Rational*a, struct Rational*interval_size);
struct RationalInterval*rational_estimate_atan2(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*y, struct Rational*x, struct Rational*interval_size);
struct FloatInterval*rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*interval_size);

struct Float*float_copy(struct Stack*output_stack, struct Float*a);
struct Float*float_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*significand, struct Integer*exponent);
bool float_equals(struct Float*a, struct Float*b);
struct Float*float_magnitude(struct Stack*output_stack, struct Float*a);
struct Float*float_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
struct Float*float_negative(struct Stack*output_stack, struct Float*a);
struct Float*float_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*minuend, struct Float*subtrahend);
struct Float*float_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
int8_t float_sign(struct Float*a);
int8_t float_compare(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b);
struct Float*float_min(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
struct Float*float_max(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b);
void float_estimate_root(struct Stack*output_stack, struct Stack*local_stack, struct Float**out_min,
    struct Float**out_max, struct Float*a, struct Rational*interval_size, struct Integer*index);
struct Rational*float_to_rational(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a);

struct RationalInterval*rational_interval_copy(struct Stack*output_stack,
    struct RationalInterval*a);
void rational_interval_expand_bounds(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalInterval*a, struct Rational*bound_candidate);
struct Rational*rational_interval_factor_interval_size(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*factor_a, struct RationalInterval*factor_b,
    struct Rational*product_interval_size);
struct RationalInterval*rational_interval_estimate_atan2(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size);
struct FloatInterval*rational_interval_to_float_interval(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct Rational*bound_interval_size);
struct FloatInterval*float_interval_copy(struct Stack*output_stack, struct FloatInterval*a);
bool float_intervals_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*a, struct FloatInterval*b);
struct FloatInterval*float_interval_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b);
bool rectangular_estimate_a_contains_b(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b);
struct RationalInterval*float_interval_to_rational_interval(struct Stack*output_stack,
    struct Stack*local_stack, struct FloatInterval*a);
bool rectangular_estimates_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b);

void pi_estimate(struct Rational*interval_size);
void pi_shrink_interval_to_one_side_of_value(struct Rational*value);

size_t polynomial_size(size_t coefficient_count);
void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count);
struct IntegerPolynomial*integer_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
void integer_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
size_t primitive_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out);
struct RationalPolynomial*integer_polynomial_to_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a);
struct RationalPolynomial*integer_polynomial_to_rational_polynomial(struct Stack*output_stack,
    struct IntegerPolynomial*a);

struct Integer*integer_residue(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b, struct Integer*characteristic);
struct Integer*modded_integer_negative(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic);
struct Integer*modded_integer_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*minuend, struct Integer*subtrahend, struct Integer*characteristic);
struct Integer*modded_integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*b, struct Integer*characteristic);
struct Integer*modded_integer_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_reduced(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*minuend, struct IntegerPolynomial*subtrahend,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
size_t squarefree_modded_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, struct IntegerPolynomial**out);

struct RationalPolynomial*rational_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Integer*b);
struct RationalPolynomial*rational_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor);
void rational_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
size_t rational_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial**out);
size_t rational_polynomial_root_count_in_rectangle(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalPolynomial*a, struct RationalInterval*real, struct RationalInterval*imaginary);
void rational_polynomial_evaluate_at_rectangular_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct RationalPolynomial*a,
    struct RectangularEstimate*argument, struct Rational*interval_size);
struct IntegerPolynomial*rational_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a);
struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a);

struct GaussianRational*gaussian_rational_copy(struct Stack*output_stack,
    struct GaussianRational*a);
bool gaussian_rational_equals(struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_add(struct Stack*output_stack, struct Stack*local_stack,
    struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_negative(struct Stack*output_stack,
    struct GaussianRational*a);
struct GaussianRational*gaussian_rational_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRational*a, struct GaussianRational*b);
struct GaussianRational*gaussian_rational_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRational*a, struct Rational*b);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*output_stack, struct Stack*local_stack, struct GaussianRationalPolynomial*a,
    struct Rational*b);

struct NestedPolynomial*nested_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend, struct NestedPolynomial*divisor);
struct RationalPolynomial*nested_polynomial_resultant(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct RationalPolynomial*number_field_element_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_reciprocal(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);
struct RationalPolynomial*number_field_element_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend, struct RationalPolynomial*divisor,
    struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomial*number_field_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial);
size_t number_field_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*generator_minimal_polynomial,
    struct NestedPolynomial**out);

void matrix_row_echelon_form(struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a,
    struct Rational**augmentation);
void matrix_diagonalize(struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a,
    struct Rational**augmentation);

struct Number*number_copy(struct Stack*output_stack, struct Number*a);
struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value);
struct Number*number_surd_initialize(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*radicand, struct Integer*index);
int8_t number_formal_compare(struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b);
struct RationalPolynomial*number_minimal_polynomial_from_annulling_polynomial(
    struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*annulling_polynomial, struct Number*a);
struct RationalPolynomial*sum_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a);
struct Number*number_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*dividend, struct Number*divisor);
struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a);
struct Number*number_integer_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Integer*exponent);
struct Number*number_rational_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Rational*exponent);
struct Number*number_root(struct Stack*output_stack, struct Stack*local_stack, struct Number*base,
    struct Integer*index);
struct Number**number_conjugates(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a);
size_t number_string(struct Stack*output_stack, struct Stack*local_stack, struct Number*a);

struct RationalPolynomial*sum_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_add(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*b);
struct Number*number_eliminate_linear_dependencies(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a);

struct RationalPolynomial*product_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_rational_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Rational*b);
struct Number*number_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*b);

struct FloatInterval*number_float_real_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_rational_real_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct FloatInterval*number_float_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_rational_imaginary_part_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
void number_rectangular_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RectangularEstimate*out, struct Number*a, struct Rational*interval_size);
void number_rectangular_estimate_halve_dimensions(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct Number*a);
struct FloatInterval*number_float_magnitude_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_rational_magnitude_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_argument_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval*number_rational_argument_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*interval_size);

struct Number**get_roots_of_unity(struct Stack*stack_a, struct Stack*stack_b,
    struct Integer*degree);

struct Stack permanent_stack;
struct Stack pi_stack_a;
struct Stack pi_stack_b;
jmp_buf memory_error_buffer;

struct Integer**primes;

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
{ number_rational_real_part_estimate, number_float_real_part_estimate };
struct EstimateGetters imaginary_estimate_getters =
{ number_rational_imaginary_part_estimate, number_float_imaginary_part_estimate };
struct EstimateGetters magnitude_estimate_getters =
{ number_rational_magnitude_estimate, number_float_magnitude_estimate };

#endif
