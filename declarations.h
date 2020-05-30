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

struct Stack
{
    void*start;
    void*end;
    void*cursor;
    void*cursor_max;
};

jmp_buf memory_error_buffer;
struct Stack permanent_stack;

void*polynomial_zero = &(size_t) { 0 };

struct ExtendedGCDInfo
{
    void*gcd;
    void*a_coefficient;
    void*b_coefficient;
    void*a_over_gcd;
    void*b_over_gcd;
};

struct Integer
{
    size_t value_count;
    int8_t sign;
    uint32_t value[];
};

struct Integer zero = { 0, 0, 0 };
struct Integer one = { 1, 1, 1 };
struct Stack prime_stack;
struct Integer**next_prime;

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

struct Float
{
    struct Integer*significand;
    struct Integer*exponent;
};

struct Float float_zero = { &zero, &zero };
struct Float float_one = { &one, &zero };

struct FloatInterval
{
    struct Float min;
    struct Float max;
};

struct Region
{
    struct FloatInterval real_interval;
    struct FloatInterval imaginary_interval;
};

struct IntegerPolynomial
{
    size_t coefficient_count;
    struct Integer*coefficients[];
};

struct IntegerPolynomial*integer_polynomial_one = &(struct IntegerPolynomial) { 1, &one };

struct IntegerPolynomialDivision
{
    struct IntegerPolynomial*quotient;
    struct IntegerPolynomial*remainder;
};

struct Rational
{
    struct Integer*numerator;
    struct Integer*denominator;
};

struct Rational rational_zero = { &zero, &one };
struct Rational rational_one = { &one, &one };

struct GaussianRational
{
    struct Rational real;
    struct Rational imaginary;
};

struct GaussianRational gaussian_rational_zero = { { &zero, &one }, { &zero, &one } };
struct GaussianRational gaussian_rational_one = { { &one, &one }, { &zero, &one } };

struct GaussianRationalPolynomial
{
    size_t coefficient_count;
    struct GaussianRational coefficients[];
};

struct GaussianRationalPolynomial gaussian_rational_polynomial_one_value =
{ 1, { { { &one, &one }, { &zero, &one } } } };
struct GaussianRationalPolynomial*gaussian_rational_polynomial_one =
    &gaussian_rational_polynomial_one_value;

struct Matrix
{
    struct Rational**rows;
    size_t width;
    size_t height;
};

struct RationalInterval
{
    struct Rational min;
    struct Rational max;
};

struct Stack pi_stack_a;
struct Stack pi_stack_b;
struct Integer*pi_sixteen_to_the_k;
struct Integer*pi_eight_k;
struct RationalInterval pi;

struct RationalPolynomial
{
    size_t coefficient_count;
    struct Rational coefficients[];
};

struct PolynomialExtendedGCDInfo
{
    void*gcd;
    struct RationalPolynomial*a_coefficient;
};

struct RationalPolynomial rational_polynomial_one_value = { 1, { &one, &one } };
struct RationalPolynomial*rational_polynomial_one = &rational_polynomial_one_value;

struct RationalPolynomialDivision
{
    struct RationalPolynomial*quotient;
    struct RationalPolynomial*remainder;
};

struct NestedPolynomial
{
    size_t coefficient_count;
    struct RationalPolynomial*coefficients[];
};

struct NestedPolynomial nested_polynomial_one_value = { 1, &rational_polynomial_one_value };
struct NestedPolynomial*nested_polynomial_one = &nested_polynomial_one_value;

struct NestedPolynomialDivision
{
    struct NestedPolynomial*quotient;
    struct NestedPolynomial*remainder;
};

struct Number
{
    char operation;
    struct RationalPolynomial*minimal_polynomial;
    union
    {
        struct Rational value;//When operation == 'r'.
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

struct Number number_one;
struct Number**roots_of_unity;

__declspec(noreturn) void crash(char*message)
{
#ifdef _DEBUG
    puts(message);
    *(int*)0 = 0;
#endif
}

#ifdef _DEBUG
#define ASSERT(condition, message) if (!(condition)) { puts(message); *(int*)0 = 0; }
#else
#define ASSERT(condition, message)
#endif

#define ALLOCATE(output_stack, type) stack_slot_allocate(output_stack, sizeof(type), _Alignof(type))

#define ARRAY_ALLOCATE(output_stack, element_count, type)\
    stack_slot_allocate(output_stack, (element_count) * sizeof(type), _Alignof(type))

#define SWAP(a, b, type) { type temp = a; a = b; b = temp; }

#define POLY(name, type, coefficient_type, coefficient_count, ...)\
union\
{\
    struct\
    {\
        size_t count;\
        coefficient_type coefficients[coefficient_count];\
    };\
    type p;\
} name = { coefficient_count, { __VA_ARGS__ } };

#include "integer/integer.h"
#include "number/number.h"
#include "stack/stack.h"

void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count,
    size_t coefficient_size, size_t coefficient_alignment)
{
    size_t*out = stack_slot_allocate(output_stack, coefficient_size * (coefficient_count + 1),
        coefficient_alignment);
    *out = coefficient_count;
    return out;
}

#define POLYNOMIAL_ALLOCATE(output_stack, coefficient_count, coefficient_type)\
polynomial_allocate(output_stack, coefficient_count, sizeof(coefficient_type),\
    _Alignof(coefficient_type))

#endif
