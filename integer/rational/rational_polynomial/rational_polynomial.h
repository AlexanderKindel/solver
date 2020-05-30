struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack,
    struct RationalPolynomial*a);
bool rational_polynomial_equals(struct RationalPolynomial*a, struct RationalPolynomial*b);
void rational_polynomial_trim_leading_zeroes(struct RationalPolynomial*a);
struct RationalPolynomial*rational_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_negate(struct Stack*restrict output_stack,
    struct RationalPolynomial*a);
struct RationalPolynomial*rational_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*minuend,
    struct RationalPolynomial*subtrahend);
struct RationalPolynomial*rational_polynomial_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Rational*b);
struct RationalPolynomial*rational_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomialDivision rational_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*dividend, struct RationalPolynomial*divisor);
struct RationalPolynomial*rational_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*base, struct Integer*exponent);
struct PolynomialExtendedGCDInfo rational_polynomial_get_extended_gcd(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b);
struct RationalPolynomial*rational_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a);
size_t rational_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial**out);
struct Rational rational_polynomial_evaluate_at_rational(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Rational*argument);
struct GaussianRational rational_polynomial_evaluate_at_gaussian_rational(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct GaussianRational*argument);
struct GaussianRationalPolynomial*rational_polynomial_evaluate_at_gaussian_rational_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct GaussianRationalPolynomial*argument);
size_t rational_polynomial_count_roots_in_rectangle(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a, struct RationalInterval*real,
    struct RationalInterval*imaginary);
size_t rational_polynomial_count_roots_in_region(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a, struct Region*region);
struct Region rational_polynomial_evaluate_at_region_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Region*argument,
    struct Rational*interval_size);
struct IntegerPolynomial*rational_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a);
struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a);

#include "integer/rational/rational_polynomial/nested_polynomial/nested_polynomial.h"
#include "integer/rational/rational_polynomial/number_field_element/number_field_element.h"