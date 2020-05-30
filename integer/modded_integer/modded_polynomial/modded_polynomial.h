struct IntegerPolynomial*modded_polynomial_modded_integer_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*a, struct Integer*b, struct Integer*characteristic);
struct Integer*modded_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic);
struct Integer*modded_integer_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_negate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend, struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomialDivision modded_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*base, struct Integer*exponent,
    struct Integer*characteristic);
struct ExtendedGCDInfo modded_polynomial_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic);
struct IntegerPolynomial*modded_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic);
size_t squarefree_modded_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial**out, struct IntegerPolynomial*a,
    struct Integer*characteristic);