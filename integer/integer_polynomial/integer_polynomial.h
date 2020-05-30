struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a);
bool integer_polynomial_equals(struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_add(struct Stack*restrict output_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_negate(struct Stack*restrict output_stack,
    struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend);
struct IntegerPolynomial*integer_polynomial_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*b);
struct IntegerPolynomial*integer_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_integer_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend, struct Integer*divisor);
struct IntegerPolynomialDivision integer_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor);
struct IntegerPolynomial*integer_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*base, struct Integer*exponent);
struct Integer*integer_polynomial_get_content(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
struct IntegerPolynomial*integer_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct PolynomialExtendedGCDInfo integer_polynomial_get_extended_gcd(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b);
struct IntegerPolynomial*integer_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
size_t integer_polynomial_squarefree_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial**out);
size_t primitive_integer_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial**out, struct IntegerPolynomial*a);
struct RationalPolynomial*integer_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a);
struct RationalPolynomial*integer_polynomial_to_rational_polynomial(struct Stack*output_stack,
    struct IntegerPolynomial*a);