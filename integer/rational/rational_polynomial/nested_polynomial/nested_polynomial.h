struct NestedPolynomial*nested_polynomial_copy(struct Stack*output_stack,
    struct NestedPolynomial*a);
bool nested_polynomial_equals(struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_negate(struct Stack*restrict output_stack,
    struct NestedPolynomial*a);
struct NestedPolynomial*nested_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*minuend,
    struct NestedPolynomial*subtrahend);
struct NestedPolynomial*nested_polynomial_rational_polynomial_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*b);
struct NestedPolynomial*nested_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_rational_polynomial_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct RationalPolynomial*divisor);
struct NestedPolynomialDivision nested_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct NestedPolynomial*divisor);
struct RationalPolynomial*nested_polynomial_get_content(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a);
struct RationalPolynomial*nested_polynomial_get_resultant(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b);
struct NestedPolynomial*nested_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a);