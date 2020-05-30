struct NestedPolynomial*number_field_polynomial_number_field_element_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomial*number_field_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomialDivision number_field_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct NestedPolynomial*divisor,
    struct RationalPolynomial*generator_minimal_polynomial);
struct NestedPolynomial*number_field_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial);
size_t number_field_polynomial_squarefree_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial**out,
    struct RationalPolynomial*generator_minimal_polynomial);
size_t number_field_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial**out, struct NestedPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial);