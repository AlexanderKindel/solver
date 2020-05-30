struct RationalPolynomial*product_get_minimal_polynomial(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a,
    struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*b);
struct Number*number_multiply(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a, struct Number*b);