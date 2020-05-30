struct RationalPolynomial*sum_get_minimal_polynomial(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a,
    struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial);
struct Number*number_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a, struct Number*b);
struct Number*number_eliminate_linear_dependencies(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);