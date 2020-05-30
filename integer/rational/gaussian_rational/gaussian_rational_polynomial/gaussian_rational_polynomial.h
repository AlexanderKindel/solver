struct GaussianRationalPolynomial*gaussian_rational_polynomial_copy(struct Stack*output_stack,
    struct GaussianRationalPolynomial*a);
bool gaussian_rational_polynomial_equals(struct GaussianRationalPolynomial*a,
    struct GaussianRationalPolynomial*b);
void gaussian_rational_polynomial_trim_leading_zeroes(struct GaussianRationalPolynomial*a);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_add(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRationalPolynomial*b);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct Rational*b);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_gaussian_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRational*b);
struct GaussianRationalPolynomial*gaussian_rational_polynomial_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRationalPolynomial*b);