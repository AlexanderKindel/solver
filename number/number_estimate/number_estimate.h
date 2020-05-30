struct FloatInterval number_get_float_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval number_get_rational_real_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct FloatInterval number_get_float_imaginary_part_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval number_get_rational_imaginary_part_estimate(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a,
    struct Rational*interval_size);
struct Region number_get_region_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
void number_halve_region_estimate_dimensions(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Region*estimate, struct Number*a);
struct FloatInterval number_get_float_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval number_get_rational_magnitude_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval number_estimate_argument(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);
struct RationalInterval number_estimate_argument_in_radians(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Rational*interval_size);