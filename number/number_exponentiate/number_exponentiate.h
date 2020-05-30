struct Number*number_surd_initialize(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*radicand, struct Integer*index);
struct Number*number_integer_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*exponent);
struct Number*number_rational_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Rational*exponent);
struct Number*number_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*exponent);
struct Number*number_take_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*index);