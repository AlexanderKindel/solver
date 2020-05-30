struct Number**number_get_conjugates(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);
struct Number*number_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a);
struct Number*number_divide(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*dividend, struct Number*divisor);