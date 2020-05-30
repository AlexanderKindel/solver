bool regions_are_disjoint(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b);
bool region_a_contains_b(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Region*a, struct Region*b);