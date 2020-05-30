void matrix_make_row_echelon_form(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Matrix*a, struct Rational*augmentation);
void matrix_diagonalize(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Matrix*a, struct Rational*augmentation);