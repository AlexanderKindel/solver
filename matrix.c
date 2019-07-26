#include "declarations.h"

//Leaves garbage allocations on local stack along with the new matrix and augmentation values.
void leaking_matrix_row_echelon_form(void*(augmentation_element_rational_multiply)(struct Stack*,
    struct Stack*, void*, struct Rational*),
    void*(augmentation_element_subtract)(struct Stack*, struct Stack*, void*, void*),
    struct Stack*output_stack, struct Stack*local_stack,
    struct Matrix*a, void**augmentation)
{
    size_t small_dimension = min(a->width, a->height);
    for (size_t i = 0; i < small_dimension; ++i)
    {
        for (size_t j = i + 1; j < a->height; ++j)
        {
            if (a->rows[i][i]->numerator->value_count == 0)
            {
                POINTER_SWAP(a->rows[i], a->rows[j]);
                POINTER_SWAP(augmentation[i], augmentation[j]);
            }
            else
            {
                for (size_t k = i + 1; k < a->height; ++k)
                {
                    struct Rational*scalar =
                        rational_divide(local_stack, output_stack, a->rows[k][i], a->rows[i][i]);
                    for (size_t l = i; l < a->width; ++l)
                    {
                        a->rows[k][l] = rational_subtract(local_stack, output_stack, a->rows[k][l],
                            rational_multiply(local_stack, output_stack, a->rows[i][l], scalar));
                    }
                    augmentation[k] = augmentation_element_subtract(local_stack, output_stack,
                        augmentation[k], augmentation_element_rational_multiply(local_stack,
                            output_stack, augmentation[i], scalar));
                }
                break;
            }
        }
    }
}

void matrix_row_echelon_form(void*(augmentation_element_rational_multiply)(struct Stack*,
    struct Stack*, void*, struct Rational*),
    void*(augmentation_element_subtract)(struct Stack*, struct Stack*, void*, void*),
    struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a, void**augmentation)
{
    void*local_stack_savepoint = local_stack->cursor;
    leaking_matrix_row_echelon_form(augmentation_element_rational_multiply,
        augmentation_element_subtract, output_stack, local_stack, a, augmentation);
    for (size_t i = 0; i < a->height; ++i)
    {
        for (size_t j = 0; j < a->width; ++j)
        {
            a->rows[i][j] = rational_copy_to_stack(output_stack, a->rows[i][j]);
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

//Leaves the matrix in an invalid state; to be used for the sake of the augmentation.
void matrix_diagonalize(struct Stack*output_stack, struct Stack*local_stack, struct Matrix*a,
    struct Rational**augmentation)
{
    leaking_matrix_row_echelon_form(rational_multiply, rational_subtract, output_stack, local_stack,
        a, augmentation);
    for (size_t i = a->height; i-- > 0;)
    {
        augmentation[i] =
            rational_divide(local_stack, output_stack, augmentation[i], a->rows[i][i]);
        for (size_t j = 0; j <= i; ++j)
        {
            a->rows[i][j] =
                rational_divide(local_stack, output_stack, a->rows[i][j], a->rows[i][i]);
        }
        for (size_t j = 0; j < i; ++j)
        {
            augmentation[j] = rational_subtract(output_stack, local_stack, augmentation[j],
                rational_multiply(local_stack, output_stack, augmentation[i], a->rows[j][i]));
        }
    }
}