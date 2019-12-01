#include "declarations.h"

//Leaves a significant amount of excess allocations on output_stack.
void matrix_make_row_echelon_form(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Matrix*a, struct Rational**augmentation)
{
    size_t small_dimension;
    if (a->width < a->height)
    {
        small_dimension = a->width;
    }
    else
    {
        small_dimension = a->height;
    }
    void*local_stack_savepoint = local_stack->cursor;
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
                        a->rows[k][l] = rational_subtract(output_stack, local_stack, a->rows[k][l],
                            rational_multiply(local_stack, output_stack, a->rows[i][l], scalar));
                    }
                    augmentation[k] = rational_subtract(output_stack, local_stack, augmentation[k],
                        rational_multiply(local_stack, output_stack, augmentation[i], scalar));
                }
                break;
            }
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

//Leaves the matrix in an invalid state; to be used for the sake of the augmentation.
void matrix_diagonalize(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Matrix*a, struct Rational**augmentation)
{
    void*local_stack_savepoint = local_stack->cursor;
    matrix_make_row_echelon_form(local_stack, output_stack, a, augmentation);
    for (size_t i = a->height; i-- > 0;)
    {
        augmentation[i] =
            rational_divide(output_stack, local_stack, augmentation[i], a->rows[i][i]);
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
    local_stack->cursor = local_stack_savepoint;
}