#include "ring_templates.c"
#include "interval_templates.c"

void instantiate_function_template(FILE*file, char*function_template,
    struct Substitution*substitution_set)
{
    char*start_of_next_block_to_write = function_template;
    char*cursor = function_template;
    while (true)
    {
        switch (*cursor)
        {
        case 0:
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), file);
            fputs("\n", file);
            return;
        case '%':
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), file);
            fputs(substitution_set[*(cursor + 1) - '0'].strings[*(cursor + 2) - '0'], file);
            start_of_next_block_to_write = cursor + 3;
            cursor += 2;
            break;
        case '#':
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), file);
            if (substitution_set[*(cursor + 1) - '0'].switches[*(cursor + 2) - '0'])
            {
                cursor += 3;
                start_of_next_block_to_write = cursor;
                while (*cursor != '#')
                {
                    if (*cursor == '%')
                    {
                        fwrite(start_of_next_block_to_write, 1,
                            (size_t)(cursor - start_of_next_block_to_write), file);
                        fputs(substitution_set[*(cursor + 1) - '0'].strings[*(cursor + 2) - '0'],
                            file);
                        cursor += 2;
                        start_of_next_block_to_write = cursor + 1;
                    }
                    ++cursor;
                }
                fwrite(start_of_next_block_to_write, 1,
                    (size_t)(cursor - start_of_next_block_to_write), file);
            }
            else
            {
                do
                {
                    ++cursor;
                } while (*cursor != '#');
            }
            start_of_next_block_to_write = cursor + 1;
        }
        ++cursor;
    }
}

void instantiate_function_template_set(FILE*file, char*function_template,
    struct Substitution**substitution_sets, unsigned short substitution_flags)
{
    for (int i = 0; i < 13; ++i)
    {
        if ((1 << i) & substitution_flags)
        {
            instantiate_function_template(file, function_template, substitution_sets[i]);
        }
    }
}

int main()
{
    FILE*template_instantiations = fopen("template_instantiations.c", "w");
    fputs("struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack, struct RationalPolynomial*a);\n\n",
        template_instantiations);
    struct Substitution*ring_substitution_sets[] =
    {
        (struct Substitution[]) { integer_substitution },
        (struct Substitution[]) { modded_integer_substitution },
        (struct Substitution[]) { rational_substitution, { 0 },
            { (char*[]){ "rational_multiply" }, 0 } },
        (struct Substitution[]) { gaussian_rational_substitution, { 0 },
            { (char*[]){ "gaussian_rational_rational_multiply" }, 0 } },
        (struct Substitution[]) { float_substitution },
        (struct Substitution[]) { integer_polynomial_substitution, integer_substitution },
        (struct Substitution[]) { modded_polynomial_substitution, modded_integer_substitution,
            integer_substitution },
        (struct Substitution[]) { rational_polynomial_substitution, rational_substitution },
        (struct Substitution[]) { gaussian_rational_polynomial_substitution,
            gaussian_rational_substitution,
            { (char*[]){ "gaussian_rational_polynomial_rational_multiply" }, 0 } },
        (struct Substitution[]) { nested_polynomial_substitution,
            rational_polynomial_substitution },
        (struct Substitution[]) { number_field_element_substitution,
            rational_polynomial_substitution, rational_polynomial_substitution },
        (struct Substitution[]) { number_field_polynomial_substitution,
            number_field_element_substitution, rational_polynomial_substitution },
        (struct Substitution[]) { number_substitution }
    };
    instantiate_function_template_set(template_instantiations,
        polynomial_copy_coefficients_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | GAUSSIAN_RATIONAL_POLYNOMIAL |
        NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_copy_template,
        ring_substitution_sets,
        INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | GAUSSIAN_RATIONAL_POLYNOMIAL |
        NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_equals_template,
        ring_substitution_sets,
        INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | GAUSSIAN_RATIONAL_POLYNOMIAL |
        NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        polynomial_trim_leading_zeroes_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | GAUSSIAN_RATIONAL_POLYNOMIAL |
        NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_add_template,
        ring_substitution_sets,
        INTEGER_POLYNOMIAL | MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL |
        GAUSSIAN_RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_negate_template,
        ring_substitution_sets,
        INTEGER_POLYNOMIAL | MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_subtract_template,
        ring_substitution_sets, INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_multiply_template,
        ring_substitution_sets,
        INTEGER_POLYNOMIAL | MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL |
        GAUSSIAN_RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL | NUMBER_FIELD_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        polynomial_multiply_by_coefficient_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL |
        GAUSSIAN_RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL | NUMBER_FIELD_POLYNOMIAL);
    ring_substitution_sets[2][0].strings[9] = "rational_unreduced_multiply";
    instantiate_function_template_set(template_instantiations,
        exponentiate_template, ring_substitution_sets,
        INTEGER | RATIONAL | FLOAT | INTEGER_POLYNOMIAL | MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL |
        NUMBER_FIELD_ELEMENT | NUMBER);
    ring_substitution_sets[2][0].strings[9] = "rational_multiply";
    instantiate_function_template_set(template_instantiations, field_polynomial_divide_template,
        ring_substitution_sets, MODDED_POLYNOMIAL | RATIONAL_POLYNOMIAL | NUMBER_FIELD_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        euclidean_domain_polynomial_divide_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        polynomial_divide_by_coefficient_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, get_gcd_template,
        ring_substitution_sets, INTEGER | MODDED_POLYNOMIAL | NUMBER_FIELD_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, get_extended_gcd_template,
        ring_substitution_sets, INTEGER | MODDED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_content_template,
        ring_substitution_sets, INTEGER_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations, polynomial_derivative_template,
        ring_substitution_sets, INTEGER_POLYNOMIAL | RATIONAL_POLYNOMIAL | NESTED_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        polynomial_squarefree_factor_template, ring_substitution_sets,
        INTEGER_POLYNOMIAL | NUMBER_FIELD_POLYNOMIAL);
    instantiate_function_template_set(template_instantiations,
        rational_polynomial_evaluate_template, ring_substitution_sets,
        RATIONAL | GAUSSIAN_RATIONAL | GAUSSIAN_RATIONAL_POLYNOMIAL);
    struct Substitution*interval_substitution_sets[] =
    {
        (struct Substitution[]) { { (char*[]){ "Rational", "rational" }, 0 } },
        (struct Substitution[]) { { (char*[]){ "Float", "float" }, 0 } }
    };
    instantiate_function_template_set(template_instantiations, interval_get_max_magnitude_template,
        interval_substitution_sets, ALL_INTERVAL_SUBSTITUTION_SETS);
    instantiate_function_template_set(template_instantiations, interval_add_template,
        interval_substitution_sets, ALL_INTERVAL_SUBSTITUTION_SETS);
    instantiate_function_template_set(template_instantiations, interval_negate_template,
        interval_substitution_sets, ALL_INTERVAL_SUBSTITUTION_SETS);
    instantiate_function_template_set(template_instantiations, interval_multiply_template,
        interval_substitution_sets, ALL_INTERVAL_SUBSTITUTION_SETS);
    fclose(template_instantiations);
    return 0;
}