#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

char*substitution_index_translations[] =
{
    "main",
    "coefficient",
    "misc"
};

#define TYPE_NAME 0
#define BASE_TYPE 1
#define BASE_TYPE_NAME 2
#define ADDITIVE_IDENTITY 3
#define MULTIPLICATIVE_IDENTITY 4
#define ADD 5
#define NEGATE 6
#define SUBTRACT 7
#define MULTIPLY 8
#define MISC_NAME 9
#define GET_QUOTIENT 10
#define INTEGER_MULTIPLY 11
#define RATIONAL_MULTIPLY 12

char*operation_index_translations[] =
{
    [TYPE_NAME] = "type_name",
    [BASE_TYPE] = "base_type",
    [BASE_TYPE_NAME] = "base_type_name",
    [ADDITIVE_IDENTITY] = "additive_identity",
    [MULTIPLICATIVE_IDENTITY] = "multiplicative_identity",
    [ADD] = "add",
    [NEGATE] = "negate",
    [SUBTRACT] = "subtract",
    [MULTIPLY] ="multiply",
    [MISC_NAME] = "misc_name",
    [GET_QUOTIENT] = "get_quotient",
    [INTEGER_MULTIPLY] = "integer_multiply",
    [RATIONAL_MULTIPLY] = "rational_multiply"
};

#define HAS_MISC_PARAMETER 0
#define ADD_AND_SUBTRACT_TAKE_MISC 1
#define ADD_TAKES_TWO_STACKS 2
#define NEGATE_TAKES_TWO_STACKS 3
#define RECIPROCAL_TAKES_TWO_STACKS 4

char*condition_index_translations[] =
{
    [HAS_MISC_PARAMETER] = "has_misc_parameter",
    [ADD_AND_SUBTRACT_TAKE_MISC] = "add_and_subtract_take_misc",
    [ADD_TAKES_TWO_STACKS] = "add_takes_two_stacks",
    [NEGATE_TAKES_TWO_STACKS] = "negate_takes_two_stacks",
    [RECIPROCAL_TAKES_TWO_STACKS] = "reciprocal_takes_two_stacks"
};

struct Substitution
{
    char**operations;
    bool*conditions;
};

struct Substitution integer_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "integer",
        [BASE_TYPE] = "Integer",
        [BASE_TYPE_NAME] = "integer",
        [ADDITIVE_IDENTITY] = "&zero",
        [MULTIPLICATIVE_IDENTITY] = "&one",
        [ADD] = "integer_add",
        [NEGATE] = "integer_negate",
        [SUBTRACT] = "integer_subtract",
        [MULTIPLY] ="integer_multiply",
        [GET_QUOTIENT] = "integer_get_quotient",
        [INTEGER_MULTIPLY] = "integer_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = false,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution modded_integer_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "modded_integer",
        [BASE_TYPE] = "Integer",
        [BASE_TYPE_NAME] = "integer",
        [ADDITIVE_IDENTITY] = "&zero",
        [MULTIPLICATIVE_IDENTITY] = "&one",
        [ADD] = "modded_integer_add",
        [NEGATE] = "modded_integer_negate",
        [SUBTRACT] = "modded_integer_subtract",
        [MULTIPLY] ="modded_integer_multiply",
        [MISC_NAME] ="characteristic"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = true,
        [ADD_AND_SUBTRACT_TAKE_MISC] = true,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = true
    }
};

struct Substitution rational_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "rational",
        [BASE_TYPE] = "Rational",
        [BASE_TYPE_NAME] = "rational",
        [ADDITIVE_IDENTITY] = "&rational_zero",
        [MULTIPLICATIVE_IDENTITY] = "&rational_one",
        [ADD] = "rational_add",
        [NEGATE] = "rational_negate",
        [SUBTRACT] = "rational_subtract",
        [MULTIPLY] ="rational_multiply",
        [GET_QUOTIENT] = "rational_divide",
        [INTEGER_MULTIPLY] = "rational_integer_multiply",
        [RATIONAL_MULTIPLY] = "rational_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false,
        [RECIPROCAL_TAKES_TWO_STACKS] = false
    }
};

struct Substitution gaussian_rational_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "gaussian_rational",
        [BASE_TYPE] = "GaussianRational",
        [BASE_TYPE_NAME] = "gaussian_rational",
        [ADDITIVE_IDENTITY] = "&gaussian_rational_zero",
        [MULTIPLICATIVE_IDENTITY] = "&gaussian_rational_one",
        [ADD] = "gaussian_rational_add",
        [NEGATE] = "gaussian_rational_negate",
        [SUBTRACT] = "gaussian_rational_subtract",
        [MULTIPLY] = "gaussian_rational_multiply",
        [RATIONAL_MULTIPLY] = "gaussian_rational_rational_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution float_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "float",
        [BASE_TYPE] = "Float",
        [BASE_TYPE_NAME] = "float",
        [ADDITIVE_IDENTITY] = "&float_zero",
        [MULTIPLICATIVE_IDENTITY] = "&float_one",
        [ADD] = "float_add",
        [NEGATE] = "float_negate",
        [SUBTRACT] = "float_subtract",
        [MULTIPLY] ="float_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution integer_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "integer_polynomial",
        [BASE_TYPE] = "IntegerPolynomial",
        [BASE_TYPE_NAME] = "integer_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&integer_polynomial_one",
        [ADD] = "integer_polynomial_add",
        [NEGATE] = "integer_polynomial_negate",
        [SUBTRACT] = "integer_polynomial_subtract",
        [MULTIPLY] ="integer_polynomial_multiply",
        [GET_QUOTIENT] = "integer_polynomial_get_quotient"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = false,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution modded_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "modded_polynomial",
        [BASE_TYPE] = "IntegerPolynomial",
        [BASE_TYPE_NAME] = "integer_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&integer_polynomial_one",
        [ADD] = "modded_polynomial_add",
        [NEGATE] = "modded_polynomial_negate",
        [SUBTRACT] = "modded_polynomial_subtract",
        [MULTIPLY] = "modded_polynomial_multiply",
        [MISC_NAME] = "characteristic",
        [GET_QUOTIENT] = "modded_polynomial_get_quotient"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = true,
        [ADD_AND_SUBTRACT_TAKE_MISC] = true,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = true
    }
};

struct Substitution rational_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "rational_polynomial",
        [BASE_TYPE] = "RationalPolynomial",
        [BASE_TYPE_NAME] = "rational_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&rational_polynomial_one",
        [ADD] = "rational_polynomial_add",
        [NEGATE] = "rational_polynomial_negate",
        [SUBTRACT] = "rational_polynomial_subtract",
        [MULTIPLY] = "rational_polynomial_multiply",
        [GET_QUOTIENT] = "rational_polynomial_get_quotient",
        [INTEGER_MULTIPLY] = "rational_polynomial_integer_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution gaussian_rational_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "gaussian_rational_polynomial",
        [BASE_TYPE] = "GaussianRationalPolynomial",
        [BASE_TYPE_NAME] = "gaussian_rational_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&gaussian_rational_polynomial_one",
        [ADD] = "gaussian_rational_polynomial_add",
        [MULTIPLY] = "gaussian_rational_polynomial_multiply",
        [RATIONAL_MULTIPLY] = "gaussian_rational_polynomial_rational_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution nested_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "nested_polynomial",
        [BASE_TYPE] = "NestedPolynomial",
        [BASE_TYPE_NAME] = "nested_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&nested_polynomial_one",
        [ADD] = "nested_polynomial_add",
        [NEGATE] = "nested_polynomial_negate",
        [SUBTRACT] = "nested_polynomial_subtract",
        [MULTIPLY] = "nested_polynomial_multiply"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false
    }
};

struct Substitution number_field_element_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "number_field_element",
        [BASE_TYPE] = "RationalPolynomial",
        [BASE_TYPE_NAME] = "rational_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&rational_polynomial_one",
        [ADD] = "rational_polynomial_add",
        [NEGATE] = "rational_polynomial_negate",
        [SUBTRACT] = "rational_polynomial_subtract",
        [MULTIPLY] = "number_field_element_multiply",
        [MISC_NAME] = "generator_minimal_polynomial",
        [GET_QUOTIENT] = "number_field_element_divide"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = true,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false,
        [RECIPROCAL_TAKES_TWO_STACKS] = true
    }
};

struct Substitution number_field_polynomial_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "number_field_polynomial",
        [BASE_TYPE] = "NestedPolynomial",
        [BASE_TYPE_NAME] = "nested_polynomial",
        [ADDITIVE_IDENTITY] = "polynomial_zero",
        [MULTIPLICATIVE_IDENTITY] = "&nested_polynomial_one",
        [ADD] = "nested_polynomial_add",
        [NEGATE] = "nested_polynomial_negate",
        [SUBTRACT] = "nested_polynomial_subtract",
        [MULTIPLY] = "number_field_polynomial_multiply",
        [MISC_NAME] = "generator_minimal_polynomial"
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = true,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true,
        [NEGATE_TAKES_TWO_STACKS] = false,
        [RECIPROCAL_TAKES_TWO_STACKS] = true
    }
};

struct Substitution number_substitution =
{
    (char*[]) {
        [TYPE_NAME] = "number",
        [BASE_TYPE] = "Number",
        [BASE_TYPE_NAME] = "number",
        [MULTIPLICATIVE_IDENTITY] = "&number_one",
        [MULTIPLY] = "number_multiply",
    },
    (bool[]) {
        [HAS_MISC_PARAMETER] = false,
        [ADD_AND_SUBTRACT_TAKE_MISC] = false,
        [ADD_TAKES_TWO_STACKS] = true
    }
};

size_t translate_index_token(FILE*output_file, char**translations, size_t translation_count,
    char*token, char*token_end)
{
    size_t token_length = (size_t)(token_end - token);
    for (size_t translation_index = 0; translation_index < translation_count; ++translation_index)
    {
        char*translation = translations[translation_index];
        size_t token_index = 0;
        while (true)
        {
            if (token_index == token_length)
            {
                if (!translation[token_index])
                {
                    return translation_index;
                }
                break;
            }
            if (!translation[token_index])
            {
                break;
            }
            if (token[token_index] != translation[token_index])
            {
                break;
            }
            ++token_index;
        }
    }
    abort();
    return translation_count;
}

#define TRANSLATE_INDEX_TOKEN(output_file, translations, token, token_end)\
translate_index_token(output_file, translations, sizeof(translations) / sizeof(translations[0]),\
    token, token_end)

char*write_token(FILE*output_file, char*token, struct Substitution*substitution_set)
{
    char*cursor = token;
    while (*cursor != '.')
    {
        ++cursor;
    }
    struct Substitution*substitution = substitution_set +
        TRANSLATE_INDEX_TOKEN(output_file, substitution_index_translations, token, cursor);
    ++cursor;
    token = cursor;
    while (*cursor != '%')
    {
        ++cursor;
    }
    fputs(substitution->operations[TRANSLATE_INDEX_TOKEN(output_file,
        operation_index_translations, token, cursor)],
        output_file);
    return cursor + 1;
}

void instantiate_function_template(FILE*output_file, char*template_string,
    struct Substitution*substitution_set)
{
    char*start_of_next_block_to_write = template_string;
    char*cursor = template_string;
    while (true)
    {
        switch (*cursor)
        {
        case 0:
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), output_file);
            fputs("\n\n", output_file);
            return;
        case '%':
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), output_file);
            ++cursor;
            cursor = write_token(output_file, cursor, substitution_set);
            start_of_next_block_to_write = cursor;
            break;
        case '#':
            fwrite(start_of_next_block_to_write, 1,
                (size_t)(cursor - start_of_next_block_to_write), output_file);
            ++cursor;
            start_of_next_block_to_write = cursor;
            while (*cursor != '.')
            {
                ++cursor;
            }
            struct Substitution*substitution = substitution_set +
                TRANSLATE_INDEX_TOKEN(output_file, substitution_index_translations,
                    start_of_next_block_to_write, cursor);
            ++cursor;
            start_of_next_block_to_write = cursor;
            while (*cursor == '_' || (*cursor >= 'a' && *cursor <= 'z'))
            {
                ++cursor;
            }
            if (substitution->conditions[TRANSLATE_INDEX_TOKEN(output_file,
                condition_index_translations, start_of_next_block_to_write, cursor)])
            {
                start_of_next_block_to_write = cursor;
                while (*cursor != '#')
                {
                    if (*cursor == '%')
                    {
                        fwrite(start_of_next_block_to_write, 1,
                            (size_t)(cursor - start_of_next_block_to_write), output_file);
                        ++cursor;
                        cursor = write_token(output_file, cursor, substitution_set);
                        start_of_next_block_to_write = cursor;
                    }
                    else
                    {
                        ++cursor;
                    }
                }
                fwrite(start_of_next_block_to_write, 1,
                    (size_t)(cursor - start_of_next_block_to_write), output_file);
            }
            else
            {
                do
                {
                    ++cursor;
                } while (*cursor != '#');
            }
            ++cursor;
            start_of_next_block_to_write = cursor;
            break;
        default:
            ++cursor;
        }
    }
}

#define INTEGER 0
#define MODDED_INTEGER 1
#define RATIONAL 2
#define GAUSSIAN_RATIONAL 3
#define FLOAT 4
#define INTEGER_POLYNOMIAL 5
#define MODDED_POLYNOMIAL 6
#define RATIONAL_POLYNOMIAL 7
#define GAUSSIAN_RATIONAL_POLYNOMIAL 8
#define NESTED_POLYNOMIAL 9
#define NUMBER_FIELD_ELEMENT 10
#define NUMBER_FIELD_POLYNOMIAL 11
#define NUMBER 12

struct Substitution**substitution_sets;

void instantiate_function_template_set(FILE*file, char*template_file_name,
    size_t*substitution_set_indices, size_t substitution_set_index_count)
{
    FILE*template_file = fopen(template_file_name, "r");
    size_t template_length = 0;
    char character;
    do
    {
        fread(&character, 1, 1, template_file);
        ++template_length;
    } while (!feof(template_file));
    char*template_string = malloc(template_length);
    rewind(template_file);
    --template_length;
    fread(template_string, template_length, 1, template_file);
    template_string[template_length] = 0;
    fclose(template_file);
    for (size_t i = 0; i < substitution_set_index_count; ++i)
    {
        instantiate_function_template(file, template_string,
            substitution_sets[substitution_set_indices[i]]);
    }
    free(template_string);
}

int main()
{
    substitution_sets = (struct Substitution*[]) {
        [INTEGER] = (struct Substitution[]) { integer_substitution },
        [MODDED_INTEGER] = (struct Substitution[]) { modded_integer_substitution },
        [RATIONAL] = (struct Substitution[]) { rational_substitution },
        [GAUSSIAN_RATIONAL] = (struct Substitution[]) { gaussian_rational_substitution },
        [FLOAT] = (struct Substitution[]) { float_substitution },
        [INTEGER_POLYNOMIAL] = (struct Substitution[]) { integer_polynomial_substitution,
            integer_substitution },
        [MODDED_POLYNOMIAL] = (struct Substitution[]) { modded_polynomial_substitution,
            modded_integer_substitution, integer_substitution },
        [RATIONAL_POLYNOMIAL] = (struct Substitution[]) { rational_polynomial_substitution,
            rational_substitution },
        [GAUSSIAN_RATIONAL_POLYNOMIAL] = (struct Substitution[]) {
            gaussian_rational_polynomial_substitution, gaussian_rational_substitution },
        [NESTED_POLYNOMIAL] = (struct Substitution[]) { nested_polynomial_substitution,
                rational_polynomial_substitution },
        [NUMBER_FIELD_ELEMENT] = (struct Substitution[]) { number_field_element_substitution,
            rational_polynomial_substitution, rational_polynomial_substitution },
        [NUMBER_FIELD_POLYNOMIAL] = (struct Substitution[]) { number_field_polynomial_substitution,
            number_field_element_substitution, rational_polynomial_substitution },
        [NUMBER] = (struct Substitution[]) { number_substitution }
    };
    FILE*template_instantiations = fopen("template_instantiations.c", "w");
    fputs("struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack, struct RationalPolynomial*a);\n\n",
        template_instantiations);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_copy_coefficients.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, GAUSSIAN_RATIONAL_POLYNOMIAL,
            NESTED_POLYNOMIAL }, 4);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_copy.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, GAUSSIAN_RATIONAL_POLYNOMIAL,
            NESTED_POLYNOMIAL }, 4);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_equals.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, GAUSSIAN_RATIONAL_POLYNOMIAL,
            NESTED_POLYNOMIAL }, 4);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_trim_leading_zeroes.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, GAUSSIAN_RATIONAL_POLYNOMIAL,
            NESTED_POLYNOMIAL }, 4);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_add.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, MODDED_POLYNOMIAL, RATIONAL_POLYNOMIAL,
            GAUSSIAN_RATIONAL_POLYNOMIAL, NESTED_POLYNOMIAL }, 5);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_negate.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, MODDED_POLYNOMIAL, RATIONAL_POLYNOMIAL,
            NESTED_POLYNOMIAL }, 4);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_subtract.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, NESTED_POLYNOMIAL }, 3);
    instantiate_function_template_set(template_instantiations, "templates/polynomial_multiply.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, MODDED_POLYNOMIAL, RATIONAL_POLYNOMIAL,
            GAUSSIAN_RATIONAL_POLYNOMIAL, NESTED_POLYNOMIAL, NUMBER_FIELD_POLYNOMIAL }, 6);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_coefficient_multiply.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, MODDED_POLYNOMIAL, RATIONAL_POLYNOMIAL,
            GAUSSIAN_RATIONAL_POLYNOMIAL, NESTED_POLYNOMIAL, NUMBER_FIELD_POLYNOMIAL }, 6);
    rational_substitution.operations[MULTIPLY] = "rational_unreduced_multiply";
    instantiate_function_template_set(template_instantiations, "templates/exponentiate.txt",
        (size_t[]) { INTEGER, RATIONAL, FLOAT, INTEGER_POLYNOMIAL, MODDED_POLYNOMIAL,
            RATIONAL_POLYNOMIAL, NUMBER_FIELD_ELEMENT, NUMBER}, 8);
    rational_substitution.operations[MULTIPLY] = "rational_multiply";
    instantiate_function_template_set(template_instantiations,
        "templates/field_polynomial_euclidean_divide.txt",
        (size_t[]) { MODDED_POLYNOMIAL, RATIONAL_POLYNOMIAL, NUMBER_FIELD_POLYNOMIAL }, 3);
    instantiate_function_template_set(template_instantiations,
        "templates/euclidean_domain_polynomial_euclidean_divide.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, NESTED_POLYNOMIAL }, 2);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_divide_by_coefficient.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, NESTED_POLYNOMIAL }, 2);
    instantiate_function_template_set(template_instantiations, "templates/get_gcd.txt",
        (size_t[]) { INTEGER, MODDED_POLYNOMIAL, NUMBER_FIELD_POLYNOMIAL }, 3);
    instantiate_function_template_set(template_instantiations, "templates/get_extended_gcd.txt",
        (size_t[]) { INTEGER, MODDED_POLYNOMIAL }, 2);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_get_content.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, NESTED_POLYNOMIAL }, 2);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_get_derivative.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, RATIONAL_POLYNOMIAL, NESTED_POLYNOMIAL }, 3);
    instantiate_function_template_set(template_instantiations,
        "templates/polynomial_squarefree_factor.txt",
        (size_t[]) { INTEGER_POLYNOMIAL, NUMBER_FIELD_POLYNOMIAL }, 2);
    instantiate_function_template_set(template_instantiations,
        "templates/rational_polynomial_evaluate.txt",
        (size_t[]) { RATIONAL, GAUSSIAN_RATIONAL, GAUSSIAN_RATIONAL_POLYNOMIAL }, 3);
    instantiate_function_template_set(template_instantiations,
        "templates/interval_get_max_magnitude.txt", (size_t[]) { RATIONAL, FLOAT }, 2);
    instantiate_function_template_set(template_instantiations, "templates/interval_add.txt",
        (size_t[]) { RATIONAL, FLOAT }, 2);
    instantiate_function_template_set(template_instantiations, "templates/interval_negate.txt",
        (size_t[]) { RATIONAL, FLOAT }, 2);
    instantiate_function_template_set(template_instantiations, "templates/interval_multiply.txt",
        (size_t[]) { RATIONAL, FLOAT }, 2);
    fclose(template_instantiations);
    return 0;
}