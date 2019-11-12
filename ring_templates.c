#include <stdbool.h>
#include <stdio.h>

struct RingStrings
{
    char*misc_name;//0
    char*type_name;//1
    char*base_type;//2
    char*base_type_name;//3
    char*additive_identity;//4
    char*multiplicative_identity;//5
    char*add;//6
    char*negative;//7
    char*subtract;//8
    char*multiply;//9
};

struct EuclideanDomainStrings
{
    struct RingStrings ring_strings;
    char*euclidean_quotient;//:
    char*euclidean_divide;//;
    char*gcd;//<
    char*coefficient_times_integer;//=
};

struct FieldStrings
{
    struct RingStrings ring_strings;
    char*divide;//:
    char*reciprocal;//;
};

struct Switches
{
    bool has_misc_parameter;//0
    bool add_and_subtract_take_misc;//1
    bool add_takes_two_stacks;//2
    bool negative_takes_two_stacks;//3
    bool reciprocal_takes_two_stacks;//4
};

struct Substitution
{
    char**strings;
    bool*switches;
};

struct Substitution integer_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            0,
            "integer",
            "Integer",
            "integer",
            "&zero",
            "&one",
            "integer_add",
            "integer_negative",
            "integer_subtract",
            "integer_multiply"
        },
        "integer_euclidean_quotient",
        "integer_euclidean_divide",
        "integer_gcd"
    },
    (bool*)&(struct Switches) { false, false, false, false }
};

struct Substitution modded_integer_substitution =
{
    (char**)&(struct FieldStrings)
    {
        {
            0,
            "modded_integer",
            "Integer",
            "integer",
            "&zero",
            "&one",
            "modded_integer_add",
            "modded_integer_negative",
            "modded_integer_subtract",
            "modded_integer_multiply"
        },
        "!not_implemented!",
        "modded_integer_reciprocal"
    },
    (bool*)&(struct Switches) { true, true, true, true, true }
};

struct Substitution rational_substitution =
{
    (char**)&(struct FieldStrings)
    {
        {
            0,
            "rational",
            "Rational",
            "rational",
            "&rational_zero",
            "&rational_one",
            "rational_add",
            "rational_negative",
            "rational_subtract",
            "rational_multiply"
        },
        "rational_divide",
        "rational_reciprocal"
    },
    (bool*)&(struct Switches) { false, false, true, false, false }
};

struct Substitution gaussian_rational_substitution =
{
    (char**)&(struct RingStrings)
    {
        0,
        "gaussian_rational",
        "GaussianRational",
        "gaussian_rational",
        "&gaussian_rational_zero",
        "&gaussian_rational_one",
        "gaussian_rational_add",
        "gaussian_rational_negative",
        "gaussian_rational_subtract",
        "gaussian_rational_multiply"
    },
    (bool*)&(struct Switches) { false, false, true, false }
};

struct Substitution float_substitution =
{
    (char**)&(struct RingStrings)
    {
        0,
        "float",
        "Float",
        "float",
        "&float_zero",
        "&float_one",
        "float_add",
        "float_negative",
        "float_subtract",
        "float_multiply"
    },
    (bool*)&(struct Switches) { false, false, true, false }
};

struct Substitution integer_polynomial_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            0,
            "integer_polynomial",
            "IntegerPolynomial",
            "integer_polynomial",
            "polynomial_zero",
            "&integer_polynomial_one",
            "integer_polynomial_add",
            "integer_polynomial_negative",
            "integer_polynomial_subtract",
            "integer_polynomial_multiply"
        },
        "integer_polynomial_euclidean_quotient",
        "integer_polynomial_euclidean_divide",
        "integer_polynomial_gcd",
        "integer_multiply"
    },
    (bool*)&(struct Switches) { false, false, false, false }
};

struct Substitution modded_polynomial_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            "characteristic",
            "modded_polynomial",
            "IntegerPolynomial",
            "integer_polynomial",
            "polynomial_zero",
            "&integer_polynomial_one",
            "modded_polynomial_add",
            "modded_polynomial_negative",
            "modded_polynomial_subtract",
            "modded_polynomial_multiply"
        },
        "modded_polynomial_euclidean_quotient",
        "modded_polynomial_euclidean_divide",
        "$not_implemented$"
    },
    (bool*)&(struct Switches) { true, true, true, true }
};

struct Substitution rational_polynomial_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            0,
            "rational_polynomial",
            "RationalPolynomial",
            "rational_polynomial",
            "polynomial_zero",
            "&rational_polynomial_one",
            "rational_polynomial_add",
            "rational_polynomial_negative",
            "rational_polynomial_subtract",
            "rational_polynomial_multiply"
        },
        "rational_polynomial_euclidean_quotient",
        "rational_polynomial_euclidean_divide",
        "rational_polynomial_gcd",
        "rational_integer_multiply"
    },
    (bool*)&(struct Switches) { false, false, true, false }
};

struct Substitution gaussian_rational_polynomial_substitution =
{
    (char**)&(struct RingStrings)
    {
        0,
        "gaussian_rational_polynomial",
        "GaussianRationalPolynomial",
        "gaussian_rational_polynomial",
        "polynomial_zero",
        "&gaussian_rational_polynomial_one",
        "gaussian_rational_polynomial_add",
        "!not_implemented!",
        "!not_implemented!",
        "gaussian_rational_polynomial_multiply"
    },
    (bool*)&(struct Switches) { false, false, true, false }
};

struct Substitution nested_polynomial_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            0,
            "nested_polynomial",
            "NestedPolynomial",
            "nested_polynomial",
            "polynomial_zero",
            "&nested_polynomial_one",
            "nested_polynomial_add",
            "nested_polynomial_negative",
            "nested_polynomial_subtract",
            "nested_polynomial_multiply"
        },
        "!not_implemented!",
        "!not_implemented!",
        "!not_implemented!",
        "rational_polynomial_integer_multiply"
    },
    (bool*)&(struct Switches) { false, false, true, false }
};

struct Substitution number_field_element_substitution =
{
    (char**)&(struct FieldStrings)
    {
        {
            "generator_minimal_polynomial",
            "number_field_element",
            "RationalPolynomial",
            "rational_polynomial",
            "polynomial_zero",
            "&rational_polynomial_one",
            "rational_polynomial_add",
            "rational_polynomial_negative",
            "rational_polynomial_subtract",
            "number_field_element_multiply"
        },
        "number_field_element_divide",
        "number_field_element_reciprocal"
    },
    (bool*)&(struct Switches) { true, false, true, false, true }
};

struct Substitution number_field_polynomial_substitution =
{
    (char**)&(struct EuclideanDomainStrings)
    {
        {
            "generator_minimal_polynomial",
            "number_field_polynomial",
            "NestedPolynomial",
            "nested_polynomial",
            "polynomial_zero",
            "&nested_polynomial_one",
            "nested_polynomial_add",
            "nested_polynomial_negative",
            "nested_polynomial_subtract",
            "number_field_polynomial_multiply"
        },
        "!not_implemented!",
        "number_field_polynomial_euclidean_divide",
        "number_field_polynomial_gcd"
    },
    (bool*)&(struct Switches) { true, false, true, false, true }
};

struct Substitution number_substitution =
{
    (char**)&(struct RingStrings)
    {
        0,
        "number",
        "Number",
        "number",
        "!not_implemented!",
        "&number_one",
        "!not_implemented!",
        "!not_implemented!",
        "!not_implemented!",
        "number_multiply"
    },
    (bool*)&(struct Switches) { false, false, true }
};

unsigned short INTEGER = 0b1;
unsigned short MODDED_INTEGER = 0b10;
unsigned short RATIONAL = 0b100;
unsigned short GAUSSIAN_RATIONAL = 0b1000;
unsigned short FLOAT = 0b10000;
unsigned short INTEGER_POLYNOMIAL = 0b100000;
unsigned short MODDED_POLYNOMIAL = 0b1000000;
unsigned short RATIONAL_POLYNOMIAL = 0b10000000;
unsigned short GAUSSIAN_RATIONAL_POLYNOMIAL = 0b100000000;
unsigned short NESTED_POLYNOMIAL = 0b1000000000;
unsigned short NUMBER_FIELD_ELEMENT = 0b10000000000;
unsigned short NUMBER_FIELD_POLYNOMIAL = 0b100000000000;
unsigned short NUMBER = 0b1000000000000;

char*polynomial_copy_coefficients_template =
"void %01_copy_coefficients(struct Stack*output_stack, struct %02*a)\n"
"{\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        a->coefficients[i] = %13_copy(output_stack, a->coefficients[i]);\n"
"    }\n"
"}\n";

char*polynomial_copy_template =
"struct %02*%01_copy(struct Stack*output_stack, struct %02*a)\n"
"{\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count);\n"
"    memcpy(out->coefficients, a->coefficients, a->coefficient_count * sizeof(void*));\n"
"    %01_copy_coefficients(output_stack, out);\n"
"    return out;\n"
"}\n";

char*polynomial_equals_template =
"bool %01_equals(struct %02*a, struct %02*b)\n"
"{\n"
"    if (a->coefficient_count != b->coefficient_count)\n"
"    {\n"
"        return false;\n"
"    }\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        if (!%13_equals(a->coefficients[i], b->coefficients[i]))\n"
"        {\n"
"            return false;\n"
"        }\n"
"    }\n"
"    return true;\n"
"}\n";

char*polynomial_trim_leading_zeroes_template =
"void %01_trim_leading_zeroes(struct %02*a)\n"
"{\n"
"    for (size_t i = a->coefficient_count; i-- > 0;)\n"
"    {\n"
"        if (%13_equals(a->coefficients[i], %14))\n"
"        {\n"
"            a->coefficient_count -= 1;\n"
"        }\n"
"        else\n"
"        {\n"
"            return;\n"
"        }\n"
"    }\n"
"}\n";

char*polynomial_add_template =
"struct %02*%01_add(struct Stack*output_stack,#02 struct Stack*local_stack,# struct %02*a, struct %02*b#01, struct %22*%00#)\n"
"{\n"
"    if (a->coefficient_count < b->coefficient_count)\n"
"    {\n"
"        POINTER_SWAP(a, b);\n"
"    }\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count);\n"
"    for (size_t i = 0; i < b->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %16(output_stack#12, local_stack#, a->coefficients[i], b->coefficients[i]#11, %00#);\n"
"    }\n"
"    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %13_copy(output_stack, a->coefficients[i]);\n"
"    }\n"
"    %03_trim_leading_zeroes(out);\n"
"    return out;\n"
"}\n";

char*polynomial_negative_template =
"struct %02*%01_negative(struct Stack*output_stack#03, struct Stack*local_stack#, struct %02*a#01, struct %22*%00#)\n"
"{\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count);\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %17(output_stack#13, local_stack#, a->coefficients[i]#01, %00#);\n"
"    }\n"
"    return out;\n"
"}\n";

char*polynomial_subtract_template =
"struct %02*%01_subtract(struct Stack*output_stack, struct Stack*local_stack, struct %02*minuend, struct %02*subtrahend)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %02*out = %06(output_stack#02, local_stack#, minuend, %07(local_stack, subtrahend));\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*polynomial_multiply_template =
"struct %02*%01_multiply(struct Stack*output_stack, struct Stack*local_stack, struct %02*a, struct %02*b#00, struct %22*%00#)\n"
"{\n"
"    if (!a->coefficient_count && !b->coefficient_count)\n"
"    {\n"
"        return %04;\n"
"    }\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count + b->coefficient_count - 1);\n"
"    for (size_t i = 0; i < out->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %14;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        for (size_t j = 0; j < b->coefficient_count; ++j)\n"
"        {\n"
"            out->coefficients[i + j] = %16(output_stack#12, local_stack#, out->coefficients[i + j],\n"
"                %19(local_stack, output_stack, a->coefficients[i], b->coefficients[j]#10, %00#)#11, %00#);\n"
"        }\n"
"    }\n"
"    %03_trim_leading_zeroes(out);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*polynomial_multiply_by_coefficient_template =
"struct %02*%01_%11_multiply(struct Stack*output_stack, struct Stack*local_stack, struct %02*a, struct %12*b#00, struct %22*%00#)\n"
"{\n"
"    if (%13_equals(b, %14))\n"
"    {\n"
"        return %04;\n"
"    }\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count);\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %19(output_stack, local_stack, a->coefficients[i], b#00, %00#);\n"
"    }\n"
"    return out;\n"
"}\n";

char*exponentiate_template =
"struct %02*%01_exponentiate(struct Stack*output_stack, struct Stack*local_stack, struct %02*base, struct Integer*exponent#00, struct %22*%00#)\n"
"{\n"
"    if (!exponent->value_count)\n"
"    {\n"
"        return %05;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %02*out = %05;\n"
"    while (true)\n"
"    {\n"
"        if (exponent->value[0] & 1)\n"
"        {\n"
"            out = %09(local_stack, output_stack, out, base#00, %00#);\n"
"        }\n"
"        exponent = integer_half(local_stack, exponent);\n"
"        if (!exponent->value_count)\n"
"        {\n"
"            out = %03_copy(output_stack, out);\n"
"            local_stack->cursor = local_stack_savepoint;\n"
"            return out;\n"
"        }\n"
"        base = %09(local_stack, output_stack, base, base#00, %00#);\n"
"    }\n"
"}\n";

char*field_polynomial_divide_template =
"void %01_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack, struct %02Division*out, struct %02*dividend, struct %02*divisor#00, struct %22*%00#)\n"
"{\n"
"    if (divisor->coefficient_count > dividend->coefficient_count)\n"
"    {\n"
"        out->quotient = %04;\n"
"        out->remainder = %03_copy(output_stack, dividend);\n"
"        return;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    out->quotient = polynomial_allocate(output_stack, 1 + dividend->coefficient_count - divisor->coefficient_count);\n"
"    out->remainder = polynomial_allocate(output_stack, dividend->coefficient_count);\n"
"    memcpy(out->remainder->coefficients, dividend->coefficients, dividend->coefficient_count * sizeof(void*));\n"
"    struct %12*leading_coefficient_reciprocal = %1;(local_stack#14, output_stack#, divisor->coefficients[divisor->coefficient_count - 1]#00, %00#);\n"
"    while (out->remainder->coefficient_count >= divisor->coefficient_count)\n"
"    {\n"
"        --out->remainder->coefficient_count;\n"
"        struct %12*quotient = %19(local_stack, output_stack, out->remainder->coefficients[out->remainder->coefficient_count], leading_coefficient_reciprocal#00, %00#);\n"
"        out->quotient->coefficients[out->remainder->coefficient_count + 1 - divisor->coefficient_count] = quotient;\n"
"        for (size_t i = 1; i < divisor->coefficient_count; ++i)\n"
"        {\n"
"            out->remainder->coefficients[out->remainder->coefficient_count - i] =\n"
"                %18(local_stack, output_stack, out->remainder->coefficients[out->remainder->coefficient_count - i],\n"
"                    %19(local_stack, output_stack, quotient, divisor->coefficients[divisor->coefficient_count - i - 1]#00, %00#)#01, %00#);\n"
"        }\n"
"    }\n"
"    %03_copy_coefficients(output_stack, out->quotient);\n"
"    %03_trim_leading_zeroes(out->remainder);\n"
"    %03_copy_coefficients(output_stack, out->remainder);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"}\n";

char*euclidean_domain_polynomial_divide_template =
"//Sets out->quotient to 0 when the coefficients of the quotient wouldn't be in the ring.\n"
"void %01_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack, struct %02Division*out, struct %02*dividend, struct %02*divisor)\n"
"{\n"
"    if (divisor->coefficient_count > dividend->coefficient_count)\n"
"    {\n"
"        out->quotient = %04;\n"
"        out->remainder = %03_copy(output_stack, dividend);\n"
"        return;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    out->quotient = polynomial_allocate(output_stack, 1 + dividend->coefficient_count - divisor->coefficient_count);\n"
"    out->remainder = polynomial_allocate(output_stack, dividend->coefficient_count);\n"
"    memcpy(out->remainder->coefficients, dividend->coefficients, dividend->coefficient_count * sizeof(void*));\n"
"    while (out->remainder->coefficient_count >= divisor->coefficient_count)\n"
"    {\n"
"        --out->remainder->coefficient_count;\n"
"        struct %12Division division;\n"
"        %1;(local_stack, output_stack, &division, out->remainder->coefficients[out->remainder->coefficient_count], divisor->coefficients[divisor->coefficient_count - 1]);\n"
"        if (%13_equals(division.remainder, %14))\n"
"        {\n"
"            out->quotient->coefficients[out->remainder->coefficient_count + 1 - divisor->coefficient_count] = division.quotient;\n"
"            for (size_t i = 1; i < divisor->coefficient_count; ++i)\n"
"            {\n"
"                out->remainder->coefficients[out->remainder->coefficient_count - i] =\n"
"                    %18(local_stack, output_stack, out->remainder->coefficients[out->remainder->coefficient_count - i],\n"
"                        %19(local_stack, output_stack, division.quotient, divisor->coefficients[divisor->coefficient_count - i - 1]));\n"
"            }\n"
"        }\n"
"        else\n"
"        {\n"
"            out->quotient = 0;\n"
"            local_stack->cursor = local_stack_savepoint;\n"
"            return;\n"
"        }\n"
"    }\n"
"    %03_copy_coefficients(output_stack, out->quotient);\n"
"    %03_trim_leading_zeroes(out->remainder);\n"
"    %03_copy_coefficients(output_stack, out->remainder);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"}\n";

char*polynomial_divide_by_coefficient_template =
"struct %02*%01_%11_divide(struct Stack*output_stack, struct Stack*local_stack, struct %02*dividend, struct %12*divisor)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %02*out = polynomial_allocate(output_stack, dividend->coefficient_count);\n"
"    for (size_t i = 0; i < dividend->coefficient_count; ++i)\n"
"    {\n"
"        out->coefficients[i] = %1:(output_stack, local_stack, dividend->coefficients[i], divisor);\n"
"    }\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*gcd_template =
"struct %02*%01_gcd(struct Stack*output_stack, struct Stack*local_stack, struct %02*a, struct %02*b#00, struct %22*%00#)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    while (!%03_equals(b, %04))\n"
"    {\n"
"        struct %02*c = b;\n"
"        struct %02Division division;\n"
"        %0;(local_stack, output_stack, &division, a, b#00, %00#);\n"
"        b = division.remainder;\n"
"        a = c;\n"
"    }\n"
"    a = %03_copy(output_stack, a);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return a;\n"
"}\n";

char*extended_gcd_template =
"void %01_extended_gcd(struct Stack*output_stack, struct Stack*local_stack, struct ExtendedGCDInfo*out, struct %02*a, struct %02*b#00, struct %22*%00#)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    out->a_coefficient = %04;\n"
"    out->b_coefficient = %05;\n"
"    out->a_over_gcd = %04;\n"
"    out->b_over_gcd = %05;\n"
"    while (!%03_equals(a, %04))\n"
"    {\n"
"        struct %02Division division;\n"
"        %0;(local_stack, output_stack, &division, b, a#00, %00#);\n"
"        struct %02*m = %08(local_stack, output_stack, out->a_coefficient,\n"
"            %09(local_stack, output_stack, out->b_over_gcd, division.quotient#00, %00#)#01, %00#);\n"
"        struct %02*n = %08(local_stack, output_stack, out->b_coefficient,\n"
"            %09(local_stack, output_stack, out->a_over_gcd, division.quotient#00, %00#)#01, %00#);\n"
"        b = a;\n"
"        a = division.remainder;\n"
"        out->a_coefficient = out->b_over_gcd;\n"
"        out->b_coefficient = out->a_over_gcd;\n"
"        out->b_over_gcd = m;\n"
"        out->a_over_gcd = n;\n"
"    }\n"
"    out->a_coefficient = %03_copy(output_stack, out->a_coefficient);\n"
"    out->a_over_gcd = %03_copy(output_stack, out->a_over_gcd);\n"
"    out->b_coefficient = %03_copy(output_stack, out->b_coefficient);\n"
"    out->b_over_gcd = %07(output_stack#01, local_stack#, out->b_over_gcd#01, %00#);\n"
"    out->gcd = %03_copy(output_stack, b);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"}\n";

char*polynomial_content_template =
"struct %12*%01_content(struct Stack*output_stack, struct Stack*local_stack, struct %02*a)\n"
"{\n"
"    if (!a->coefficient_count)\n"
"    {\n"
"        return %15;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %12*out = a->coefficients[0];\n"
"    for (size_t i = 1; i < a->coefficient_count; ++i)\n"
"    {\n"
"        out = %1<(local_stack, output_stack, out, a->coefficients[i]);\n"
"    }\n"
"    out = %13_copy(output_stack, out);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*polynomial_derivative_template =
"struct %02*%01_derivative(struct Stack*output_stack, struct Stack*local_stack, struct %02*a)\n"
"{\n"
"    if (!a->coefficient_count)\n"
"    {\n"
"        return a;\n"
"    }\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %02*out = polynomial_allocate(output_stack, a->coefficient_count - 1);\n"
"    struct Integer*multiplier = &zero;\n"
"    for (size_t i = 1; i < a->coefficient_count; ++i)\n"
"    {\n"
"        multiplier = integer_add(local_stack, multiplier, &one);\n"
"        out->coefficients[i - 1] = %0=(output_stack, local_stack, a->coefficients[i], multiplier);\n"
"    }\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*polynomial_squarefree_factor_template =
"size_t %01_squarefree_factor(struct Stack*output_stack, struct Stack*local_stack, struct %02*a, struct %02**out#00, struct %12*%00#)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    size_t factor_count = 0;\n"
"    struct %02*b = a;\n"
"    struct %02*c = %03_derivative(local_stack, output_stack, a);\n"
"    a = %0<(local_stack, output_stack, b, c#00, %00#);\n"
"    struct %02Division division;\n"
"    %0;(local_stack, output_stack, &division, b, a#00, %00#);\n"
"    b = division.quotient;\n"
"    do\n"
"    {\n"
"        %0;(local_stack, output_stack, &division, c, a#00, %00#);\n"
"        c = %08(local_stack, output_stack, division.quotient, %03_derivative(local_stack, output_stack, b));\n"
"        a = %0<(output_stack, local_stack, b, c#00, %00#);\n"
"        if (a->coefficient_count > 1)\n"
"        {\n"
"            out[factor_count] = a;\n"
"            ++factor_count;\n"
"        }\n"
"        %0;(local_stack, output_stack, &division, b, a#00, %00#);\n"
"        b = division.quotient;\n"
"    } while (b->coefficient_count > 1);\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return factor_count;\n"
"}\n";

char*rational_polynomial_evaluate_template =
"struct %02*rational_polynomial_evaluate_at_%01(struct Stack*output_stack, struct Stack*local_stack, struct RationalPolynomial*a, struct %02*argument)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    struct %02*out = %04;\n"
"    struct %02*argument_power = %04;\n"
"    for (size_t i = 0; i < a->coefficient_count; ++i)\n"
"    {\n"
"        out = %06(local_stack, output_stack, out, %20(local_stack, output_stack, argument_power, a->coefficients[i]));\n"
"        argument_power = %09(local_stack, output_stack, argument_power, argument);\n"
"    }\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";