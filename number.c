#include "declarations.h"

struct Number*number_copy(struct Stack*output_stack, struct Number*a)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = a->operation;
    out->minimal_polynomial = rational_polynomial_copy(output_stack, a->minimal_polynomial);
    switch (a->operation)
    {
    case 'r':
        out->value.numerator = integer_copy(output_stack, a->value.numerator);
        out->value.denominator = integer_copy(output_stack, a->value.denominator);
        return out;
    case '^':
    case '*':
        out->left = number_copy(output_stack, a->left);
        out->right = number_copy(output_stack, a->right);
        return out;
    case '+':
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, a->term_count, struct RationalPolynomial*);
        for (size_t i = 0; i < a->term_count; ++i)
        {
            out->terms_in_terms_of_generator[i] =
                rational_polynomial_copy(output_stack, a->terms_in_terms_of_generator[i]);
        }
    }
    out->term_count = a->term_count;
    out->terms = ARRAY_ALLOCATE(output_stack, out->term_count, struct Number*);
    for (size_t i = 0; i < a->term_count; ++i)
    {
        out->terms[i] = number_copy(output_stack, a->terms[i]);
    }
    return out;
}

struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = 'r';
    out->value.numerator = integer_copy(output_stack, value->numerator);
    out->value.denominator = integer_copy(output_stack, value->denominator);
    out->minimal_polynomial = rational_minimal_polynomial(output_stack, value);
    return out;
}

//Assigns radicand and exponent directly to the left and right fields of the return value without
//copying them to output_stack.
struct Number*number_surd_initialize(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*radicand, struct Number*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '^';
    out->left = radicand;
    out->right = exponent;
    size_t surd_index = integer_to_size_t(exponent->value.denominator);
    struct RationalPolynomial*annulling_polynomial = polynomial_allocate(local_stack,
        surd_index * (radicand->minimal_polynomial->coefficient_count - 1) + 1);
    annulling_polynomial->coefficients[0] = radicand->minimal_polynomial->coefficients[0];
    for (size_t i = 1; i < radicand->minimal_polynomial->coefficient_count; ++i)
    {
        size_t annulling_polynomial_coefficient_index = i * surd_index;
        annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
            radicand->minimal_polynomial->coefficients[i];
        for (size_t j = 1; j < surd_index; ++j)
        {
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                &rational_zero;
        }
    }
    out->minimal_polynomial = number_minimal_polynomial_from_annulling_polynomial(output_stack,
        local_stack, annulling_polynomial, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_minimal_polynomial_from_annulling_polynomial(
    struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*annulling_polynomial, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial**candidates = ARRAY_ALLOCATE(local_stack,
        annulling_polynomial->coefficient_count, struct RationalPolynomial*);
    size_t candidate_count =
        rational_polynomial_factor(local_stack, output_stack, annulling_polynomial, candidates);
    if (candidate_count == 1)
    {
        struct RationalPolynomial*out = rational_polynomial_copy(output_stack, candidates[0]);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct FloatInterval zero_estimate = { &float_zero, &float_zero };
    struct RectangularEstimate a_estimate;
    number_rectangular_estimate(local_stack, output_stack, &a_estimate, a, &rational_one);
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct RectangularEstimate evaluation_estimate;
            rational_polynomial_evaluate_at_rectangular_estimate(local_stack, output_stack,
                &evaluation_estimate, candidates[i], &a_estimate);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &evaluation_estimate,
                &(struct RectangularEstimate){ &zero_estimate, &zero_estimate }))
            {
                candidate_count -= 1;
                candidates[i] = candidates[candidate_count];
                if (candidate_count == 1)
                {
                    struct RationalPolynomial*out =
                        rational_polynomial_copy(output_stack, candidates[0]);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                ++i;
            }
        }
        number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate, a);
    }
}

struct RationalPolynomial*sum_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*power = &nested_polynomial_one;
    struct NestedPolynomial d = { 2, &(struct RationalPolynomial){2, &rational_zero, &rational_one},
        &(struct RationalPolynomial){1, &(struct Rational){&INT(1, -), &one}} };
    struct NestedPolynomial*t = (struct NestedPolynomial*)&polynomial_zero;
    for (size_t i = 0; i < left_minimal_polynomial->coefficient_count; ++i)
    {
        t = nested_polynomial_add(local_stack, output_stack, t,
            nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, power,
                &(struct RationalPolynomial){ 1, left_minimal_polynomial->coefficients[i] }));
        power = nested_polynomial_multiply(local_stack, output_stack, power, &d);
    }
    struct RationalPolynomial*out =
        number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
            nested_polynomial_resultant(local_stack, output_stack, t,
                rational_polynomial_to_nested_polynomial(local_stack, right_minimal_polynomial)),
            a);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out;
    switch (a->operation)
    {
    case 'r':
    {
        if (!a->value.numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return number_divide_by_zero_error;
        }
        out = number_rational_initialize(output_stack, rational_reciprocal(local_stack, &a->value));
        break;
    }
    case '^':
    {
        out = number_divide(output_stack, local_stack,
            number_exponentiate(local_stack, output_stack, a->left,
                &(struct Rational){ integer_add(local_stack, a->right->value.denominator,
                    &INT(1, -)), a->right->value.denominator }), a->left);
        break;
    }
    case '*':
    {
        out = number_multiply(output_stack, local_stack,
            number_reciprocal(local_stack, output_stack, a->left),
            number_reciprocal(local_stack, output_stack, a->right));
        break;
    }
    case '+':
    {
        struct Number**conjugates = number_conjugates(local_stack, output_stack, a);
        out = &number_one;
        for (size_t i = 0; i < a->minimal_polynomial->coefficient_count - 1; ++i)
        {
            out = number_multiply(local_stack, output_stack, out, conjugates[i]);
        }
        struct Rational*coefficient =
            rational_reciprocal(local_stack, a->minimal_polynomial->coefficients[0]);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out = number_rational_multiply(output_stack, local_stack, out, coefficient);
        break;
    }
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal = number_reciprocal(local_stack, output_stack, divisor);
    if (reciprocal == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    return number_multiply(output_stack, local_stack, dividend, reciprocal);
}

struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return rational_copy(output_stack, &a->value);
    case '^':
        return &rational_one;
    case '*':
        return number_rational_factor(output_stack, local_stack, a->left);
    case '+':
    case 'g':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational*out = &rational_one;
        for (size_t i = 0; i < a->term_count; ++i)
        {
            struct Rational*term_factor =
                number_rational_factor(local_stack, output_stack, a->terms[i]);
            out = rational_reduced(local_stack, output_stack,
                integer_gcd(local_stack, output_stack, out->numerator, term_factor->numerator),
                integer_lcm(local_stack, output_stack, out->denominator, term_factor->denominator));
        }
        out = rational_copy(output_stack, out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Rational*exponent)
{
    if (!exponent->numerator->value_count)
    {
        return &number_one;
    }
    if (exponent->numerator->sign < 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        base = number_reciprocal(local_stack, output_stack, base);
        if (base == number_divide_by_zero_error)
        {
            local_stack->cursor = local_stack_savepoint;
            return number_divide_by_zero_error;
        }
        struct Number*out = number_exponentiate(output_stack, local_stack, base,
            rational_negative(local_stack, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        if (base->value.numerator->value_count == 0)
        {
            return number_copy(output_stack, base);
        }
        void*local_stack_savepoint = local_stack->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*out = number_rational_initialize(output_stack,
                rational_exponentiate(local_stack, output_stack, &base->value,
                    exponent->numerator));
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*out = number_rational_multiply(output_stack, local_stack,
                number_exponentiate(local_stack, output_stack,
                    number_rational_initialize(local_stack,
                        &(struct Rational){integer_multiply(local_stack, output_stack,
                            base->value.numerator,
                            integer_exponentiate(local_stack, output_stack, base->value.denominator,
                                integer_add(local_stack, exponent->denominator, &INT(1, -)))),
                        &one}), exponent), &(struct Rational){&one, base->value.denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(exponent->denominator, &one))
        {
            struct Integer*radicand = integer_magnitude(local_stack, base->value.numerator);
            struct Factor*factors;
            size_t factor_count = integer_factor(local_stack, output_stack, &factors, radicand);
            struct Integer*coefficient = &one;
            struct Integer*multiplicity_gcd = exponent->denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct IntegerDivision multiplicity_reduction;
                integer_euclidean_divide(local_stack, output_stack, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(local_stack, output_stack, coefficient,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = integer_gcd(local_stack, output_stack, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            struct Integer*new_degree;
            if (base->value.numerator->sign > 0)
            {
                new_degree = integer_euclidean_quotient(local_stack, output_stack,
                    exponent->denominator, multiplicity_gcd);
                if (integer_equals(new_degree, &one))
                {
                    struct Rational*rational_radicand = ALLOCATE(output_stack, struct Rational);
                    rational_radicand->numerator = integer_copy(output_stack, radicand);
                    rational_radicand->denominator = &one;
                    struct Number*out = number_rational_initialize(output_stack, rational_radicand);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                new_degree = exponent->denominator;
            }
            radicand = integer_initialize(local_stack, 1, base->value.numerator->sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                    output_stack, factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(local_stack, output_stack,
                    factors[factor_index].value, reduced_multiplicity);
                radicand = integer_multiply(local_stack, output_stack, radicand, exponentiation);
            }
            struct Number*out = number_rational_multiply(output_stack, local_stack, 
                number_surd_initialize(local_stack, output_stack,
                    number_rational_initialize(local_stack, &(struct Rational){radicand, &one}),
                    number_rational_initialize(local_stack, &(struct Rational){&one, new_degree})),
                &(struct Rational){coefficient, &one});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        return base;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_exponentiate(output_stack, local_stack, base->left,
            rational_multiply(local_stack, output_stack, exponent, &base->right->value));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_multiply(output_stack, local_stack,
            number_exponentiate(local_stack, output_stack, base->left, exponent),
            number_exponentiate(local_stack, output_stack, base->right, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*out = number_exponentiate(output_stack, local_stack,
                generic_exponentiate(&(struct RingOperations){number_copy, 0, 0, &number_one, 0,
                    0, number_generic_multiply},
                    local_stack, output_stack, base, exponent->numerator, 0),
                &(struct Rational){&one, exponent->denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Rational*base_rational_factor =
            number_rational_factor(local_stack, output_stack, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(local_stack, output_stack, &factors, base_rational_factor->numerator);
        struct Rational base_cancelling_rational_factor =
            { integer_exponentiate(local_stack, output_stack, base_rational_factor->denominator,
                exponent->denominator), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor->denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                output_stack, factors[factor_index].multiplicity, exponent->denominator);
            base_cancelling_rational_factor.denominator =
                integer_multiply(local_stack, output_stack,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        integer_multiply(local_stack, output_stack, reduced_multiplicity,
                            exponent->denominator)), base_cancelling_rational_factor.denominator);
            product_rational_factor.numerator =
                integer_multiply(local_stack, output_stack, product_rational_factor.numerator,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        reduced_multiplicity));
        }
        struct Number*out = number_rational_multiply(output_stack, local_stack,
            number_surd_initialize(local_stack, output_stack,
                number_rational_multiply(local_stack, output_stack, base,
                    &base_cancelling_rational_factor), 
                number_rational_initialize(local_stack, exponent)), &product_rational_factor);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number**number_sum_or_product_conjugates(struct Number*(operation)(struct Stack*,
    struct Stack*, struct Number*, struct Number*), struct Stack*output_stack,
    struct Stack*local_stack, struct Number**a_conjugates, size_t a_conjugate_count,
    struct Number**b_conjugates, size_t b_conjugate_count,
    struct RationalPolynomial*minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number**out =
        ARRAY_ALLOCATE(output_stack, minimal_polynomial->coefficient_count - 1, struct Number*);
    size_t conjugate_count = 0;
    for (size_t i = 0; i < a_conjugate_count; ++i)
    {
        for (size_t j = 0; j < b_conjugate_count; ++j)
        {
            struct Number*candidate =
                operation(local_stack, output_stack, a_conjugates[i], b_conjugates[j]);
            if (rational_polynomial_equals(minimal_polynomial, candidate->minimal_polynomial))
            {
                out[conjugate_count] = number_copy(output_stack, candidate);
                ++conjugate_count;
                if (conjugate_count == minimal_polynomial->coefficient_count - 1)
                {
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
        }
    }
    crash("Not enough conjugates found.");
}

struct Number**number_conjugates(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        struct Number**out = ALLOCATE(output_stack, struct Number*);
        out[0] = number_copy(output_stack, a);
        return out;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**radicand_conjugates = number_conjugates(local_stack, output_stack, a->left);
        for (size_t i = 0; i < a->left->minimal_polynomial->coefficient_count - 1; ++i)
        {
            radicand_conjugates[i] = number_exponentiate(local_stack, output_stack,
                radicand_conjugates[i], &a->right->value);
        }
        size_t roots_of_unity_count = integer_to_size_t(a->right->value.denominator);
        struct Number**roots_of_unity =
            get_roots_of_unity(output_stack, local_stack, a->right->value.denominator);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, radicand_conjugates, a->left->minimal_polynomial->coefficient_count - 1,
            roots_of_unity, roots_of_unity_count, a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**left_conjugates = number_conjugates(local_stack, output_stack, a->left);
        struct Number**right_conjugates = number_conjugates(local_stack, output_stack, a->right);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, left_conjugates, a->left->minimal_polynomial->coefficient_count - 1,
            right_conjugates, a->right->minimal_polynomial->coefficient_count - 1,
            a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**out = number_conjugates(local_stack, output_stack, a->terms[0]);
        struct RationalPolynomial*left_minimal_polynomial = a->terms[0]->minimal_polynomial;
        size_t original_term_count = a->term_count;
        a->term_count = 2;
        while (a->term_count < original_term_count)
        {
            struct Number*right_term = a->terms[a->term_count - 1];
            struct RationalPolynomial*new_left_minimal_polynomial =
                sum_minimal_polynomial(local_stack, output_stack, a, left_minimal_polynomial,
                    right_term->minimal_polynomial);
            out = number_sum_or_product_conjugates(number_add, local_stack, output_stack, out,
                left_minimal_polynomial->coefficient_count - 1,
                number_conjugates(local_stack, output_stack, right_term),
                right_term->minimal_polynomial->coefficient_count - 1, new_left_minimal_polynomial);
            left_minimal_polynomial = new_left_minimal_polynomial;
            ++a->term_count;
        }
        out = number_sum_or_product_conjugates(number_add, output_stack, local_stack, out,
            left_minimal_polynomial->coefficient_count - 1,
            number_conjugates(local_stack, output_stack, a->terms[a->term_count - 1]),
            a->terms[a->term_count - 1]->minimal_polynomial->coefficient_count - 1,
            a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_evaluate(struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    if (a->operation == 'r')
    {
        a->minimal_polynomial = rational_minimal_polynomial(output_stack, &a->value);
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    a->left = number_evaluate(local_stack, output_stack, a->left);
    if (a->left == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    a->right = number_evaluate(local_stack, output_stack, a->right);
    if (a->right == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    struct Number*out;
    switch (a->operation)
    {
    case '+':
        out = number_add(output_stack, local_stack, a->left, a->right);
        break;
    case '-':
    {
        out = number_add(output_stack, local_stack, a->left,
            number_rational_multiply(local_stack, output_stack, a->right,
                &(struct Rational){ integer_initialize(local_stack, 1, -1), &one }));
        break;
    }
    case '*':
        out = number_multiply(output_stack, local_stack, a->left, a->right);
        break;
    case '/':
        out = number_divide(output_stack, local_stack, a->left, a->right);
        break;
    case '^':
        if (a->right->operation != 'r')
        {
            puts("The input expression contains an exponentiation whose exponent is not both real "
                "and rational; this program doesn't handle transcendental numbers.");
            return 0;
        }
        out = number_exponentiate(output_stack, local_stack, a->left, &a->right->value);
        break;
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t number_string(struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        size_t char_count = integer_string(output_stack, local_stack, a->value.numerator);
        if (!integer_equals(a->value.denominator, &one))
        {
            *(char*)ALLOCATE(output_stack, char) = '/';
            char_count += 1 + integer_string(output_stack, local_stack, a->value.denominator);
        }
        return char_count;
    }
    case '^':
    {
        size_t char_count;
        if (a->left->operation != 'r' || a->left->value.numerator->sign < 0)
        {
            *(char*)ALLOCATE(output_stack, char) = '(';
            char_count = 2 + number_string(output_stack, local_stack, a->left);
            *(char*)ALLOCATE(output_stack, char) = ')';
        }
        else
        {
            char_count = number_string(output_stack, local_stack, a->left);
        }
        char*exponent = ARRAY_ALLOCATE(output_stack, 2, char);
        exponent[0] = '^';
        exponent[1] = '(';
        char_count += 3 + number_string(output_stack, local_stack, a->right);
        *(char*)ALLOCATE(output_stack, char) = ')';
        return char_count;
    }
    case '*':
    {
        size_t char_count;
        if (a->left->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            if (integer_equals(a->left->value.numerator,
                integer_initialize(local_stack, 1, -1)))
            {
                *(char*)ALLOCATE(output_stack, char) = '-';
                char_count = 1 + number_string(output_stack, local_stack, a->right);
            }
            else if (!integer_equals(a->left->value.numerator, &one))
            {
                char_count =
                    integer_string(output_stack, local_stack, a->left->value.numerator);
                char*right_string = output_stack->cursor;
                size_t right_char_count = number_string(output_stack, local_stack, a->right);
                if (right_string[0] != '(')
                {
                    ALLOCATE(output_stack, char);
                    memcpy(right_string + 1, right_string, right_char_count);
                    *right_string = '*';
                    char_count += 1 + right_char_count;
                }
            }
            else
            {
                char_count = number_string(output_stack, local_stack, a->right);
            }
            if (!integer_equals(a->left->value.denominator, &one))
            {
                *(char*)ALLOCATE(output_stack, char) = '/';
                char_count +=
                    1 + integer_string(output_stack, local_stack, a->left->value.denominator);
            }
            local_stack->cursor = local_stack_savepoint;
        }
        else
        {
            char_count = number_string(output_stack, local_stack, a->left) + 1;
            *(char*)ALLOCATE(output_stack, char) = '*';
            char_count += number_string(output_stack, local_stack, a->right);
        }
        return char_count;
    }
    case '+':
    {
        size_t char_count = number_string(output_stack, local_stack, a->terms[0]);
        for (size_t i = 1; i < a->term_count; ++i)
        {
            *(char*)ALLOCATE(output_stack, char) = '+';
            char_count += number_string(output_stack, local_stack, a->terms[i]) + 1;
        }
        return char_count;
    }
    default:
        crash("Number operation not recognized.");
    }
}