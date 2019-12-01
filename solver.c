#include "declarations.h"
#include "template_instantiations.c"
#include "stack.c"
#include "integer.c"
#include "rational.c"
#include "float.c"
#include "interval.c"
#include "pi.c"
#include "integer_polynomial.c"
#include "modded_polynomial.c"
#include "rational_polynomial.c"
#include "gaussian_rational.c"
#include "nested_polynomial.c"
#include "matrix.c"
#include "number.c"
#include "number_add.c"
#include "number_multiply.c"
#include "number_estimate.c"
#include "roots_of_unity.c"

struct Token*quit_program = (struct Token*)1;

struct Token*get_input(struct Stack*output_number_stack, struct Stack*output_token_stack)
{
    char next_char = getchar();
    switch (next_char)
    {
    case '\n':
        return 0;
    case 'q':
        return quit_program;
    }
    struct Token*out = ALLOCATE(output_token_stack, struct Token);
    struct Token*token = out;
    token->previous = 0;
    token->next = 0;
    while (true)
    {
        if (isdigit(next_char))
        {
            token->is_parsed = true;
            token->value = ALLOCATE(output_number_stack, struct UnevaluatedNumber);
            token->value->operation = 'r';
            void*token_stack_savepoint = output_token_stack->cursor;
            token->value->value = char_to_integer(output_token_stack, next_char);
            next_char = getchar();
            while (isdigit(next_char))
            {
                token->value->value = integer_add(output_token_stack,
                    integer_multiply(output_token_stack, output_number_stack, token->value->value,
						INT(10, 1)),
					char_to_integer(output_token_stack, next_char));
                next_char = getchar();
            }
            token->value->value = integer_copy(output_number_stack, token->value->value);
            output_token_stack->cursor = token_stack_savepoint;
        }
        else
        {
            switch (next_char)
            {
            case '+':
            case '-':
            case '*':
            case '/':
            case '^':
            case '(':
            case ')':
                break;
            default:
                printf("'%c' is an invalid character.\n", next_char);
                while (getchar() != '\n')
                {}
                return 0;
            }
            token->is_parsed = false;
            token->value = ALLOCATE(output_number_stack, struct UnevaluatedNumber);
            token->value->left = 0;
            token->value->right = 0;
            token->value->operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            return out;
        }
        struct Token*next_token = ALLOCATE(output_token_stack, struct Token);
        next_token->previous = token;
        next_token->next = 0;
        token->next = next_token;
        token = next_token;
    }
}

bool is_unparsed_operation(struct Token*a, char operation)
{
    return !a->is_parsed && a->value->operation == operation;
}

bool parse_binary_operation(struct Token*a, char operation)
{
    if (!is_unparsed_operation(a, operation))
    {
        return true;
    }
    if (!a->previous || !a->previous->is_parsed)
    {
        printf("%c missing left operand.\n", operation);
        return false;
    }
    if (a->next && a->next->is_parsed)
    {
        a->is_parsed = true;
        a->value->left = a->previous->value;
        a->value->right = a->next->value;
        a->previous = a->previous->previous;
        if (a->previous)
        {
            a->previous->next = a;
        }
        a->next = a->next->next;
        if (a->next)
        {
            a->next->previous = a;
        }
    }
    else
    {
        printf("%c missing right operand.\n", operation);
        return false;
    }
    return true;
}

void rewind_to_first(struct Token**a)
{
    while ((*a)->previous)
    {
        *a = (*a)->previous;
    }
}

bool parse_binary_operation_pair(struct Token**a, char operation_a, char operation_b)
{
    rewind_to_first(a);
    while (true)
    {
        if (!parse_binary_operation(*a, operation_a))
        {
            return false;
        }
        else if (!parse_binary_operation(*a, operation_b))
        {
            return false;
        }
        if (!(*a)->next)
        {
            return true;
        }
        *a = (*a)->next;
    }
}

struct UnevaluatedNumber*token_parse(struct Stack*output_stack, struct Stack*local_stack,
    struct Token*a)
{
    if (!a)
    {
        puts("Empty expression.");
        return 0;
    }
    while (true)
    {
        if (a->value->operation == ')')
        {
            puts("Unmatched ).");
            return 0;
        }
        if (a->value->operation == '(')
        {
            int unmatched_paren_count = 1;
            struct Token*nested_espression = a;
            while (unmatched_paren_count > 0)
            {
                nested_espression = nested_espression->next;
                if (!nested_espression)
                {
                    puts("Unmatched (.");
                    return 0;
                }
                if (nested_espression->value->operation == '(')
                {
                    unmatched_paren_count += 1;
                }
                else if (nested_espression->value->operation == ')')
                {
                    unmatched_paren_count -= 1;
                }
            }
            struct Token*previous = a->previous;
            struct Token*next = nested_espression->next;
            a->next->previous = 0;
            nested_espression->previous->next = 0;
            struct Token*parsed_nested_expression = ALLOCATE(local_stack, struct Token);
            parsed_nested_expression->is_parsed = true;
            parsed_nested_expression->value = token_parse(output_stack, local_stack, a->next);
            if (parsed_nested_expression)
            {
                parsed_nested_expression->previous = previous;
                if (previous)
                {
                    previous->next = parsed_nested_expression;
                }
                parsed_nested_expression->next = next;
                if (next)
                {
                    next->previous = parsed_nested_expression;
                }
                a = parsed_nested_expression;
            }
            else
            {
                return 0;
            }
        }
        if (!a->next)
        {
            break;
        }
        a = a->next;
    }
    rewind_to_first(&a);
    while (true)
    {
        struct Token*next = a->next;
        if (is_unparsed_operation(a, '+') && (!a->previous || !a->previous->is_parsed))
        {
            if (!next || (!next->is_parsed && next->value->operation != '+' &&
                next->value->operation != '-'))
            {
                puts("+ missing right operand.");
                return 0;
            }
            next->previous = a->previous;
            if (a->previous)
            {
                a->previous->next = next;
            }
        }
        if (!next)
        {
            break;
        }
        a = next;
    }
    rewind_to_first(&a);
    while (true)
    {
        if (!parse_binary_operation(a, '^'))
        {
            return 0;
        }
        if (!a->next)
        {
            break;
        }
        a = a->next;
    }
    rewind_to_first(&a);
    while (true)
    {
        if (is_unparsed_operation(a, '-') && (!a->previous || !a->previous->is_parsed))
        {
            if (!a->next)
            {
                puts("- missing right operand.");
                return 0;
            }
            if (a->next->is_parsed)
            {
                a->is_parsed = true;
                a->value->operation = 'r';
                a->value->value = integer_initialize(output_stack, 1, -1);
                struct Token*times = ALLOCATE(local_stack, struct Token);
                times->is_parsed = true;
                times->value = ALLOCATE(output_stack, struct UnevaluatedNumber);
                times->value->operation = '*';
                times->value->left = a->value;
                times->value->right = a->next->value;
                times->next = a->next->next;
                times->previous = a->previous;
                if (a->previous)
                {
                    a->previous->next = times;
                }
                a = times;
                if (a->next)
                {
                    a->next->previous = a;
                }
                else
                {
                    break;
                }
                a = a->next;
            }
            else if (is_unparsed_operation(a->next, '-'))
            {
                struct Token*new_next = a->next->next;
                if (!new_next)
                {
                    puts("- missing right operand.");
                    return 0;
                }
                new_next->previous = a->previous;
                if (a->previous)
                {
                    a->previous->next = new_next;
                }
                a = new_next;
            }
            else
            {
                puts("- missing left operand.");
                return 0;
            }
        }
        else
        {
            if (!a->next)
            {
                break;
            }
            a = a->next;
        }
    }
    if (!parse_binary_operation_pair(&a, '*', '/'))
    {
        return 0;
    }
    if (!parse_binary_operation_pair(&a, '+', '-'))
    {
        return 0;
    }
    return a->value;
}

struct Number*unevaluated_number_evaluate(struct Stack*output_stack, struct Stack*local_stack,
    struct UnevaluatedNumber*a)
{
    if (a->operation == 'r')
    {
        return number_rational_initialize(output_stack, &(struct Rational){ a->value, &one });
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*left = unevaluated_number_evaluate(local_stack, output_stack, a->left);
    if (!left)
    {
        return 0;
    }
    struct Number*right = unevaluated_number_evaluate(local_stack, output_stack, a->right);
    if (!right)
    {
        return 0;
    }
    struct Number*out;
    switch (a->operation)
    {
    case '+':
        out = number_add(output_stack, local_stack, left, right);
        break;
    case '-':
    {
        out = number_add(output_stack, local_stack, left,
            number_rational_multiply(local_stack, output_stack, right,
                &(struct Rational){ integer_initialize(local_stack, 1, -1), &one }));
        break;
    }
    case '*':
        out = number_multiply(output_stack, local_stack, left, right);
        break;
    case '/':
        out = number_divide(output_stack, local_stack, left, right);
        break;
    case '^':
        if (right->operation != 'r')
        {
            puts("The input expression contains an exponentiation whose exponent is not both real "
                "and rational; this program doesn't handle transcendental numbers.");
            return 0;
        }
        out = number_rational_exponentiate(output_stack, local_stack, left, right->value);
        break;
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void init_permanent_memory()
{
    number_one.operation = 'r';
    number_one.value = &rational_one;
    number_one.minimal_polynomial = polynomial_allocate(&permanent_stack, 2);
    number_one.minimal_polynomial->coefficients[0] = ALLOCATE(&permanent_stack, struct Rational);
    number_one.minimal_polynomial->coefficients[0]->numerator =
        integer_initialize(&permanent_stack, 1, -1);
    number_one.minimal_polynomial->coefficients[0]->denominator = &one;
    number_one.minimal_polynomial->coefficients[1] = &rational_one;
    roots_of_unity[0] = &number_one;
    roots_of_unity[1] = &number_one;
    roots_of_unity[2] = number_rational_initialize(&permanent_stack,
        &(struct Rational){ integer_initialize(&permanent_stack, 1, -1), &one });
    primes[0] = integer_initialize(&permanent_stack, 2, 1);
    primes[1] = integer_initialize(&permanent_stack, 3, 1);
    pi.min = ALLOCATE(&pi_stack_a, struct Rational);
    pi.min->numerator = integer_initialize(&pi_stack_a, 47, 1);
    pi.min->denominator = integer_initialize(&pi_stack_a, 15, 1);
    pi.max = ALLOCATE(&pi_stack_a, struct Rational);
    pi.max->numerator = integer_initialize(&pi_stack_a, 40189, 1);
    pi.max->denominator = integer_initialize(&pi_stack_a, 12285, 1);
    pi_sixteen_to_the_k = integer_initialize(&pi_stack_a, 16, 1);
    pi_eight_k = integer_initialize(&pi_stack_a, 8, 1);
}

void init(struct Stack*stack_a, struct Stack*stack_b)
{
    SET_PAGE_SIZE();
    void*base_address = RESERVE_MEMORY(UINT32_MAX);
    stack_initialize(&pi_stack_a, base_address, page_size);
    stack_initialize(&pi_stack_b, pi_stack_a.end, page_size);
    primes = (struct Integer**)pi_stack_b.end;
    COMMIT_PAGE(primes);
    roots_of_unity = (struct Number**)((uint8_t*)pi_stack_b.end + page_size);
    COMMIT_PAGE(roots_of_unity);
    size_t arena_size = page_size * ((UINT32_MAX - 4 * page_size) / (3 * page_size));
    void*permanent_stack_start = (uint8_t*)pi_stack_b.end + 2 * page_size;
    stack_initialize(&permanent_stack, permanent_stack_start, arena_size);
    stack_initialize(stack_a, permanent_stack.end, arena_size);
    stack_initialize(stack_b, stack_a->end, arena_size);
    init_permanent_memory();
}

void reset_permanent_memory()
{
    stack_reset(&pi_stack_a);
    stack_reset(&pi_stack_b);
    stack_reset(&permanent_stack);
    memset(primes, 0, page_size);
    memset(roots_of_unity, 0, page_size);
    init_permanent_memory();
}

int main()
{
    struct Stack stack_a;
    struct Stack stack_b;
    init(&stack_a, &stack_b);
    while (true)
    {
        if (setjmp(memory_error_buffer))
        {
            while (getchar() != '\n')
            {}
            reset_permanent_memory();
        }
        else
        {
            struct Token*input = get_input(&stack_a, &stack_b);
            if (input == quit_program)
            {
                return 0;
            }
            if (setjmp(memory_error_buffer))
            {
                reset_permanent_memory();
            }
            else if (input)
            {
                struct Token*token = input;
                while (token->next)
                {
                    if ((token->value->operation == 'r' && token->next->value->operation == '(') ||
                        (token->value->operation == ')' &&
                        (token->next->value->operation == 'r' ||
                            token->next->value->operation == '(')))
                    {
                        struct Token*times = ALLOCATE(&stack_b, struct Token);
                        times->value = ALLOCATE(&stack_a, struct UnevaluatedNumber);
                        times->value->operation = '*';
                        times->previous = token;
                        times->next = token->next;
                        times->is_parsed = false;
                        token->next->previous = times;
                        token->next = times;
                    }
                    token = token->next;
                }
                struct UnevaluatedNumber*parsed_input = token_parse(&stack_a, &stack_b, input);
                stack_b.cursor = (void*)stack_b.start;
                if (parsed_input)
                {
                    struct Number*evaluation =
                        unevaluated_number_evaluate(&stack_b, &stack_a, parsed_input);
                    if (evaluation)
                    {
                        evaluation =
                            number_eliminate_linear_dependencies(&stack_a, &stack_b, evaluation);
                        stack_b.cursor = (void*)stack_b.start;
                        char*string = ARRAY_ALLOCATE(&stack_b, 2, char);
                        string[0] = '=';
                        string[1] = '\n';
                        number_to_string(&stack_b, &stack_a, evaluation);
                        *(char*)ALLOCATE(&stack_b, char) = 0;
                        puts(string);
                    }
                }
            }
        }
        puts("");
        stack_reset(&stack_a);
        stack_reset(&stack_b);
    }
}
