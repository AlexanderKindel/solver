#include "declarations.h"
#include "stack.c"
#include "generic.c"
#include "integer.c"
#include "rational.c"
#include "float.c"
#include "interval.c"
#include "pi.c"
#include "integer_polynomial.c"
#include "modded_polynomial.c"
#include "rational_polynomial.c"
#include "number_field_element.c"
#include "nested_polynomial.c"
#include "number_field_polynomial.c"
#include "matrix.c"
#include "algebraic_number.c"
#include "number.c"
#include "number_estimate.c"
#include "roots_of_unity.c"

struct Number*get_input(struct Stack*output_stack, struct Stack*local_stack)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    struct Number*number = out;
    while (true)
    {
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.numerator = integer_from_char(local_stack, next_char);
            number->value.denominator = &one;
            next_char = getchar();
            while (isdigit(next_char))
            {
                number->value.numerator = integer_add(local_stack, integer_multiply(local_stack,
                    output_stack, number->value.numerator, &INT(10, +)),
                    integer_from_char(local_stack, next_char));
                next_char = getchar();
            }
            number->value.numerator = integer_copy(output_stack, number->value.numerator);
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
                putc(next_char, stdout);
                puts(" is an invalid character.");
                while (getchar() != '\n')
                {}
                return 0;
            }
            number->operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Number*next_number = ALLOCATE(output_stack, struct Number);
        next_number->previous = number;
        number->next = next_number;
        number = next_number;
    }
}

bool is_numeric(struct Number*a)
{
    if (!a)
    {
        return false;
    }
    return a->operation == 'r' || a->left;
}

bool is_unparsed_operation(struct Number*a, char operation)
{
    return a->operation == operation && !a->left;
}

bool parse_binary_operation(struct Number*a, char operation)
{
    if (!is_unparsed_operation(a, operation))
    {
        return true;
    }
    if (!a->previous || !is_numeric(a->previous))
    {
        putc(operation, stdout);
        puts(" missing left operand.");
        return false;
    }
    if (a->next && is_numeric(a->next))
    {
        a->left = a->previous;
        a->right = a->next;
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
        putc(operation, stdout);
        puts(" missing right operand.");
        return false;
    }
    return true;
}

void rewind_to_first(struct Number**a)
{
    while ((*a)->previous)
    {
        *a = (*a)->previous;
    }
}

bool parse_binary_operation_pair(struct Number**number, char operation_a, char operation_b)
{
    rewind_to_first(number);
    while (true)
    {
        if (!parse_binary_operation(*number, operation_a))
        {
            return false;
        }
        else if (!parse_binary_operation(*number, operation_b))
        {
            return false;
        }
        if (!(*number)->next)
        {
            return true;
        }
        *number = (*number)->next;
    }
}

struct Number*parse_input(struct Stack*output_stack, struct Number*input)
{
    if (!input)
    {
        puts("Empty expression.");
        return 0;
    }
    while (true)
    {
        if (input->operation == ')')
        {
            puts("Unmatched ).");
            return 0;
        }
        if (input->operation == '(')
        {
            int unmatched_paren_count = 1;
            struct Number*nested_number = input;
            while (unmatched_paren_count > 0)
            {
                nested_number = nested_number->next;
                if (!nested_number)
                {
                    puts("Unmatched (.");
                    return 0;
                }
                if (nested_number->operation == '(')
                {
                    unmatched_paren_count += 1;
                }
                else if (nested_number->operation == ')')
                {
                    unmatched_paren_count -= 1;
                }
            }
            struct Number*previous = input->previous;
            struct Number*next = nested_number->next;
            input->next->previous = 0;
            nested_number->previous->next = 0;
            struct Number*parsed_nested_expression = parse_input(output_stack, input->next);
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
                input = parsed_nested_expression;
            }
            else
            {
                return 0;
            }
        }
        if (!input->next)
        {
            break;
        }
        input = input->next;
    }
    rewind_to_first(&input);
    while (true)
    {
        struct Number*next_number = input->next;
        if (is_unparsed_operation(input, '+') && !is_numeric(input->previous))
        {
            if (!next_number || (!is_numeric(next_number) &&
                !is_unparsed_operation(next_number, '+') &&
                !is_unparsed_operation(next_number, '-')))
            {
                puts("+ missing right operand.");
                return 0;
            }
            next_number->previous = input->previous;
            if (input->previous)
            {
                input->previous->next = next_number;
            }
        }
        if (!next_number)
        {
            break;
        }
        input = next_number;
    }
    rewind_to_first(&input);
    while (true)
    {
        if (!parse_binary_operation(input, '^'))
        {
            return 0;
        }
        if (!input->next)
        {
            break;
        }
        input = input->next;
    }
    rewind_to_first(&input);
    while (true)
    {
        if (is_unparsed_operation(input, '-') && !is_numeric(input->previous))
        {
            if (!input->next)
            {
                puts("- missing right operand.");
                return 0;
            }
            if (is_numeric(input->next))
            {
                input->operation = 'r';
                input->value.numerator = integer_initialize(output_stack, 1, -1);
                input->value.denominator = &one;
                struct Number*times = ALLOCATE(output_stack, struct Number);
                times->operation = '*';
                times->left = input;
                times->right = input->next;
                times->next = input->next->next;
                if (input->previous)
                {
                    input->previous->next = times;
                    times->previous = input;
                }
                input = times;
                if (input->next)
                {
                    input->next->previous = input;
                }
                else
                {
                    break;
                }
                input = input->next;
            }
            else if (is_unparsed_operation(input->next, '-'))
            {
                struct Number*new_next = input->next->next;
                if (!new_next)
                {
                    puts("- missing right operand.");
                    return 0;
                }
                new_next->previous = input->previous;
                if (input->previous)
                {
                    input->previous->next = new_next;
                }
                input = new_next;
            }
            else
            {
                puts("- missing left operand.");
                return 0;
            }
        }
        else
        {
            if (!input->next)
            {
                break;
            }
            input = input->next;
        }
    }
    if (!parse_binary_operation_pair(&input, '*', '/'))
    {
        return 0;
    }
    if (!parse_binary_operation_pair(&input, '+', '-'))
    {
        return 0;
    }
    return input;
}

void init(struct Stack*stack_a, struct Stack*stack_b)
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = max(system_info.dwAllocationGranularity, system_info.dwPageSize);
    stack_initialize(&pi_stack_a, (size_t)system_info.lpMinimumApplicationAddress,
        (size_t)system_info.lpMinimumApplicationAddress + page_size);
    stack_initialize(&pi_stack_b, pi_stack_a.end, pi_stack_a.end + page_size);
    primes =
        VirtualAlloc((void*)pi_stack_b.end, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    roots_of_unity = VirtualAlloc((void*)(pi_stack_b.end + page_size), page_size,
        MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    size_t permanent_stack_start = pi_stack_b.end + 2 * page_size;
    size_t arena_size = page_size * (((size_t)system_info.lpMaximumApplicationAddress -
        permanent_stack_start) / (3 * page_size));
    stack_initialize(&permanent_stack, permanent_stack_start, permanent_stack_start + arena_size);
    stack_initialize(stack_a, permanent_stack.end, permanent_stack.end + arena_size);
    stack_initialize(stack_b, stack_a->end, stack_a->end + arena_size);
    roots_of_unity[0] = &number_one;
    roots_of_unity[1] = &number_one;
    roots_of_unity[2] = number_rational_initialize(&permanent_stack,
        &(struct Rational){integer_initialize(&permanent_stack, 1, -1), &one});
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
    number_one.operation = 'r';
    number_one.value = rational_one;
    number_one.minimal_polynomial = polynomial_allocate(&permanent_stack, 2);
    number_one.minimal_polynomial->coefficients[0] = ALLOCATE(&permanent_stack, struct Rational);
    number_one.minimal_polynomial->coefficients[0]->numerator =
        integer_initialize(&permanent_stack, 1, -1);
    number_one.minimal_polynomial->coefficients[0]->denominator = &one;
    number_one.minimal_polynomial->coefficients[1] = &rational_one;
}

int main()
{
    struct Stack stack_a;
    struct Stack stack_b;
    init(&stack_a, &stack_b);
    while (true)
    {
        struct Number*input = get_input(&stack_a, &stack_b);
        if (input)
        {
            struct Number*number = input;
            while (number->next)
            {
                if ((number->operation == 'r' && number->next->operation == '(') ||
                    (number->operation == ')' &&
                    (number->next->operation == 'r' || number->next->operation == '(')))
                {
                    struct Number*times = ALLOCATE(&stack_a, struct Number);
                    times->operation = '*';
                    times->previous = number;
                    times->next = number->next;
                    number->next->previous = times;
                    number->next = times;
                }
                number = number->next;
            }
            input = parse_input(&stack_a, (struct Number*)input);
            if (input)
            {
                struct Number*evaluation = number_evaluate(&stack_b, &stack_a, input);
                if (evaluation != number_divide_by_zero_error)
                {
                    char*string = ARRAY_ALLOCATE(&stack_a, 2, char);
                    string[0] = '=';
                    string[1] = '\n';
                    number_string(&stack_a, &stack_b, evaluation);
                    *(char*)ALLOCATE(&stack_a, char) = 0;
                    puts(string);
                }
            }
        }
        puts("");
        stack_free(&stack_a);
        stack_free(&stack_b);
    }
    return 0;
}