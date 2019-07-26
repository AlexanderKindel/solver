#include "declarations.h"
#include "arenas.c"
#include "generic.c"
#include "integer.c"
#include "rational.c"
#include "integer_polynomial.c"
#include "modded_polynomial.c"
#include "rational_polynomial.c"
#include "number_field_element.c"
#include "nested_polynomial.c"
#include "number_field_polynomial.c"
#include "matrix.c"
#include "algebraic_number.c"
#include "float.c"
#include "number.c"
#include "number_operations.c"
#include "number_estimate.c"
#include "number_conjugates.c"
#include "roots_of_unity.c"

struct Number*get_input(struct PoolSet*transient_pool_set, struct Stack*stack_a,
    struct Stack*stack_b)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    void*stack_a_savepoint = stack_a->cursor;
    struct Number*out = number_allocate(transient_pool_set);
    struct Number*number = out;
    while (true)
    {
        if (isdigit(next_char))
        {
            number->operation = 'r';
            number->value.denominator = pool_integer_initialize(transient_pool_set, 1, 1);
            number->value.numerator = integer_from_char(stack_a, next_char);
            next_char = getchar();
            while (isdigit(next_char))
            {
                number->value.numerator =
                    integer_multiply(stack_a, stack_b, number->value.numerator, &INT(10, +));
                number->value.numerator = integer_add(stack_a, number->value.numerator,
                    integer_from_char(stack_a, next_char));
                next_char = getchar();
            }
            integer_move_to_pool(transient_pool_set, &number->value.numerator);
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
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        struct Number*next_number = number_allocate(transient_pool_set);
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

struct Number*parse_input(struct PoolSet*transient_pool_set, struct Number*input)
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
            struct Number*parsed_nested_expression = parse_input(transient_pool_set, input->next);
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
                number_parse_node_free(transient_pool_set, input);
                number_parse_node_free(transient_pool_set, nested_number);
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
            number_parse_node_free(transient_pool_set, input);
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
                input->value.numerator = pool_integer_initialize(transient_pool_set, 1, -1);
                input->value.denominator = pool_integer_initialize(transient_pool_set, 1, 1);
                struct Number*times = number_allocate(transient_pool_set);
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
                number_parse_node_free(transient_pool_set, input->next);
                number_parse_node_free(transient_pool_set, input);
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

void init(struct PoolSet*transient_pool_set, struct Stack*stack_a, struct Stack*stack_b)
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = max(system_info.dwAllocationGranularity, system_info.dwPageSize);
    size_t arena_size = page_size * (((size_t)system_info.lpMaximumApplicationAddress -
        (size_t)system_info.lpMinimumApplicationAddress) / (5 * page_size));
    primes = VirtualAlloc(system_info.lpMinimumApplicationAddress, page_size,
        MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    stack_initialize(&permanent_stack, (size_t)system_info.lpMinimumApplicationAddress + page_size,
        (size_t)system_info.lpMinimumApplicationAddress + arena_size);
    roots_of_unity = VirtualAlloc((void*)permanent_stack.end, page_size, MEM_RESERVE | MEM_COMMIT,
        PAGE_READWRITE);
    pool_memory_start = permanent_stack.end + page_size;
    pool_set_initialize(&permanent_pool_set, pool_memory_start, permanent_stack.end + arena_size);
    pool_memory_end = permanent_pool_set.slot_page_stack.end + arena_size;
    pool_set_initialize(transient_pool_set, permanent_pool_set.slot_page_stack.end,
        pool_memory_end);
    stack_initialize(stack_a, pool_memory_end, pool_memory_end + arena_size);
    stack_initialize(stack_b, stack_a->end, stack_a->end + arena_size);
    roots_of_unity[0] = pool_value_allocate(&permanent_pool_set, sizeof(struct Number));
    roots_of_unity[0]->operation = 'r';
    roots_of_unity[0]->value = rational_one;
    roots_of_unity[0]->magnitude_estimate =
        pool_value_allocate(&permanent_pool_set, sizeof(struct FloatInterval));
    roots_of_unity[0]->magnitude_estimate->min = &float_one;
    roots_of_unity[0]->magnitude_estimate->max = &float_one;
    roots_of_unity[1] = roots_of_unity[0];
    roots_of_unity[2] = pool_value_allocate(&permanent_pool_set, sizeof(struct Number));
    roots_of_unity[2]->operation = 'r';
    roots_of_unity[2]->value.numerator = pool_integer_initialize(&permanent_pool_set, 1, -1);
    roots_of_unity[2]->value.denominator = &one;
    roots_of_unity[2]->magnitude_estimate = roots_of_unity[0]->magnitude_estimate;
    primes[0] = stack_integer_initialize(&permanent_stack, 2, 1);
    primes[1] = stack_integer_initialize(&permanent_stack, 3, 1);
    pi_estimate_min.numerator = pool_integer_initialize(&permanent_pool_set, 47, 1);
    pi_estimate_min.denominator = pool_integer_initialize(&permanent_pool_set, 15, 1);
    pi_interval_size.numerator = pool_integer_initialize(&permanent_pool_set, 1696, 1);
    pi_interval_size.denominator = pool_integer_initialize(&permanent_pool_set, 12285, 1);
    pi_sixteen_to_the_k = pool_integer_initialize(&permanent_pool_set, 16, 1);
    pi_eight_k = pool_integer_initialize(&permanent_pool_set, 8, 1);
}

int main()
{
    struct PoolSet transient_pool_set;
    struct Stack stack_a;
    struct Stack stack_b;
    init(&transient_pool_set, &stack_a, &stack_b);
    while (true)
    {
        struct Number*input = get_input(&transient_pool_set, &stack_a, &stack_b);
        if (input)
        {
            struct Number*number = input;
            while (number->next)
            {
                if ((number->operation == 'r' && number->next->operation == '(') ||
                    (number->operation == ')' &&
                    (number->next->operation == 'r' || number->next->operation == '(')))
                {
                    struct Number*times = number_allocate(&transient_pool_set);
                    times->operation = '*';
                    times->previous = number;
                    times->next = number->next;
                    number->next->previous = times;
                    number->next = times;
                }
                number = number->next;
            }
            input = parse_input(&transient_pool_set, (struct Number*)input);
            if (input)
            {
                struct Number*evaluation =
                    number_evaluate(&transient_pool_set, &stack_a, &stack_b, input);
                if (evaluation != number_divide_by_zero_error)
                {
                    char*string = stack_slot_allocate(&stack_a, 2 * sizeof(char), _Alignof(char));
                    string[0] = '=';
                    string[1] = '\n';
                    number_string(&stack_a, &stack_b, evaluation);
                    *(char*)STACK_SLOT_ALLOCATE(&stack_a, char) = 0;
                    puts(string);
                }
            }
        }
        puts("");
        stack_free(&transient_pool_set.slot_page_stack);
        memset(transient_pool_set.pool_page, 0, page_size);
        stack_free(&stack_a);
        stack_free(&stack_b);
    }
    return 0;
}