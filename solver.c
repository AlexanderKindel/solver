#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

DWORD page_size;

size_t next_page_end(size_t address)
{
    size_t distance_from_page_start = address % page_size;
    if (distance_from_page_start)
    {
        return address + page_size - distance_from_page_start;
    }
    return address;
}

void advance_arena_cursor(void**cursor, size_t step_size, size_t step_count)
{
    size_t page_end = next_page_end(*cursor);
    *cursor = (size_t)*cursor + step_size * step_count;
    if (*cursor > page_end)
    {
        VirtualAlloc(page_end, (size_t)*cursor - page_end, MEM_COMMIT, PAGE_READWRITE);
    }
}

void rewind_arena_cursor(void**cursor, void*cursor_target)
{
    size_t page_end = next_page_end(cursor_target);
    if (*cursor > page_end)
    {
        VirtualFree(page_end, (size_t)*cursor - page_end, MEM_DECOMMIT);
    }
    *cursor = cursor_target;
}

struct Integer
{
    size_t value_count;
    uint32_t*values;
    int8_t sign;
};

struct Integer int_new(uint32_t**int_arena_cursor, int32_t value)
{
    if (value == 0)
    {
        struct Integer zero = { 0, 0, 0 };
        return zero;
    }
    struct Integer integer;
    integer.value_count = 1;
    integer.values = *int_arena_cursor;
    advance_arena_cursor(int_arena_cursor, sizeof(uint32_t), 1);
    if (value > 0)
    {
        integer.values[0] = value;
        integer.sign = 1;
    }
    else
    {
        integer.values[0] = -value;
        integer.sign = -1;
    }
    return integer;
}

struct Integer int_copy(uint32_t**int_arena_cursor, struct Integer*a)
{
    struct Integer copy = { a->value_count, *int_arena_cursor, a->sign };
    advance_arena_cursor(int_arena_cursor, sizeof(uint32_t), a->value_count);
    memcpy(copy.values, a->values, a->value_count * sizeof(uint32_t));
    return copy;
}

void int_move_values_to_cursor(uint32_t**int_arena_cursor, struct Integer*a)
{
    memcpy(*int_arena_cursor, a->values, a->value_count * sizeof(uint32_t));
    a->values = *int_arena_cursor;
    *int_arena_cursor += a->value_count;
}

void calculate_sum_values(struct Integer*short_integer, struct Integer*sum)
{
    for (int i = 0; i < short_integer->value_count; ++i)
    {
        uint64_t remainder = short_integer->values[i];
        for (int j = i; j < sum->value_count; ++j)
        {
            uint64_t sum_value = sum->values[j] + remainder;
            sum->values[j] = (uint32_t)sum_value;
            remainder = (sum_value & 0xffffffff00000000) >> 32;
            if (remainder == 0)
            {
                break;
            }
        }
    }
}

void twos_complement(struct Integer*integer)
{
    for (int i = 0; i < integer->value_count; ++i)
    {
        integer->values[i] = ~integer->values[i];
    }
    for (int i = 0; i < integer->value_count; ++i)
    {
        uint32_t power = 1;
        for (int j = 0; j < 32; ++j)
        {
            integer->values[i] ^= power;
            if ((integer->values[i] & power) != 0)
            {
                return;
            }
            power = power << 1;
        }
    }
}

void trim_leading_zeroes(struct Integer*a)
{
    for (int i = a->value_count - 1; i >= 0; --i)
    {
        if (a->values[i] != 0)
        {
            a->value_count = i + 1;
            return;
        }
    }
    a->value_count = 0;
    a->sign = 0;
}

struct Integer int_add(uint32_t**int_arena_cursor, struct Integer*a, struct Integer*b)
{
    if (a->sign == 0)
    {
        return int_copy(int_arena_cursor, b);
    }
    if (b->sign == 0)
    {
        return int_copy(int_arena_cursor, a);
    }
    struct Integer*long_integer = a;
    struct Integer*short_integer = b;
    if (a->value_count < b->value_count)
    {
        long_integer = b;
        short_integer = a;
    }
    struct Integer sum;
    sum.sign = short_integer->sign;
    sum.value_count = long_integer->value_count + 1;
    sum.values = *int_arena_cursor;
    advance_arena_cursor(int_arena_cursor, sizeof(uint32_t), sum.value_count);
    memcpy(sum.values, long_integer->values, long_integer->value_count * sizeof(uint32_t));
    sum.values[long_integer->value_count] = 0;
    if (short_integer->sign == long_integer->sign)
    {
        calculate_sum_values(short_integer, &sum);
    }
    else
    {
        twos_complement(&sum);
        calculate_sum_values(short_integer, &sum);
        if (sum.values[long_integer->value_count] != 0)
        {
            twos_complement(&sum);
            sum.sign *= -1;
        }
    }
    trim_leading_zeroes(&sum);
    rewind_arena_cursor(int_arena_cursor, sum.values + sum.value_count);
    return sum;
}

void int_add_to_a_in_place(struct Integer*a, struct Integer*b)
{
    if (b->sign == 0)
    {
        return;
    }
    if (a->sign == 0)
    {
        a->value_count = b->value_count;
        a->sign = b->sign;
        memcpy(a->values, b->values, b->value_count * sizeof(uint32_t));
        return;
    }
    if (a->value_count < b->value_count)
    {
        memset(a->values + a->value_count, 0, (b->value_count - a->value_count) * sizeof(uint32_t));
        a->value_count = b->value_count + 1;
    }
    else
    {
        a->value_count += 1;
    }
    if (a->sign == b->sign)
    {
        calculate_sum_values(b, a);
    }
    else
    {
        twos_complement(a);
        calculate_sum_values(b, a);
        if (a->values[a->value_count - 1] != 0)
        {
            twos_complement(a);
            a->sign *= -1;
        }
    }
    int i = a->value_count - 1;
    while (true)
    {
        if (a->values[i] != 0)
        {
            break;
        }
        if (i == 0)
        {
            a->sign = 0;
            break;
        }
        --i;
    }
    a->value_count = i + 1;
}

struct Integer int_negative(uint32_t**int_arena_cursor, struct Integer*a)
{
    struct Integer copy = int_copy(int_arena_cursor, a);
    copy.sign = -copy.sign;
    return copy;
}

struct Integer int_subtract(uint32_t**int_arena_cursor, struct Integer*minuend,
    struct Integer*subtrahend)
{
    uint32_t*local_int_arena_cursor = *int_arena_cursor;
    size_t ints_to_advance = 1;
    ints_to_advance += max(minuend->value_count, subtrahend->value_count);
    advance_arena_cursor(&local_int_arena_cursor, sizeof(uint32_t), ints_to_advance);
    struct Integer negative = int_negative(&local_int_arena_cursor, subtrahend);
    struct Integer difference = int_add(int_arena_cursor, minuend, &negative);
    rewind_arena_cursor(&local_int_arena_cursor, *int_arena_cursor);
    return difference;
}

struct Integer int_multiply(uint32_t**int_arena_cursor, struct Integer*a, struct Integer*b)
{
    struct Integer product = { 0, *int_arena_cursor, 0 };
    uint32_t*local_int_arena_cursor = *int_arena_cursor;
    advance_arena_cursor(&local_int_arena_cursor, sizeof(uint32_t),
        a->value_count + b->value_count);
    for (int i = 0; i < a->value_count; ++i)
    {
        for (int j = 0; j < b->value_count; ++j)
        {
            uint64_t product_component = (uint64_t)a->values[i] * b->values[j];
            size_t shift = i + j;
            struct Integer integer_component;
            integer_component.values = local_int_arena_cursor;
            integer_component.sign = 0;
            memset(local_int_arena_cursor, 0, shift * sizeof(uint32_t));
            if (product_component > 0)
            {
                integer_component.values[shift] = product_component;
                integer_component.value_count = shift + 1;
                integer_component.sign = 1;
            }
            uint32_t high_bytes = (product_component & 0xffffffff00000000) >> 32;
            if (high_bytes > 0)
            {
                integer_component.values[shift + 1] = high_bytes;
                integer_component.value_count = shift + 2;
                integer_component.sign = 1;
            }
            int_add_to_a_in_place(&product, &integer_component);
        }
    }
    product.sign = a->sign * b->sign;
    *int_arena_cursor += product.value_count;
    rewind_arena_cursor(&local_int_arena_cursor, int_arena_cursor);
    return product;
}

struct Division
{
    struct Integer quotient;
    struct Integer remainder;
};

int leading_digit_place(struct Integer*a)
{
    uint32_t divisor_leading_digit = 0x80000000;
    uint32_t last_value = a->values[a->value_count - 1];
    for (int i = 31; i >= 0; --i)
    {
        if ((last_value & divisor_leading_digit) != 0)
        {
            return i;
        }
        divisor_leading_digit = divisor_leading_digit >> 1;
    }
    return 0;
}

void calculate_division_values(uint32_t*int_arena_cursor, struct Integer*positive_divisor,
    struct Division*division, size_t quotient_value_index, uint32_t quotient_digit)
{
    memset(division->quotient.values, 0, division->quotient.value_count * sizeof(uint32_t));
    while (true)
    {
        for (int i = 32; i > 0; --i)
        {
            struct Integer difference =
                int_subtract(&int_arena_cursor, &division->remainder, positive_divisor);
            if (difference.sign >= 0)
            {
                division->quotient.values[quotient_value_index] |= quotient_digit;
                division->remainder = difference;
            }
            if (quotient_digit == 1)
            {
                if (quotient_value_index == 0)
                {
                    return;
                }
                quotient_digit = 0x80000000;
                quotient_value_index -= 1;
            }
            else
            {
                quotient_digit = quotient_digit >> 1;
            }
            for (int j = 0; j < positive_divisor->value_count - 1; ++j)
            {
                positive_divisor->values[j] = positive_divisor->values[j] >> 1 |
                    positive_divisor->values[j + 1] << 31;
            }
            positive_divisor->values[positive_divisor->value_count - 1] =
                positive_divisor->values[positive_divisor->value_count - 1] >> 1;
        }
    }
}

struct Division int_euclidean_divide(uint32_t**int_arena_cursor, struct Integer*dividend,
    struct Integer*divisor)
{
    int8_t quotient_sign = dividend->sign * divisor->sign;
    int dividend_leading_digit_place = leading_digit_place(dividend);
    int divisor_leading_digit_place = leading_digit_place(divisor);
    if (dividend->value_count > divisor->value_count ||
        (dividend->value_count == divisor->value_count &&
        dividend_leading_digit_place >= divisor_leading_digit_place))
    {
        uint32_t*local_int_arena_cursor = *int_arena_cursor;
        advance_arena_cursor(&local_int_arena_cursor, sizeof(uint32_t), dividend->value_count);
        struct Integer positive_divisor = { dividend->value_count, local_int_arena_cursor, 1 };
        advance_arena_cursor(&local_int_arena_cursor, sizeof(uint32_t), dividend->value_count);
        size_t quotient_value_index = dividend->value_count - divisor->value_count;
        memset(positive_divisor.values, 0, quotient_value_index * sizeof(uint32_t));
        memcpy(positive_divisor.values + quotient_value_index, divisor->values,
            divisor->value_count * sizeof(uint32_t));
        int shift = dividend_leading_digit_place - divisor_leading_digit_place;
        uint32_t quotient_digit;
        if (shift >= 0)
        {
            positive_divisor.values[quotient_value_index] =
                positive_divisor.values[quotient_value_index] << shift;
            for (int i = quotient_value_index + 1; i < positive_divisor.value_count; ++i)
            {
                positive_divisor.values[i] = positive_divisor.values[i] << shift |
                    positive_divisor.values[i - 1] >> 32 - shift;
            }
            quotient_digit = 1 << shift;
        }
        else
        {
            shift *= -1;
            for (int i = quotient_value_index; i < positive_divisor.value_count; ++i)
            {
                positive_divisor.values[i - 1] = positive_divisor.values[i] << 32 - shift |
                    positive_divisor.values[i - 1] >> shift;
            }
            positive_divisor.values[positive_divisor.value_count - 1] =
                positive_divisor.values[positive_divisor.value_count - 1] >> shift;
            quotient_digit = 1 << 32 - shift;
            quotient_value_index -= 1;
        }
        int8_t dividend_sign = dividend->sign;
        dividend->sign = 1;
        struct Division division =
            { {dividend->value_count, *int_arena_cursor, quotient_sign}, *dividend };
        calculate_division_values(local_int_arena_cursor, &positive_divisor, &division,
            quotient_value_index, quotient_digit);
        trim_leading_zeroes(&division.quotient);
        *int_arena_cursor += division.quotient.value_count;
        int_move_values_to_cursor(int_arena_cursor, &division.remainder);
        rewind_arena_cursor(&local_int_arena_cursor, *int_arena_cursor);
        if (division.remainder.sign != 0)
        {
            division.remainder.sign = dividend->sign;
        }
        dividend->sign = dividend_sign;
        return division;
    }
    struct Division division =
        { {0, 0, 0}, { dividend->value_count, *int_arena_cursor, quotient_sign } };
    advance_arena_cursor(int_arena_cursor, sizeof(uint32_t), dividend->value_count);
    memcpy(division.remainder.values, dividend->values,
        dividend->value_count * sizeof(uint32_t));
    return division;
}

void push_char(char**char_arena_cursor, char c)
{
    char*destination = *char_arena_cursor;
    advance_arena_cursor(char_arena_cursor, sizeof(char), 1);
    *destination = c;
}

void int_string(char**char_arena_cursor, uint32_t*int_arena_cursor, struct Integer*a)
{
    if (a->sign == 0)
    {
        push_char(char_arena_cursor, '0');
        return;
    }
    if (a->sign < 0)
    {
        push_char(char_arena_cursor, '-');
    }
    char*buffer_start = *char_arena_cursor;
    advance_arena_cursor(char_arena_cursor, sizeof(char), 10 * a->value_count + 1);
    char*next_char = *char_arena_cursor;
    struct Integer quotient = *a;
    struct Integer power = int_new(&int_arena_cursor, 10);
    while (quotient.sign != 0)
    {
        struct Division division = int_euclidean_divide(&int_arena_cursor, &quotient, &power);
        quotient = division.quotient;
        next_char -= 1;
        if (division.remainder.sign != 0)
        {
            *next_char = division.remainder.values[0] + '0';
        }
        else
        {
            *next_char = '0';
        }
    }
    size_t char_count = *char_arena_cursor - next_char;
    memcpy(buffer_start, next_char, char_count);
    rewind_arena_cursor(char_arena_cursor, buffer_start + char_count);
}

struct Number
{
    union
    {
        struct Integer value;
        struct Node
        {
            struct Number*left;
            struct Number*right;
        };
    };
    char operation;
};

struct NumberSlot
{
    struct Number number;
    struct NumberSlot*previous;
    struct NumberSlot*next;
};

struct NumberSlot*number_slot_new(struct NumberSlot**arena_cursor, struct NumberSlot**free_list)
{
    struct NumberSlot*allocation = *free_list;
    if (allocation)
    {
        *free_list = allocation->next;
    }
    else
    {
        allocation = *arena_cursor;
        advance_arena_cursor(arena_cursor, sizeof(struct NumberSlot), 1);
        allocation->next = 0;
    }
    return allocation;
}

void free_number_slot(struct NumberSlot**number_slot_free_list, struct NumberSlot*slot)
{
    slot->next = *number_slot_free_list;
    *number_slot_free_list = slot;
}

bool get_input(struct NumberSlot**number_slot_arena_cursor, struct NumberSlot**number_free_list,
    uint32_t**int_arena_cursor)
{
    char next_char = getchar();
    if (next_char == '\n')
    {
        return false;
    }
    struct NumberSlot*previous_slot = 0;
    while (true)
    {
        struct NumberSlot*slot = number_slot_new(number_slot_arena_cursor, number_free_list);
        slot->previous = previous_slot;
        slot->next = *number_slot_arena_cursor;
        if (isdigit(next_char))
        {
            uint32_t*int_arena_savepoint = *int_arena_cursor;
            slot->number.operation = 'i';
            slot->number.value = int_new(int_arena_cursor, next_char - '0');
            struct Integer ten = int_new(int_arena_cursor, 10);
            next_char = getchar();
            while (isdigit(next_char))
            {
                slot->number.value = int_multiply(int_arena_cursor, &slot->number.value, &ten);
                struct Integer digit = int_new(int_arena_cursor, next_char - '0');
                slot->number.value = int_add(int_arena_cursor, &slot->number.value, &digit);
                next_char = getchar();
            }
            int_move_values_to_cursor(&int_arena_savepoint, &slot->number.value);
            rewind_arena_cursor(int_arena_cursor, int_arena_savepoint);
        }
        else
        {
            switch (next_char)
            {
            case '+':
                break;
            case '-':
                break;
            case '*':
                break;
            case '/':
                break;
            case '^':
                break;
            case '(':
                break;
            case ')':
                break;
            default:
                printf("\"%c\"%s", next_char, " is an invalid character.");
                while (getchar() != '\n')
                {}
                return false;
            }
            slot->number.operation = next_char;
            next_char = getchar();
        }
        if (next_char == '\n')
        {
            slot->next = 0;
            return true;
        }
        previous_slot = slot;
    }
}

bool is_numeric(struct NumberSlot*slot)
{
    if (!slot)
    {
        return false;
    }
    return slot->number.operation == 'i' || slot->number.left;
}

bool is_unparsed_operation(struct NumberSlot*slot, char operation)
{
    return slot->number.operation == operation && !slot->number.left;
}

bool parse_binary_operation(struct NumberSlot*slot, char operation)
{
    if (!is_unparsed_operation(slot, operation))
    {
        return true;
    }
    if (!slot->previous || !is_numeric(slot->previous))
    {
        printf("%c missing left operand.", operation);
        return false;
    }
    if (slot->next && is_numeric(slot->next))
    {
        slot->number.left = &slot->previous->number;
        slot->number.right = &slot->next->number;
        slot->previous = slot->previous->previous;
        if (slot->previous)
        {
            slot->previous->next = slot;
        }
        slot->next = slot->next->next;
        if (slot->next)
        {
            slot->next->previous = slot;
        }
    }
    else
    {
        printf("%c missing right operand.", operation);
        return false;
    }
    return true;
}

void rewind_to_first_slot(struct NumberSlot **slot)
{
    while ((*slot)->previous)
    {
        *slot = (*slot)->previous;
    }
}

bool parse_binary_operation_pair(struct NumberSlot**slot, char operation_a, char operation_b)
{
    rewind_to_first_slot(slot);
    while (true)
    {
        if (is_unparsed_operation(*slot, operation_a) &&
            !parse_binary_operation(*slot, operation_a))
        {
            return false;
        }
        else if (is_unparsed_operation(*slot, operation_b) &&
            !parse_binary_operation(*slot, operation_b))
        {
            return false;
        }
        if (!(*slot)->next)
        {
            return true;
        }
        *slot = (*slot)->next;
    }
}

struct NumberSlot*parse_input(struct NumberSlot**number_slot_arena_cursor,
    struct NumberSlot**number_slot_free_list, uint32_t**int_arena_cursor,
    struct NumberSlot*number_slot)
{
    if (!number_slot)
    {
        printf("Empty expression.");
        return 0;
    }
    while (true)
    {
        if (number_slot->number.operation == ')')
        {
            printf("Unmatched ).");
            return 0;
        }
        if (number_slot->number.operation == '(')
        {
            int unmatched_paren_count = 1;
            struct NumberSlot*nested_number_slot = number_slot;
            while (unmatched_paren_count > 0)
            {
                nested_number_slot = nested_number_slot->next;
                if (!nested_number_slot)
                {
                    printf("Unmatched (.");
                    return 0;
                }
                if (nested_number_slot->number.operation == '(')
                {
                    unmatched_paren_count += 1;
                }
                else if (nested_number_slot->number.operation == ')')
                {
                    unmatched_paren_count -= 1;
                }
            }
            struct NumberSlot*previous = number_slot->previous;
            struct NumberSlot*next = nested_number_slot->next;
            number_slot->next->previous = 0;
            nested_number_slot->previous->next = 0;
            struct NumberSlot*nested_expression = parse_input(number_slot_arena_cursor,
                number_slot_free_list, int_arena_cursor, number_slot->next);
            if (nested_expression)
            {
                nested_expression->previous = previous;
                if (previous)
                {
                    previous->next = nested_expression;
                }
                nested_expression->next = next;
                if (next)
                {
                    next->previous = nested_expression;
                }
                free_number_slot(number_slot_free_list, number_slot);
                free_number_slot(number_slot_free_list, nested_number_slot);
                number_slot = nested_expression;
            }
            else
            {
                return 0;
            }
        }
        if (!number_slot->next)
        {
            break;
        }
        number_slot = number_slot->next;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        struct NumberSlot*next_slot = number_slot->next;
        if (is_unparsed_operation(number_slot, '+') && !is_numeric(number_slot->previous))
        {
            if (!next_slot || (!is_numeric(next_slot) && !is_unparsed_operation(next_slot, '+') &&
                !is_unparsed_operation(next_slot, '-')))
            {
                printf("+ missing right operand.");
                return 0;
            }
            next_slot->previous = number_slot->previous;
            if (number_slot->previous)
            {
                number_slot->previous->next = next_slot;
            }
            free_number_slot(number_slot_free_list, number_slot);
        }
        if (!next_slot)
        {
            break;
        }
        number_slot = next_slot;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        if (!parse_binary_operation(number_slot, '^'))
        {
            return 0;
        }
        if (!number_slot->next)
        {
            break;
        }
        number_slot = number_slot->next;
    }
    rewind_to_first_slot(&number_slot);
    while (true)
    {
        if (is_unparsed_operation(number_slot, '-') && !is_numeric(number_slot->previous))
        {
            if (!number_slot->next)
            {
                printf("- missing right operand.");
                return 0;
            }
            if (is_numeric(number_slot->next))
            {
                number_slot->number.operation = 'i';
                number_slot->number.value = int_new(int_arena_cursor, -1);
                struct NumberSlot*times =
                    number_slot_new(number_slot_arena_cursor, number_slot_free_list);
                times->number.operation = '*';
                times->number.left = &number_slot->number;
                times->number.right = &number_slot->next->number;
                times->next = number_slot->next->next;
                if (number_slot->previous)
                {
                    number_slot->previous->next = times;
                }
                number_slot = times;
                if (number_slot->next)
                {
                    number_slot->next->previous = number_slot;
                }
                else
                {
                    break;
                }
                number_slot = number_slot->next;
            }
            else if (is_unparsed_operation(number_slot->next, '-'))
            {
                struct NumberSlot*new_next = number_slot->next->next;
                if (!new_next)
                {
                    printf("- missing right operand.");
                    return 0;
                }
                new_next->previous = number_slot->previous;
                if (number_slot->previous)
                {
                    number_slot->previous->next = new_next;
                }
                free_number_slot(number_slot_free_list, number_slot->next);
                free_number_slot(number_slot_free_list, number_slot);
                number_slot = new_next;
            }
            else
            {
                printf("- missing left operand.");
                return 0;
            }
        }
        else
        {
            if (!number_slot->next)
            {
                break;
            }
            number_slot = number_slot->next;
        }
    }
    if (!parse_binary_operation_pair(&number_slot, '*', '/'))
    {
        return 0;
    }
    if (!parse_binary_operation_pair(&number_slot, '+', '-'))
    {
        return 0;
    }
    return number_slot;
}

bool int_equals(struct Integer*a, struct Integer*b)
{
    if (a->sign != b->sign)
    {
        return false;
    }
    if (a->value_count != b->value_count)
    {
        return false;
    }
    for (int i = 0; i < a->value_count; ++a)
    {
        if (a->values[i] != b->values[i])
        {
            return false;
        }
    }
    return true;
}

struct ExtendedGCDInfo
{
    struct Integer gcd;
    struct Integer a_coefficient;
    struct Integer b_coefficient;
    struct Integer a_over_gcd;
    struct Integer b_over_gcd;
};

struct ExtendedGCDInfo extended_gcd(uint32_t**int_arena_cursor, struct Integer a, struct Integer b)
{
    uint32_t*int_arena_savepoint = *int_arena_cursor;
    advance_arena_cursor(int_arena_cursor, sizeof(uint32_t), 5 * max(a.value_count, b.value_count));
    struct ExtendedGCDInfo info;
    struct Integer zero = { 0, 0, 0 };
    struct Integer one = int_new(int_arena_cursor, 1);
    info.a_coefficient = zero;
    info.b_coefficient = one;
    info.b_over_gcd = one;
    info.a_over_gcd = zero;
    while (!int_equals(&a, &zero))
    {
        struct Division division = int_euclidean_divide(int_arena_cursor, &b, &a);
        struct Integer product =
            int_multiply(int_arena_cursor, &info.b_over_gcd, &division.quotient);
        struct Integer m = int_subtract(int_arena_cursor, &info.a_coefficient, &product);
        product = int_multiply(int_arena_cursor, &info.a_over_gcd, &division.quotient);
        struct Integer n = int_subtract(int_arena_cursor, &info.b_coefficient, &product);
        b = a;
        a = division.remainder;
        info.a_coefficient = info.b_over_gcd;
        info.b_coefficient = info.a_over_gcd;
        info.b_over_gcd = m;
        info.a_over_gcd = n;
    }
    info.gcd = b;
    info.b_over_gcd.sign = -info.b_over_gcd.sign;
    int_move_values_to_cursor(&int_arena_savepoint, &info.gcd);
    int_move_values_to_cursor(&int_arena_savepoint, &info.a_coefficient);
    int_move_values_to_cursor(&int_arena_savepoint, &info.b_coefficient);
    int_move_values_to_cursor(&int_arena_savepoint, &info.a_over_gcd);
    int_move_values_to_cursor(&int_arena_savepoint, &info.b_over_gcd);
    rewind_arena_cursor(int_arena_cursor, int_arena_savepoint);
    return info;
}

bool evaluate(struct NumberSlot**number_slot_arena_cursor, struct NumberSlot**number_slot_free_list,
    uint32_t**int_arena_cursor, struct Number*a)
{
    switch (a->operation)
    {
    case 'i':
        return true;
    case '-':
    {
        struct NumberSlot*negative =
            number_slot_new(number_slot_arena_cursor, number_slot_free_list);
        negative->number.operation = 'i';
        negative->number.value = int_new(int_arena_cursor, -1);
        struct NumberSlot*times =
            number_slot_new(number_slot_arena_cursor, number_slot_free_list);
        times->number.operation = '*';
        times->number.left = negative;
        times->number.right = a->right;
        a->operation = '+';
        a->right = times;
    }
    }
    if (!evaluate(number_slot_arena_cursor, number_slot_free_list, int_arena_cursor, a->left))
    {
        return false;
    }
    if (!evaluate(number_slot_arena_cursor, number_slot_free_list, int_arena_cursor, a->right))
    {
        return false;
    }
    if (a->left->operation == 'i')
    {
        if (a->right->operation == 'i')
        {
            switch (a->operation)
            {
            case '+':
            {
                struct Integer sum = int_add(int_arena_cursor, &a->left->value, &a->right->value);
                a->operation = 'i';
                free_number_slot(number_slot_free_list, a->left);
                free_number_slot(number_slot_free_list, a->right);
                a->value = sum;
                return true;
            }
            case '*':
            {
                struct Integer product =
                    int_multiply(int_arena_cursor, &a->left->value, &a->right->value);
                a->operation = 'i';
                free_number_slot(number_slot_free_list, a->left);
                free_number_slot(number_slot_free_list, a->right);
                a->value = product;
                return true;
            }
            case '/':
            {
                if (a->right->value.sign == 0)
                {
                    printf("Tried to divide by 0.");
                    return false;
                }
                if (a->left->value.sign == 0)
                {
                    a->operation = 'i';
                    free_number_slot(number_slot_free_list, a->left);
                    free_number_slot(number_slot_free_list, a->right);
                    a->value.value_count = 0;
                    a->value.sign = 0;
                    return true;
                }
                uint32_t*int_arena_savepoint = *int_arena_cursor;
                advance_arena_cursor(int_arena_cursor, sizeof(uint32_t),
                    a->left->value.value_count + a->right->value.value_count);
                struct ExtendedGCDInfo gcd_info =
                    extended_gcd(int_arena_cursor, a->left->value, a->right->value);
                a->left->value = gcd_info.a_over_gcd;
                a->right->value = gcd_info.b_over_gcd;
                if (a->right->value.sign < 0)
                {
                    a->left->value.sign = -a->left->value.sign;
                    a->right->value.sign = -a->right->value.sign;
                }
                int_move_values_to_cursor(&int_arena_savepoint, &a->left->value);
                struct Integer one = int_new(int_arena_cursor, 1);
                if (int_equals(&a->right->value, &one))
                {
                    struct Number*numerator = a->left;
                    a->operation = 'i';
                    a->value = numerator->value;
                    free_number_slot(number_slot_free_list, numerator);
                    free_number_slot(number_slot_free_list, a->right);
                }
                else
                {
                    int_move_values_to_cursor(&int_arena_savepoint, &a->right->value);
                }
                rewind_arena_cursor(int_arena_cursor, int_arena_savepoint);
                return true;
            }
            case '^':
            {
                if (a->right->value.sign < 0)
                {
                    a->right->value.sign = 1;
                    struct NumberSlot*one =
                        number_slot_new(number_slot_arena_cursor, number_slot_free_list);
                    one->number.operation = 'i';
                    one->number.value = int_new(int_arena_cursor, 1);
                    struct NumberSlot*over =
                        number_slot_new(number_slot_arena_cursor, number_slot_free_list);
                    over->number.operation = '/';
                    over->number.left = one;
                    over->number.right = a->left;
                    a->left = over;
                    return evaluate(number_slot_arena_cursor, number_slot_free_list,
                        int_arena_cursor, a);
                }
                uint32_t*node_int_arena_savepoint = *int_arena_cursor;
                struct Integer exponentiation = int_new(int_arena_cursor, 1);
                struct Integer base_to_a_power_of_two = a->left->value;
                struct Integer exponent = a->right->value;
                struct Integer two = int_new(int_arena_cursor, 2);
                while (exponent.sign > 0)
                {
                    struct Division division =
                        int_euclidean_divide(int_arena_cursor, &exponent, &two);
                    if (division.remainder.sign > 0)
                    {
                        exponentiation = int_multiply(int_arena_cursor, &exponentiation,
                            &base_to_a_power_of_two);
                    }
                    base_to_a_power_of_two = int_multiply(int_arena_cursor, &base_to_a_power_of_two,
                        &base_to_a_power_of_two);
                    exponent = division.quotient;
                }
                memcpy(node_int_arena_savepoint, exponentiation.values,
                    exponentiation.value_count * sizeof(uint32_t));
                exponentiation.values = node_int_arena_savepoint;
                rewind_arena_cursor(int_arena_cursor,
                    node_int_arena_savepoint + exponentiation.value_count);
                a->operation = 'i';
                free_number_slot(number_slot_free_list, a->left);
                free_number_slot(number_slot_free_list, a->right);
                a->value = exponentiation;
                return true;
            }
            }
        }
    }
}

void print_number(struct uint32_t*int_arena_cursor, char*char_arena_cursor, struct Number*number)
{
    if (number->operation == 'i')
    {
        char*char_arena_savepoint = char_arena_cursor;
        int_string(&char_arena_cursor, int_arena_cursor, &number->value);
        push_char(&char_arena_cursor, 0);
        printf("%s", char_arena_savepoint);
        rewind_arena_cursor(&char_arena_cursor, char_arena_savepoint);
    }
    else
    {
        printf("(");
        print_number(int_arena_cursor, char_arena_cursor, number->left);
        printf(")");
        printf("%c", number->operation);
        printf("(");
        print_number(int_arena_cursor, char_arena_cursor, number->right);
        printf(")");
    }
}

int main()
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    page_size = system_info.dwPageSize;
    struct NumberSlot*number_slot_arena = VirtualAlloc(0, 1073741824, MEM_RESERVE, PAGE_READWRITE);
    struct NumberSlot*number_slot_arena_cursor = number_slot_arena;
    struct NumberSlot*number_slot_free_list = 0;
    struct uint32_t*int_arena = VirtualAlloc(0, 1073741824, MEM_RESERVE, PAGE_READWRITE);
    struct uint32_t*int_arena_cursor = int_arena;
    char*char_arena = VirtualAlloc(0, 1073741824, MEM_RESERVE, PAGE_READWRITE);
    char*char_arena_cursor = char_arena;
    while (true)
    {
        if (get_input(&number_slot_arena_cursor, &number_slot_free_list, &int_arena_cursor,
            &int_arena_cursor))
        {
            struct NumberSlot*number_slot = number_slot_arena;
            while (number_slot->next)
            {
                if (number_slot->next->number.operation == '(' &&
                    (number_slot->number.operation == 'i' || number_slot->number.operation == ')'))
                {
                    struct NumberSlot*times =
                        number_slot_new(&number_slot_arena_cursor, &number_slot_free_list);
                    times->number.operation = '*';
                    times->previous = number_slot;
                    times->next = number_slot->next;
                    number_slot->next->previous = times;
                    number_slot->next = times;
                }
                number_slot = number_slot->next;
            }
            struct NumberSlot*input = parse_input(&number_slot_arena_cursor, &number_slot_free_list,
                &int_arena_cursor, number_slot_arena);
            if (input && evaluate(&number_slot_arena_cursor, &number_slot_free_list,
                &int_arena_cursor, &input->number))
            {
                printf("=\n");
                print_number(int_arena_cursor, char_arena_cursor, &input->number);
            }
        }
        printf("\n\n");
        rewind_arena_cursor(&number_slot_arena_cursor, number_slot_arena);
        number_slot_free_list = 0;
        rewind_arena_cursor(&int_arena_cursor, int_arena);
        rewind_arena_cursor(&char_arena_cursor, char_arena);
    }
    return 0;
}
