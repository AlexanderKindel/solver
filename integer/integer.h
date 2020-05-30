struct SmallInteger
{
    size_t value_count;
    int8_t sign;
    uint32_t value;
};

#define INT(value, sign) (struct Integer*)&(struct SmallInteger) { 1, sign, value }

struct Integer*integer_allocate(struct Stack*output_stack, size_t value_count);
struct Integer*integer_copy(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_initialize(struct Stack*output_stack, uint32_t value, int8_t sign);
struct Integer*char_to_integer(struct Stack*output_stack, char a);
struct Integer*size_t_to_integer(struct Stack*output_stack, size_t a);
size_t integer_to_size_t(struct Integer*a);
bool integer_equals(struct Integer*a, struct Integer*b);
struct Integer*integer_get_magnitude(struct Stack*output_stack, struct Integer*a);
void integer_add_to_a_in_place(struct Integer*restrict a, struct Integer*restrict b);
struct Integer*integer_add(struct Stack*output_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_negate(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend);
struct Integer*integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
struct IntegerDivision integer_euclidean_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*dividend, struct Integer*divisor);
struct Integer*integer_double(struct Stack*output_stack, struct Integer*a);
void integer_halve_in_place(struct Integer*a);
struct Integer*integer_halve(struct Stack*output_stack, struct Integer*a);
struct Integer*integer_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*base, struct Integer*exponent);
struct ExtendedGCDInfo integer_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
int8_t integer_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Integer*a, struct Integer*b);
struct Integer*integer_get_lcm(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b);
struct Integer*integer_take_square_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a);
struct Integer*get_next_prime(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b);
size_t size_t_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    size_t**out, size_t a);
size_t integer_factor(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Factor**out, struct Integer*a);
size_t integer_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*a);

#include "integer/float/float.h"
#include "integer/integer_polynomial/integer_polynomial.h"
#include "integer/modded_integer/modded_integer.h"
#include "integer/rational/rational.h"