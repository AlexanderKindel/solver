﻿struct IntervalStrings
{
    char*bound_type;//0
    char*bound_type_name;//1
};

unsigned short ALL_INTERVAL_SUBSTITUTION_SETS = 0b11;

char*interval_max_magnitude_template =
"struct %00*%01_interval_max_magnitude(struct Stack*output_stack, struct Stack*local_stack, struct %00Interval*a)\n"
"{\n"
"    void*local_stack_savepoint = local_stack->cursor;\n"
"    void*out = %01_max(output_stack, local_stack, %01_magnitude(local_stack, a->min), %01_magnitude(local_stack, a->max));\n"
"    local_stack->cursor = local_stack_savepoint;\n"
"    return out;\n"
"}\n";

char*interval_add_template =
"struct %00Interval*%01_interval_add(struct Stack*output_stack, struct Stack*local_stack, struct %00Interval*a, struct %00Interval*b)\n"
"{\n"
"    struct %00Interval*out = ALLOCATE(output_stack, struct %00Interval);\n"
"    out->min = %01_add(output_stack, local_stack, a->min, b->min);\n"
"    out->max = %01_add(output_stack, local_stack, a->max, b->max);\n"
"    return out;\n"
"}\n";

char*interval_negative_template =
"struct %00Interval*%01_interval_negative(struct Stack*output_stack, struct Stack*local_stack, struct %00Interval*a)\n"
"{\n"
"    struct %00Interval*out = ALLOCATE(output_stack, struct %00Interval);\n"
"    if (%01_sign(a->min) != -%01_sign(a->max))\n"
"    {\n"
"        out->min = %01_negative(output_stack, a->max);\n"
"        out->max = %01_negative(output_stack, a->min);\n"
"    }\n"
"    else\n"
"    {\n"
"        void*local_stack_savepoint = local_stack->cursor;\n"
"        out->min = %01_copy(output_stack, %01_min(output_stack, local_stack, %01_negative(local_stack, a->max), a->min));\n"
"        out->max = %01_copy(output_stack, %01_max(output_stack, local_stack, %01_negative(local_stack, a->min), a->max));\n"
"        local_stack->cursor = local_stack_savepoint;\n"
"    }\n"
"    return out;\n"
"}\n";

char*interval_multiply_template =
"struct %00Interval*%01_interval_multiply(struct Stack*output_stack, struct Stack*local_stack, struct %00Interval*a, struct %00Interval*b)\n"
"{\n"
"    struct %00Interval*out = ALLOCATE(output_stack, struct %00Interval);\n"
"    if (%01_sign(a->min) >= 0)\n"
"    {\n"
"        if (%01_sign(b->min) >= 0)\n"
"        {\n"
"            out->min = %01_multiply(output_stack, local_stack, a->min, b->min);\n"
"            out->max = %01_multiply(output_stack, local_stack, a->max, b->max);\n"
"        }\n"
"        else\n"
"        {\n"
"            out->min = %01_multiply(output_stack, local_stack, a->max, b->min);\n"
"            if (%01_sign(b->max) <= 0)\n"
"            {\n"
"                out->max = %01_multiply(output_stack, local_stack, a->min, b->max);\n"
"            }\n"
"            else\n"
"            {\n"
"                out->max = %01_multiply(output_stack, local_stack, a->max, b->max);\n"
"            }\n"
"        }\n"
"    }\n"
"    else if (%01_sign(a->max) <= 0)\n"
"    {\n"
"        if (%01_sign(b->min) >= 0)\n"
"        {\n"
"            out->min = %01_multiply(output_stack, local_stack, a->min, b->max);\n"
"            out->max = %01_multiply(output_stack, local_stack, a->max, b->min);\n"
"        }\n"
"        else\n"
"        {\n"
"            if (%01_sign(b->max) <= 0)\n"
"            {\n"
"                out->min = %01_multiply(output_stack, local_stack, a->max, b->max);\n"
"            }\n"
"            else\n"
"            {\n"
"                out->min = %01_multiply(output_stack, local_stack, a->min, b->max);\n"
"            }\n"
"            out->max = %01_multiply(output_stack, local_stack, a->min, b->min);\n"
"        }\n"
"    }\n"
"    else\n"
"    {\n"
"        if (%01_sign(b->min) >= 0)\n"
"        {\n"
"            out->min = %01_multiply(output_stack, local_stack, a->min, b->max);\n"
"            out->max = %01_multiply(output_stack, local_stack, a->max, b->max);\n"
"        }\n"
"        else if (%01_sign(b->max) <= 0)\n"
"        {\n"
"            out->min = %01_multiply(output_stack, local_stack, a->max, b->min);\n"
"            out->max = %01_multiply(output_stack, local_stack, a->min, b->min);\n"
"        }\n"
"        else\n"
"        {\n"
"            void*local_stack_savepoint = local_stack->cursor;\n"
"            out->min = %01_copy(output_stack, %01_min(output_stack, local_stack,\n"
"                %01_multiply(local_stack, output_stack, a->min, b->max),\n"
"                %01_multiply(local_stack, output_stack, a->max, b->min)));\n"
"            out->max = %01_copy(output_stack, %01_max(output_stack, local_stack,\n"
"                %01_multiply(local_stack, output_stack, a->min, b->min),\n"
"                %01_multiply(local_stack, output_stack, a->max, b->max)));\n"
"            local_stack->cursor = local_stack_savepoint;\n"
"        }\n"
"    }\n"
"    return out;\n"
"}\n";