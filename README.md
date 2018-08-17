Solver is a symbolic computation program that converts radical expressions to as close to a canonical form as possible. The allowed input characters are '+', '-', '\*', '/', '^', '(', ')', '0', '1', '2', '3', '4', '5', '6', '7', '8', and '9'; no white space is allowed, nor are decimal points. Exponentiation by anything the program isn't able to reduce to a rational number is an error, but otherwise, any string of these characters that is grammatical by the usual rules is grammatical here. Note that towers of exponentiations without parentheses, like '2^2^3', follow the convention of evaluating from bottom to top. Entering a string that begins with the character 'q' exits the program.

For each input expression, the program returns an expression equal to that input. The context in which each operation is allowed is limited in output expressions; which expression the program will choose for the output can be understood in terms of these limitations. 

### Addition and Subtraction

The terms of any sums are linearly indepedent over the rationals:

2^(1/2)+3^(1/2)+(5+2\*6^(1/2))^(1/2)=2\*2^(1/2)+2\*3^(1/2)

Note that the program doesn't attempt to make a canonical choice of which terms to express a sum in terms of so, for example, changing the order of the terms in the sum above results in a different choice:

(5+2\*6^(1/2))^(1/2)+2^(1/2)+3^(1/2)=2(5+2\*6^(1/2))^(1/2)

### Multiplication

The factors of any product consist of one or more surds, or a rational number and one or more surds. Thus all products of rational numbers are consolidated, and all products involving one or more sums are distributed:

(1/2+2^(1/2))(1/3+3^(1/2))=1/6+6^(1/2)+2^(1/2)/3+3^(1/2)/2

The program attempts to consolidate the factors of products of surds as much as possible. It is always possible to consolidate a product of real-valued surds into a single surd:

2^(1/3)\*3^(1/2)=108^(1/6)

For complex-valued surds, however, the exponentiation identity used for such consolidation is true only up to multiplication by a root of unity. Roughly speaking, roots of unity other than 1 come into play in cases where the sum of the complex arguments of two terms in a multiplication is greater than or equal to 2pi. In general, performing the consolidation then multiplying by the necessary root of unity causes infinite looping problems. These problems don't arise in products of surds of the same index, so the program consolidates all surds of the same index, then consolidates the resulting distinct-index factors in pairs until no two factors have arguments whose sum is less than 2pi.

### Division

The denominator of any division is a positive integer. Thus all denominators are rationalized: 

1/(1+2^(1/3))=1/3-2^(1/3)/3+4^(1/3)/3

and the signs of negative denominators are transferred to the numerators:

1/-2=-1/2

In addition, no division is allowed inside radicands:

(1/2+2^(1/2))^(1/2)=(2+4\*2^(1/2))^(1/2)/2

no numerator is a sum:

(1+2^(1/2))/2=1/2+2^(1/2)/2

and no numerator has an integer factor that is divisible by the denominator:

6\*2^(1/2)/4=3\*2^(1/2)/2

### Exponentiation

Any exponent is a fraction with numerator 1:

2^(-4/3)=4^(1/3)/4

any radicand is either an integer or a sum:

(2\*2^(1/2))^(1/3)=8^(1/6)

and no radicand of a degree n surd has an integer factor that is an nth power:

(8+12\*2^(1/2))^(1/2)=2(2+3\*2^(1/2))^(1/2)
