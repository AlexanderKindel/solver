Solver is a symbolic computation program that converts radical expressions to as close to a canonical form as possible. The allowed input characters are '+', '-', '\*', '/', '^', '(', ')', '0', '1', '2', '3', '4', '5', '6', '7', '8', and '9'; no white space is allowed, nor are decimal points. Exponentiation by anything the program isn't able to reduce to a rational number is an error, but otherwise, any string of these characters that is grammatical by the usual rules is grammatical here. Note that towers of exponentiations without parentheses, like '2^2^3', follow the convention of evaluating from bottom to top. Entering a string that begins with the character 'q' exits the program.

For each input expression, the program returns an expression equal to it that obeys the following structural rules, in addition to a few stylistic ones:

1.) The terms of any sum are linearly independent over the rationals: 

1/2+2^(1/2)+2\*2^(1/2)+2/3=7/6+3\*2^(1/2)

2.) Any product consists of one rational number and one surd:

2\*3^(1/2)\*2^(1/3)\*2/3=4\*108^(1/6)/3

3.) The denominator of any division is a positive integer: 

1/(1+2^(1/2))=-1/2+2^(1/2)/2

4.) No numerator of a division has an integer factor that is divisible by the denominator:

6\*2^(1/2)/4=3\*2^(1/2)/2

5.) No numerator of a division is a sum:

(1+2^(1/2))/2=1/2+2^(1/2)/2

6.) Any exponent is a fraction with numerator 1:

2^(-4/3)=4^(1/3)/4

7.) There are no divisions within radicands:

(1/2+2^(1/2))^(1/2)=(2+4\*2^(1/2))^(1/2)/2

8.) Any radicand is either an integer or a sum:

(2\*2^(1/2))^(1/3)=8^(1/6)

9.) No radicand of a degree n surd has an integer factor that is an nth power:

(8+12\*2^(1/2))^(1/2)=2(2+3\*2^(1/2))^(1/2)
