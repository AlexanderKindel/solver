using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleSolver
{
    class Math
    {
        public abstract class Number : IComparable
        {
            protected abstract Number add(Number number);
            protected abstract Number multiply(Number number);
            public abstract Number negative();
            public abstract Number reciprocal();
            public abstract int CompareTo(Object obj);
            public static Number operator +(Number a, Number b)
            {
                return a.add(b);
            }
            public static Number operator -(Number a, Number b)
            {
                return a.add(b.negative());
            }
            public static Number operator *(Number a, Number b)
            {
                return a.multiply(b);
            }
            public static Number operator /(Number a, Number b)
            {
                return a.multiply(b.reciprocal());
            }
            public override bool Equals(object obj)
            {
                if (obj.GetType() != GetType())
                    return false;
                return CompareTo(obj) == 0;
            }
        }
        public abstract class Term : Number
        {
            protected override abstract Number add(Number number);
            protected override abstract Number multiply(Number number);
            public override abstract Number negative();
            public override abstract Number reciprocal();
            public override abstract int CompareTo(Object obj);
            public static Number operator +(Term a, Number b)
            {
                return a.add(b);
            }
            public static Number operator -(Term a, Number b)
            {
                return a.add(b.negative());
            }
            public static Number operator *(Term a, Number b)
            {
                return a.multiply(b);
            }
            public static Number operator /(Term a, Number b)
            {
                return a.multiply(b.reciprocal());
            }
            public override bool Equals(object obj)
            {
                if (obj.GetType() != GetType())
                    return false;
                return CompareTo(obj) == 0;
            }
        }
        public abstract class Rational : Term
        {
            protected override abstract Number add(Number number);
            protected override abstract Number multiply(Number number);
            public override abstract Number negative();
            public override abstract Number reciprocal();
            public abstract bool isIntegral();
            public override abstract int CompareTo(Object obj);
            public static Number operator +(Rational a, Number b)
            {
                return a.add(b);
            }
            public static Number operator -(Rational a, Term b)
            {
                return a.add(b.negative());
            }
            public static Number operator *(Rational a, Term b)
            {
                return a.multiply(b);
            }
            public static Number operator /(Rational a, Term b)
            {
                return a.multiply(b.reciprocal());
            }
        }
        public class Fraction : Rational
        {
            public int Numerator { get; }
            public int Denominator { get; }
            public Fraction(int numerator, int denominator)
            {
                if (denominator == 0)
                    throw new DivideByZeroException();
                if (numerator == 0)
                    denominator = 1;
                else
                {
                    int GCD = getGCD(numerator, denominator);
                    numerator /= GCD;
                    denominator /= GCD;
                    if (denominator < 0)
                    {
                        numerator *= -1;
                        denominator *= -1;
                    }
                }
                Numerator = numerator;
                Denominator = denominator;
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Fraction(Numerator * fraction.Denominator +
                        fraction.Numerator * Denominator, Denominator * fraction.Denominator);
                }
                return number + this;
            }
            public override Number negative()
            {
                return new Fraction(-Numerator, Denominator);
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Fraction(Numerator * fraction.Numerator,
                        Denominator * fraction.Denominator);
                }
                return number * this;
            }
            public override Number reciprocal()
            {
                return new Fraction(Denominator, Numerator);
            }
            public override bool isIntegral()
            {
                return Denominator == 1;
            }
            public override int CompareTo(Object obj)
            {
                if (obj is Fraction)
                {
                    Fraction fraction = (Fraction)obj;
                    return Numerator * fraction.Denominator - Denominator * fraction.Numerator;
                }
                ComplexNumber complexNumber = (ComplexNumber)obj;
                return complexNumber.CompareTo(this);
            }
            public override int GetHashCode()
            {
                return Numerator | Denominator;
            }
            public override string ToString()
            {
                if (Denominator == 1)
                    return Numerator.ToString();
                return Numerator.ToString() + "/" + Denominator.ToString();
            }
        }
        public class ComplexNumber : Rational
        {
            public Fraction Real { get; }
            public Fraction Imaginary { get; }
            public ComplexNumber(Fraction real, Fraction imaginary)
            {
                Real = real;
                Imaginary = imaginary;
            }
            public ComplexNumber(int real, int imaginary)
            {
                Real = new Fraction(real, 1);
                Imaginary = new Fraction(imaginary, 1);
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new ComplexNumber((Fraction)(Real + fraction), Imaginary);
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    Fraction real = (Fraction)(Real + complexNumber.Real);
                    Fraction imaginary = (Fraction)(Imaginary + complexNumber.Imaginary);
                    if (imaginary.Numerator == 0)
                        return real;
                    return new ComplexNumber(real, imaginary);
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    Fraction real = (Fraction)(Real * fraction);
                    Fraction imaginary = (Fraction)(Imaginary * fraction);
                    if (imaginary.Numerator == 0)
                        return real;
                    return new ComplexNumber((Fraction)(Real * fraction),
                        (Fraction)(Imaginary * fraction));
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    Fraction real = (Fraction)(Real * complexNumber.Real -
                        Imaginary * complexNumber.Imaginary);
                    Fraction imaginary = (Fraction)(Real * complexNumber.Imaginary +
                        Imaginary * complexNumber.Real);
                    if (imaginary.Numerator == 0)
                        return real;
                    return new ComplexNumber(real, imaginary);
                }
                return number * this;
            }
            public override Number negative()
            {
                return new ComplexNumber((Fraction)Real.negative(), (Fraction)Imaginary.negative());
            }
            public override Number reciprocal()
            {
                Fraction denominator = (Fraction)(Real * Real + Imaginary * Imaginary);
                return new ComplexNumber((Fraction)(Real / denominator),
                    (Fraction)(Imaginary.negative() / denominator));
            }
            public override bool isIntegral()
            {
                return Real.Denominator == 1 && Imaginary.Denominator == 1;
            }
            public override int CompareTo(Object obj)
            {
                int comparison;
                if (obj is Fraction)
                {
                    Fraction fraction = (Fraction)obj;
                    comparison = Real.CompareTo(fraction);
                    if (comparison != 0)
                        return comparison;
                    return Imaginary.Numerator;
                }
                ComplexNumber complexNumber = (ComplexNumber)obj;
                comparison = Real.CompareTo(complexNumber.Real);
                if (comparison != 0)
                    return comparison;
                return Imaginary.CompareTo(complexNumber.Imaginary);
            }
            public override int GetHashCode()
            {
                return Real.GetHashCode() | Imaginary.GetHashCode();
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                if (Real.Numerator != 0)
                {
                    output.Append(Real.ToString());
                    if (Imaginary.Numerator > 0)
                        output.Append('+');
                }
                if (Imaginary.Numerator != 1)
                    if (Imaginary.Numerator == -1)
                        output.Append('-');
                    else
                        output.Append(Imaginary.Numerator);
                output.Append('i');
                if (Imaginary.Denominator != 1)
                    output.Append('/' + Imaginary.Denominator.ToString());
                return output.ToString();
            }
        }
        public class Exponentiation : Term
        {
            public Number Base { get; }
            public Number Exponent { get; }
            public Exponentiation(Number expBase, Number exponent)
            {
                Base = expBase;
                Exponent = exponent;
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return this;
                    return new Expression(new List<Term> { this, fraction });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Expression(new List<Term> { this, complexNumber });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Exponent.CompareTo(exponentiation.Exponent) == 0)
                        return new Exponentiation(Base + exponentiation.Base, Exponent);
                    return new Expression(new List<Term> { this, exponentiation });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == fraction.Denominator)
                        return this;
                    return new Product(fraction, new List<Exponentiation> { this });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Product(complexNumber, new List<Exponentiation> { this });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Exponent.Equals(exponentiation.Exponent))
                        return new Exponentiation(Base * exponentiation.Base, Exponent);
                    if (Base.Equals(exponentiation.Base))
                        return new Exponentiation(Base, Exponent + exponentiation.Exponent);
                    return new Product(new Fraction(1, 1),
                        new List<Exponentiation> { this, exponentiation });
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Exponentiation(Base.negative(), Exponent);
            }
            public override Number reciprocal()
            {
                return new Exponentiation(Base, Exponent.negative());
            }
            public override int CompareTo(Object obj)
            {
                Exponentiation exponentiation = (Exponentiation)obj;
                int comparison = Base.CompareTo(exponentiation.Base);
                if (comparison != 0)
                    return comparison;
                return Exponent.CompareTo(exponentiation.Exponent);
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                if (Base is Fraction)
                {
                    Fraction baseFraction = (Fraction)Base;
                    if (baseFraction.Denominator == 1 && baseFraction.Numerator >= 0)
                        output.Append(baseFraction.Numerator);
                    else
                        output.Append('(' + baseFraction.ToString() + ')');
                }
                else if (Base is ComplexNumber)
                {
                    ComplexNumber complexBase = (ComplexNumber)Base;
                    if (complexBase.Real.Numerator == 0 &&
                        complexBase.Imaginary.Numerator == complexBase.Imaginary.Denominator)
                        output.Append('i');
                    else
                        output.Append('(' + complexBase.ToString() + ')');
                }
                else
                    output.Append('(' + Base.ToString() + ')');
                output.Append('^');
                if (Exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)Exponent;
                    if (exponentFraction.Denominator == 1 && exponentFraction.Numerator >= 0)
                        return output.Append(exponentFraction.Numerator).ToString();
                }
                if (Exponent is ComplexNumber)
                {
                    ComplexNumber complexExponent = (ComplexNumber)Exponent;
                    if (complexExponent.Real.Numerator == 0 &&
                        complexExponent.Imaginary.Numerator ==
                        complexExponent.Imaginary.Denominator)
                        return output.Append('i').ToString();
                }
                return output.Append('(' + Exponent.ToString() + ')').ToString();
            }
        }
        public class Product : Term
        {
            public Rational Scalar { get; }
            public List<Exponentiation> Exponentiations { get; }
            public Product(Rational rational, List<Exponentiation> exponentiations)
            {
                Scalar = rational;
                Exponentiations = exponentiations;
                Exponentiations.Sort();
                for (int i = 0; i < Exponentiations.Count - 1;)
                {
                    if (Exponentiations[i].Base.CompareTo(Exponentiations[i + 1].Base) == 0)
                    {
                        Exponentiations[i] = new Exponentiation(Exponentiations[i].Base,
                            Exponentiations[i].Exponent + Exponentiations[i + 1].Exponent);
                        Exponentiations.RemoveAt(i + 1);
                        continue;
                    }
                    ++i;
                }
            }
            public Product(Exponentiation exponentiation)
            {
                Scalar = new Fraction(1, 1);
                Exponentiations = new List<Exponentiation> { exponentiation };
            }
            public Product(Rational rational)
            {
                Scalar = rational;
                Exponentiations = new List<Exponentiation>();
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Expression(new List<Term> { this, fraction });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Expression(new List<Term> { this, complexNumber });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    return new Expression(new List<Term> { this, exponentiation });
                }
                if (number is Product)
                {
                    Product product = (Product)number;
                    return new Expression(new List<Term> { this, product });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Product((Rational)(Scalar * fraction), Exponentiations);
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Product((Rational)(Scalar * complexNumber), Exponentiations);
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    List<Exponentiation> exponentiations = Exponentiations;
                    for (int i = 0; i < exponentiations.Count; ++i)
                    {
                        Number product = exponentiations[i] * exponentiation;
                        if (product is Exponentiation)
                        {
                            exponentiations[i] = (Exponentiation)product;
                            return new Product(Scalar, exponentiations);
                        }
                    }
                    exponentiations.Add(exponentiation);
                    return new Product(Scalar, exponentiations);
                }
                if (number is Product)
                {
                    Product product = (Product)number;
                    List<Exponentiation> exponentiations = Exponentiations;
                    for (int i = 0; i < product.Exponentiations.Count; ++i)
                    {
                        bool addExponentiation = true;
                        for (int j = 0; j < exponentiations.Count; ++j)
                        {
                            Number term = product.Exponentiations[i] * exponentiations[j];
                            if (term is Exponentiation)
                            {
                                exponentiations[i] = (Exponentiation)term;
                                addExponentiation = false;
                                break;
                            }
                        }
                        if (addExponentiation)
                            exponentiations.Add(product.Exponentiations[i]);
                    }
                    return new Product((Rational)(Scalar * product.Scalar), exponentiations);
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Product((Rational)Scalar.negative(), Exponentiations);
            }
            public override Number reciprocal()
            {
                List<Exponentiation> exponentiations = new List<Exponentiation>();
                foreach (Exponentiation exponentation in Exponentiations)
                    exponentiations.Add((Exponentiation)exponentation.reciprocal());
                return new Product((Rational)Scalar.reciprocal(), exponentiations);
            }
            public override int CompareTo(object obj)
            {
                Product product = (Product)obj;
                int comparison = Scalar.CompareTo(product.Scalar);
                if (comparison != 0)
                    return comparison;
                int countComparison = Exponentiations.Count - product.Exponentiations.Count;
                if (countComparison != 0)
                    return countComparison;
                for (int i = 0; i < Exponentiations.Count; ++i)
                {
                    comparison = Exponentiations[i].CompareTo(product.Exponentiations[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder(Exponentiations[0].ToString());
                for (int i = 1; i < Exponentiations.Count; ++i)
                {
                    String exponentiationString = Exponentiations[i].ToString();
                    if (output[output.Length - 1] != ')' && exponentiationString[0] != '(')
                        output.Append('*');
                    output.Append(exponentiationString);
                }
                StringBuilder scalarString = new StringBuilder();
                if (Scalar is Fraction)
                {
                    Fraction scalarFraction = (Fraction)Scalar;
                    if (scalarFraction.Numerator != 1)
                    {
                        if (scalarFraction.Numerator == -1)
                            scalarString.Append('-');
                        else
                            scalarString.Append(scalarFraction.Numerator);
                    }
                    if (output[0] != '(' && scalarString.Length > 0)
                        scalarString.Append('*');
                    if (scalarFraction.Denominator != 1)
                        output.Append('/' + scalarFraction.Denominator.ToString());
                }
                else
                {
                    ComplexNumber complexScalar = (ComplexNumber)Scalar;
                    if (complexScalar.Real.Numerator == 0)
                    {
                        if (complexScalar.Imaginary.Numerator == 1)
                            scalarString.Append('i');
                        else
                        {
                            if (complexScalar.Imaginary.Numerator == -1)
                                scalarString.Append("-i");
                            else
                                scalarString.Append(
                                    complexScalar.Imaginary.Numerator.ToString() + 'i');
                        }
                        if (output[0] != '(' && scalarString.Length > 0)
                            scalarString.Append('*');
                        if (complexScalar.Imaginary.Denominator != 1)
                            output.Append('/' + complexScalar.Imaginary.Denominator.ToString());
                    }
                    else
                        output.Insert(0, '(' + Scalar.ToString() + ')');
                }
                return output.Insert(0, scalarString).ToString();
            }
        }
        public class Expression : Number
        {
            public List<Term> m_terms;
            public Expression(List<Term> terms)
            {
                m_terms = terms;
                m_terms.Sort();
            }
            protected override Number add(Number number)
            {
                List<Term> terms = m_terms;
                if (number is Expression)
                {
                    Expression expression = (Expression)number;
                    terms.AddRange(expression.m_terms);
                }
                else
                    terms.Add((Term)number);
                return new Expression(terms);
            }
            protected override Number multiply(Number number)
            {
                List<Term> terms = new List<Term>();
                if (number is Expression)
                {
                    Expression multiplier = (Expression)number;
                    foreach (Term multiplicandProduct in m_terms)
                        foreach (Term multiplierProduct in multiplier.m_terms)
                            terms.Add((Term)(multiplicandProduct * multiplierProduct));
                }
                else
                    foreach (Term term in m_terms)
                        terms.Add((Term)(term * number));
                return new Expression(terms);
            }
            public override Number negative()
            {
                List<Term> terms = new List<Term>();
                foreach (Term term in m_terms)
                    terms.Add((Term)term.negative());
                return new Expression(terms);
            }
            public override Number reciprocal()
            {
                return new Exponentiation(this, new Fraction(-1, 1));
            }
            public override int CompareTo(object obj)
            {
                Expression expression = (Expression)obj;
                int countComparison = m_terms.Count - expression.m_terms.Count;
                if (countComparison != 0)
                    return countComparison;
                for (int i = 0; i < m_terms.Count; ++i)
                {
                    int comparison = m_terms[i].CompareTo(expression.m_terms[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override String ToString()
            {
                StringBuilder output = new StringBuilder(m_terms[0].ToString());
                for (int i = 1; i < m_terms.Count; ++i)
                    output.Append('+' + m_terms[i].ToString());
                return output.ToString();
            }
        }
        public class NumberList : Number
        {
            List<Number> Numbers { get; }
            public NumberList(List<Number> numbers)
            {
                Numbers = numbers;
                for (int i = 0; i < Numbers.Count; ++i)
                    for (int j = i + 1; j < Numbers.Count; ++j)
                        if (Numbers[i].Equals(Numbers[j]))
                            Numbers.RemoveAt(j);
            }
            protected override Number add(Number number)
            {
                List<Number> output = new List<Number>();
                if (number is NumberList)
                {
                    NumberList list = (NumberList)number;
                    foreach (Number n in Numbers)
                        foreach (Number m in list.Numbers)
                            output.Add(n + m);
                }
                else
                    foreach (Number n in Numbers)
                        output.Add(n + number);
                return new NumberList(output);
            }
            public override Number negative()
            {
                List<Number> output = new List<Number>();
                foreach (Number n in Numbers)
                    output.Add(n.negative());
                return new NumberList(output);
            }
            protected override Number multiply(Number number)
            {
                List<Number> output = new List<Number>();
                if (number is NumberList)
                {
                    NumberList list = (NumberList)number;
                    foreach (Number n in Numbers)
                        foreach (Number m in list.Numbers)
                            output.Add(n * m);
                }
                else
                    foreach (Number n in Numbers)
                        output.Add(n * number);
                return new NumberList(output);
            }
            public override Number reciprocal()
            {
                List<Number> output = new List<Number>();
                foreach (Number n in Numbers)
                    output.Add(n.reciprocal());
                return new NumberList(output);
            }
            public override int CompareTo(object obj)
            {
                NumberList list = (NumberList)obj;
                return Numbers.Count - list.Numbers.Count;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder('{' + Numbers[0].ToString());
                for (int i = 1; i < Numbers.Count; ++i)
                    output.Append(',' + Numbers[i].ToString());
                output.Append('}');
                return output.ToString();
            }
        }
        static public Dictionary<int, int> getFactorization(int x)
        {//The keys respresent the factors, and the values represent their multiplicities.
            Dictionary<int, int> factors = new Dictionary<int, int>();
            if (x < 0)
                x *= -1;
            int factor = 2;
            while (factor <= x / 2)
                if (x % factor == 0)
                {
                    if (factors.ContainsKey(factor))
                        factors[factor] += 1;
                    else
                        factors.Add(factor, 1);
                    x /= factor;
                }
                else
                    ++factor;
            if (factors.ContainsKey(x))
                factors[x] += 1;
            else if (x != 1)
                factors.Add(x, 1);
            return factors;
        }
        static public Dictionary<Rational, int> getFactorization(Rational gaussian)
        {
            Dictionary<Rational, int> factors = new Dictionary<Rational, int>();
            int sumOfRealAndImaginary = 2;
            while (true)
            {
                int realPart = sumOfRealAndImaginary / 2;
                int imaginaryPart = sumOfRealAndImaginary - realPart;
                for (int i = 0; realPart - i >= 0;)
                {
                    int gaussianMagnitudeSquared;
                    if (gaussian is Fraction)
                    {
                        Fraction fractionGaussian = (Fraction)gaussian;
                        gaussianMagnitudeSquared = fractionGaussian.Numerator *
                            fractionGaussian.Numerator;
                    }
                    else
                    {
                        ComplexNumber complexGaussian = (ComplexNumber)gaussian;
                        gaussianMagnitudeSquared = complexGaussian.Real.Numerator *
                            complexGaussian.Real.Numerator + complexGaussian.Imaginary.Numerator *
                            complexGaussian.Imaginary.Numerator;
                    }
                    if ((realPart - i) * (realPart - i) + (imaginaryPart + i) *
                    (imaginaryPart + i) > gaussianMagnitudeSquared)
                    {
                        if (factors.ContainsKey(gaussian))
                            factors[gaussian] += 1;
                        else
                            factors.Add(gaussian, 1);
                        return factors;
                    }
                    Rational factor;
                    bool testFactor()
                    {
                        Rational quotient = (Rational)(gaussian / factor);
                        if (quotient.isIntegral())
                        {
                            if (factors.ContainsKey(factor))
                                factors[factor] += 1;
                            else
                                factors.Add(factor, 1);
                            gaussian = quotient;
                            return true;
                        }
                        return false;
                    }
                    if (imaginaryPart + i == 0)
                    {
                        factor = new Fraction(realPart - i, 1);
                        if (testFactor())
                            continue;
                    }
                    else
                    {
                        factor = new ComplexNumber(realPart - i, imaginaryPart + i);
                        if (testFactor())
                            continue;
                        factor = new ComplexNumber(realPart - i, -imaginaryPart - i);
                        if (testFactor())
                            continue;
                    }
                    ++i;
                }
                ++sumOfRealAndImaginary;
            }
        }
        static public int getGCD(int x, int y)
        {
            if (x < 0)
                x *= -1;
            if (y < 0)
                y *= -1;
            while (true)
                if (x > y)
                {
                    if (x % y == 0)
                        return y;
                    x %= y;
                }
                else
                {
                    if (y % x == 0)
                        return x;
                    y %= x;
                }
        }
        static public int getLCM(int x, int y)
        {
            return x / getGCD(x, y) * y;
        }
        static public Number exponentiate(Number expBase, Number exponent)
        {
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                if (expBase is Fraction)
                {
                    Fraction baseFraction = (Fraction)expBase;
                    if (exponentFraction.Numerator < 0)
                    {
                        exponentFraction = (Fraction)(exponentFraction.negative());
                        baseFraction = (Fraction)(baseFraction.reciprocal());
                    }
                    Dictionary<int, int> rationalizerFactors =
                        getFactorization(baseFraction.Denominator);
                    Dictionary<int, int> radicandFactors = getFactorization(baseFraction.Numerator);
                    foreach (int factor in rationalizerFactors.Keys)
                        if (radicandFactors.ContainsKey(factor))
                            radicandFactors[factor] += (exponentFraction.Denominator - 1) *
                                rationalizerFactors[factor];
                        else
                            radicandFactors.Add(factor, (exponentFraction.Denominator - 1) *
                                rationalizerFactors[factor]);
                    int root = 1;
                    int radicand = 1;
                    foreach (int factor in radicandFactors.Keys)
                    {
                        for (int i = 0; i < radicandFactors[factor] /
                            exponentFraction.Denominator; ++i)
                            root *= factor;
                        for (int i = 0; i < radicandFactors[factor] %
                            exponentFraction.Denominator; ++i)
                            radicand *= factor;
                    }
                    Fraction rootToPower = new Fraction(1, 1);
                    int radicandToPower = 1;
                    for (int i = 0; i < exponentFraction.Numerator %
                        exponentFraction.Denominator; ++i)
                    {
                        rootToPower = (Fraction)(rootToPower *
                            new Fraction(root, baseFraction.Denominator));
                        radicandToPower *= radicand;
                    }
                    for (int i = 0; i < exponentFraction.Numerator /
                        exponentFraction.Denominator; ++i)
                        rootToPower = (Fraction)(rootToPower * baseFraction);
                    if (radicandToPower == 1)
                        return rootToPower;
                    if (radicandToPower == -1)
                        return rootToPower * new ComplexNumber(0, -1);
                    return new Exponentiation(new Fraction(radicandToPower, 1),
                        new Fraction(1, exponentFraction.Denominator)) * rootToPower;
                }
                if (expBase is ComplexNumber)
                {
                    ComplexNumber baseComplexNumber = (ComplexNumber)expBase;
                    if (exponentFraction.Numerator < 0)
                    {
                        exponentFraction = (Fraction)(exponentFraction.negative());
                        baseComplexNumber = (ComplexNumber)(baseComplexNumber.reciprocal());
                    }
                    int rationalizer = getLCM(baseComplexNumber.Real.Denominator,
                        baseComplexNumber.Imaginary.Denominator);
                    baseComplexNumber = (ComplexNumber)(baseComplexNumber *
                        new Fraction(rationalizer, 1));
                    Dictionary<Rational, int> radicandFactors = getFactorization(baseComplexNumber);
                    Rational root = new Fraction(1, 1);
                    Rational radicand = new Fraction(1, 1);
                    foreach (Rational factor in radicandFactors.Keys)
                    {
                        for (int i = 0; i < radicandFactors[factor] /
                            exponentFraction.Denominator; ++i)
                            root = (Rational)(root * factor);
                        for (int i = 0; i < radicandFactors[factor] %
                            exponentFraction.Denominator; ++i)
                            radicand = (Rational)(radicand * factor);
                    }
                    Rational rootToPower = new Fraction(1, 1);
                    Rational radicandToPower = new Fraction(1, 1);
                    for (int i = 0; i < exponentFraction.Numerator %
                        exponentFraction.Denominator; ++i)
                    {
                        rootToPower = (Rational)(rootToPower * root *
                            new Fraction(1, rationalizer));
                        radicandToPower = (Rational)(radicandToPower * radicand);
                    }
                    for (int i = 0; i < exponentFraction.Numerator /
                        exponentFraction.Denominator; ++i)
                        rootToPower = (Rational)(rootToPower * baseComplexNumber *
                            new Fraction(1, rationalizer));
                    if (radicand is Fraction)
                    {
                        Fraction radicandFraction = (Fraction)radicand;
                        if (radicandFraction.Numerator == radicandFraction.Denominator)
                            return rootToPower;
                        if (radicandFraction.Numerator == -radicandFraction.Denominator)
                            return rootToPower * new ComplexNumber(0, -1);
                    }
                    return new Exponentiation(radicandToPower,
                        new Fraction(1, exponentFraction.Denominator)) * rootToPower;
                }
            }
            return new Exponentiation(expBase, exponent);
        }
    }
    class Solver
    {
        internal static Math.Number evaluateInput(List<char> operations, List<Math.Number> numbers)
        {
            for (int i = 0; i < operations.Count;)
            {
                if (operations[i] == '(')
                {
                    int numOfUnmatchedParens = 1;
                    int matchingParenIndex = i;
                    while (numOfUnmatchedParens > 0)
                    {
                        matchingParenIndex += 1;
                        if (operations[matchingParenIndex] == '(')
                            numOfUnmatchedParens += 1;
                        else if (operations[matchingParenIndex] == ')')
                            numOfUnmatchedParens -= 1;
                    }
                    numbers[i] = evaluateInput(
                        operations.GetRange(i + 1, matchingParenIndex - i - 1),
                        numbers.GetRange(i + 1, matchingParenIndex - i - 1));
                    operations[i] = ' ';
                    numbers.RemoveRange(i + 1, matchingParenIndex - i);
                    operations.RemoveRange(i + 1, matchingParenIndex - i);
                }
                else
                    ++i;
            }
            for (int i = 0; i < operations.Count;)
            {
                if (i - 1 < 0 || numbers[i - 1] == null)
                {
                    if (operations[i] == '+')
                    {
                        numbers.RemoveAt(i);
                        operations.RemoveAt(i);
                    }
                    else if (operations[i] == '-')
                    {
                        numbers[i] = new Math.Fraction(-1, 1);
                        numbers.Insert(i + 1, null);
                        operations[i] = ' ';
                        operations.Insert(i + 1, '*');
                        ++i;
                    }
                    else
                        ++i;
                }
                else
                    ++i;
            }
            for (int i = 0; i < operations.Count;)
            {
                if (operations[i] == '^')
                {
                    numbers[i - 1] = Math.exponentiate(numbers[i - 1], numbers[i + 1]);
                    numbers.RemoveRange(i, 2);
                    operations.RemoveRange(i, 2);
                }
                else
                    ++i;
            }
            for (int i = 0; i < operations.Count;)
            {
                if (operations[i] == '*')
                {
                    numbers[i - 1] *= numbers[i + 1];
                    numbers.RemoveRange(i, 2);
                    operations.RemoveRange(i, 2);
                }
                else if (operations[i] == '/')
                {
                    numbers[i - 1] /= numbers[i + 1];
                    numbers.RemoveRange(i, 2);
                    operations.RemoveRange(i, 2);
                }
                else
                    ++i;
            }
            for (int i = 0; i < operations.Count;)
            {
                if (operations[i] == '+')
                {
                    numbers[i - 1] += numbers[i + 1];
                    numbers.RemoveRange(i, 2);
                    operations.RemoveRange(i, 2);
                }
                else if (operations[i] == '-')
                {
                    numbers[i - 1] -= numbers[i + 1];
                    numbers.RemoveRange(i, 2);
                    operations.RemoveRange(i, 2);
                }
                else
                    ++i;
            }
            return numbers[0];
        }
        static void Main(string[] args)
        {
            while (true)
            {
                string input = Console.ReadLine();
                if (input[0] == 'q')
                    return;
                List<char> operations = new List<char>();
                List<Math.Number> numbers = new List<Math.Number>();
                StringBuilder numberCollector = new StringBuilder();
                bool lastCharWasDigit = false;
                foreach (char c in input)
                {
                    if (lastCharWasDigit)
                    {
                        if (char.IsDigit(c))
                            numberCollector.Append(c);
                        else
                        {
                            if (c == 'i')
                            {
                                Int32 number = Int32.Parse(numberCollector.ToString());
                                if (number == 0)
                                    numbers.Add(new Math.Fraction(0, 1));
                                else
                                    numbers.Add(new Math.ComplexNumber(0, number));
                                operations.Add(' ');
                            }
                            else
                            {
                                numbers.Add(new Math.Fraction(
                                    Int32.Parse(numberCollector.ToString()), 1));
                                operations.Add(' ');
                                operations.Add(c);
                                numbers.Add(null);
                            }
                            lastCharWasDigit = false;
                        }
                    }
                    else if (char.IsDigit(c))
                    {
                        numberCollector = new StringBuilder(c.ToString());
                        lastCharWasDigit = true;
                    }
                    else
                    {
                        if (c == 'i')
                        {
                            numbers.Add(new Math.ComplexNumber(0, 1));
                            operations.Add(' ');
                        }
                        else
                        {
                            operations.Add(c);
                            numbers.Add(null);
                        }
                    }
                }
                if (lastCharWasDigit)
                {
                    operations.Add(' ');
                    numbers.Add(new Math.Fraction(Int32.Parse(numberCollector.ToString()), 1));
                }
                for (int i = 0; i < operations.Count;)
                {
                    if (operations[i] == '(' && i > 0 && numbers[i - 1] != null)
                    {
                        operations.Insert(i, '*');
                        numbers.Insert(i, null);
                        i += 2;
                    }
                    else if ((operations[i] == ')' || numbers[i] is Math.ComplexNumber) &&
                        i + 1 < operations.Count && (numbers[i + 1] != null ||
                        operations[i + 1] == '('))
                    {
                        operations.Insert(i + 1, '*');
                        numbers.Insert(i + 1, null);
                        i += 2;
                    }
                    else
                        ++i;
                }
                try
                {
                    Console.WriteLine('=' + evaluateInput(operations, numbers).ToString());
                }
                catch (DivideByZeroException)
                {
                    Console.WriteLine("Division by 0 error.");
                }
            }
        }
    }
}
