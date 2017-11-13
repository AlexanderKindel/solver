using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleSolver
{
    class Math
    {
        public abstract class Number : IComparable, IEquatable<Number>
        {
            protected abstract Number add(Number number);
            protected abstract Number multiply(Number number);
            public abstract Number negative();
            public abstract Number reciprocal();
            public abstract Number exponentiate(Number number);
            public Number intExponentiate(int exponent)
            {
                Number output = new Fraction(1, 1);
                if (exponent < 0)
                {
                    for (int i = 0; i > exponent; --i)
                        output *= this;
                    return output.reciprocal();
                }
                for (int i = 0; i < exponent; ++i)
                    output *= this;
                return output;
            }
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
            public abstract int CompareTo(Object obj);
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public bool Equals(Number number)
            {
                return CompareTo(number) == 0;
            }
            public abstract override int GetHashCode();
        }
        public abstract class Term : Number
        {
            public static Term createTerm(Rational coefficient, List<Factor> factors)
            {//To be used instead of the Product constructor if it's unkown whether the given input
             //could be represented by a simpler type.
                if (factors.Count == 0)
                    return coefficient;
                if (coefficient is Fraction)
                {
                    Fraction coefficientFraction = (Fraction)coefficient;
                    if (coefficientFraction.Numerator == 0)
                        return coefficientFraction;
                    if (factors.Count == 1 && coefficientFraction.Numerator ==
                        coefficientFraction.Denominator)
                        return factors[0];
                }
                return new Product(coefficient, factors);
            }
        }
        public abstract class Factor : Term
        {
        }
        public abstract class Rational : Factor
        {
            static public Number exponentiateRational(Rational expBase, Fraction exponent,
                int rationalizer)
            {
                for (int i = 0; i < exponent.Denominator; ++i)
                    expBase = (Rational)(expBase * new Fraction(rationalizer, 1));
                Dictionary<Rational, int> radicandFactors = getFactorization(expBase);
                List<int> exponentDivisors = new List<int>();
                int divisor = 1;
                while (divisor <= exponent.Denominator / 2)
                {
                    if (exponent.Denominator % divisor == 0)
                        exponentDivisors.Add(divisor);
                    ++divisor;
                }
                exponentDivisors.Add(exponent.Denominator);
                Dictionary<int, Rational> termComponents = new Dictionary<int, Rational>();
                termComponents.Add(1, new Fraction(1, rationalizer));
                termComponents.Add(exponent.Denominator, new Fraction(1, 1));
                foreach (Rational factor in radicandFactors.Keys)
                    for (int i = exponentDivisors.Count - 1; i >= 0; --i)
                        if (exponentDivisors[i] <= radicandFactors[factor])
                        {
                            int index = exponent.Denominator / exponentDivisors[i];
                            if (!termComponents.ContainsKey(index))
                                termComponents.Add(index, new Fraction(1, 1));
                            for (int j = 0;
                                j < radicandFactors[factor] / exponentDivisors[i]; ++j)
                                termComponents[index] =
                                    (Rational)(termComponents[index] * factor);
                            for (int j = 0;
                                j < radicandFactors[factor] % exponentDivisors[i]; ++j)
                                termComponents[exponent.Denominator] =
                                    (Rational)(termComponents[exponent.Denominator] * factor);
                            break;
                        }
                Dictionary<int, Rational> termComponentstoPower =
                    new Dictionary<int, Rational>();
                foreach (int index in termComponents.Keys)
                {
                    termComponentstoPower.Add(index, new Fraction(1, 1));
                    for (int i = 0; i < exponent.Numerator; ++i)
                        termComponentstoPower[index] =
                            (Rational)(termComponentstoPower[index] * termComponents[index]);
                }
                Rational coefficient = termComponentstoPower[1];
                termComponentstoPower.Remove(1);
                if (termComponentstoPower[exponent.Denominator] is Fraction)
                {
                    Fraction radicandFraction =
                        (Fraction)termComponentstoPower[exponent.Denominator];
                    if (radicandFraction.Numerator == radicandFraction.Denominator)
                        termComponentstoPower.Remove(exponent.Denominator);
                    else if (radicandFraction.Numerator == -radicandFraction.Denominator)
                    {
                        if (exponent.Denominator % 2 == 1)
                        {
                            coefficient = (Rational)coefficient.negative();
                            termComponentstoPower.Remove(exponent.Denominator);
                        }
                        else if (exponent.Denominator % 4 == 0)
                        {
                            if (termComponentstoPower.ContainsKey(2))
                                termComponentstoPower[2] = (Rational)(termComponentstoPower[2] *
                                    new ComplexNumber(0, 1));
                            else
                                termComponentstoPower.Add(2, new ComplexNumber(0, 1));
                            termComponentstoPower.Remove(exponent.Denominator);
                        }
                        else if (exponent.Denominator % 2 == 0)
                        {
                            coefficient = (Rational)(coefficient * new ComplexNumber(0, 1));
                            termComponentstoPower.Remove(exponent.Denominator);
                        }
                    }
                }
                else
                {
                    ComplexNumber complexRadicand =
                        (ComplexNumber)termComponentstoPower[exponent.Denominator];
                    if (complexRadicand.Real.Numerator == 0)
                    {
                        if (complexRadicand.Imaginary.Numerator ==
                            complexRadicand.Imaginary.Denominator)
                        {
                            if (exponent.Denominator % 4 == 3)
                            {
                                coefficient =
                                    (Rational)(coefficient * new ComplexNumber(0, -1));
                                termComponentstoPower.Remove(exponent.Denominator);
                            }
                            else if (exponent.Denominator % 4 == 1)
                            {
                                coefficient = (Rational)(coefficient * new ComplexNumber(0, 1));
                                termComponentstoPower.Remove(exponent.Denominator);
                            }
                        }
                        else if (complexRadicand.Imaginary.Numerator ==
                            -complexRadicand.Imaginary.Denominator)
                        {
                            if (exponent.Denominator % 4 == 3)
                            {
                                coefficient = (Rational)(coefficient * new ComplexNumber(0, 1));
                                termComponentstoPower.Remove(exponent.Denominator);
                            }
                            else if (exponent.Denominator % 4 == 1)
                            {
                                coefficient =
                                    (Rational)(coefficient * new ComplexNumber(0, -1));
                                termComponentstoPower.Remove(exponent.Denominator);
                            }
                        }
                    }
                }
                List<Factor> factors = new List<Factor>();
                int largestIndex = 1;
                foreach (int index in termComponentstoPower.Keys)
                {
                    factors.Add(new Exponentiation(termComponentstoPower[index],
                        new Fraction(1, index)));
                    if (index > largestIndex)
                        largestIndex = index;
                }
                if (returnAllRoots && largestIndex < exponent.Denominator)
                {
                    returnAllRoots = false;
                    Number cosine = cos(new Fraction(1, exponent.Denominator));
                    Number rootOfUnity = cosine + (new Fraction(1, 1) -
                        cosine * cosine).exponentiate(new Fraction(1, 2)) *
                        new ComplexNumber(0, -1);
                    List<Number> roots =
                        new List<Number> { createTerm(coefficient, factors) };
                    for (int i = 1; i < exponent.Denominator / largestIndex; ++i)
                        roots.Add(roots[i - 1] * rootOfUnity);
                    returnAllRoots = true;
                    return new NumberList(roots);
                }
                return createTerm(coefficient, factors);
            }
            public abstract bool isGaussian();
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
            public override Number negative()
            {
                return new Fraction(-Numerator, Denominator);
            }
            public override Number reciprocal()
            {
                return new Fraction(Denominator, Numerator);
            }
            public override Number exponentiate(Number number)
            {
                if (Numerator == Denominator || Numerator == 0)
                    return this;
                if (number is Fraction)
                {
                    Fraction exponent = (Fraction)number;
                    if (exponent.Denominator == 1)
                        return intExponentiate(exponent.Numerator);
                    if (exponent.Numerator < 0)
                        return exponentiateRational((Rational)reciprocal(),
                            (Fraction)exponent.negative(), Numerator);
                    return exponentiateRational(this, exponent, Denominator);
                }
                return new Exponentiation(this, number);
            }
            public override bool isGaussian()
            {
                return Denominator == 1;
            }
            public override int CompareTo(Object obj)
            {
                int comparison = getTypeIndex(this) - getTypeIndex((Number)obj);
                if (comparison != 0)
                    return comparison;
                Fraction fraction = (Fraction)obj;
                return Numerator * fraction.Denominator - Denominator * fraction.Numerator;
            }
            public override int GetHashCode()
            {
                if (Numerator == 0)
                    return 0;
                if (Numerator > 0)
                    return (Numerator - 2 + Denominator) * (Numerator - 1 + Denominator) / 2 +
                        Numerator;
                return (Numerator + 2 - Denominator) * (Numerator + 1 - Denominator) / 2 +
                    Numerator;
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
            public override Number exponentiate(Number number)
            {
                if (number is Fraction)
                {
                    Fraction exponent = (Fraction)number;
                    ComplexNumber expBase = this;
                    if (exponent.Denominator == 1)
                        return intExponentiate(exponent.Numerator);
                    if (exponent.Numerator < 0)
                    {
                        exponent = (Fraction)(exponent.negative());
                        expBase = (ComplexNumber)expBase.reciprocal();
                    }
                    int rationalizer = getLCM(Real.Denominator, Imaginary.Denominator);
                    return exponentiateRational(expBase, exponent, rationalizer);
                }
                return new Exponentiation(this, number);
            }
            public override bool isGaussian()
            {
                return Real.Denominator == 1 && Imaginary.Denominator == 1;
            }
            public override int CompareTo(Object obj)
            {
                int comparison = getTypeIndex(this) - getTypeIndex((Number)obj);
                if (comparison != 0)
                    return comparison;
                ComplexNumber complexNumber = (ComplexNumber)obj;
                comparison = Real.CompareTo(complexNumber.Real);
                if (comparison != 0)
                    return comparison;
                return Imaginary.CompareTo(complexNumber.Imaginary);
            }
            public override int GetHashCode()
            {
                return Real.GetHashCode() ^ Imaginary.GetHashCode();
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
        public enum Constant
        {
            TAU
        }
        public class Transcendental : Factor
        {
            public Constant Value { get; }
            public Transcendental(Constant value)
            {
                Value = value;
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return this;
                    return new Sum(new List<Term> { fraction, this });
                }
                if (number is ComplexNumber)
                    return new Sum(new List<Term> { (ComplexNumber)number, this });
                if (number is Transcendental)
                {
                    Transcendental transcendental = (Transcendental)number;
                    if (transcendental.Value == Value)
                        return new Product(new Fraction(2, 1), new List<Factor> { this });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return fraction;
                    if (fraction.Numerator == fraction.Denominator)
                        return this;
                    return new Product(fraction, new List<Factor> { this });
                }
                if (number is ComplexNumber)
                    return new Product((ComplexNumber)number, new List<Factor> { this });
                if (number is Transcendental)
                {
                    Transcendental transcendental = (Transcendental)number;
                    if (transcendental.Value == Value)
                        return new Exponentiation(this, new Fraction(2, 1));
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Product(new Fraction(-1, 1), new List<Factor> { this });
            }
            public override Number reciprocal()
            {
                return new Exponentiation(this, new Fraction(-1, 1));
            }
            public override Number exponentiate(Number number)
            {
                if (number is Fraction)
                {
                    Fraction exponent = (Fraction)number;
                    if (exponent.Denominator == 1)
                        return intExponentiate(exponent.Numerator);
                }
                return new Exponentiation(this, number);
            }
            public override int CompareTo(object obj)
            {
                if (obj is Transcendental)
                {
                    Transcendental transcendental = (Transcendental)obj;
                    return Value - transcendental.Value;
                }
                return getTypeIndex(this) - getTypeIndex((Number)obj);
            }
            public override int GetHashCode()
            {
                return (int)Value;
            }
            public override string ToString()
            {
                return "tau";
            }
        }
        public class Cosine : Factor
        {
            public Number Argument { get; }
            public Cosine(Number argument)
            {
                Argument = argument;
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return this;
                }
                if (number is Cosine)
                {
                    Cosine cosine = (Cosine)number;
                    if (cosine.Argument == Argument)
                        return new Product(new Fraction(2, 1), new List<Factor> { this });
                }
                if (number is Term)
                    return new Sum(new List<Term> { this, (Term)number });
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return fraction;
                    if (fraction.Numerator == fraction.Denominator)
                        return this;
                    return new Product(fraction, new List<Factor> { this });
                }
                if (number is ComplexNumber)
                    return new Product((ComplexNumber)number, new List<Factor> { this });
                if (number is Transcendental)
                    return new Product(new Fraction(1, 1),
                        new List<Factor> { (Transcendental)number, this });
                if (number is Cosine)
                {
                    Cosine cosine = (Cosine)number;
                    if (cosine.Argument.Equals(Argument))
                        return new Exponentiation(this, new Fraction(2, 1));
                    return new Product(new Fraction(1, 1), new List<Factor> { this, cosine });
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Product(new Fraction(-1, 1), new List<Factor> { this });
            }
            public override Number reciprocal()
            {
                return new Exponentiation(this, new Fraction(-1, 1));
            }
            public override Number exponentiate(Number number)
            {
                if (number is Fraction)
                {
                    Fraction exponent = (Fraction)number;
                    if (exponent.Denominator == 1)
                        return intExponentiate(exponent.Numerator);
                }
                return new Exponentiation(this, number);
            }
            public override int CompareTo(object obj)
            {
                if (obj is Cosine)
                {
                    Cosine cosine = (Cosine)obj;
                    return Argument.CompareTo(cosine.Argument);
                }
                return getTypeIndex(this) - getTypeIndex((Number)obj);
            }
            public override int GetHashCode()
            {
                return Argument.GetHashCode();
            }
            public override string ToString()
            {
                return "cos(" + Argument.ToString() + ')';
            }
        }
        public class Exponentiation : Factor
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
                    return new Sum(new List<Term> { this, fraction });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (CompareTo(exponentiation) == 0)
                        return new Product(new Fraction(2, 1), new List<Factor> { this });
                    if (CompareTo(exponentiation.negative()) == 0)
                        return new Fraction(0, 1);
                    return new Sum(new List<Term> { this, exponentiation });
                }
                if (number is Factor)
                    return new Sum(new List<Term> { this, (Factor)number });
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == fraction.Denominator)
                        return this;
                    return new Product(fraction, new List<Factor> { this });
                }
                if (number is ComplexNumber)
                    return new Product((ComplexNumber)number, new List<Factor> { this });
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Base.Equals(exponentiation.Base))
                        return Base.exponentiate(Exponent + exponentiation.Exponent);
                    if (Exponent.Equals(exponentiation.Exponent))
                    {
                        Number outputBase = Base * exponentiation.Base;
                        if (outputBase is Exponentiation)
                        {
                            Exponentiation outputExponentiation = (Exponentiation)outputBase;
                            return new Exponentiation(outputExponentiation.Base,
                                Exponent * Exponent);
                        }
                        return (Base * exponentiation.Base).exponentiate(Exponent);
                    }
                    return new Product(new Fraction(1, 1),
                        new List<Factor> { this, exponentiation });
                }
                if (number is Factor)
                    return new Product(new Fraction(1, 1),
                        new List<Factor> { (Factor)number, this });
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
            public override Number exponentiate(Number number)
            {
                return Base.exponentiate(Exponent * number);
            }
            public override int CompareTo(Object obj)
            {
                int comparison = getTypeIndex(this) - getTypeIndex((Number)obj);
                if (comparison != 0)
                    return comparison;
                Exponentiation exponentiation = (Exponentiation)obj;
                comparison = Base.CompareTo(exponentiation.Base);
                if (comparison != 0)
                    return comparison;
                return Exponent.CompareTo(exponentiation.Exponent);
            }
            public override int GetHashCode()
            {
                return Base.GetHashCode() ^ Exponent.GetHashCode();
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                if (Base is Fraction)
                {
                    Fraction baseFraction = (Fraction)Base;
                    if (baseFraction.Denominator == 1 && baseFraction.Numerator > 0)
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
                    if (exponentFraction.Denominator == 1 && exponentFraction.Numerator > 0)
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
            public Rational Coefficient { get; }
            public List<Factor> Factors { get; }
            public Product(Rational coefficient, List<Factor> factors)
            {
                Coefficient = coefficient;
                Factors = factors;
                Factors.Sort();
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 0)
                        return this;
                }
                if (number is Factor)
                    return new Sum(new List<Term> { this, (Factor)number });
                if (number is Product)
                {
                    Product product = (Product)number;
                    if (Factors.Count == product.Factors.Count)
                    {
                        for (int i = 0; i < Factors.Count; ++i)
                            if (Factors[i].CompareTo(product.Factors[i]) != 0)
                                return new Sum(new List<Term> { this, product });
                        return createTerm((Rational)(Coefficient + product.Coefficient), Factors);
                    }
                    return new Sum(new List<Term> { this, product });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Rational)
                    return createTerm((Rational)(Coefficient * number), Factors);
                if (number is Factor)
                {
                    List<Factor> factors = new List<Factor>(Factors);
                    List<Factor> productFactors = new List<Factor>(Factors);
                    for (int i = 0; i < factors.Count; ++i)
                    {
                        Number product = number * factors[i];
                        productFactors.Add((Factor)number);
                        if (product != new Product(Coefficient, productFactors)) 
                        {
                            factors.RemoveAt(i);
                            return createTerm(Coefficient, factors) * product;
                        }
                    }
                    factors.Add((Factor)number);
                    return createTerm(Coefficient, factors);
                }
                if (number is Product)
                {
                    Product product = (Product)number;
                    Number output = new Product(Coefficient, Factors);
                    output = output * product.Coefficient;
                    for (int i = 0; i < product.Factors.Count; ++i) 
                        output = output * product.Factors[i];
                    return output;
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Product((Rational)Coefficient.negative(), Factors);
            }
            public override Number reciprocal()
            {
                List<Factor> factors = new List<Factor>();
                foreach (Factor factor in Factors)
                    factors.Add((Factor)factor.reciprocal());
                return new Product((Rational)Coefficient.reciprocal(), factors);
            }
            public override Number exponentiate(Number number)
            {
                Number output = Coefficient;
                foreach (Factor factor in Factors)
                    output *= factor.exponentiate(number);
                return output;
            }
            public override int CompareTo(object obj)
            {
                int comparison = getTypeIndex(this) - getTypeIndex((Number)obj);
                if (comparison != 0)
                    return comparison;
                Product product = (Product)obj;
                comparison = Coefficient.CompareTo(product.Coefficient);
                if (comparison != 0)
                    return comparison;
                comparison = Factors.Count - product.Factors.Count;
                if (comparison != 0)
                    return comparison;
                for (int i = 0; i < Factors.Count; ++i)
                {
                    comparison = Factors[i].CompareTo(product.Factors[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override int GetHashCode()
            {
                int output = Coefficient.GetHashCode();
                foreach (Factor factor in Factors)
                    output = output | factor.GetHashCode();
                return output;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder(Factors[0].ToString());
                for (int i = 1; i < Factors.Count; ++i)
                {
                    String exponentiationString = Factors[i].ToString();
                    if (output[output.Length - 1] != ')' && exponentiationString[0] != '(')
                        output.Append('*');
                    output.Append(exponentiationString);
                }
                StringBuilder coefficientString = new StringBuilder();
                if (Coefficient is Fraction)
                {
                    Fraction coefficientFraction = (Fraction)Coefficient;
                    if (coefficientFraction.Numerator != 1)
                    {
                        if (coefficientFraction.Numerator == -1)
                            coefficientString.Append('-');
                        else
                        {
                            coefficientString.Append(coefficientFraction.Numerator);
                            if (output[0] != '(')
                                coefficientString.Append('*');
                        }
                    }
                    if (coefficientFraction.Denominator != 1)
                        output.Append('/' + coefficientFraction.Denominator.ToString());
                }
                else
                {
                    ComplexNumber complexCoefficient = (ComplexNumber)Coefficient;
                    if (complexCoefficient.Real.Numerator == 0)
                    {
                        if (complexCoefficient.Imaginary.Numerator == 1)
                            coefficientString.Append('i');
                        else
                        {
                            if (complexCoefficient.Imaginary.Numerator == -1)
                                coefficientString.Append("-i");
                            else
                                coefficientString.Append(
                                    complexCoefficient.Imaginary.Numerator.ToString() + 'i');
                        }
                        if (output[0] != '(' && coefficientString.Length > 0)
                            coefficientString.Append('*');
                        if (complexCoefficient.Imaginary.Denominator != 1)
                            output.Append('/' + complexCoefficient.Imaginary.Denominator.ToString());
                    }
                    else
                        output.Insert(0, '(' + Coefficient.ToString() + ')');
                }
                return output.Insert(0, coefficientString).ToString();
            }
        }
        public class Sum : Number
        {
            public List<Term> Terms { get; }
            public Sum(List<Term> terms)
            {
                Terms = terms;
                Terms.Sort();
            }
            protected override Number add(Number number)
            {
                List<Term> terms = new List<Term>(Terms);
                if (number is Sum)
                {
                    Sum sum = (Sum)number;
                    Number output = new Sum(Terms);
                    foreach (Term term in sum.Terms)
                        output += term;
                }
                else
                {
                    for (int i = 0; i < Terms.Count; ++i)
                    {
                        Number sum = Terms[i] + number;
                        if (!(sum is Sum))
                        {
                            Terms.RemoveAt(i);
                            if (Terms.Count == 0)
                                return sum;
                            return new Sum(Terms) + sum;
                        }
                    }
                    terms.Add((Term)number);
                }
                return new Sum(terms);
            }
            protected override Number multiply(Number number)
            {
                Number output = new Fraction(0, 1);
                List<Term> terms = new List<Term>();
                if (number is Sum)
                {
                    Sum multiplier = (Sum)number;
                    foreach (Term multiplicandTerm in Terms)
                        foreach (Term multiplierTerm in multiplier.Terms)
                            output = output + multiplicandTerm * multiplierTerm;
                }
                else
                    foreach (Term term in Terms)
                        output = output + term * number;
                return output;
            }
            public override Number negative()
            {
                List<Term> terms = new List<Term>();
                foreach (Term term in Terms)
                    terms.Add((Term)term.negative());
                return new Sum(terms);
            }
            public override Number reciprocal()
            {
                return new Exponentiation(this, new Fraction(-1, 1));
            }
            public override Number exponentiate(Number number)
            {
                if (number is Fraction)
                {
                    Fraction exponent = (Fraction)number;
                    if (exponent.Denominator == 1)
                        return intExponentiate(exponent.Numerator);
                }
                return new Exponentiation(this, number);
            }
            public override int CompareTo(object obj)
            {
                int comparison = getTypeIndex(this) - getTypeIndex((Number)obj);
                if (comparison != 0)
                    return comparison;
                Sum expression = (Sum)obj;
                comparison = Terms.Count - expression.Terms.Count;
                if (comparison != 0)
                    return comparison;
                for (int i = 0; i < Terms.Count; ++i)
                {
                    comparison = Terms[i].CompareTo(expression.Terms[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override int GetHashCode()
            {
                int output = 0;
                foreach (Term term in Terms)
                    output = output ^ term.GetHashCode();
                return output;
            }
            public override String ToString()
            {
                StringBuilder output = new StringBuilder(Terms[0].ToString());
                for (int i = 1; i < Terms.Count; ++i)
                {
                    String termString = Terms[i].ToString();
                    if (termString[0] != '-')
                        output.Append('+');
                    output.Append(termString);
                }
                return output.ToString();
            }
        }
        public class NumberList : Term
        {
            public List<Number> Numbers { get; }
            public NumberList(List<Number> numbers)
            {
                Numbers = numbers;
                numbers.Sort();
            }
            protected override Number add(Number number)
            {
                List<Number> outputList = new List<Number>();
                if (number is NumberList)
                {
                    NumberList list = (NumberList)number;
                    foreach (Number n in Numbers)
                        foreach (Number m in list.Numbers)
                        {
                            Number sum = n + m;
                            if (!outputList.Contains(sum))
                                outputList.Add(sum);
                        }
                }
                else
                    foreach (Number n in Numbers)
                        outputList.Add(n + number);
                return new NumberList(outputList);
            }
            protected override Number multiply(Number number)
            {
                List<Number> outputList = new List<Number>();
                if (number is NumberList)
                {
                    NumberList list = (NumberList)number;
                    foreach (Number n in Numbers)
                        foreach (Number m in list.Numbers)
                        {
                            Number product = n * m;
                            if (!outputList.Contains(product)) 
                                outputList.Add(product);
                        }
                }
                else
                {
                    if (number is Fraction)
                    {
                        Fraction fraction = (Fraction)number;
                        if (fraction.Numerator == 0)
                            return new Fraction(0, 1);
                    }
                    foreach (Number n in Numbers)
                    {
                        Number product = n * number;
                        if (outputList.Contains(product))
                            outputList.Add(product);
                    }
                }
                return new NumberList(outputList);
            }
            public override Number negative()
            {
                List<Number> output = new List<Number>();
                foreach (Number n in Numbers)
                    output.Add(n.negative());
                return new NumberList(output);
            }
            public override Number reciprocal()
            {
                List<Number> output = new List<Number>();
                foreach (Number n in Numbers)
                    output.Add(n.reciprocal());
                return new NumberList(output);
            }
            public override Number exponentiate(Number number)
            {
                List<Number> output = new List<Number>();
                foreach (Number n in Numbers)
                    output.Add(n.exponentiate(number));
                return new NumberList(output);
            }
            public override int CompareTo(object obj)
            {
                NumberList list = (NumberList)obj;
                int comparison = Numbers.Count - list.Numbers.Count;
                if (comparison != 0)
                    return comparison;
                for (int i = 0; i < Numbers.Count; ++i) 
                {
                    comparison = Numbers[i].CompareTo(list.Numbers[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override int GetHashCode()
            {
                int output = 0;
                foreach (Number number in Numbers)
                    output = output ^ number.GetHashCode();
                return output;
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
        protected static int getTypeIndex(Number number)
        {
            if (number is Fraction)
                return 1;
            if (number is ComplexNumber)
                return 2;
            if (number is Transcendental)
                return 3;
            if (number is Cosine)
                return 4;
            if (number is Exponentiation)
                return 5;
            if (number is Product)
                return 6;
            return 7;
        }
        static public Dictionary<Rational, int> getFactorization(Rational gaussian)
        {//The keys represent the factors and the values, their multiplicities.
            Dictionary<Rational, int> factors = new Dictionary<Rational, int>();
            int sumOfRealAndImaginary = 2;
            while (true)
            {
                int imaginaryPart = sumOfRealAndImaginary / 2;
                int realPart = sumOfRealAndImaginary - imaginaryPart;
                for (int i = 0; imaginaryPart - i >= 0;)
                {
                    int gaussianMagnitudeSquared;
                    Fraction fractionGaussian = null;
                    ComplexNumber complexGaussian = null;
                    if (gaussian is Fraction)
                    {
                        fractionGaussian = (Fraction)gaussian;
                        gaussianMagnitudeSquared = fractionGaussian.Numerator *
                            fractionGaussian.Numerator;
                    }
                    else
                    {
                        complexGaussian = (ComplexNumber)gaussian;
                        gaussianMagnitudeSquared = complexGaussian.Real.Numerator *
                            complexGaussian.Real.Numerator + complexGaussian.Imaginary.Numerator *
                            complexGaussian.Imaginary.Numerator;
                    }
                    if ((realPart + i) * (realPart + i) + (imaginaryPart - i) *
                        (imaginaryPart - i) > gaussianMagnitudeSquared)
                    {
                        if (gaussian is Fraction && fractionGaussian.Numerator ==
                            fractionGaussian.Denominator)
                            return factors;
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
                        if (quotient.isGaussian())
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
                    if (realPart + i == 0)
                    {
                        factor = new Fraction(0, imaginaryPart - i);
                        if (testFactor())
                            continue;
                    }
                    else
                    {
                        factor = new ComplexNumber(realPart + i, imaginaryPart - i);
                        if (testFactor())
                            continue;
                        factor = new ComplexNumber(-realPart - i, imaginaryPart - i);
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
        static public Number cos(Number argument)
        {//The argument is expressed in turns; 90 degrees is represented as 1/4.
            if (argument is Fraction)
            {
                Fraction fraction = (Fraction)argument;
                int denominator = fraction.Denominator;
                int numerator = fraction.Numerator % denominator;
                int multiplicityOfTwo = 0;
                while (denominator % 2 == 0)
                {
                    denominator /= 2;
                    ++multiplicityOfTwo;
                }
                Number cosine;
                switch (denominator)
                {
                    case 1:
                        cosine = new Fraction(1, 1);
                        break;
                    case 3:
                        cosine = new Fraction(-1, 2);
                        break;
                    case 5:
                        cosine = new Exponentiation(new Fraction(5, 1), new Fraction(1, 2)) *
                            new Fraction(1, 4) - new Fraction(1, 4);
                        break;
                    case 15:
                        cosine = new Sum(new List<Term> { new Fraction(1, 8), new Product(
                            new Fraction(1, 8), new List<Factor> { new Exponentiation(
                            new Fraction(5, 1), new Fraction(1, 2)) }), new Product(
                            new Fraction(1, 8), new List<Factor> { new Exponentiation(
                            new Fraction(30, 1) - new Exponentiation(new Fraction(5, 1),
                            new Fraction(1, 2)) * new Fraction(6, 1), new Fraction(1, 2)) }) });
                        break;
                    case 17:
                        cosine = (new Fraction(-1, 1) + new Exponentiation(new Fraction(17, 1),
                            new Fraction(1, 2)) + new Exponentiation(new Fraction(34, 1) -
                            new Exponentiation(new Fraction(17, 1), new Fraction(1, 2)) *
                            new Fraction(2, 1), new Fraction(1, 2)) + new Exponentiation(
                            new Fraction(17, 1) + new Exponentiation(new Fraction(17, 1),
                            new Fraction(1, 2)) * new Fraction(3, 1) - new Exponentiation(
                            new Fraction(34, 1) - new Exponentiation(new Fraction(17, 1),
                            new Fraction(1, 2)) * new Fraction(2, 1), new Fraction(1, 2)) -
                            new Exponentiation(new Fraction(34, 1) + new Exponentiation(
                            new Fraction(17, 1), new Fraction(1, 2)) * new Fraction(2, 1),
                            new Fraction(1, 2)) * new Fraction(2, 1), new Fraction(1, 2)) *
                            new Fraction(2, 1)) * new Fraction(1, 16);
                        break;
                    default:
                        return new Cosine(argument * new Transcendental(Constant.TAU));
                }
                for (int i = 0; i < multiplicityOfTwo; ++i)
                {
                    denominator *= 2;
                    if (4 < 3 * denominator && denominator < 4)
                        cosine = (cosine * new Fraction(1, 2) + new Fraction(1, 2)).exponentiate(
                        new Fraction(1, 2)).negative();
                    else
                        cosine = (cosine * new Fraction(1, 2) + new Fraction(1, 2)).exponentiate(
                        new Fraction(1, 2));
                }
                Number s = new Fraction(1, 1);
                Number t = cosine;
                Number angleMultipleCosine = cosine;
                for (int i = 1; i < numerator; ++i)
                {
                    angleMultipleCosine = new Fraction(2, 1) * cosine * t - s;
                    s = t;
                    t = angleMultipleCosine;
                }
                return angleMultipleCosine;
            }
            return new Cosine(argument * new Transcendental(Constant.TAU));
        }
        static protected bool returnAllRoots = true;
    }
    class Solver
    {
        internal static Math.Number evaluateExpression(List<char> operations,
            List<Math.Number> numbers)
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
                    numbers[i] = evaluateExpression(
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
                    numbers[i - 1] = numbers[i - 1].exponentiate(numbers[i + 1]);
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
                    Console.Write("=\n" +
                        checked(evaluateExpression(operations, numbers).ToString() + "\n\n"));
                }
                catch (DivideByZeroException)
                {
                    Console.WriteLine("Division by 0 error.");
                }
            }
        }
    }
}
