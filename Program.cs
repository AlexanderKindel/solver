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
        }
        public class Fraction : Number
        {
            public int Numerator { get; private set; }
            public int Denominator { get; private set; }
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
            public override int CompareTo(Object obj)
            {
                Fraction fraction = (Fraction)obj;
                return Numerator * fraction.Denominator - Denominator * fraction.Numerator;
            }
            public override string ToString()
            {
                if (Denominator == 1)
                    return Numerator.ToString();
                return Numerator.ToString() + "/" + Denominator.ToString();
            }
        }
        public class ComplexNumber : Number
        {
            public Fraction Real { get; private set; }
            public Fraction Imaginary { get; private set; }
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
                    return new ComplexNumber((Fraction)(Real + complexNumber.Real),
                        (Fraction)(Imaginary + complexNumber.Imaginary));
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new ComplexNumber((Fraction)(Real * fraction),
                        (Fraction)(Imaginary * fraction));
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new ComplexNumber((Fraction)(Real * complexNumber.Real -
                        Imaginary * complexNumber.Imaginary), (Fraction)
                        (Real * complexNumber.Imaginary + Imaginary * complexNumber.Real));
                }
                return number * this;
            }
            public override Number negative()
            {
                return new ComplexNumber((Fraction)Real.negative(),
                    (Fraction)Imaginary.negative());
            }
            public override Number reciprocal()
            {
                Fraction denominator = (Fraction)(Real * Real + Imaginary * Imaginary);
                return new ComplexNumber((Fraction)(Real / denominator),
                    (Fraction)(Imaginary.negative() / denominator));
            }
            public override int CompareTo(Object obj)
            {
                ComplexNumber complexNumber = (ComplexNumber)obj;
                int comparison = Real.CompareTo(complexNumber.Real);
                if (comparison != 0) 
                    return comparison;
                return Imaginary.CompareTo(complexNumber.Imaginary);
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                if (Imaginary.Numerator != 0)
                {
                    if (Imaginary.Numerator > 0)
                    {
                        if (Real.Numerator != 0)
                            output.Append(Real.ToString() + '+');
                        if (Imaginary.Numerator != 1)
                            output.Append(Imaginary.Numerator);
                    }
                    else if (Imaginary.Numerator < 0)
                    {
                        if (Real.Numerator != 0)
                            output.Append(Real.ToString());
                        output.Append('-');
                        if (Imaginary.Numerator != -1)
                            output.Append(-Imaginary.Numerator);
                    }
                    output.Append('i');
                    if (Imaginary.Denominator != 1)
                    {
                        output.Append('/' + Imaginary.Denominator.ToString());
                    }
                }
                else if (Real.Numerator != 0)
                    output.Append(Real.ToString());
                else
                    return "0";
                return output.ToString();
            }
        }
        public class Exponentiation : Number
        {
            public Number Base { get; private set; }
            public Number Exponent { get; private set; }
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
                    return new Expression(new List<Term> { new Term(this),
                        new Term(new ComplexNumber(fraction, new Fraction(0, 1))) });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    if (complexNumber.Real.Numerator == 0 &&
                        complexNumber.Imaginary.Numerator == 0)
                        return this;
                    return new Expression(new List<Term> { new Term(this),
                        new Term(complexNumber) });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Exponent.Equals(exponentiation.Exponent))
                        return new Exponentiation(Base + exponentiation.Base, Exponent);
                    return new Expression(new List<Term> { new Term(this),
                        new Term(exponentiation) });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator == 1 && fraction.Denominator == 1)
                        return this;
                    return new Term(new ComplexNumber(fraction, new Fraction(0, 1)),
                        new List<Exponentiation> { this });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    if (complexNumber.Real.Numerator == complexNumber.Real.Denominator &&
                        complexNumber.Imaginary.Numerator == 0)
                        return this;
                    return new Term(complexNumber, new List<Exponentiation> { this });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Exponent.Equals(exponentiation.Exponent))
                        return new Exponentiation(Base * exponentiation.Base, Exponent);
                    if (Base.Equals(exponentiation.Base))
                        return new Exponentiation(Base, Exponent + exponentiation.Exponent);
                    return new Term(new ComplexNumber(1, 0),
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
                if (Base is Fraction)
                {
                    Fraction baseFraction = (Fraction)Base;
                    if (Exponent is Fraction)
                    {
                        Fraction exponentFraction = (Fraction)Exponent;
                        if (baseFraction.Denominator == 1)
                            if (exponentFraction.Denominator == 1)
                                return Base.ToString() + "^" + Exponent;
                            else
                                return Base.ToString() + "^(" + Exponent + ')';
                        else
                        {
                            if (exponentFraction.Denominator == 1)
                                return '(' + Base.ToString() + ")^" + Exponent;
                            else
                                return '(' + Base.ToString() + ")^(" + Exponent + ')';
                        }
                    }
                }
                return '(' + Base.ToString() + ")^(" + Exponent + ')';
            }
        }
        public class Term : Number
        {
            public ComplexNumber Scalar { get; private set; }
            public List<Exponentiation> Exponentiations { get; private set; }
            public Term(ComplexNumber complexNumber, List<Exponentiation> exponentiations)
            {
                Scalar = complexNumber;
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
            public Term(Exponentiation exponentiation)
            {
                Scalar = new ComplexNumber(1, 0);
                Exponentiations = new List<Exponentiation> { exponentiation };
            }
            public Term(ComplexNumber complexNumber)
            {
                Scalar = complexNumber;
                Exponentiations = new List<Exponentiation>();
            }
            protected override Number add(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Expression(new List<Term> { this,
                        new Term(new ComplexNumber(fraction, new Fraction(0, 1))) });
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Expression(new List<Term> { this, new Term(complexNumber) });
                }
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    return new Expression(new List<Term> { this, new Term(exponentiation) });
                }
                if (number is Term)
                {
                    Term term = (Term)number;
                    return new Expression(new List<Term> { this, term });
                }
                return number + this;
            }
            protected override Number multiply(Number number)
            {
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return new Term((ComplexNumber)(Scalar * fraction), Exponentiations);
                }
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return new Term((ComplexNumber)(Scalar * complexNumber), Exponentiations);
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
                            return new Term(Scalar, exponentiations);
                        }
                    }
                    exponentiations.Add(exponentiation);
                    return new Term(Scalar, exponentiations);
                }
                if (number is Term)
                {
                    Term term = (Term)number;
                    List<Exponentiation> exponentiations = Exponentiations;
                    for (int i = 0; i < term.Exponentiations.Count; ++i)
                    {
                        bool addExponentiation = true;
                        for (int j = 0; j < exponentiations.Count; ++j)
                        {
                            Number product = term.Exponentiations[i] * exponentiations[j];
                            if (product is Exponentiation)
                            {
                                exponentiations[i] = (Exponentiation)product;
                                addExponentiation = false;
                                break;
                            }
                        }
                        if (addExponentiation)
                            exponentiations.Add(term.Exponentiations[i]);
                    }
                    return new Term((ComplexNumber)(Scalar * term.Scalar), exponentiations);
                }
                return number * this;
            }
            public override Number negative()
            {
                return new Term((ComplexNumber)Scalar.negative(), Exponentiations);
            }
            public override Number reciprocal()
            {
                List<Exponentiation> exponentiations = new List<Exponentiation>();
                foreach (Exponentiation exponentation in Exponentiations)
                    exponentiations.Add((Exponentiation)exponentation.reciprocal());
                return new Term((ComplexNumber)Scalar.reciprocal(), exponentiations);
            }
            public override int CompareTo(object obj)
            {
                Term term = (Term)obj;
                int comparison = Scalar.CompareTo(term.Scalar);
                if (comparison != 0)
                    return comparison;
                int countComparison = Exponentiations.Count - term.Exponentiations.Count;
                if (countComparison != 0)
                    return countComparison;
                for (int i = 0; i < Exponentiations.Count; ++i)
                {
                    comparison = Exponentiations[i].CompareTo(term.Exponentiations[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                if (Scalar.Real.Numerator != Scalar.Real.Denominator ||
                    Scalar.Imaginary.Numerator != 0)
                    output.Append(Scalar);
                foreach (Exponentiation exponentiation in Exponentiations)
                    output.Append('(' + exponentiation.ToString() + ')');
                return output.ToString();
            }
        }
        public class Expression : Number
        {
            public List<Term> m_terms;
            public Expression(List<Term> terms)
            {
                m_terms = terms;
                m_terms.Sort();
                for (int i=0;i<m_terms.Count-1;)
                {
                    if (m_terms[i].Exponentiations.Count != m_terms[i + 1].Exponentiations.Count) 
                    {
                        ++i;
                        continue;
                    }
                    bool termsAreAddable = true;
                    for (int j = 0; j < m_terms[i].Exponentiations.Count; ++j)
                        if (m_terms[i].Exponentiations[j].CompareTo(
                            m_terms[i + 1].Exponentiations[j]) != 0) 
                        {
                            termsAreAddable = false;
                            break;
                        }
                    if (termsAreAddable)
                    {
                        m_terms[i] = new Term((ComplexNumber)(m_terms[i].Scalar +
                            m_terms[i + 1].Scalar), m_terms[i].Exponentiations);
                        m_terms.RemoveAt(i + 1);
                        continue;
                    }
                    ++i;
                }
            }
            protected override Number add(Number number)
            {
                List<Term> terms = m_terms;
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    if (fraction.Numerator != 0)
                        terms.Add(new Term(new ComplexNumber(fraction, new Fraction(0, 1))));
                }
                else if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    if (complexNumber.Real.Numerator != 0 ||
                        complexNumber.Imaginary.Numerator != 0)
                        terms.Add(new Term(complexNumber));
                }
                else if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    terms.Add(new Term(exponentiation));
                }
                else if (number is Term)
                {
                    Term term = (Term)number;
                    terms.Add(term);
                }
                else
                {
                    Expression expression = (Expression)number;
                    terms.AddRange(expression.m_terms);
                }
                return new Expression(terms);
            }
            protected override Number multiply(Number number)
            {
                List<Term> terms = new List<Term>();
                if (number is Fraction)
                {
                    Fraction multiplier = (Fraction)number;
                    foreach (Term term in m_terms)
                        terms.Add((Term)(term * multiplier));
                }
                else if (number is ComplexNumber)
                {
                    ComplexNumber multiplier = (ComplexNumber)number;
                    foreach (Term term in m_terms)
                        terms.Add((Term)(term * multiplier));
                }
                else if (number is Exponentiation)
                {
                    Exponentiation multiplier = (Exponentiation)number;
                    foreach (Term term in m_terms)
                        terms.Add((Term)(term * multiplier));
                }
                else if (number is Term)
                {
                    Term multiplier = (Term)number;
                    foreach (Term term in m_terms)
                        terms.Add((Term)(term * multiplier));
                }
                else
                {
                    Expression multiplier = (Expression)number;
                    foreach (Term multiplicandTerm in m_terms)
                        foreach (Term multiplierTerm in multiplier.m_terms)
                            terms.Add((Term)(multiplicandTerm * multiplierTerm));
                }
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
        static public Dictionary<ComplexNumber, int> getFactorization(ComplexNumber gaussian)
        {
            Dictionary<ComplexNumber, int> factors = new Dictionary<ComplexNumber, int>();
            int sumOfRealAndImaginary = 2;
            while (true)
            {
                int realPart = sumOfRealAndImaginary / 2;
                int imaginaryPart = sumOfRealAndImaginary - realPart;
                for (int i = 0; realPart - i >= 0;)
                {
                    if ((realPart - i) * (realPart - i) + (imaginaryPart + i) *
                        (imaginaryPart + i) > gaussian.Real.Numerator *
                        gaussian.Real.Numerator + gaussian.Imaginary.Numerator *
                        gaussian.Imaginary.Numerator)
                    {
                        if (factors.ContainsKey(gaussian))
                            factors[gaussian] += 1;
                        else
                            factors.Add(gaussian, 1);
                        return factors;
                    }
                    ComplexNumber factor =
                        new ComplexNumber(realPart - i, imaginaryPart + i);
                    ComplexNumber quotient = (ComplexNumber)(gaussian / factor);
                    if (quotient.Real.Denominator == 1 && quotient.Imaginary.Denominator == 1)
                    {
                        if (factors.ContainsKey(factor))
                            factors[factor] += 1;
                        else
                            factors.Add(factor, 1);
                        gaussian = (ComplexNumber)(gaussian / factor);
                        continue;
                    }
                    factor = new ComplexNumber(realPart - i, -imaginaryPart - i);
                    quotient = (ComplexNumber)(gaussian / factor);
                    if (quotient.Real.Denominator == 1 &&
                        quotient.Imaginary.Denominator == 1)
                    {
                        if (factors.ContainsKey(factor))
                            factors[factor] += 1;
                        else
                            factors.Add(factor, 1);
                        gaussian = (ComplexNumber)(gaussian / factor);
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
                    Dictionary<ComplexNumber, int> radicandFactors =
                        getFactorization(baseComplexNumber);
                    ComplexNumber root = new ComplexNumber(1, 0);
                    ComplexNumber radicand = new ComplexNumber(1, 0);
                    foreach (ComplexNumber factor in radicandFactors.Keys)
                    {
                        for (int i = 0; i < radicandFactors[factor] /
                            exponentFraction.Denominator; ++i)
                            root = (ComplexNumber)(root * factor);
                        for (int i = 0; i < radicandFactors[factor] %
                            exponentFraction.Denominator; ++i)
                            radicand = (ComplexNumber)(radicand * factor);
                    }
                    ComplexNumber rootToPower = new ComplexNumber(1, 0);
                    ComplexNumber radicandToPower = new ComplexNumber(1, 0);
                    for (int i = 0; i < exponentFraction.Numerator %
                        exponentFraction.Denominator; ++i)
                    {
                        rootToPower = (ComplexNumber)(rootToPower * root *
                            new Fraction(1, rationalizer));
                        radicandToPower = (ComplexNumber)(radicandToPower * radicand);
                    }
                    for (int i = 0; i < exponentFraction.Numerator /
                        exponentFraction.Denominator; ++i)
                        rootToPower = (ComplexNumber)(rootToPower * baseComplexNumber *
                            new Fraction(1, rationalizer));
                    if (radicand.Equals(new ComplexNumber(1, 0)))
                        return rootToPower;
                    if (radicand.Equals(new ComplexNumber(-1, 0)))
                        return rootToPower * new ComplexNumber(0, -1);
                    return new Exponentiation(radicandToPower,
                        new Fraction(1, exponentFraction.Denominator)) * rootToPower;
                }
            }
            if (exponent is ComplexNumber)
            {
                ComplexNumber exponentComplexNumber = (ComplexNumber)exponent;
                if (exponentComplexNumber.Imaginary.Numerator == 0)
                    return exponentiate(expBase, exponentComplexNumber.Real);
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
                        numbers[i + 1] = numbers[i + 1].negative();
                        numbers.RemoveAt(i);
                        operations.RemoveAt(i);
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
                                numbers.Add(new Math.ComplexNumber(0,
                                    Int32.Parse(numberCollector.ToString())));
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
                    Console.WriteLine(evaluateInput(operations, numbers));
                }
                catch (DivideByZeroException)
                {
                    Console.WriteLine("Division by 0 error.");
                }
            }
        }
    }
}