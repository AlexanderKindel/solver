using System;
using System.Collections.Generic;
using System.Text;

namespace Solver
{
    class Solver
    {
        static List<List<T>> generateCartesianProduct<T>(List<List<T>> sets)
        {
            List<List<T>> elements = new List<List<T>>();
            void generateElements(int setIndex, List<T> element)
            {
                if (setIndex < sets.Count)
                    foreach (T t in sets[setIndex])
                    {
                        List<T> enlargedElement = new List<T>(element);
                        enlargedElement.Add(t);
                        generateElements(setIndex + 1, enlargedElement);
                    }
                else
                    elements.Add(element);
            }
            generateElements(0, new List<T>());
            return elements;
        }
        abstract class Number : IComparable, IEquatable<Number>
        {
            protected abstract Number add(Number number);
            public static Number operator +(Number a, Number b)
            {
                return a.add(b);
            }
            public abstract Number negative();
            public static Number operator -(Number a, Number b)
            {
                return a.add(b.negative());
            }
            protected abstract Number multiply(Number number);
            public static Number operator *(Number a, Number b)
            {
                return a.multiply(b);
            }
            public abstract Number reciprocal();
            public static Number operator /(Number a, Number b)
            {
                return a.multiply(b.reciprocal());
            }
            Polynomial minimalPolynomial = null;
            public Polynomial MinimalPolynomial
            {
                get
                {
                    if (minimalPolynomial == null)
                        minimalPolynomial = calculateMinimalPolynomial();
                    return minimalPolynomial;
                }
            }
            protected abstract Polynomial calculateMinimalPolynomial();
            public abstract List<Number> getConjugates();
            protected List<Number> testConjugateCandidates(List<Number> conjugateCandidates)
            {
                List<Number> conjugates = new List<Number>();
                foreach (Number conjugate in conjugateCandidates)
                    if (conjugate.MinimalPolynomial.Equals(MinimalPolynomial))
                        conjugates.Add(conjugate);
                return conjugates;
            }
            public virtual Integer getDenominatorLCM()
            {
                return One;
            }
            public virtual Integer getGreatestIntegerFactor()
            {
                return One;
            }
            public virtual Number exponentiate(Number exponent)
            {
                if (exponent is Integer)
                {
                    Integer intExponent = (Integer)exponent;
                    if (intExponent.Sign < 0)
                        return exponentiate(exponent).reciprocal();
                    Number output = One;
                    Number baseToAPowerOfTwo = this;
                    while (intExponent.Sign > 0)
                    {
                        IntegerDivision division = intExponent.euclideanDivideBy(new Integer(2)); 
                        if (division.remainder == One)
                            output *= baseToAPowerOfTwo;
                        baseToAPowerOfTwo *= baseToAPowerOfTwo;
                        intExponent = division.quotient;
                    }
                    return output;
                }
                Integer exponentIntegerFactor = exponent.getGreatestIntegerFactor();
                if (exponentIntegerFactor > One)
                    return exponentiate(exponentIntegerFactor).exponentiate(
                        exponent / exponentIntegerFactor);
                if (exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)exponent;
                    Integer denominatorLCM = getDenominatorLCM();
                    if (denominatorLCM != One)
                    {
                        IntegerDivision division =
                            exponentFraction.Numerator.euclideanDivideBy(exponentFraction.Denominator);
                        return (this * denominatorLCM).exponentiate(exponent) *
                            denominatorLCM.exponentiate(Fraction.create(exponentFraction.Denominator -
                            division.remainder, exponentFraction.Denominator)) /
                            denominatorLCM.exponentiate(division.quotient + One);
                    }
                    Integer numeratorGCD = getGreatestIntegerFactor();
                    if (numeratorGCD != One)
                    {
                        Number factor = numeratorGCD.exponentiate(exponent);
                        if (!(factor is Exponentiation))
                            return (this / numeratorGCD).exponentiate(exponent) * factor;
                    }
                }
                return Exponentiation.create(this, exponent);
            }          
            protected int getTypeIndex()
            {
                if (this is Integer)
                    return 1;
                if (this is Fraction)
                    return 2;
                if (this is Exponentiation)
                    return 3;
                if (this is ComplexExponential)
                    return 4;
                if (this is Product)
                    return 5;
                if (this is Sum)
                    return 6;
                if (this is ComplexNumber)
                    return 7;
                return 8;
            }
            public virtual int CompareTo(object obj)
            {
                return getTypeIndex() - ((Number)obj).getTypeIndex();
            }
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public bool Equals(Number number)
            {
                return CompareTo(number) == 0;
            }
            public abstract override int GetHashCode();

            //Treats str as the string representation of a numerical constant, and returns the
            //string representation of the reference object multiplied by that constant.
            public abstract string insertString(string str);
            public abstract override string ToString();
        }
        abstract class Term : Number
        { }
        abstract class Factor : Term
        { }
        abstract class Rational : Factor
        {
            protected override Polynomial calculateMinimalPolynomial()
            {
                return new Polynomial(new List<Rational> { (Rational)negative(), One });
            }
            public override List<Number> getConjugates()
            {
                return new List<Number> { this };
            }
        }
        struct IntegerDivision
        {
            public Integer quotient;
            public Integer remainder;
        }
        static Integer Zero = new Integer(0);
        static Integer One = new Integer(1);
        class Integer : Rational
        {
            uint[] Values;
            public sbyte Sign { get; private set; }
            Integer(uint[] values, sbyte sign)
            {
                int lastNonzeroValueIndex = 0;
                for (int i = values.Length - 1; i >= 0; --i)
                    if (values[i] != 0)
                    {
                        lastNonzeroValueIndex = i;
                        break;
                    }
                Values = new uint[lastNonzeroValueIndex + 1];
                for (int i = 0; i <= lastNonzeroValueIndex; ++i)
                    Values[i] = values[i];
                if (Values.Length == 1 && Values[0] == 0)
                    Sign = 0;
                else
                    Sign = sign;
            }
            public Integer(int value)
            {
                if (value == 0)
                {
                    Sign = 0;
                    Values = new uint[] { (uint)value };
                }
                else if (value > 0)
                {
                    Sign = 1;
                    Values = new uint[] { (uint)value };
                }
                else
                {
                    Sign = -1;
                    Values = new uint[] { (uint)-value };
                }
            }
            public static Integer parse(string decimalString)
            {
                Integer output = Zero;
                Integer multiplier = One;
                Integer multiplierStepSize = new Integer(1000000000);
                int startIndex = decimalString.Length - 9;
                while (startIndex >= 0)
                {
                    output = output + multiplier *
                        new Integer(int.Parse(decimalString.Substring(startIndex, 9)));
                    multiplier = multiplier * multiplierStepSize;
                    startIndex -= 9;
                }
                if (startIndex > -9)
                    output = output + multiplier *
                        new Integer(int.Parse(decimalString.Substring(0, startIndex + 9)));
                return output;
            }
            public static Integer operator +(Integer a, Integer b)
            {
                uint[] aValues;
                uint[] bValues;
                uint[] takeTwosComplement(uint[] values)
                {
                    uint[] outputValues = new uint[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                        outputValues[i] = ~values[i];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        uint power = 1;
                        for (int j = 0; j < 32; ++j)
                        {
                            outputValues[i] ^= power;
                            if ((values[i] & power) != 0)
                                return outputValues;
                            power = power << 1;
                        }
                    }
                    return outputValues;
                }
                void padValueArrays(Integer shortInteger, Integer longInteger)
                {
                    aValues = new uint[longInteger.Values.Length];
                    shortInteger.Values.CopyTo(aValues, 0);
                    if (shortInteger.Sign < 0)
                        aValues = takeTwosComplement(aValues);
                    if (longInteger.Sign < 0)
                        bValues = takeTwosComplement(longInteger.Values);
                    else
                    bValues = longInteger.Values;
                }
                if (a.Values.Length < b.Values.Length)
                    padValueArrays(a, b);
                else
                    padValueArrays(b, a);
                uint[] sum = new uint[aValues.Length];
                bool remainder = false;
                for (int i = 0; i < aValues.Length; ++i)
                {
                    uint power = 1;
                    for (int j = 0; j < 32; ++j)
                    {
                        bool aDigitIsZero = (aValues[i] & power) == 0;
                        bool bDigitIsZero = (bValues[i] & power) == 0;
                        if (remainder)
                        {
                            if (aDigitIsZero)
                            {
                                if (bDigitIsZero)
                                {
                                    sum[i] |= power;
                                    remainder = false;
                                }
                            }
                            else if (!bDigitIsZero)
                                sum[i] |= power;
                        }
                        else if (!aDigitIsZero)
                        {
                            if (!bDigitIsZero)
                                remainder = true;
                            else
                                sum[i] |= power;
                        }
                        else if (!bDigitIsZero)
                            sum[i] |= power;
                        power = power << 1;
                    }
                }
                if (a.Sign < 0)
                {
                    if (b.Sign < 0)
                    {
                        if (!remainder)
                        {
                            uint[] sumWithRemainder = new uint[sum.Length + 1];
                            takeTwosComplement(sum).CopyTo(sumWithRemainder, 0);
                            sumWithRemainder[sum.Length] = 1;
                            return new Integer(sumWithRemainder, -1);
                        }
                        return new Integer(takeTwosComplement(sum), -1);
                    }
                    if (remainder)
                        return new Integer(sum, 1);
                    return new Integer(takeTwosComplement(sum), -1);
                }
                if (b.Sign < 0)
                {
                    if (remainder)
                        return new Integer(sum, 1);
                    return new Integer(takeTwosComplement(sum), -1);
                }
                if (remainder)
                {
                    uint[] sumWithRemainder = new uint[sum.Length + 1];
                    sum.CopyTo(sumWithRemainder, 0);
                    sumWithRemainder[sum.Length] = 1;                
                    return new Integer(sumWithRemainder, 1);
                }
                return new Integer(sum, 1);
            }               
            public static Integer operator ++(Integer a)
            {
                return a + One;
            }
            protected override Number add(Number number)
            {
                if (number is Integer)
                    return (Integer)number + this;
                return number + this;
            }
            public static Integer operator -(Integer a)
            {
                return new Integer(a.Values, (sbyte)-a.Sign);
            }
            public static Integer operator -(Integer a, Integer b)
            {
                return a + -b;
            }
            public static Integer operator --(Integer a)
            {
                return a + new Integer(-1);
            }
            public override Number negative()
            {
                return new Integer(Values, (sbyte)-Sign);
            }        
            static uint[] shiftValuesLeft(uint[] values, int valuePlaces, int digitPlaces)
            {
                uint[] shiftedValues;
                if (digitPlaces == 0)
                {
                    shiftedValues = new uint[values.Length + valuePlaces];
                    values.CopyTo(shiftedValues, valuePlaces);
                }
                else if (digitPlaces > 0)
                {
                    shiftedValues = new uint[values.Length + valuePlaces + 1];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        shiftedValues[valuePlaces + i] += values[i] << digitPlaces;
                        shiftedValues[valuePlaces + i + 1] += values[i] >> 32 - digitPlaces;
                    }
                }
                else
                {
                    shiftedValues = new uint[values.Length + valuePlaces];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        if (valuePlaces + i - 1 >= 0)
                            shiftedValues[valuePlaces + i - 1] += values[i] << 32 + digitPlaces;
                        shiftedValues[valuePlaces + i] += values[i] >> -digitPlaces;
                    }
                }
                return shiftedValues;
            }
            Integer shiftLeft(int valuePlaces, int digitPlaces)
            {
                return new Integer(shiftValuesLeft(Values, valuePlaces, digitPlaces), Sign);
            }
            public static Integer operator *(Integer a, Integer b)
            {
                Integer product = Zero;
                for (int i = 0; i < a.Values.Length; ++i)
                {
                    uint power = 1;
                    for (int j = 0; j < 32; ++j)
                    {
                        if ((a.Values[i] & power) != 0)
                            product = product + b.shiftLeft(i, j);
                        power = power << 1;
                    }
                }
                product.Sign = (sbyte)(b.Sign * a.Sign);
                return product;
            }
            protected override Number multiply(Number number)
            {
                if (number is Integer)
                    return (Integer)number * this;
                return number * this;
            }
            public override Number reciprocal()
            {
                return Fraction.create(One, this);
            }
            public IntegerDivision euclideanDivideBy(Integer divisor)
            {
                if (divisor.Sign == 0)
                    throw new DivideByZeroException();            
                int divisorLeadingDigitPlace = 0;
                uint power = 0b10000000000000000000000000000000;
                for (int i = 32; i > 0; --i)
                {
                    if ((divisor.Values[divisor.Values.Length - 1] & power) != 0)
                    {
                        divisorLeadingDigitPlace = i;
                        break;
                    }
                    power = power >> 1;
                }
                Integer positiveDivisor = new Integer(divisor.Values, 1);
                Integer remainder = new Integer(Values, 1);
                uint[] quotient = new uint[Values.Length];
                void calculateValue(int valuePlace, int stoppingDigitPlace)
                {
                    power = 0b10000000000000000000000000000000;
                    for (int i = 32; i >= stoppingDigitPlace; --i)
                    {
                        Integer shiftedDivisor = positiveDivisor.shiftLeft(valuePlace -
                            positiveDivisor.Values.Length + 1, i - divisorLeadingDigitPlace);                    
                        Integer difference = remainder - shiftedDivisor;
                        if (difference.Sign >= 0)
                        {
                            quotient[valuePlace] |= power;
                            remainder = difference;
                        }
                        power = power >> 1;
                    }
                }
                for (int i = Values.Length - 1; i >= positiveDivisor.Values.Length; --i)
                    calculateValue(i, 1);
                calculateValue(positiveDivisor.Values.Length - 1, divisorLeadingDigitPlace);
                IntegerDivision division = new IntegerDivision();
                division.quotient = new Integer(quotient,
                    (sbyte)(Sign * divisor.Sign)).shiftLeft(0, 1 - divisorLeadingDigitPlace);
                division.remainder = remainder;
                return division;
            }
            public override Number exponentiate(Number exponent)
            {
                if (Sign == 0)
                    return Zero;
                Integer exponentIntegerFactor = exponent.getGreatestIntegerFactor();
                if (exponentIntegerFactor > One)
                    return base.exponentiate(exponentIntegerFactor).exponentiate(
                        exponent / exponentIntegerFactor);
                if (exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)exponent;
                    Integer radicand = this;
                    Dictionary<Integer, Integer> radicandFactors =
                        new Dictionary<Integer, Integer> { };
                    if (radicand.Sign < 0)
                    {
                        radicand = -radicand;
                        radicandFactors.Add(new Integer(-1), One);
                    }
                    Integer factor;
                    for (Integer i = new Integer(2);
                        i <= radicand.euclideanDivideBy(new Integer(2)).quotient;)
                    {
                        IntegerDivision division = radicand.euclideanDivideBy(i);
                        if (division.remainder.Sign == 0)
                        {
                            factor = i;
                            if (radicandFactors.ContainsKey(factor))
                                ++radicandFactors[factor];
                            else
                                radicandFactors.Add(factor, One);
                            radicand = division.quotient;
                        }
                        else
                            ++i;
                    }
                    factor = radicand;
                    if (radicandFactors.ContainsKey(factor))
                        ++radicandFactors[factor];
                    else
                        radicandFactors.Add(factor, One);
                    List<Integer> exponentDivisors = exponentFraction.Denominator.getDivisors();
                    Dictionary<Integer, Number> termComponents = new Dictionary<Integer, Number>();
                    termComponents.Add(One, One);
                    termComponents.Add(exponentFraction.Denominator, One);
                    foreach (Integer n in radicandFactors.Keys)
                        for (int i = exponentDivisors.Count - 1; i >= 0; --i)
                            if (exponentDivisors[i] <= radicandFactors[n])
                            {
                                Integer index =
                                    (Integer)(exponentFraction.Denominator / exponentDivisors[i]);
                                if (!termComponents.ContainsKey(index))
                                    termComponents.Add(index, One);
                                IntegerDivision division =
                                    radicandFactors[n].euclideanDivideBy(exponentDivisors[i]);
                                for (Integer j = Zero; j < division.quotient; ++j)
                                    termComponents[index] = termComponents[index] * n;
                                for (Integer j = Zero; j < division.remainder; ++j)
                                    termComponents[exponentFraction.Denominator] =
                                        termComponents[exponentFraction.Denominator] * n;
                                break;
                            }
                    Number coefficient = termComponents[One];
                    termComponents.Remove(One);
                    Integer highestIndexComponent =
                        (Integer)termComponents[exponentFraction.Denominator];
                    if (highestIndexComponent == One)
                        termComponents.Remove(exponentFraction.Denominator);
                    else if (highestIndexComponent.Sign < 0)
                    {
                        if (highestIndexComponent == new Integer(-1))
                            termComponents.Remove(exponentFraction.Denominator);
                        else
                            termComponents[exponentFraction.Denominator] =
                                termComponents[exponentFraction.Denominator].negative();
                        Integer two = new Integer(2);
                        if (exponentFraction.Denominator.euclideanDivideBy(two).remainder.Sign == 0)
                            coefficient *= ComplexExponential.create(new Fraction(One,
                                two * exponentFraction.Denominator));
                        else
                            coefficient = coefficient.negative();
                    }
                    List<Factor> factors = new List<Factor>();
                    foreach (Integer index in termComponents.Keys)               
                        factors.Add(new Surd(termComponents[index], index));
                    Number output = Product.create(One, factors) * coefficient;
                    return output;
                }
                return base.exponentiate(exponent);
            }
            public override Integer getGreatestIntegerFactor()
            {
                return this;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                return (this - (Integer)obj).Sign;
            }
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public static bool operator ==(Integer a, Integer b)
            {
                if (a.Sign != b.Sign || a.Values.Length != b.Values.Length)
                    return false;
                for (int i = 0; i < a.Values.Length; ++i)
                    if (a.Values[i] != b.Values[i])
                        return false;
                return true;
            }
            public static bool operator !=(Integer a, Integer b)
            {
                return !(a == b);
            }
            public static bool operator <(Integer a, Integer b)
            {
                if ((a - b).Sign < 0)
                    return true;
                return false;
            }
            public static bool operator >(Integer a, Integer b)
            {
                if ((a - b).Sign > 0)
                    return true;
                return false;
            }
            public static bool operator <=(Integer a, Integer b)
            {
                return !(a > b);
            }
            public static bool operator >=(Integer a, Integer b)
            {
                return !(a < b);
            }        
            public override int GetHashCode()
            {
                uint output = 0;
                foreach (uint value in Values)
                    output = output ^ value;
                return (int)(output * Sign);
            }
            public override String insertString(String str)
            {
                if (this != One)
                    if (this == new Integer(-1))
                        return '-' + str;
                    else
                        return ToString() + str;
                return str;
            }
            public override string ToString()
            {
                if (Sign == 0)
                    return "0";            
                StringBuilder output = new StringBuilder();
                Integer quotient = this;
                Integer power = new Integer(10);
                while (quotient.Sign != 0)
                {
                    IntegerDivision division = quotient.euclideanDivideBy(power);
                    quotient = division.quotient;
                    output.Insert(0, division.remainder.Values[0]);
                }
                if (Sign < 0)
                    output.Insert(0, '-');
                return output.ToString();
            }
            public static Integer getGCD(Integer a, Integer b)
            {
                Integer c;
                while (b.Sign != 0)
                {
                    c = b;
                    b = a.euclideanDivideBy(b).remainder;
                    a = c;
                }
                if (a.Sign < 0)
                    a = -a;
                return a;
            }
            public static Integer getLCM(Integer a, Integer b)
            {
                return (Integer)(a / getGCD(a, b) * b);
            }
            public List<Integer> getDivisors()
            {
                Integer x = this;
                List<Integer> divisors = new List<Integer>();
                Integer divisor = One;
                Integer half = x.euclideanDivideBy(new Integer(2)).quotient;
                while (divisor <= half)
                {
                    if (x.euclideanDivideBy(divisor).remainder.Sign == 0)
                        divisors.Add(divisor);
                    ++divisor;
                }
                divisors.Add(x);
                return divisors;
            }
        }
        class Fraction : Rational
        {
            public Integer Numerator { get; }
            public Integer Denominator { get; }
            public Fraction(Integer numerator, Integer denominator)
            {
                Numerator = numerator;
                Denominator = denominator;
            }
            public static Rational create(Integer numerator, Integer denominator)
            {
                if (denominator.Sign == 0)
                    throw new DivideByZeroException();
                if (numerator.Sign == 0)
                    return numerator;
                Integer GCD = Integer.getGCD(numerator, denominator);
                numerator = numerator.euclideanDivideBy(GCD).quotient;
                denominator = denominator.euclideanDivideBy(GCD).quotient;
                if (denominator.Sign < 0)
                {
                    numerator = -numerator;
                    denominator = -denominator;
                }
                if (denominator == One)
                    return numerator;
                return new Fraction(numerator, denominator);
            }
            protected override Number add(Number number)
            {
                if (number is Integer)
                    return create((Integer)(Numerator + number * Denominator), Denominator);
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return create(Numerator * fraction.Denominator + fraction.Numerator * Denominator,
                        Denominator * fraction.Denominator);
                }
                return number + this;
            }
            public override Number negative()
            {
                return new Fraction(-Numerator, Denominator);
            }
            protected override Number multiply(Number number)
            {
                if (number is Integer)
                    return create((Integer)(Numerator * number), Denominator);
                if (number is Fraction)
                {
                    Fraction fraction = (Fraction)number;
                    return create(Numerator * fraction.Numerator, Denominator * fraction.Denominator);
                }
                return number * this;
            }
            public override Number reciprocal()
            {
                return create(Denominator, Numerator);
            }
            public override Integer getDenominatorLCM()
            {
                return Denominator;
            }
            public override Integer getGreatestIntegerFactor()
            {
                return Numerator;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                Fraction fraction = (Fraction)obj;
                comparison = (Numerator - fraction.Numerator).Sign;
                if (comparison != 0)
                    return comparison;
                return (Denominator - fraction.Denominator).Sign;
            }
            public override int GetHashCode()
            {
                return Numerator.GetHashCode() ^ Denominator.GetHashCode();
            }
            public override String insertString(String str)
            {
                return Numerator.insertString(str) + "/" + Denominator;
            }
            public override String ToString()
            {
                return Numerator + "/" + Denominator;
            }
        }
        abstract class Exponentiation : Factor
        {
            public Number Base { get; }
            public Number Exponent { get; }
            protected Exponentiation(Number expBase, Number exponent)
            {
                Base = expBase;
                Exponent = exponent;
            }
            public static Factor create(Number expBase, Number exponent)
            {
                if (exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)exponent;
                    if (exponentFraction.Numerator == One)
                        return new Surd(expBase, exponentFraction.Denominator);
                }
                return new Transcendental(expBase, exponent);
            }
            protected override Number add(Number number)
            {
                if (number is Integer)
                {
                    Integer integer = (Integer)number;
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Fraction)
                    return new Sum(new List<Term> { this, (Fraction)number });
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (CompareTo(exponentiation) == 0)
                        return Product.create(new Integer(2), new List<Factor> { this });
                    if (CompareTo(exponentiation.negative()) == 0)
                        return Zero;
                    return new Sum(new List<Term> { this, exponentiation });
                }
                return number + this;
            }
            public override Number negative()
            {
                return Product.create(new Integer(-1), new List<Factor> { this });
            }
            protected override Number multiply(Number number)
            {
                if (number is Rational)
                    return Product.create((Rational)number, new List<Factor> { this });
                if (number is Exponentiation)
                {
                    Exponentiation exponentiation = (Exponentiation)number;
                    if (Base.Equals(exponentiation.Base))
                        return Base.exponentiate(Exponent + exponentiation.Exponent);
                    if (Exponent.Equals(exponentiation.Exponent))
                    {
                        Number outputBase = Base * exponentiation.Base;
                        if (outputBase is Exponentiation)
                            return create(((Exponentiation)outputBase).Base, Exponent * Exponent);
                        return (outputBase).exponentiate(Exponent);
                    }
                    return Product.create(One, new List<Factor> { this, exponentiation });
                }
                return number * this;
            }
            public override Number exponentiate(Number exponent)
            {
                return Base.exponentiate(Exponent * exponent);
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                Exponentiation exponentiation = (Exponentiation)obj;
                comparison = Exponent.CompareTo(exponentiation.Exponent);
                if (comparison != 0)
                    return comparison;
                return Base.CompareTo(exponentiation.Base);
            }
            public override int GetHashCode()
            {
                return Base.GetHashCode() ^ Exponent.GetHashCode();
            }
            public override String insertString(String str)
            {
                return ToString() + '*' + str;
            }
            public override String ToString()
            {
                String encloseNumber(Number number)
                {
                    string numberString = number.ToString();
                    if ((number is Integer && ((Integer)number).Sign >= 0) || numberString == "i")
                        return numberString;
                    return '(' + numberString + ')';
                }
                return encloseNumber(Base) + '^' + encloseNumber(Exponent);
            }
        }
        class Surd : Exponentiation
        {
            public Integer Index { get; }
            public Surd(Number radicand, Integer index) :
                base(radicand, new Fraction(One, index))
            {
                Index = index;
            }
            protected override Number multiply(Number number)
            {
                if (number is Surd)
                {
                    Surd surd = (Surd)number;
                    Integer indexLCM = Integer.getLCM(Index, surd.Index);
                    return (Base.exponentiate(indexLCM / Index) * surd.Base.exponentiate(indexLCM /
                        surd.Index)).exponentiate(Fraction.create(One, indexLCM));
                }
                return base.multiply(number);
            }
            public override Number reciprocal()
            {
                return Base.exponentiate(Fraction.create(Index - One, Index)) / Base;
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                Polynomial baseMinimalPolynomial = Base.MinimalPolynomial;
                if (baseMinimalPolynomial == null)
                    return null;
                List<Rational> coefficients =
                    new List<Rational> { baseMinimalPolynomial.Coefficients[0] };
                for (int i = 1; i < baseMinimalPolynomial.Coefficients.Count; ++i)
                {
                    for (Integer j = One; j < Index; ++j)
                        coefficients.Add(Zero);
                    coefficients.Add(baseMinimalPolynomial.Coefficients[i]);
                }
                return new Polynomial(coefficients);
            }
            public override List<Number> getConjugates()
            {
                List<Number> baseConjugates = Base.getConjugates();
                if (baseConjugates == null)
                    return null;
                List<Number> rootsOfUnity = new List<Number>();
                for (Integer i = Zero; i < Index; ++i)
                    rootsOfUnity.Add(ComplexExponential.create(Fraction.create(i, Index)));
                List<Number> conjugateCandidates = new List<Number>();
                foreach (Number conjugate in baseConjugates)
                    foreach (Number root in rootsOfUnity)
                        conjugateCandidates.Add(new Surd(conjugate, Index) * root);
                return testConjugateCandidates(conjugateCandidates);
            }
        }
        class Transcendental : Exponentiation
        {
            public Transcendental(Number expBase, Number exponent) : base(expBase, exponent)
            { }
            public override Number reciprocal()
            {
                return create(Base, Exponent.negative());
            }            
            protected override Polynomial calculateMinimalPolynomial()
            {
                return null;
            }
            public override List<Number> getConjugates()
            {
                return null;
            }
        }
        class ComplexExponential : Factor
        {//Represents e^(tau*i*Exponent).
            public Number Exponent { get; }
            ComplexExponential(Number exponent)
            {
                if (exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)exponent;
                    Exponent = exponentFraction - exponentFraction.Numerator.euclideanDivideBy(
                        exponentFraction.Denominator).quotient;
                    if (exponentFraction.Numerator.Sign < 0)
                        Exponent += One;
                }
                else
                    Exponent = exponent;
            }
            public static Number create(Number exponent)
            {
                if (exponent is Integer)
                    return One;
                if (exponent is Fraction)
                {
                    Fraction exponentFraction = (Fraction)exponent;
                    Integer denominator = exponentFraction.Denominator;
                    Integer numerator =
                        exponentFraction.Numerator.euclideanDivideBy(denominator).remainder;                
                    Integer two = new Integer(2);
                    IntegerDivision division = denominator.euclideanDivideBy(two);
                    int multiplicityOfTwo = 0;
                    while (division.remainder.Sign == 0)
                    {
                        denominator = division.quotient;
                        division = denominator.euclideanDivideBy(two);
                        ++multiplicityOfTwo;
                    }
                    Number cosine;
                    if (denominator == One)
                        cosine = One;
                    else if (denominator == new Integer(3)) 
                        cosine = new Fraction(new Integer(-1), two);
                    else if(denominator == new Integer(5))
                        cosine = (new Surd(new Integer(5),
                            new Integer(2)) - One) / new Integer(4);
                    else if (denominator == new Integer(15))
                    {
                        Surd a = new Surd(new Integer(5), two);
                        cosine = (One + a + new Surd(new Integer(30) -
                            new Integer(6) * a, two)) / new Integer(8);         
                    }
                    else if (denominator == new Integer(17))
                    {
                        Integer seventeen = new Integer(17);
                        Surd a = new Surd(seventeen, two);
                        Surd b = new Surd(seventeen * two - a * two, two);
                        cosine = (new Integer(-1) + a + b + new Surd(seventeen + a * new Integer(3) -
                            b - new Surd((seventeen + a) * two, two) * two, two) * two) *
                            new Fraction(One, new Integer(16));
                    }
                    else
                        return new ComplexExponential(exponent);
                    Rational half = new Fraction(One, two);
                    Integer three = new Integer(3); 
                    Integer four = new Integer(4);
                    for (int i = 0; i < multiplicityOfTwo; ++i)
                    {
                        denominator = denominator * two;
                        cosine = (cosine * half + half).exponentiate(half);
                        if (four < three * denominator && denominator < four)
                            cosine = cosine.negative();
                    }
                    Number r = One;
                    Number s = cosine;
                    Number angleMultipleCosine = cosine;
                    for (Integer i = One; i < numerator; ++i)
                    {
                        angleMultipleCosine = two * cosine * s - r;
                        r = s;
                        s = angleMultipleCosine;
                    }
                    Number output = ComplexNumber.create(angleMultipleCosine,
                        (One - angleMultipleCosine * angleMultipleCosine).exponentiate(half));
                    Integer t = four * numerator;
                    if (denominator < t && t < two * denominator ||
                        three * denominator < t && t < four * denominator)
                        return output.negative();
                    return output;
                }
                return new ComplexExponential(exponent);
            }
            protected override Number add(Number number)
            {
                if (number is Integer)
                {
                    Integer integer = (Integer)number;
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Fraction || number is Exponentiation)
                    return new Sum(new List<Term> { this, (Term)number });
                if (number is ComplexExponential)
                {
                    Number exponent = ((ComplexExponential)number).Exponent;
                    if (Exponent.CompareTo(exponent) == 0)
                        return Product.create(new Integer(2), new List<Factor> { this });
                }
                return number + this;
            }
            public override Number negative()
            {
                return Product.create(new Integer(-1), new List<Factor> { this });
            }
            protected override Number multiply(Number number)
            {
                if (number is Integer)
                    return Product.create((Integer)number, new List<Factor> { this });
                if (number is Fraction)
                    return Product.create((Fraction)number, new List<Factor> { this, });
                if (number is Exponentiation)
                    return Product.create(One, new List<Factor> { this, (Factor)number });
                if (number is ComplexExponential)
                    return create(Exponent + ((ComplexExponential)number).Exponent);
                return number * this;
            }
            public override Number reciprocal()
            {
                return create(Exponent.negative());
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                if (!(Exponent is Fraction))
                    return null;            
                List<Rational> coefficients = new List<Rational> { new Integer(-1) };
                Integer index = ((Fraction)Exponent).Denominator;
                for (Integer i = One; i < index; ++i)
                    coefficients.Add(Zero);
                coefficients.Add(One);
                return new Polynomial(coefficients);
            }
            public override List<Number> getConjugates()
            {
                if (!(Exponent is Fraction))
                    return null;
                List<Number> conjugates = new List<Number>();
                Integer index = ((Fraction)Exponent).Denominator;
                for (Integer i = Zero; i < index; ++i)
                    conjugates.Add(create(Fraction.create(i, index)));
                return conjugates;
            }
            public override Number exponentiate(Number exponent)
            {
                return create(exponent * exponent);
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                return Exponent.CompareTo(((ComplexExponential)obj).Exponent);
            }
            public override int GetHashCode()
            {
                return Exponent.GetHashCode();
            }
            public override string insertString(string str)
            {
                return ToString() + '*' + str;
            }
            public override string ToString()
            {
                Number exponentWith_i = Exponent * ComplexNumber.create(Zero, One);
                return "e^(" + exponentWith_i.insertString("tau") + ')';
            }
        }
        class Product : Term
        {
            public Rational Coefficient { get; }
            public List<Factor> Factors { get; }
            Product(Rational coefficient, List<Factor> factors)
            {
                Coefficient = coefficient;
                Factors = factors;
                Factors.Sort();
            }
            public static Term create(Rational coefficient, List<Factor> factors)
            {
                if (coefficient is Integer)
                {
                    Integer coefficientInteger = (Integer)coefficient;
                    if (coefficientInteger.Sign == 0)
                        return Zero;
                    if (factors.Count == 1 && coefficientInteger == One)
                        return factors[0];
                }
                if (factors.Count == 0)
                    return coefficient;
                return new Product(coefficient, factors);
            }
            protected override Number add(Number number)
            {
                if (number is Integer)
                {
                    Integer integer = (Integer)number;
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Factor)
                {
                    if (Factors.Count == 1 && number.Equals(Factors[0]))
                        return create((Rational)(Coefficient + One), Factors);
                    return new Sum(new List<Term> { this, (Term)number });
                }
                if (number is Product)
                {
                    Product product = (Product)number;
                    if (Factors.Count == product.Factors.Count)
                    {
                        for (int i = 0; i < Factors.Count; ++i)
                            if (Factors[i].CompareTo(product.Factors[i]) != 0)
                                return new Sum(new List<Term> { this, product });
                        return create((Rational)(Coefficient + product.Coefficient), Factors);
                    }
                    return new Sum(new List<Term> { this, product });
                }
                return number + this;
            }
            public override Number negative()
            {
                return create((Rational)Coefficient.negative(), Factors);
            }
            protected override Number multiply(Number number)
            {
                if (number is Rational)
                    return create((Rational)(Coefficient * number), Factors);
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
                            return create(Coefficient, factors) * product;
                        }
                    }
                    factors.Add((Factor)number);
                    return create(Coefficient, factors);
                }
                if (number is Product)
                {
                    Product product = (Product)number;
                    Number output =
                        new Product((Rational)(Coefficient * product.Coefficient), Factors);
                    for (int i = 0; i < product.Factors.Count; ++i)
                        output = output * product.Factors[i];
                    return output;
                }
                return number * this;
            }
            public override Number reciprocal()
            {
                List<Factor> factors = new List<Factor>();
                foreach (Factor factor in Factors)
                    factors.Add((Factor)factor.reciprocal());
                return create((Rational)Coefficient.reciprocal(), factors);
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                Polynomial[] minimalPolynomials = new Polynomial[Factors.Count];
                for (int i = 0; i < Factors.Count; ++i)
                {
                    Polynomial polynomial = Factors[i].MinimalPolynomial;
                    if (polynomial == null)
                        return null;
                    minimalPolynomials[i] = polynomial;
                }
                MultivariatePolynomial variableForm = new MultivariatePolynomial(minimalPolynomials);
                int[] indices = new int[minimalPolynomials.Length];
                for (int i = 0; i < indices.Length; ++i)
                    indices[i] = 1;
                variableForm.setCoefficient(Coefficient, indices);
                return variableForm.getMinimalPolynomial();            
            }
            public override List<Number> getConjugates()
            {
                List<List<Number>> factorConjugates = new List<List<Number>>();
                foreach (Factor factor in Factors)
                {
                    List<Number> conjugates = factor.getConjugates();
                    if (conjugates == null)
                        return null;
                    factorConjugates.Add(conjugates);
                }
                List<List<Number>> conjugateCombinations = generateCartesianProduct(factorConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = One;
                    foreach (Number number in combination)
                        conjugate *= number;
                    conjugateCandidates.Add(conjugate);
                }
                return testConjugateCandidates(conjugateCandidates);
            }
            public override Integer getDenominatorLCM()
            {
                return Coefficient.getDenominatorLCM();
            }
            public override Integer getGreatestIntegerFactor()
            {
                return Coefficient.getGreatestIntegerFactor();
            }
            public override Number exponentiate(Number exponent)
            {
                Number output = Coefficient.exponentiate(exponent);
                foreach (Factor factor in Factors)
                    output *= factor.exponentiate(exponent);
                return output;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                Product product = (Product)obj;
                comparison = Coefficient.CompareTo(Coefficient);
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
                    output = output ^ factor.GetHashCode();
                return output;
            }
            StringBuilder formatFactors()
            {
                StringBuilder output = new StringBuilder(Factors[0].ToString());
                for (int i = 1; i < Factors.Count; ++i)
                    output.Append('*' + Factors[i].ToString());
                return output;
            }
            string addCoefficient(StringBuilder factors)
            {
                Integer numerator;
                if (Coefficient is Integer)
                    numerator = (Integer)Coefficient;
                else
                {
                    Fraction coefficientFraction = (Fraction)Coefficient;
                    factors.Append("/" + coefficientFraction.Denominator);
                    numerator = coefficientFraction.Numerator;
                }
                if (numerator != One)
                {
                    if (numerator == new Integer(-1))
                        return factors.Insert(0, '-').ToString();
                    else
                    {
                        if (char.IsDigit(factors[0]))
                            return factors.Insert(0, numerator + "*").ToString();
                        return factors.Insert(0, numerator).ToString();
                    }
                }
                return factors.ToString();
            }
            public override string insertString(string str)
            {
                StringBuilder output = formatFactors().Append('*' + str);
                return addCoefficient(output);
            }
            public override string ToString()
            {
                StringBuilder output = formatFactors();
                return addCoefficient(output);
            }
        }
        class Sum : Number
        {
            public List<Term> Terms { get; }
            public Sum(List<Term> terms)
            {
                terms.Sort();
                Terms = terms;
            }
            protected override Number add(Number number)
            {
                if (number is Term)
                {
                    List<Term> terms = new List<Term>(Terms);
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
                    return new Sum(terms);
                }
                else if (number is Sum)
                {
                    Sum sum = (Sum)number;
                    Number output = this;
                    foreach (Term term in sum.Terms)
                        output += term;
                    return output;
                }
                return number + this;
            }
            public override Number negative()
            {
                List<Term> terms = new List<Term>();
                foreach (Term term in Terms)
                    terms.Add((Term)term.negative());
                return new Sum(terms);
            }
            protected override Number multiply(Number number)
            {
                List<Term> terms = new List<Term>();
                if (number is Term)
                {
                    Number output = Zero;
                    foreach (Term term in Terms)
                        output += term * number;
                    return output;
                }
                if (number is Sum)
                {
                    Sum multiplier = (Sum)number;
                    Number output = Zero;
                    foreach (Term multiplicandTerm in Terms)
                        foreach (Term multiplierTerm in multiplier.Terms)
                            output += multiplicandTerm * multiplierTerm;
                    return output;
                }
                return number * this;
            }
            public override Number reciprocal()
            {
                List<Number> conjugates = getConjugates();
                if (conjugates == null) 
                    return new Transcendental(this, new Integer(-1));
                conjugates.Remove(this);
                Number numerator = One;
                foreach(Number conjugate in conjugates)                
                    numerator *= conjugate;
                return numerator / (numerator * this);
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                Polynomial[] minimalPolynomials = new Polynomial[Terms.Count];
                Rational constant;
                int startIndex;
                if (Terms[0] is Rational) 
                {
                    constant = (Rational)Terms[0];
                    minimalPolynomials = new Polynomial[Terms.Count - 1];
                    startIndex = 1;
                }
                else
                {
                    constant = Zero;
                    minimalPolynomials = new Polynomial[Terms.Count];
                    startIndex = 0;
                }
                for (int i = startIndex; i < Terms.Count; ++i)
                {
                    Polynomial polynomial = Terms[i].MinimalPolynomial;
                    if (polynomial == null)
                        return null;
                    minimalPolynomials[i - startIndex] = polynomial;
                }
                MultivariatePolynomial variableForm = new MultivariatePolynomial(minimalPolynomials);
                int[] indices = new int[minimalPolynomials.Length];
                variableForm.setCoefficient(constant, indices);
                for (int i = 0; i < indices.Length; ++i)
                {
                    indices[i] = 1;
                    variableForm.setCoefficient(One, indices);
                    indices[i] = 0;
                }
                return variableForm.getMinimalPolynomial();
            }
            public override List<Number> getConjugates()
            {
                List<List<Number>> termConjugates = new List<List<Number>>();
                foreach (Term term in Terms)
                {
                    List<Number> conjugates = term.getConjugates();
                    if (conjugates == null)
                        return null;
                    termConjugates.Add(term.getConjugates());
                }
                List<List<Number>> conjugateCombinations = generateCartesianProduct(termConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = Zero;
                    foreach (Number number in combination)
                        conjugate += number;
                    conjugateCandidates.Add(conjugate);
                }
                return testConjugateCandidates(conjugateCandidates);
            }
            public override Integer getDenominatorLCM()
            {
                Integer LCM = One;
                foreach (Term term in Terms)
                    LCM = Integer.getLCM(LCM, term.getDenominatorLCM());
                return LCM;
            }
            public override Integer getGreatestIntegerFactor()
            {
                Integer GCD = Terms[0].getGreatestIntegerFactor();
                for (int i = 1; i < Terms.Count; ++i)
                    GCD = Integer.getGCD(GCD, Terms[i].getGreatestIntegerFactor());
                return GCD;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
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
            public override string insertString(string str)
            {
                return '(' + ToString() + ')' + str;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder(Terms[0].ToString());
                for (int i = 1; i < Terms.Count; ++i)
                {
                    string termString = Terms[i].ToString();
                    if (termString[0] != '-')
                        output.Append('+');
                    output.Append(termString);
                }
                return output.ToString();
            }
        }
        class ComplexNumber : Number
        {
            public Number Real { get; }
            public Number Imaginary { get; }
            protected ComplexNumber(Number real, Number imaginary)
            {
                Real = real;
                Imaginary = imaginary;
            }
            public static Number create(Number real, Number imaginary)
            {
                if (imaginary is Integer)
                {
                    Integer imaginaryInteger = (Integer)imaginary;
                    if (imaginaryInteger.Sign == 0)
                        return real;
                }
                return new ComplexNumber(real, imaginary);
            }
            protected override Number add(Number number)
            {
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return create(Real + complexNumber.Real, Imaginary + complexNumber.Imaginary);
                }
                return create(Real + number, Imaginary);
            }
            public override Number negative()
            {
                return create(Real.negative(), Imaginary.negative());
            }
            protected override Number multiply(Number number)
            {
                if (number is ComplexNumber)
                {
                    ComplexNumber complexNumber = (ComplexNumber)number;
                    return create(Real * complexNumber.Real - Imaginary * complexNumber.Imaginary,
                        Real * complexNumber.Imaginary + Imaginary * complexNumber.Real);
                }
                return create(Real * number, Imaginary * number);
            }
            public override Number reciprocal()
            {
                Number denominator = Real * Real + Imaginary * Imaginary;
                return create(Real / denominator, Imaginary.negative() / denominator);
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                MultivariatePolynomial variableForm = new MultivariatePolynomial(
                    new Polynomial[] { Real.MinimalPolynomial, Imaginary.MinimalPolynomial,
                        new Polynomial(new List<Rational> { One, Zero, One })});
                variableForm.setCoefficient(One, new int[] { 1, 0, 0 });
                variableForm.setCoefficient(One, new int[] { 0, 1, 1 });
                return variableForm.getMinimalPolynomial();
            }
            public override List<Number> getConjugates()
            {
                List<Number> realPartConjugates = Real.getConjugates();
                List<Number> imaginaryPartConjugates = Imaginary.getConjugates();
                List<Number> conjugateCandidates = new List<Number>();
                foreach (Number a in realPartConjugates)
                    foreach (Number b in imaginaryPartConjugates)
                        conjugateCandidates.Add(a + create(Zero, One) * b);
                return testConjugateCandidates(conjugateCandidates);
            }
            public override Integer getDenominatorLCM()
            {
                return Integer.getLCM(Real.getDenominatorLCM(), Imaginary.getDenominatorLCM());
            }
            public override Integer getGreatestIntegerFactor()
            {
                return Integer.getGCD(Real.getGreatestIntegerFactor(),
                    Imaginary.getGreatestIntegerFactor());
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                ComplexNumber gaussianInteger = (ComplexNumber)obj;
                comparison = Real.CompareTo(gaussianInteger.Real);
                if (comparison != 0)
                    return comparison;
                return Imaginary.CompareTo(gaussianInteger.Imaginary);
            }
            public override int GetHashCode()
            {
                return Real.GetHashCode() ^ Imaginary.GetHashCode();
            }
            public override string insertString(string str)
            {
                if (Real is Integer && ((Integer)Real).Sign == 0)
                    return Imaginary.insertString("i*" + str);
                return '(' + ToString() + ')' + str;
            }
            public override string ToString()
            {
                StringBuilder output =
                    new StringBuilder(Imaginary.insertString("i").ToString());
                if (!(Real is Integer && ((Integer)Real).Sign == 0))
                {
                    if (output[0] != '-')
                        output.Insert(0, '+');
                    output.Insert(0, Real.ToString());
                }
                return output.ToString();
            }
        }
        class Polynomial : IComparable
        {
            //The index of the coefficient represents the degree of its term.
            public List<Rational> Coefficients { get; }
            public Polynomial(List<Rational> coefficients)
            {
                while (coefficients.Count != 0 &&
                    coefficients[coefficients.Count - 1].Equals(Zero))
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
            }
            public Polynomial getIntegerPolynomial()
            {
                if (Coefficients.Count == 0)
                    return this;
                Integer LCM = One;
                foreach (Rational coefficient in Coefficients)
                    LCM = Integer.getLCM(LCM, coefficient.getDenominatorLCM());
                Integer integerLCM = LCM;
                List<Rational> integerCoefficients = new List<Rational>();
                foreach (Rational coefficient in Coefficients)
                    integerCoefficients.Add((Rational)(coefficient * integerLCM));
                Integer content = (Integer)integerCoefficients[0];
                for (int i = 1; i < integerCoefficients.Count; ++i)
                    content = Integer.getGCD(content, (Integer)integerCoefficients[i]);
                List<Rational> reducedCoefficients = new List<Rational>();
                foreach (Rational coefficient in integerCoefficients)
                    reducedCoefficients.Add((Integer)(coefficient / content));
                return new Polynomial(reducedCoefficients);
            }
            public static Polynomial operator +(Polynomial a, Polynomial b)
            {
                if (a.Coefficients.Count < b.Coefficients.Count)
                    return b + a;
                List<Rational> output = a.Coefficients;
                for (int i = 0; i < b.Coefficients.Count; ++i)
                    output[i] = (Rational)(output[i] + b.Coefficients[i]);
                return new Polynomial(output);
            }
            public Polynomial negative()
            {
                List<Rational> output = new List<Rational>();
                foreach (Number coefficient in Coefficients)
                    output.Add((Rational)coefficient.negative());
                return new Polynomial(output);
            }
            public static Polynomial operator -(Polynomial a, Polynomial b)
            {
                return a + b.negative();
            }
            public static Polynomial operator *(Polynomial a, Polynomial b)
            {
                List<Rational> output = new List<Rational>();
                for (int i = 0; i < a.Coefficients.Count + b.Coefficients.Count - 1; ++i)
                    output.Add(Zero);
                for (int i = 0; i < a.Coefficients.Count; ++i)
                    for (int j = 0; j < b.Coefficients.Count; ++j)
                        output[i + j] =
                            (Rational)(output[i + j] + a.Coefficients[i] * b.Coefficients[j]);
                return new Polynomial(output);
            }
            public static Polynomial operator *(Rational a, Polynomial b)
            {
                return new Polynomial(new List<Rational> { a }) * b;
            }
            public static Polynomial operator *(Polynomial a, Rational b)
            {
                return b * a;
            }
            struct PolynomialDivision
            {
                public Polynomial Quotient;
                public Polynomial Remainder;
            }
            PolynomialDivision divide(Polynomial divisor)
            {
                List<Rational> quotient = new List<Rational>();
                List<Rational> remainder = new List<Rational>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, (Rational)(remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1]));
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            remainder[remainder.Count - j] =
                                (Rational)(remainder[remainder.Count - j] - quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1]);
                    }
                PolynomialDivision division;
                division.Quotient = new Polynomial(quotient);
                division.Remainder = new Polynomial(remainder);
                return division;
            }
            public static Polynomial operator /(Polynomial a, Polynomial b)
            {
                return a.divide(b).Quotient;
            }
            public static Polynomial operator %(Polynomial a, Polynomial b)
            {
                return a.divide(b).Remainder;
            }
            public Rational evaluateAt(Integer input)
            {
                Rational output = Coefficients[0];
                for (int i = 1; i < Coefficients.Count; ++i)
                    output = (Rational)(output + Coefficients[i] * input.exponentiate(new Integer(i)));
                return output;
            }
            public Polynomial getDerivative()
            {
                List<Rational> derivativeCoefficients = new List<Rational>();
                for (int i = 1; i < Coefficients.Count; ++i)
                    derivativeCoefficients.Add((Rational)(Coefficients[i] * new Integer(i)));
                return new Polynomial(derivativeCoefficients);
            }
            public virtual int CompareTo(object obj)
            {
                Polynomial polynomial = (Polynomial)obj;
                int comparison = Coefficients.Count - polynomial.Coefficients.Count;
                if (comparison != 0)
                    return comparison;
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    comparison = Coefficients[i].CompareTo(polynomial.Coefficients[i]);
                    if (comparison != 0)
                        return comparison;
                }
                return 0;
            }
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public override int GetHashCode()
            {
                int output = 0;
                foreach (Rational coefficient in Coefficients)
                    output ^= coefficient.GetHashCode();
                return output;
            }
            public bool Equals(Polynomial polynomial)
            {
                return CompareTo(polynomial) == 0;
            }
            public static Polynomial getGCD(Polynomial a, Polynomial b)
            {
                while (b.Coefficients.Count != 0)
                {
                    Polynomial c = a % b;
                    a = b;
                    b = c;
                }
                return a;
            }
            public List<Polynomial> getFactors()
            {
                Polynomial integerPolynomial = getIntegerPolynomial();          
                List<Polynomial> factors = new List<Polynomial>();
                Polynomial derivative = integerPolynomial.getDerivative();
                Polynomial a = getGCD(integerPolynomial, derivative);
                Polynomial b = integerPolynomial / a;
                Polynomial c = derivative / a - b.getDerivative();
                while (!(b.Coefficients.Count == 1 && (b.Coefficients[0].Equals(One)
                    || b.Coefficients[0].Equals(new Integer(-1)))))
                {
                    a = getGCD(b, c);
                    if (!factors.Contains(a))
                        factors.Add(a.getIntegerPolynomial());
                    b = b / a;
                    c = c / a - b.getDerivative();
                }
                List<Polynomial> irreducibleFactors = new List<Polynomial>();
                void attemptFactorization()
                {
                    foreach (Polynomial factor in factors)
                    {
                        for (int i = 1; i < (factor.Coefficients.Count + 1) / 2; ++i)
                        {
                            List<Integer> outputValues = new List<Integer>();
                            List<List<Integer>> outputValueDivisors = new List<List<Integer>>();
                            outputValues.Add((Integer)integerPolynomial.evaluateAt(Zero));
                            outputValueDivisors.Add(outputValues[0].getDivisors());
                            for (int j = 1; j <= i; ++j)
                            {
                                outputValues.Add(
                                    (Integer)integerPolynomial.evaluateAt(new Integer(j)));
                                outputValueDivisors.Add(outputValues[j].getDivisors());
                                int numberOfPositiveDivisors = outputValueDivisors[j].Count;
                                for (int k = 0; k < numberOfPositiveDivisors; ++k)
                                    outputValueDivisors[j].Add(-outputValueDivisors[j][k]);
                            }
                            List<List<Integer>> outputCombinations =
                                generateCartesianProduct(outputValueDivisors);
                            foreach (List<Integer> combination in outputCombinations)
                            {
                                Polynomial factorCandidate = new Polynomial(new List<Rational>());
                                for (int j = 0; j < outputValueDivisors.Count; ++j)
                                {
                                    Polynomial numerator =
                                        new Polynomial(new List<Rational> { outputValues[j] });
                                    int denominator = 1;
                                    for (int k = 0; k < outputValueDivisors.Count; ++k)
                                        if (k != j)
                                        {
                                            numerator *= new Polynomial(
                                                new List<Rational> { new Integer(-k), One });
                                            denominator *= j - k;
                                        }
                                    List<Rational> basisPolynomialCoefficients = new List<Rational>();
                                    foreach (Rational coefficient in numerator.Coefficients)
                                        basisPolynomialCoefficients.Add((Rational)(coefficient /
                                            new Integer(denominator)));
                                    factorCandidate += new Polynomial(basisPolynomialCoefficients);
                                }
                                if ((factor % factorCandidate).Coefficients.Count == 0) 
                                {
                                    if (!irreducibleFactors.Contains(factorCandidate))
                                        irreducibleFactors.Add(factorCandidate);
                                    factors.Remove(factor);
                                    Polynomial reducedFactor = factor / factorCandidate;
                                    if (!factors.Contains(reducedFactor))
                                        factors.Add(reducedFactor);
                                    return;
                                }
                            }
                        }
                        if (!irreducibleFactors.Contains(factor))
                            irreducibleFactors.Add(factor);
                        factors.Remove(factor);
                        return;
                    }
                }
                while (factors.Count > 0)
                    attemptFactorization();
                irreducibleFactors.Sort();
                return irreducibleFactors;
            }
#if DEBUG
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                foreach (Rational coefficient in Coefficients)
                    output.Append(coefficient.ToString() + ',');
                return output.ToString();
            }
#endif
        }
        class MultivariatePolynomial
        {//Represents a polynomial whose variables each represent a specific number.
            Polynomial[] MinimalPolynomials;        
            public Rational[] Coefficients { get; }
            public MultivariatePolynomial(Polynomial[] minimalPolynomials)
            {//The nth minimalPolynomial entry is the minimal polynomial of the number the nth variable
             //represents.
                MinimalPolynomials = minimalPolynomials;
                int numberOfCoefficients = 1;
                for (int i = 0; i < MinimalPolynomials.Length; ++i)
                    numberOfCoefficients *= MinimalPolynomials[i].Coefficients.Count - 1;
                Coefficients = new Rational[numberOfCoefficients];
                for (int i = 0; i < Coefficients.Length; ++i)
                    Coefficients[i] = Zero;
            }
            int getInternalIndex(params int[] externalIndices)
            {
                int index = 0;
                int multiplier = 1;
                for (int i = 0; i < externalIndices.Length; ++i)
                {
                    index += externalIndices[i] * multiplier;
                    multiplier *= MinimalPolynomials[i].Coefficients.Count - 1;
                }
                return index;
            }
            int[] getExternalIndices(int internalIndex)
            {
                int[] indices = new int[MinimalPolynomials.Length];
                for (int i = 0; i < MinimalPolynomials.Length; ++i)
                {
                    indices[i] = internalIndex % (MinimalPolynomials[i].Coefficients.Count - 1);
                    internalIndex /= MinimalPolynomials[i].Coefficients.Count - 1;
                }
                return indices;
            }
            public void setCoefficient(Rational coefficient, int[] indices)
            {//The nth index is the degree with respect to the nth variable of the term being set.
                Coefficients[getInternalIndex(indices)] = coefficient;
            }
            public static MultivariatePolynomial operator +(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial sum = new MultivariatePolynomial(a.MinimalPolynomials);
                for (int i = 0; i < a.Coefficients.Length; ++i)
                    sum.Coefficients[i] = (Rational)(a.Coefficients[i] + b.Coefficients[i]);
                return sum;
            }
            static MultivariatePolynomial reductionlessMultiply(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial product =
                    new MultivariatePolynomial(a.MinimalPolynomials);
                for (int i = 0; i < a.Coefficients.Length; ++i)
                    for (int j = 0; j < b.Coefficients.Length; ++j)
                    {
                        Rational coefficient = (Rational)(a.Coefficients[i] * b.Coefficients[j]);
                        if (!coefficient.Equals(Zero)) 
                            product.Coefficients[i + j] =
                                (Rational)(a.Coefficients[i] * b.Coefficients[j]);
                    }
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(a.MinimalPolynomials);
                for (int i = 0; i < a.Coefficients.Length; ++i)
                    for (int j = 0; j < b.Coefficients.Length; ++j)
                    {                    
                        Rational coefficient = (Rational)(a.Coefficients[i] * b.Coefficients[j]);
                        if (coefficient.Equals(Zero))
                            continue;
                        MultivariatePolynomial productComponent =
                            new MultivariatePolynomial(a.MinimalPolynomials);
                        productComponent.Coefficients[0] = coefficient;
                        int[] indices = a.getExternalIndices(i);
                        int[] indicesB = a.getExternalIndices(j);
                        int[] termComponentIndices = new int[indices.Length];
                        for (int k = 0; k < indices.Length; ++k)
                        {
                            indices[k] += indicesB[k];
                            MultivariatePolynomial termComponent =
                                new MultivariatePolynomial(a.MinimalPolynomials);                        
                            if (indices[k] < a.MinimalPolynomials[k].Coefficients.Count - 1)
                            {
                                termComponentIndices[k] = indices[k];
                                termComponent.setCoefficient(One, termComponentIndices);
                            }
                            else
                            {
                                int remainder =
                                    indices[k] - a.MinimalPolynomials[k].Coefficients.Count + 1;
                                if (remainder != 0)
                                {
                                    termComponentIndices[k] = remainder;
                                    MultivariatePolynomial remainderComponent =
                                        new MultivariatePolynomial(a.MinimalPolynomials);
                                    remainderComponent.setCoefficient(One,
                                        termComponentIndices);
                                    productComponent =
                                        reductionlessMultiply(productComponent, remainderComponent);
                                }
                                for (int l = 0; l < a.MinimalPolynomials[k].Coefficients.Count - 1;
                                    ++l)
                                {
                                    termComponentIndices[k] = l;
                                    termComponent.setCoefficient((Rational)a.MinimalPolynomials[k].
                                        Coefficients[l].negative(), termComponentIndices);
                                }                            
                            }
                            productComponent = reductionlessMultiply(productComponent, termComponent);
                            termComponentIndices[k] = 0;
                        }
                        product += productComponent;
                    }
                return product;
            }
            public static MultivariatePolynomial operator *(Rational a, MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(b.MinimalPolynomials);
                for (int i = 0; i < b.Coefficients.Length; ++i)
                    product.Coefficients[i] = (Rational)(a * b.Coefficients[i]);
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a, Rational b)
            {
                return b * a;
            }
            public Polynomial getMinimalPolynomial()
            {            
                MultivariatePolynomial[] powers = new MultivariatePolynomial[Coefficients.Length];
                MultivariatePolynomial power = new MultivariatePolynomial(MinimalPolynomials);
                power.Coefficients[0] = One;
                MultivariatePolynomial[] matrix = new MultivariatePolynomial[Coefficients.Length];
                Polynomial[] augmentation = new Polynomial[Coefficients.Length];
                List<Rational> augmentationRow = new List<Rational> { One };
                for (int i = 0; i < Coefficients.Length; ++i)
                {
                    power *= this;                
                    powers[i] = power;
                    matrix[i] = power;                
                    augmentationRow.Insert(0, Zero);
                    augmentation[i] = new Polynomial(new List<Rational>(augmentationRow));                
                }
                for (int i = 1; i < matrix.Length; ++i)
                    for (int j = i; j < matrix.Length; ++j)
                    {
                        if (matrix[i - 1].Coefficients[Coefficients.Length - i].Equals(Zero))
                        {
                            MultivariatePolynomial tempRow = matrix[j];
                            Polynomial tempAugmentationRow = augmentation[j];
                            matrix[j] = matrix[i - 1];
                            augmentation[j] = augmentation[i - 1];
                            matrix[i - 1] = tempRow;
                            augmentation[i - 1] = tempAugmentationRow;
                        }
                        else
                        {
                            for (int k = i; k < matrix.Length; ++k)
                            {
                                Rational scalar =
                                    (Rational)(matrix[k].Coefficients[Coefficients.Length - i] /
                                    matrix[i - 1].Coefficients[Coefficients.Length - i].negative());
                                matrix[k] += matrix[i - 1] * scalar;
                                augmentation[k] += augmentation[i - 1] * scalar;
                            }
                            break;
                        }
                    }
                List<Rational> annullingPolynomialCoefficients =
                    augmentation[matrix.Length - 1].Coefficients;
                annullingPolynomialCoefficients[0] = (Rational)(
                    annullingPolynomialCoefficients[0] - matrix[matrix.Length - 1].Coefficients[0]);
                Polynomial minimalPolynomial = new Polynomial(annullingPolynomialCoefficients);
                List<Polynomial> factors = minimalPolynomial.getFactors();            
                foreach (Polynomial factor in factors)
                {
                    MultivariatePolynomial sum = new MultivariatePolynomial(MinimalPolynomials);
                    sum.Coefficients[0] = factor.Coefficients[0];
                    for (int i = 1; i < factor.Coefficients.Count; ++i)
                        sum += factor.Coefficients[i] * powers[i - 1];
                    bool sumIsZero = true;
                    foreach (Rational coefficient in sum.Coefficients)
                        if (!coefficient.Equals(Zero)) 
                        {
                            sumIsZero = false;
                            break;
                        }
                    if(sumIsZero)
                    {
                        minimalPolynomial = factor;
                        break;
                    }
                }
                return minimalPolynomial * (Rational)minimalPolynomial.Coefficients[
                    minimalPolynomial.Coefficients.Count - 1].reciprocal();
            }
#if DEBUG
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                for (int i = 0; i < Coefficients.Length; ++i)
                {
                    output.Append(Coefficients[i].ToString());
                    int[] externalCoefficients = getExternalIndices(i);
                    for (int j = 0; j < externalCoefficients.Length; ++j)
                        output.Append("(" + externalCoefficients[j] + ")");
                    output.Append("+");
                }
                return output.ToString();
            }
#endif
        }    
        static Number evaluateExpression(List<char> operations, List<Number> numbers)
        {
            for (int i = 0; i < operations.Count;)
            {
                if (operations[i] == ')')
                    throw new InvalidUserInput("Unmatched parentheses.");
                if (operations[i] == '(')
                {
                    int numOfUnmatchedParens = 1;
                    int matchingParenIndex = i;
                    while (numOfUnmatchedParens > 0)
                    {
                        matchingParenIndex += 1;
                        try
                        {
                            if (operations[matchingParenIndex] == '(')
                                numOfUnmatchedParens += 1;
                            else if (operations[matchingParenIndex] == ')')
                                numOfUnmatchedParens -= 1;
                        }
                        catch (ArgumentOutOfRangeException)
                        {
                            throw new InvalidUserInput("Unmatched parentheses.");
                        }
                    }
                    numbers[i] = evaluateExpression(operations.GetRange(i + 1, matchingParenIndex -
                        i - 1), numbers.GetRange(i + 1, matchingParenIndex - i - 1));
                    operations[i] = ' ';
                    numbers.RemoveRange(i + 1, matchingParenIndex - i);
                    operations.RemoveRange(i + 1, matchingParenIndex - i);
                }
                else
                    ++i;
            }
            for (int i = 0; i < operations.Count;)
            {
                if (i == 0 || numbers[i - 1] == null)
                {
                    if (operations[i] == '+')
                    {
                        numbers[i] = One;
                        numbers.Insert(i + 1, null);
                        operations[i] = ' ';
                        operations.Insert(i + 1, '*');
                    }
                    else if (operations[i] == '-')
                    {
                        numbers[i] = new Integer(-1);
                        numbers.Insert(i + 1, null);
                        operations[i] = ' ';
                        operations.Insert(i + 1, '*');
                    }
                    ++i;
                }
                else
                    ++i;
            }
            try
            {
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
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new InvalidUserInput("Operator missing operand.");
            }
            if (numbers.Count == 0)
                return null;
            return numbers[0];
        }
        public class InvalidUserInput : Exception
        {
            public InvalidUserInput(string message) : base(message)
            { }
        }
        static void Main(string[] args)
        {
            while (true)
            {
                try
                {
#if DEBUG
                    string input = "1/(1+2^(1/3))";
#else
                    string input = Console.ReadLine();
                    if (input[0] == 'q')
                        return;
#endif
                    List<char> operations = new List<char>();
                    List<Number> numbers = new List<Number>();
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
                                    numbers.Add(ComplexNumber.create(Zero,
                                        Integer.parse(numberCollector.ToString())));
                                    operations.Add(' ');
                                }
                                else if (!"()+-*/^".Contains(c.ToString()))
                                    throw new InvalidUserInput(c + " is an invalid character.");
                                else
                                {
                                    numbers.Add(Integer.parse(numberCollector.ToString()));
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
                        else if (c == 'i')
                        {
                            numbers.Add(ComplexNumber.create(Zero, One));
                            operations.Add(' ');
                        }
                        else if (!"()+-*/^".Contains(c.ToString()))
                            throw new InvalidUserInput(c + " is an invalid character.");
                        else
                        {
                            operations.Add(c);
                            numbers.Add(null);
                        }
                    }
                    if (lastCharWasDigit)
                    {
                        operations.Add(' ');
                        numbers.Add(Integer.parse(numberCollector.ToString()));
                    }
                    for (int i = 0; i < operations.Count;)
                    {
                        if (operations[i] == '(' && i > 0 && numbers[i - 1] != null)
                        {
                            operations.Insert(i, '*');
                            numbers.Insert(i, null);
                            i += 2;
                        }
                        else if ((operations[i] == ')' || numbers[i] is ComplexNumber) &&
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
                    Number expression = evaluateExpression(operations, numbers);
                    if (expression == null)
                        Console.WriteLine("=\n");
                    else
                        Console.Write("=\n" + expression.ToString() + "\n\n");
                }
                catch (InvalidUserInput e)
                {
                    Console.WriteLine(e.Message);
                }
                catch (DivideByZeroException e)
                {
                    Console.WriteLine(e.Message);
                }
#if DEBUG
                return;
#endif
            }
        }
    }
}
