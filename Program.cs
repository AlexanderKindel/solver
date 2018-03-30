using System;
using System.Collections.Generic;
using System.Text;

namespace Solver
{
    class Solver
    {
        delegate T multiplier<T>(T a, T b);
        static T integerExponentiate<T>(multiplier<T> multiply, T unit, T expBase, Integer exponent)
        {
            T output = unit;
            T baseToAPowerOfTwo = expBase;
            while (exponent.Sign > 0)
            {
                IntegerDivision division = exponent.euclideanDivideBy(new Integer(2));
                if (division.remainder == One)
                    output = multiply(output, baseToAPowerOfTwo);
                baseToAPowerOfTwo = multiply(baseToAPowerOfTwo, baseToAPowerOfTwo);
                exponent = division.quotient;
            }
            return output;
        }
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
            protected void removeNonConjugates(List<Number> conjugateCandidates)
            {
                for (int i = 0;
                    conjugateCandidates.Count > MinimalPolynomial.Coefficients.Count - 1;)
                    if (!conjugateCandidates[i].MinimalPolynomial.Equals(MinimalPolynomial))
                        conjugateCandidates.RemoveAt(i);
                    else
                        ++i;
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
                    Integer integerExponent = (Integer)exponent;
                    if (integerExponent.Sign < 0)
                        return exponentiate(exponent.negative()).reciprocal();
                    return integerExponentiate(delegate (Number a, Number b) { return a * b; },
                        One, this, integerExponent);
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
            public Integer magnitude()
            {
                return new Integer(Values, 1);
            }
            static uint[] shiftValuesLeft(uint[] values, int valuePlaces, int digitPlaces)
            {//Negative valuePlaces and digitPlaces values shift to the right. Left shifts preserve
             //all digits by adding extra uints, while right shifts truncate the rightmost digits.
                uint[] shiftedValues;
                int smallestValuePlace;
                if (valuePlaces < 0)
                    smallestValuePlace = -valuePlaces;
                else
                    smallestValuePlace = 0;
                if (digitPlaces == 0)
                {
                    shiftedValues = new uint[values.Length + valuePlaces];
                    for (int i = smallestValuePlace; i < values.Length; ++i)
                        shiftedValues[valuePlaces + i] = values[i];
                }
                else if (digitPlaces > 0)
                {
                    shiftedValues = new uint[values.Length + valuePlaces + 1];
                    for (int i = smallestValuePlace; i < values.Length; ++i)
                    {
                        shiftedValues[valuePlaces + i] += values[i] << digitPlaces;
                        shiftedValues[valuePlaces + i + 1] += values[i] >> 32 - digitPlaces;
                    }
                }
                else
                {
                    shiftedValues = new uint[values.Length + valuePlaces];
                    for (int i = smallestValuePlace; i < values.Length; ++i)
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
                IntegerDivision division = new IntegerDivision();
                division.remainder = new Integer(Values, 1);
                if (divisor.Values.Length > Values.Length ||
                    (divisor.Values.Length == Values.Length &&
                    divisor.Values[divisor.Values.Length - 1] > Values[Values.Length - 1])) 
                {
                    division.quotient = Zero;
                    return division;
                }
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
                uint[] quotient = new uint[Values.Length];
                void calculateValue(int valuePlace, int stoppingDigitPlace)
                {
                    power = 0b10000000000000000000000000000000;
                    for (int i = 32; i >= stoppingDigitPlace; --i)
                    {
                        Integer shiftedDivisor = positiveDivisor.shiftLeft(valuePlace -
                            positiveDivisor.Values.Length + 1, i - divisorLeadingDigitPlace);
                        Integer difference = division.remainder - shiftedDivisor;
                        if (difference.Sign >= 0)
                        {
                            quotient[valuePlace] |= power;
                            division.remainder = difference;
                        }
                        power = power >> 1;
                    }
                }
                for (int i = Values.Length - 1; i >= positiveDivisor.Values.Length; --i)
                    calculateValue(i, 1);
                calculateValue(positiveDivisor.Values.Length - 1, divisorLeadingDigitPlace);                
                division.quotient = new Integer(quotient, (sbyte)(Sign * divisor.Sign)).shiftLeft(
                    1 - positiveDivisor.Values.Length, 1 - divisorLeadingDigitPlace);
                if (Sign < 0 && division.remainder.Sign != 0) 
                    division.remainder = positiveDivisor - division.remainder;
                return division;
            }
            protected override Polynomial calculateMinimalPolynomial()
            {
                return new Polynomial(new List<Integer> { -this, One });
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
                return (a - b).Sign < 0;
            }
            public static bool operator >(Integer a, Integer b)
            {
                return (a - b).Sign > 0;
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
                return a.magnitude();
            }
            public static Integer getGCD(List<Integer> list)
            {
                Integer GCD = list[0];
                for (int i = 1; i < list.Count; ++i)
                    GCD = getGCD(GCD, list[i]);
                return GCD;
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
            public Integer thisChooseK(Integer k)
            {
                Integer numerator = One;
                for (Integer i = this - k + One; i <= this; ++i)
                    numerator *= i;
                Integer denominator = One;
                for (Integer i = new Integer(2); i <= k; ++i)
                    denominator *= i;
                return (Integer)(numerator / denominator);
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
            protected override Polynomial calculateMinimalPolynomial()
            {
                return new Polynomial(new List<Integer> { -Numerator, Denominator });
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
                List<Integer> coefficients =
                    new List<Integer> { baseMinimalPolynomial.Coefficients[0] };
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
                removeNonConjugates(conjugateCandidates);
                return conjugateCandidates;
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
            public ComplexExponential(Number exponent)
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
                    else if (denominator == new Integer(5))
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
                List<Integer> coefficients = new List<Integer> { new Integer(-1), One };
                List<Polynomial> dividends = new List<Polynomial>();
                Integer index = ((Fraction)Exponent).Denominator;
                for (Integer i = One; i < index; ++i)
                {
                    dividends.Add(new Polynomial(new List<Integer>(coefficients)));
                    coefficients.Insert(1, Zero);
                }
                Polynomial annullingPolynomial = dividends[dividends.Count - 1];
                dividends.RemoveAt(dividends.Count - 1);
                List<Polynomial> factors = annullingPolynomial.getFactors();
                foreach (Polynomial factor in factors)
                {
                    bool isDivisor = false;
                    foreach (Polynomial dividend in dividends)
                        if ((dividend % factor).Coefficients.Count == 0)
                        {
                            isDivisor = true;
                            break;
                        }
                    if (!isDivisor)
                        return factor;
                }
                throw new Exception("This function should have returned during the loop above.");
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
                    for (int i = 0; i < factors.Count; ++i)
                    {
                        Number product = number * factors[i];
                        if (!product.Equals(new Product(Coefficient,
                            new List<Factor> { (Factor)number, factors[i] })))
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
                MultivariatePolynomial variableForm =
                    new MultivariatePolynomial(minimalPolynomials);
                int[] indices = new int[minimalPolynomials.Length];
                for (int i = 0; i < indices.Length; ++i)
                    indices[i] = 1;
                variableForm.setCoefficient(indices, Coefficient);
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
                List<List<Number>> conjugateCombinations =
                    generateCartesianProduct(factorConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = One;
                    foreach (Number number in combination)
                        conjugate *= number;
                    conjugateCandidates.Add(conjugate);
                }
                removeNonConjugates(conjugateCandidates);
                return conjugateCandidates;
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
                foreach (Number conjugate in conjugates)
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
                MultivariatePolynomial variableForm =
                    new MultivariatePolynomial(minimalPolynomials);
                variableForm.setCoefficient(new int[minimalPolynomials.Length], constant);
                for (int i = 0; i < minimalPolynomials.Length; ++i)
                {
                    int[] indices = new int[minimalPolynomials.Length];
                    indices[i] = 1;
                    variableForm.setCoefficient(indices, One);
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
                    termConjugates.Add(conjugates);
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
                removeNonConjugates(conjugateCandidates);
                return conjugateCandidates;
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
                        new Polynomial(new List<Integer> { One, Zero, One })});
                variableForm.setCoefficient(new int[] { 1, 0, 0 }, One);
                variableForm.setCoefficient(new int[] { 0, 1, 1 }, One);
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
                removeNonConjugates(conjugateCandidates);
                return conjugateCandidates;
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
        static List<Integer> primes = new List<Integer> { new Integer(2), new Integer(3) };
        class Polynomial : IComparable
        {
            //The index of the coefficient represents the degree of its term.
            public List<Integer> Coefficients { get; }
            Integer Characteristic = Zero;
            public Polynomial(List<Integer> coefficients)
            {
                while (coefficients.Count != 0 && coefficients[coefficients.Count - 1].Equals(Zero))
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
            }
            public Polynomial(List<Integer> coefficients, Integer characteristic)
            {
                Coefficients = new List<Integer>();
                Characteristic = characteristic;
                if (!Characteristic.Equals(Zero))
                    foreach (Integer coefficient in coefficients)
                        Coefficients.Add(coefficient.euclideanDivideBy(characteristic).remainder);
                else
                    Coefficients = coefficients;
                while (Coefficients.Count != 0 && Coefficients[Coefficients.Count - 1].Equals(Zero))
                    Coefficients.RemoveAt(Coefficients.Count - 1);
            }
            public static Polynomial operator +(Polynomial a, Polynomial b)
            {
                if (a.Coefficients.Count < b.Coefficients.Count)
                    return b + a;
                List<Integer> output = a.Coefficients;
                for (int i = 0; i < b.Coefficients.Count; ++i)
                    output[i] = output[i] + b.Coefficients[i];
                return new Polynomial(output, a.Characteristic);
            }
            public Polynomial negative()
            {
                List<Integer> output = new List<Integer>();
                foreach (Integer coefficient in Coefficients)
                    output.Add(-coefficient);
                return new Polynomial(output, Characteristic);
            }
            public static Polynomial operator -(Polynomial a, Polynomial b)
            {
                return a + b.negative();
            }
            public static Polynomial operator *(Polynomial a, Polynomial b)
            {
                List<Integer> output = new List<Integer>();
                for (int i = 0; i < a.Coefficients.Count + b.Coefficients.Count - 1; ++i)
                    output.Add(Zero);
                for (int i = 0; i < a.Coefficients.Count; ++i)
                    for (int j = 0; j < b.Coefficients.Count; ++j)
                        output[i + j] += a.Coefficients[i] * b.Coefficients[j];
                return new Polynomial(output, a.Characteristic);
            }            
            public static Polynomial operator *(Polynomial a, Integer b)
            {
                return a * new Polynomial(new List<Integer> { b });
            }
            public static Polynomial operator *(Integer a, Polynomial b)
            {
                return b * a;
            }
            struct PolynomialDivision
            {
                public Polynomial Quotient;
                public Polynomial Remainder;
            }
            delegate Number Divider(Integer a, Integer b);
            PolynomialDivision divideBy(Polynomial divisor)
            {
                Divider divide;
                if (Characteristic == Zero)
                    divide = delegate (Integer a, Integer b) { return a / b; };
                else
                    divide = delegate (Integer a, Integer b) {
                        return a * b.exponentiate(Characteristic - new Integer(2)); };
                PolynomialDivision division;
                List<Integer> quotient = new List<Integer>();
                List<Integer> remainder = new List<Integer>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        Number quotientCoefficient = divide(remainder[remainder.Count - 1],
                            divisor.Coefficients[divisor.Coefficients.Count - 1]);
                        if (!(quotientCoefficient is Integer))
                        {
                            division.Quotient = null;
                            division.Remainder = null;
                            return division;
                        }
                        quotient.Insert(0, (Integer)quotientCoefficient);
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            remainder[remainder.Count - j] =
                                remainder[remainder.Count - j] - quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                    }
                division.Quotient = new Polynomial(quotient, Characteristic);
                division.Remainder = new Polynomial(remainder, Characteristic);
                return division;
            }
            public static Polynomial operator /(Polynomial a, Polynomial b)
            {
                return a.divideBy(b).Quotient;
            }
            public static Polynomial operator %(Polynomial a, Polynomial b)
            {
                return a.divideBy(b).Remainder;
            }
            public Polynomial exponentiate(Integer exponent)
            {
                return integerExponentiate(delegate (Polynomial A, Polynomial B) { return A * B; },
                    new Polynomial(new List<Integer> { One }, Characteristic), this, exponent);
            }
            public Integer evaluateAt(Integer input)
            {
                Integer output = Coefficients[0];
                for (int i = 1; i < Coefficients.Count; ++i)
                    output =
                        (Integer)(output + Coefficients[i] * input.exponentiate(new Integer(i)));
                return output;
            }
            public Polynomial getDerivative()
            {
                List<Integer> derivativeCoefficients = new List<Integer>();
                for (int i = 1; i < Coefficients.Count; ++i)
                    derivativeCoefficients.Add(Coefficients[i] * new Integer(i));
                return new Polynomial(derivativeCoefficients, Characteristic);
            }
            Integer reduceToPrimitivePart()
            {//Returns the content used in the reduction.
                Integer content = Integer.getGCD(Coefficients);
                for (int i = 0; i < Coefficients.Count; ++i)
                    Coefficients[i] = (Integer)(Coefficients[i] / content);
                return content;
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
            static Polynomial getGCD(Polynomial a, Polynomial b)
            {
                Polynomial A;
                Polynomial B;
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    A = new Polynomial(new List<Integer>(b.Coefficients), b.Characteristic);
                    B = new Polynomial(new List<Integer>(a.Coefficients), a.Characteristic);
                }
                else
                {
                    A = new Polynomial(new List<Integer>(a.Coefficients), a.Characteristic);
                    B = new Polynomial(new List<Integer>(b.Coefficients), b.Characteristic);
                }
                if (b.Coefficients.Count == 0)
                    return a;
                Integer d = Integer.getGCD(A.reduceToPrimitivePart(), B.reduceToPrimitivePart());
                Integer g = One;
                Number h = One;
                Integer degree = new Integer(A.Coefficients.Count - B.Coefficients.Count);
                Polynomial remainder = ((Integer)B.Coefficients[B.Coefficients.Count - 1].
                    exponentiate(degree + One) * A).divideBy(B).Remainder;
                while (remainder.Coefficients.Count > 1)
                {
                    A = B;
                    Number divisor = (g * h).exponentiate(degree);
                    B = remainder;
                    for (int i = 0; i < B.Coefficients.Count; ++i)
                        B.Coefficients[i] = (Integer)(B.Coefficients[i] / divisor);
                    g = A.Coefficients[A.Coefficients.Count - 1];
                    h = h.exponentiate(One - degree) * g.exponentiate(degree);
                    degree = new Integer(A.Coefficients.Count - B.Coefficients.Count);
                    remainder = ((Integer)B.Coefficients[B.Coefficients.Count - 1].exponentiate(
                        degree + One) * A).divideBy(B).Remainder;
                }
                if (remainder.Coefficients.Count == 1)
                    return new Polynomial(new List<Integer> { d }, a.Characteristic);
                B.reduceToPrimitivePart();
                return d * B;
            }
            public List<Polynomial> getFactors()
            {
                List<Polynomial> irreducibleFactors = new List<Polynomial>();
                if (Characteristic != Zero)
                {//Assumes *this is square-free, which is currently admissible because the only time
                 //getFactors() is called on a Polynomial over a field of nonzero characteristic is
                 //by itself, when called on Polynomial over the integers, and the Polynomial
                 //produced in that context is guaranteed to be square-free.
                    Dictionary<int, Polynomial> distinctDegreeFactors =
                        new Dictionary<int, Polynomial>();
                    Polynomial V = this;
                    Polynomial W = new Polynomial(new List<Integer> { Zero, One }, Characteristic);
                    int d = 0;
                    while (d < (Coefficients.Count + 1) / 2)
                    {
                        ++d;
                        W = W.exponentiate(Characteristic).divideBy(this).Remainder;
                        Polynomial degreeDFactorProduct =
                            getGCD(W - new Polynomial(new List<Integer> { Zero, One }), V);
                        if (degreeDFactorProduct.Coefficients.Count > 1)
                        {
                            distinctDegreeFactors.Add(d, degreeDFactorProduct);
                            V = V.divideBy(degreeDFactorProduct).Quotient;
                            W = W.divideBy(V).Remainder;
                        }
                    }
                    void CZSplit(Polynomial factorProduct, int degree)
                    {
                        if ((factorProduct.Coefficients.Count - 1) / degree == 1)
                        {
                            irreducibleFactors.Add(factorProduct);
                            return;
                        }
                        Polynomial B;
                        if (Characteristic == new Integer(2))
                        {
                            Polynomial T =
                                new Polynomial(new List<Integer> { Zero, One }, Characteristic);
                            Polynomial C =
                                new Polynomial(new List<Integer>(T.Coefficients), Characteristic);
                            for (int j = 1; j < degree; ++j)
                                C = T + (C * C).divideBy(factorProduct).Remainder;
                            B = getGCD(factorProduct, C);
                            while (B.Coefficients.Count < 2 ||
                                B.Coefficients.Count == factorProduct.Coefficients.Count)
                            {
                                T *= new Polynomial(new List<Integer> { Zero, Zero, One },
                                    Characteristic);
                                C = new Polynomial(new List<Integer>(T.Coefficients),
                                    Characteristic);
                                for (int j = 1; j < degree; ++j)
                                    C = T + (C * C).divideBy(factorProduct).Remainder;
                                B = getGCD(factorProduct, C);
                            }
                        }
                        else
                        {
                            List<Integer> coefficients = new List<Integer>();
                            Integer p = Characteristic - One;
                            for (int j = 1; j < 2 * degree; ++j)
                                coefficients.Add(p);
                            coefficients.Add(One);
                            Integer power = (Integer)((Characteristic.exponentiate(
                                new Integer(degree)) - One) / new Integer(2));
                            B = getGCD(factorProduct, new Polynomial(new List<Integer>(
                                coefficients), Characteristic).exponentiate(power) -
                                new Polynomial(new List<Integer> { One }, Characteristic));
                            while (B.Coefficients.Count < 2 ||
                                B.Coefficients.Count == factorProduct.Coefficients.Count)
                            {
                                int j = coefficients.Count - 2;
                                while (coefficients[j] == Zero)
                                {
                                    coefficients[j] = p;
                                    --j;
                                }
                                --coefficients[j];
                                B = getGCD(factorProduct, new Polynomial(new List<Integer>(
                                    coefficients), Characteristic).exponentiate(power) -
                                    new Polynomial(new List<Integer> { One }, Characteristic));
                            }
                        }
                        CZSplit(B, degree);
                        CZSplit(factorProduct.divideBy(B).Quotient, degree);
                    }
                    foreach (int degree in distinctDegreeFactors.Keys)
                        CZSplit(distinctDegreeFactors[degree], degree);
                    return irreducibleFactors;
                }
                List<Polynomial> squarefreeFactors = new List<Polynomial>();
                Polynomial derivative = getDerivative();
                Polynomial a = getGCD(this, derivative);
                Polynomial b = this / a;
                Polynomial c = derivative / a - b.getDerivative();
                while (!(b.Coefficients.Count == 1 && (b.Coefficients[0].Equals(One)
                    || b.Coefficients[0].Equals(new Integer(-1)))))
                {
                    a = getGCD(b, c);
                    if (!squarefreeFactors.Contains(a))
                        squarefreeFactors.Add(a);
                    b = b / a;
                    c = c / a - b.getDerivative();
                }
                List<Polynomial> splitFactor(Polynomial factor)
                {
                    Polynomial moddedFactor;
                    int i = 0;
                    while (true)
                    {
                        if (i == primes.Count)
                        {
                            Integer primeCandidate = primes[primes.Count - 1] + new Integer(2);
                            while (true)
                            {
                                bool isDivisible = false;
                                for (int j = 0; primes[j] <=
                                    primeCandidate.euclideanDivideBy(new Integer(2)).quotient; ++j)
                                    if (primeCandidate.euclideanDivideBy(
                                        primes[j]).remainder.Equals(Zero))
                                    {
                                        isDivisible = true;
                                        break;
                                    }
                                if (!isDivisible)
                                    break;
                                primeCandidate += new Integer(2);
                            }
                            primes.Add(primeCandidate);
                        }
                        moddedFactor = new Polynomial(factor.Coefficients, primes[i]);
                        Polynomial GCD = getGCD(factor, factor.getDerivative());
                        if (GCD.Coefficients.Count == 1 && GCD.Coefficients[0] == One)
                        {
                            Integer leadingCoefficientInverse = (Integer)
                                moddedFactor.Coefficients[moddedFactor.Coefficients.Count -
                                1].exponentiate(moddedFactor.Characteristic - new Integer(2));
                            for (int j = 0; j < moddedFactor.Coefficients.Count; ++j)
                                moddedFactor.Coefficients[j] *= leadingCoefficientInverse;
                            break;
                        }
                        ++i;
                    }
                    List<Polynomial> irreducibleModdedFactors = moddedFactor.getFactors();
                    Integer bound = Zero;
                    foreach (Integer coefficient in factor.Coefficients)
                        bound += coefficient * coefficient;
                    Integer squareRoot = One;
                    while (squareRoot * squareRoot < bound)
                        ++squareRoot;
                    Integer factorDegree = new Integer(factor.Coefficients.Count - 1);
                    Integer k = factorDegree.euclideanDivideBy(new Integer(4)).quotient;
                    bound = (squareRoot * factorDegree.thisChooseK(k) +
                        factor.Coefficients[factor.Coefficients.Count - 1].magnitude() *
                        factorDegree.thisChooseK(k - One)) * new Integer(2) *
                        factor.Coefficients[factor.Coefficients.Count - 1].magnitude();
                    Integer e = One;
                    Integer characteristicPower = moddedFactor.Characteristic;
                    while (characteristicPower < bound)
                    {
                        ++e;
                        characteristicPower *= moddedFactor.Characteristic;
                    }
                    List<Polynomial> liftedFactors = new List<Polynomial>();
                    Polynomial productOfUnliftedFactors =
                        new Polynomial(factor.Coefficients, moddedFactor.Characteristic);
                    for (int j = 0; j < irreducibleModdedFactors.Count - 1; ++j)
                    {
                        productOfUnliftedFactors /= irreducibleModdedFactors[j];
                        Integer characteristicExponent = Zero;
                        Integer characteristicToPower = One;
                        Polynomial A = irreducibleModdedFactors[j];
                        Polynomial B = productOfUnliftedFactors;
                        A.Characteristic = Zero;
                        B.Characteristic = Zero;
                        while (characteristicExponent < e)
                        {                            
                            ++characteristicExponent;
                            characteristicToPower *= moddedFactor.Characteristic;
                            Polynomial ACoefficient = new Polynomial(new List<Integer> { One },
                                moddedFactor.Characteristic);
                            Polynomial GCD = A;
                            Polynomial BCoefficient =
                                new Polynomial(new List<Integer>(), moddedFactor.Characteristic);
                            Polynomial Y = B;
                            PolynomialDivision division;
                            while (Y.Coefficients.Count > 0)
                            {
                                division = GCD.divideBy(Y);
                                Polynomial T = ACoefficient - BCoefficient * division.Quotient;
                                ACoefficient = BCoefficient;
                                GCD = Y;
                                BCoefficient = T;
                                Y = division.Remainder;
                            }
                            BCoefficient = (GCD - A * ACoefficient) / productOfUnliftedFactors;
                            List<Integer> fCoefficients = (factor - A * B).Coefficients;
                            for (int l = 0; l < fCoefficients.Count; ++l)
                                fCoefficients[l] =
                                    (Integer)(fCoefficients[l] / characteristicToPower);
                            Polynomial f =
                                new Polynomial(fCoefficients, moddedFactor.Characteristic);
                            division = (BCoefficient * f).divideBy(A);
                            division.Quotient.Characteristic = Zero;
                            division.Remainder.Characteristic = Zero;
                            A += characteristicToPower * division.Remainder;
                            B += characteristicToPower * (ACoefficient * f + B * division.Quotient);
                        }
                        liftedFactors.Add(A);
                    }
                    liftedFactors.Add(productOfUnliftedFactors);
                    List<Polynomial> factorSplit = new List<Polynomial>();
                    int d = 1;
                    List<int> testFactorCombination(int elementsToAdd,
                        int indexOfPreviousElement, List<int> combinationIndices)
                    {
                        if (elementsToAdd == 0) 
                        {
                            List<int> finalizedCombinationIndices;
                            Polynomial V = new Polynomial(new List<Integer> {
                                factor.Coefficients[factor.Coefficients.Count - 1] });
                            if (2 * combinationIndices.Count >= Coefficients.Count)
                            {
                                finalizedCombinationIndices = new List<int>();
                                for (int j = liftedFactors.Count - 1; j >= 0; --j)
                                    if (j == combinationIndices[combinationIndices.Count - 1])
                                        combinationIndices.RemoveAt(combinationIndices.Count - 1);
                                    else
                                        finalizedCombinationIndices.Add(j);
                            }
                            else
                                finalizedCombinationIndices = combinationIndices;
                            foreach (int index in finalizedCombinationIndices)
                                V *= liftedFactors[index];
                            Integer two = new Integer(2);
                            for (int j = 0; j < V.Coefficients.Count; ++j)
                            {
                                Integer remainder = V.Coefficients[j].euclideanDivideBy(
                                    characteristicPower).remainder;
                                if (remainder * two <= characteristicPower)
                                    V.Coefficients[j] = remainder;
                                else
                                    V.Coefficients[j] = remainder - characteristicPower;
                            }
                            PolynomialDivision division = factor.divideBy(V);
                            if (division.Remainder.Coefficients.Count == 0)
                            {
                                while (division.Remainder.Coefficients.Count == 0) 
                                {
                                    factor = division.Quotient;
                                    division = factor.divideBy(V);
                                }
                                V.reduceToPrimitivePart();
                                factorSplit.Add(V);
                                return finalizedCombinationIndices;
                            }
                            return new List<int>();
                        }
                        for (int j = indexOfPreviousElement + 1;
                            j < liftedFactors.Count - elementsToAdd; ++j)
                        {
                            List<int> enlargedCombinationIndices =
                                new List<int>(combinationIndices);
                            enlargedCombinationIndices.Add(j);
                            enlargedCombinationIndices = testFactorCombination(elementsToAdd - 1,
                                j, enlargedCombinationIndices);
                            if (enlargedCombinationIndices.Count != 0)
                                return enlargedCombinationIndices;
                        }
                        return new List<int>();
                    }
                    List<int> factorsToRemove;
                    while (2 * d < liftedFactors.Count)   
                    {
                        factorsToRemove = testFactorCombination(d, -1, new List<int>());
                        while (factorsToRemove.Count != 0)
                        {
                            foreach (int index in factorsToRemove)
                                liftedFactors.RemoveAt(index);
                            factorsToRemove = testFactorCombination(d, -1, new List<int>());
                        }
                        ++d;
                    }
                    if (2 * d == liftedFactors.Count)
                    {
                        factorsToRemove = testFactorCombination(d, 0, new List<int> { 0 });
                        foreach (int index in factorsToRemove)
                            liftedFactors.RemoveAt(index);
                    }
                    Polynomial finalFactor = new Polynomial(new List<Integer> { One });
                    foreach (Polynomial liftedFactor in liftedFactors)
                        finalFactor *= liftedFactor;
                    if (finalFactor.Coefficients.Count > 1)
                    {
                        finalFactor.reduceToPrimitivePart();
                        factorSplit.Add(finalFactor);
                    }
                    return factorSplit;
                }
                foreach (Polynomial factor in squarefreeFactors)
                {
                    Polynomial X = new Polynomial(new List<Integer> { Zero, One });
                    if (factor.Coefficients[0].Equals(Zero))
                    {
                        irreducibleFactors.Add(X);
                        factor.Coefficients.RemoveAt(0);
                        while (factor.Coefficients[0].Equals(Zero))
                            factor.Coefficients.RemoveAt(0);
                    }
                    if (factor.Coefficients[0].magnitude() <
                        factor.Coefficients[factor.Coefficients.Count - 1].magnitude())
                    {
                        factor.Coefficients.Reverse();
                        List<Polynomial> splitFactors = splitFactor(factor);
                        for (int i = 0; i < splitFactors.Count; ++i)
                        {
                            splitFactors[i].Coefficients.Reverse();
                            irreducibleFactors.Add(splitFactors[i]);
                        }
                    }
                    else
                        irreducibleFactors.AddRange(splitFactor(factor));
                }
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
            public Dictionary<int[], Rational> Coefficients { get; }
            public MultivariatePolynomial(Polynomial[] minimalPolynomials)
            {//The nth minimalPolynomial entry is the minimal polynomial of the number the nth
             //variable represents.
                MinimalPolynomials = minimalPolynomials;
                Coefficients = new Dictionary<int[], Rational>();
            }
            public void setCoefficient(int[] indices, Rational coefficient)
            {//The nth index is the degree with respect to the nth variable of the term being set.
                foreach (int[] keyIndices in Coefficients.Keys)
                    if (areEqualByValue(indices, keyIndices))
                    {
                        Coefficients[indices] = coefficient;
                        return;
                    }
                Coefficients.Add(indices, coefficient);
            }
            public static MultivariatePolynomial operator +(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial sum = new MultivariatePolynomial(a.MinimalPolynomials);
                foreach (int[] indices in a.Coefficients.Keys)
                    sum.Coefficients.Add(indices, a.Coefficients[indices]);
                foreach (int[] indicesB in b.Coefficients.Keys)
                {
                    bool indicesFound = false;
                    foreach (int[] indicesA in sum.Coefficients.Keys)
                        if (areEqualByValue(indicesA, indicesB))
                        {
                            Rational termSum =
                                (Rational)(sum.Coefficients[indicesA] + b.Coefficients[indicesB]);
                            if (termSum.Equals(Zero))
                                sum.Coefficients.Remove(indicesA);
                            else
                                sum.Coefficients[indicesA] = termSum;
                            indicesFound = true;
                            break;
                        }
                    if (!indicesFound)
                        sum.Coefficients.Add(indicesB, b.Coefficients[indicesB]);
                }
                return sum;
            }
            static MultivariatePolynomial reductionlessMultiply(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial product =
                    new MultivariatePolynomial(a.MinimalPolynomials);
                foreach (int[] indicesA in a.Coefficients.Keys)
                    foreach (int[] indicesB in b.Coefficients.Keys)
                    {
                        int[] productIndices = new int[indicesA.Length];
                        for (int i = 0; i < productIndices.Length; ++i)
                            productIndices[i] += indicesB[i] + indicesA[i];
                        bool indicesFound = false;
                        foreach (int[] keyIndices in product.Coefficients.Keys)
                            if (areEqualByValue(productIndices, keyIndices))
                            {
                                product.Coefficients[productIndices] =
                                    (Rational)(product.Coefficients[productIndices] +
                                    a.Coefficients[indicesA] * b.Coefficients[indicesB]);
                                indicesFound = true;
                                break;
                            }
                        if (!indicesFound)
                            product.Coefficients.Add(productIndices,
                                (Rational)(a.Coefficients[indicesA] * b.Coefficients[indicesB]));
                    }
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(a.MinimalPolynomials);
                foreach (int[] indicesA in a.Coefficients.Keys)
                    foreach (int[] indicesB in b.Coefficients.Keys)
                    {
                        MultivariatePolynomial productComponent =
                            new MultivariatePolynomial(a.MinimalPolynomials);
                        productComponent.Coefficients.Add(new int[indicesA.Length],
                            (Rational)(a.Coefficients[indicesA] * b.Coefficients[indicesB]));
                        for (int i = 0; i < indicesA.Length; ++i)
                        {
                            MultivariatePolynomial termComponent =
                                new MultivariatePolynomial(a.MinimalPolynomials);
                            int index = indicesA[i] + indicesB[i];
                            if (index < a.MinimalPolynomials[i].Coefficients.Count - 1)
                            {
                                int[] termComponentIndices = new int[indicesA.Length];
                                termComponentIndices[i] = index;
                                termComponent.Coefficients.Add(termComponentIndices, One);
                            }
                            else
                            {
                                index -= a.MinimalPolynomials[i].Coefficients.Count - 1;
                                if (index != 0)
                                {
                                    int[] termComponentIndices = new int[indicesA.Length];
                                    termComponentIndices[i] = index;
                                    MultivariatePolynomial reducedDegreeFactor =
                                        new MultivariatePolynomial(a.MinimalPolynomials);
                                    reducedDegreeFactor.setCoefficient(termComponentIndices, One);
                                    productComponent = reductionlessMultiply(productComponent,
                                        reducedDegreeFactor);
                                }
                                for (int j = 0; j < a.MinimalPolynomials[i].Coefficients.Count - 1;
                                    ++j)
                                {
                                    if (a.MinimalPolynomials[i].Coefficients[j].Equals(Zero))
                                        continue;
                                    int[] termComponentIndices = new int[indicesA.Length];
                                    termComponentIndices[i] = j;
                                    termComponent.Coefficients[termComponentIndices] =
                                        (Rational)(a.MinimalPolynomials[i].Coefficients[j] /
                                        a.MinimalPolynomials[i].Coefficients[
                                        a.MinimalPolynomials[i].Coefficients.Count - 1].negative());
                                }
                            }
                            productComponent =
                                reductionlessMultiply(productComponent, termComponent);
                        }
                        product += productComponent;
                    }
                return product;
            }
            public static MultivariatePolynomial operator *(Rational a, MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(b.MinimalPolynomials);
                if (a.Equals(Zero))
                    return product;
                foreach (int[] indices in b.Coefficients.Keys)
                    product.Coefficients.Add(indices, (Rational)(a * b.Coefficients[indices]));
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a, Rational b)
            {
                return b * a;
            }
            public Polynomial getMinimalPolynomial()
            {
                List<MultivariatePolynomial> powers = new List<MultivariatePolynomial>();
                MultivariatePolynomial power = new MultivariatePolynomial(MinimalPolynomials);
                power.Coefficients[new int[MinimalPolynomials.Length]] = One;
                List<List<Rational>> matrix = new List<List<Rational>>();
                List<int[]> termsPresent = new List<int[]> { new int[MinimalPolynomials.Length] };
                List<Polynomial> augmentation = new List<Polynomial>();
                List<Integer> augmentationRow = new List<Integer> { One };
                bool constantIsPresent = false;
                while (!constantIsPresent || powers.Count < termsPresent.Count)
                {
                    power *= this;
                    powers.Add(power);
                    List<Rational> matrixRow = new List<Rational>();
                    for (int i = 0; i < termsPresent.Count; ++i)
                        matrixRow.Add(Zero);
                    foreach (int[] indices in power.Coefficients.Keys)
                    {
                        bool indicesFound = false;
                        for (int i = 0; i < termsPresent.Count; ++i)
                            if (areEqualByValue(indices, termsPresent[i]))
                            {
                                matrixRow[i] = power.Coefficients[indices];
                                indicesFound = true;
                                break;
                            }
                        if (!indicesFound)
                        {
                            termsPresent.Add(indices);
                            for (int i = 0; i < matrix.Count; ++i)
                                matrix[i].Add(Zero);
                            matrixRow.Add(power.Coefficients[indices]);
                        }
                    }
                    matrix.Add(matrixRow);
                    if (!matrixRow[0].Equals(Zero))
                        constantIsPresent = true;
                    augmentationRow.Insert(0, Zero);
                    augmentation.Add(new Polynomial(new List<Integer>(augmentationRow)));
                }
                for (int i = 1; i < matrix.Count; ++i)
                    for (int j = i; j < matrix.Count; ++j)
                    {
                        if (matrix[i - 1][matrix.Count - i].Equals(Zero))
                        {
                            List<Rational> tempRow = matrix[j];
                            Polynomial tempAugmentationRow = augmentation[j];
                            matrix[j] = matrix[i - 1];
                            augmentation[j] = augmentation[i - 1];
                            matrix[i - 1] = tempRow;
                            augmentation[i - 1] = tempAugmentationRow;
                        }
                        else
                        {
                            for (int k = i; k < matrix.Count; ++k)
                            {
                                Integer scalarA =
                                    matrix[k][matrix.Count - i].getDenominatorLCM() *
                                    matrix[i - 1][matrix.Count - i].getGreatestIntegerFactor();
                                Integer scalarB =
                                    matrix[k][matrix.Count - i].getGreatestIntegerFactor() *
                                    matrix[i - 1][matrix.Count - i].getDenominatorLCM();
                                for (int l = 0; l < matrix.Count; ++l)
                                    matrix[k][l] = (Rational)(
                                        matrix[k][l] * scalarA - matrix[i - 1][l] * scalarB);
                                augmentation[k] =
                                    augmentation[k] * scalarA - augmentation[i - 1] * scalarB;
                            }
                            break;
                        }
                    }
                List<Integer> annullingPolynomialCoefficients = augmentation[matrix.Count - 1].Coefficients;
                for (int i = 0; i < annullingPolynomialCoefficients.Count; ++i)
                    annullingPolynomialCoefficients[i] *= matrix[matrix.Count - 1][0].getDenominatorLCM();
                annullingPolynomialCoefficients[0] -= matrix[matrix.Count - 1][0].getGreatestIntegerFactor();
                Polynomial minimalPolynomial = new Polynomial(annullingPolynomialCoefficients);
                List<Polynomial> factors = minimalPolynomial.getFactors();
                foreach (Polynomial factor in factors)
                {
                    MultivariatePolynomial sum = new MultivariatePolynomial(MinimalPolynomials);
                    sum.Coefficients.Add(new int[MinimalPolynomials.Length],
                        factor.Coefficients[0]);
                    for (int i = 1; i < factor.Coefficients.Count; ++i)
                        sum += factor.Coefficients[i] * powers[i - 1];
                    if (sum.Coefficients.Keys.Count == 0)
                    {
                        minimalPolynomial = factor;
                        break;
                    }
                }
                return minimalPolynomial;
            }
            static bool areEqualByValue(int[] a, int[] b)
            {
                for (int i = 0; i < a.Length; ++i)
                    if (a[i] != b[i])
                        return false;
                return true;
            }
#if DEBUG
            public override string ToString()
            {
                StringBuilder output = new StringBuilder();
                foreach (int[] indices in Coefficients.Keys)
                {
                    output.Append("{");
                    foreach (int index in indices)
                        output.Append(index + ",");
                    output.Append("}(" + Coefficients[indices].ToString() + ")+");
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
                    string input = "1/(1+2^(1/2))";
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
