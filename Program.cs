using System;
using System.Collections.Generic;
using System.Text;

namespace Solver
{
    class Solver
    {
        static List<List<T>> GenerateCartesianProduct<T>(List<List<T>> sets)
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
        static T Exponentiate<T>(T expBase, Integer exponent) where T : IArithmetic<T>
        {
            T output = expBase.GetMultiplicativeIdentity();
            T baseToAPowerOfTwo = expBase;
            while (exponent.Sign > 0)
            {
                Division<Integer> division = exponent.EuclideanDivideBy(new Integer(2));
                if (division.Remainder == One)
                    output = output.Times(baseToAPowerOfTwo);
                baseToAPowerOfTwo = baseToAPowerOfTwo.Times(baseToAPowerOfTwo);
                exponent = division.Quotient;
            }
            return output;
        }
        struct Division<T>
        {
            public T Quotient;
            public T Remainder;
        }
        struct ExtendedGCDInfo<T>
        {
            public T GCD;
            public T ACoefficient;
            public T BCoefficient;
            public T AOverGCD;
            public T BOverGCD;
        }
        static ExtendedGCDInfo<T> ExtendedGCD<T>(T a, T b) where T :
            IArithmetic<T>, IDivisible<T>
        {
            ExtendedGCDInfo<T> output = new ExtendedGCDInfo<T>();
            output.ACoefficient = a.GetAdditiveIdentity();
            output.BCoefficient = a.GetMultiplicativeIdentity();                       
            output.BOverGCD = a.GetMultiplicativeIdentity();
            output.AOverGCD = a.GetAdditiveIdentity();
            while (!a.Equals(a.GetAdditiveIdentity())) 
            {
                Division<T> division = b.EuclideanDivideBy(a);
                T m = output.ACoefficient.Minus(output.BOverGCD.Times(division.Quotient));
                T n = output.BCoefficient.Minus(output.AOverGCD.Times(division.Quotient));
                b = a;
                a = division.Remainder;
                output.ACoefficient = output.BOverGCD;
                output.BCoefficient = output.AOverGCD;
                output.BOverGCD = m;
                output.AOverGCD = n;
            }
            output.GCD = b;
            output.BOverGCD = a.GetAdditiveIdentity().Minus(output.BOverGCD);
            return output;
        }
        interface IArithmetic<T>
        {
            T GetAdditiveIdentity();
            T GetMultiplicativeIdentity();
            T Minus(T n);
            T Times(T n);
        }
        abstract class Number : IArithmetic<Number>, IComparable, IEquatable<Number>
        {            
            public Number GetAdditiveIdentity()
            {
                return Zero;
            }
            protected abstract Number Add(Number number);
            public static Number operator +(Number a, Number b)
            {
                return a.Add(b);
            }
            public Number Minus(Number number)
            {
                return this - number;
            }
            public abstract Number Negative();
            public static Number operator -(Number a, Number b)
            {
                return a.Add(b.Negative());
            }
            public Number GetMultiplicativeIdentity()
            {
                return One;
            }
            public abstract Number Times(Number number);
            public static Number operator *(Number a, Number b)
            {
                return a.Times(b);
            }
            public abstract Number Reciprocal();
            public static Number operator /(Number a, Number b)
            {
                return a.Times(b.Reciprocal());
            }
            Polynomial minimalPolynomial = null;
            public Polynomial MinimalPolynomial
            {
                get
                {
                    if (minimalPolynomial == null)
                        minimalPolynomial = CalculateMinimalPolynomial();
                    return minimalPolynomial;
                }
            }
            protected abstract Polynomial CalculateMinimalPolynomial();
            public abstract List<Number> GetConjugates();
            protected void RemoveNonConjugates(List<Number> conjugateCandidates)
            {
                for (int i = 0;
                    conjugateCandidates.Count > MinimalPolynomial.Coefficients.Count - 1;)
                    if (!conjugateCandidates[i].MinimalPolynomial.Equals(MinimalPolynomial))
                        conjugateCandidates.RemoveAt(i);
                    else
                        ++i;
            }
            public virtual Integer GetDenominatorLCM()
            {
                return One;
            }
            public virtual Integer GetGreatestIntegerFactor()
            {
                return One;
            }
            public virtual Number Exponentiate(Number exponent)
            {
                if (exponent is Integer integerExponent)
                {
                    if (integerExponent.Sign < 0)
                        return Exponentiate(exponent.Negative()).Reciprocal();
                    return Solver.Exponentiate(this, integerExponent);
                }
                Integer exponentIntegerFactor = exponent.GetGreatestIntegerFactor();
                if (exponentIntegerFactor != One) 
                    return Exponentiate(exponentIntegerFactor).Exponentiate(
                        exponent / exponentIntegerFactor);
                if (exponent is Fraction exponentFraction)
                {
                    Integer denominatorLCM = GetDenominatorLCM();
                    if (denominatorLCM != One)                    
                        return (this * denominatorLCM).Exponentiate(exponent) *
                            denominatorLCM.Exponentiate(new Fraction(exponentFraction.Denominator - One,
                            exponentFraction.Denominator)) / denominatorLCM;                    
                    Integer numeratorGCD = GetGreatestIntegerFactor();
                    if (numeratorGCD != One)
                    {
                        Number factor = numeratorGCD.Exponentiate(exponent);
                        if (!(factor is Exponentiation))
                            return (this / numeratorGCD).Exponentiate(exponent) * factor;
                    }
                }
                return Exponentiation.Create(this, exponent);
            }
            protected int GetTypeIndex()
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
                return 0;
            }
            public virtual int CompareTo(object obj)
            {
                Number number = obj as Number;
                if (ReferenceEquals(number, null))
                    return 1;
                return GetTypeIndex() - number.GetTypeIndex();
            }
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public bool Equals(Number number)
            {
                return CompareTo(number) == 0;
            }
            public static bool operator ==(Number a, Number b)
            {
                if (ReferenceEquals(a, null)) 
                {
                    if (ReferenceEquals(b, null))
                        return true;
                    return false;
                }
                return a.CompareTo(b) == 0;
            }
            public static bool operator !=(Number a, Number b)
            {
                return !(a == b);
            }
            public abstract override int GetHashCode();

            //Treats str as the string representation of a numerical constant, and returns the
            //string representation of the reference object multiplied by that constant.
            public abstract string InsertString(string str);
            public abstract override string ToString();
        }
        abstract class Term : Number
        { }
        abstract class Factor : Term
        { }
        abstract class Rational : Factor
        {            
            public override List<Number> GetConjugates()
            {
                return new List<Number> { this };
            }
        }
        interface IDivisible<T> : IArithmetic<T> 
        {            
            Division<T> EuclideanDivideBy(T divisor);
        }
        static Integer Zero = new Integer(0);
        static Integer One = new Integer(1);
        class Integer : Rational, IDivisible<Integer>
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
            public static Integer Parse(string decimalString)
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
            new public Integer GetAdditiveIdentity()
            {
                return Zero;
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
            public Integer Add(Integer a)
            {
                return this + a;
            }
            protected override Number Add(Number number)
            {
                if (number is Integer integer)
                    return integer + this;
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
            public Integer Minus(Integer a)
            {
                return this - a;
            }
            public override Number Negative()
            {
                return new Integer(Values, (sbyte)-Sign);
            }
            public Integer Magnitude()
            {
                if (Sign < 0)
                    return new Integer(Values, 1);
                else
                    return this;
            }
            static uint[] ShiftValuesLeft(uint[] values, int valuePlaces, int digitPlaces)
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
            public Integer ShiftLeft(int valuePlaces, int digitPlaces)
            {
                return new Integer(ShiftValuesLeft(Values, valuePlaces, digitPlaces), Sign);
            }
            new public Integer GetMultiplicativeIdentity()
            {
                return One;
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
                            product = product + b.ShiftLeft(i, j);
                        power = power << 1;
                    }
                }
                product.Sign = (sbyte)(b.Sign * a.Sign);
                return product;
            }
            public Integer Times(Integer a)
            {
                return this * a;
            }
            public override Number Times(Number number)
            {
                if (number is Integer integer)
                    return integer * this;
                return number * this;
            }
            public override Number Reciprocal()
            {
                if (Sign == 0)
                    throw new DivideByZeroException();
                if (Magnitude() == One)
                    return this;
                return new Fraction(One, this);
            }
            public Integer Reciprocal(Integer characteristic)
            {
                return (Integer)Exponentiate(characteristic - new Integer(2));
            }
            public Division<Integer> EuclideanDivideBy(Integer divisor)
            {
                if (divisor.Sign == 0)
                    throw new DivideByZeroException();
                Division<Integer> division = new Division<Integer>();
                if (divisor.Values.Length > Values.Length ||
                    (divisor.Values.Length == Values.Length &&
                    divisor.Values[divisor.Values.Length - 1] > Values[Values.Length - 1]))
                {
                    division.Quotient = Zero;
                    division.Remainder = this;
                    return division;
                }
                division.Remainder = new Integer(Values, 1);
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
                        Integer shiftedDivisor = positiveDivisor.ShiftLeft(valuePlace -
                            positiveDivisor.Values.Length + 1, i - divisorLeadingDigitPlace);
                        Integer difference = division.Remainder - shiftedDivisor;
                        if (difference.Sign >= 0)
                        {
                            quotient[valuePlace] |= power;
                            division.Remainder = difference;
                        }
                        power = power >> 1;
                    }
                }
                for (int i = Values.Length - 1; i >= positiveDivisor.Values.Length; --i)
                    calculateValue(i, 1);
                calculateValue(positiveDivisor.Values.Length - 1, divisorLeadingDigitPlace);
                division.Quotient = new Integer(quotient, (sbyte)(Sign * divisor.Sign)).ShiftLeft(
                    1 - positiveDivisor.Values.Length, 1 - divisorLeadingDigitPlace);
                if (division.Remainder.Sign != 0)
                    division.Remainder.Sign = Sign;
                return division;
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                return new Polynomial(new List<Integer> { -this, One });
            }
            public override Number Exponentiate(Number exponent)
            {
                if (Sign == 0)
                    return Zero;
                if (exponent == Zero) 
                    return One;
                Integer exponentIntegerFactor = exponent.GetGreatestIntegerFactor();
                if (exponentIntegerFactor != One) 
                    return base.Exponentiate(exponentIntegerFactor).Exponentiate(
                        exponent / exponentIntegerFactor);
                if (exponent is Fraction exponentFraction)
                {
                    Integer radicand = this;
                    Dictionary<Integer, Integer> radicandFactors =
                        new Dictionary<Integer, Integer>();
                    if (radicand.Sign < 0)
                    {
                        radicand = -radicand;
                        radicandFactors.Add(new Integer(-1), One);
                    }
                    Integer factor;
                    for (Integer i = new Integer(2);
                        i <= radicand.EuclideanDivideBy(new Integer(2)).Quotient;)
                    {
                        Division<Integer> division = radicand.EuclideanDivideBy(i);
                        if (division.Remainder.Sign == 0)
                        {
                            factor = i;
                            if (radicandFactors.ContainsKey(factor))
                                ++radicandFactors[factor];
                            else
                                radicandFactors.Add(factor, One);
                            radicand = division.Quotient;
                        }
                        else
                            ++i;
                    }
                    factor = radicand;
                    if (radicandFactors.ContainsKey(factor))
                        ++radicandFactors[factor];
                    else
                        radicandFactors.Add(factor, One);
                    List<Integer> exponentDivisors = exponentFraction.Denominator.GetDivisors();
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
                                Division<Integer> division =
                                    radicandFactors[n].EuclideanDivideBy(exponentDivisors[i]);
                                for (Integer j = Zero; j < division.Quotient; ++j)
                                    termComponents[index] = termComponents[index] * n;
                                for (Integer j = Zero; j < division.Remainder; ++j)
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
                                termComponents[exponentFraction.Denominator].Negative();
                        Integer two = new Integer(2);
                        if (exponentFraction.Denominator.EuclideanDivideBy(two).Remainder.Sign == 0)
                            coefficient *= ComplexExponential.Create(new Fraction(One,
                                two * exponentFraction.Denominator));
                        else
                            coefficient = coefficient.Negative();
                    }
                    List<Factor> factors = new List<Factor>();
                    foreach (Integer index in termComponents.Keys)
                        factors.Add(new Surd(termComponents[index], index));
                    Number output = Product.Create(One, factors) * coefficient;
                    return output;
                }
                return base.Exponentiate(exponent);
            }
            public override Integer GetGreatestIntegerFactor()
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
            public override String InsertString(String str)
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
                    Division<Integer> division = quotient.EuclideanDivideBy(power);
                    quotient = division.Quotient;
                    output.Insert(0, division.Remainder.Values[0]);
                }
                if (Sign < 0)
                    output.Insert(0, '-');
                return output.ToString();
            }
            public static Integer GetGCD(Integer a, Integer b)
            {
                Integer c;
                while (b.Sign != 0)
                {
                    c = b;
                    b = a.EuclideanDivideBy(b).Remainder;
                    a = c;
                }
                return a.Magnitude();
            }
            public static Integer GetGCD(List<Integer> list)
            {
                Integer GCD = list[0];
                for (int i = 1; i < list.Count; ++i)
                    GCD = GetGCD(GCD, list[i]);
                return GCD;
            }
            public static Integer GetLCM(Integer a, Integer b)
            {
                return (Integer)(a / GetGCD(a, b) * b);
            }
            public List<Integer> GetDivisors()
            {
                Integer x = this;
                List<Integer> divisors = new List<Integer>();
                Integer divisor = One;
                Integer half = x.EuclideanDivideBy(new Integer(2)).Quotient;
                while (divisor <= half)
                {
                    if (x.EuclideanDivideBy(divisor).Remainder.Sign == 0)
                        divisors.Add(divisor);
                    ++divisor;
                }
                divisors.Add(x);
                return divisors;
            }
            public Integer ThisChooseK(Integer k)
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
            public static Rational Create(Integer numerator, Integer denominator)
            {
                if (denominator.Sign == 0)
                    throw new DivideByZeroException();
                if (numerator.Sign == 0)
                    return numerator;
                ExtendedGCDInfo<Integer> extendedGCD = ExtendedGCD(numerator, denominator);
                numerator = extendedGCD.AOverGCD;
                denominator = extendedGCD.BOverGCD;
                if (denominator.Sign < 0)
                {
                    numerator = -numerator;
                    denominator = -denominator;
                }
                if (denominator == One)
                    return numerator;
                return new Fraction(numerator, denominator);
            }
            protected override Number Add(Number number)
            {
                if (number is Integer integer)
                    return Create(Numerator + integer * Denominator, Denominator);
                if (number is Fraction fraction)                
                    return Create(Numerator * fraction.Denominator +
                        fraction.Numerator * Denominator, Denominator * fraction.Denominator);                
                return number + this;
            }
            public override Number Negative()
            {
                return new Fraction(-Numerator, Denominator);
            }
            public override Number Times(Number number)
            {
                if (number is Integer integer)
                    return Create(Numerator * integer, Denominator);
                if (number is Fraction fraction)                
                    return Create(Numerator * fraction.Numerator,
                        Denominator * fraction.Denominator);                
                return number * this;
            }
            public override Number Reciprocal()
            {
                return Create(Denominator, Numerator);
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                return new Polynomial(new List<Integer> { -Numerator, Denominator });
            }
            public override Integer GetDenominatorLCM()
            {
                return Denominator;
            }
            public override Integer GetGreatestIntegerFactor()
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
            public override String InsertString(String str)
            {
                return Numerator.InsertString(str) + "/" + Denominator;
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
            public static Factor Create(Number expBase, Number exponent)
            {
                if (exponent is Fraction exponentFraction)                
                    if (exponentFraction.Numerator == One)
                        return new Surd(expBase, exponentFraction.Denominator);                
                return new Transcendental(expBase, exponent);
            }
            protected override Number Add(Number number)
            {
                if (number is Integer integer)
                {
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Fraction fraction)
                    return new Sum(new List<Term> { this, fraction });
                if (number is Exponentiation exponentiation)
                {
                    if (CompareTo(exponentiation) == 0)
                        return Product.Create(new Integer(2), new List<Factor> { this });
                    if (CompareTo(exponentiation.Negative()) == 0)
                        return Zero;
                    return new Sum(new List<Term> { this, exponentiation });
                }
                return number + this;
            }
            public override Number Negative()
            {
                return Product.Create(new Integer(-1), new List<Factor> { this });
            }
            public override Number Times(Number number)
            {
                if (number is Rational rational)
                    return Product.Create(rational, new List<Factor> { this });
                if (number is Exponentiation exponentiation)
                {
                    if (Base == exponentiation.Base) 
                        return Base.Exponentiate(Exponent + exponentiation.Exponent);
                    if (Exponent == exponentiation.Exponent)
                    {
                        Number outputBase = Base * exponentiation.Base;
                        if (outputBase is Exponentiation baseExponentiation)
                            return Create(baseExponentiation.Base, Exponent * Exponent);
                        return (outputBase).Exponentiate(Exponent);
                    }
                    return Product.Create(One, new List<Factor> { this, exponentiation });
                }
                return number * this;
            }
            public override Number Exponentiate(Number exponent)
            {
                return Base.Exponentiate(Exponent * exponent);
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
            public override String InsertString(String str)
            {
                return ToString() + '*' + str;
            }
            public override String ToString()
            {
                String encloseNumber(Number number)
                {
                    string numberString = number.ToString();
                    if ((number is Integer integer && integer.Sign >= 0) || numberString == "i")
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
            public override Number Times(Number number)
            {
                if (number is Surd surd)
                {
                    Integer indexLCM = Integer.GetLCM(Index, surd.Index);
                    return (Base.Exponentiate(indexLCM / Index) * surd.Base.Exponentiate(indexLCM /
                        surd.Index)).Exponentiate(Fraction.Create(One, indexLCM));
                }
                return base.Times(number);
            }
            public override Number Reciprocal()
            {
                return Base.Exponentiate(Fraction.Create(Index - One, Index)) / Base;
            }
            protected override Polynomial CalculateMinimalPolynomial()
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
            public override List<Number> GetConjugates()
            {
                List<Number> baseConjugates = Base.GetConjugates();
                if (baseConjugates == null)
                    return null;
                List<Number> rootsOfUnity = new List<Number>();
                for (Integer i = Zero; i < Index; ++i)
                    rootsOfUnity.Add(ComplexExponential.Create(Fraction.Create(i, Index)));
                List<Number> conjugateCandidates = new List<Number>();
                foreach (Number conjugate in baseConjugates)
                    foreach (Number root in rootsOfUnity)
                        conjugateCandidates.Add(new Surd(conjugate, Index) * root);
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
        }
        class Transcendental : Exponentiation
        {
            public Transcendental(Number expBase, Number exponent) : base(expBase, exponent)
            { }
            public override Number Reciprocal()
            {
                return Create(Base, Exponent.Negative());
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                return null;
            }
            public override List<Number> GetConjugates()
            {
                return null;
            }
        }
        class ComplexExponential : Factor
        {//Represents e^(tau*i*Exponent).
            public Number Exponent { get; }
            public ComplexExponential(Number exponent)
            {
                if (exponent is Fraction fraction)
                {
                    Exponent = fraction -
                        fraction.Numerator.EuclideanDivideBy(fraction.Denominator).Quotient;
                    if (fraction.Numerator.Sign < 0)
                        Exponent += One;
                }
                else
                    Exponent = exponent;
            }
            public static Number Create(Number exponent)
            {
                if (exponent is Integer)
                    return One;
                if (exponent is Fraction exponentFraction)
                {
                    Integer denominator = exponentFraction.Denominator;
                    Integer numerator =
                        exponentFraction.Numerator.EuclideanDivideBy(denominator).Remainder;
                    Integer two = new Integer(2);
                    Division<Integer> division = denominator.EuclideanDivideBy(two);
                    int multiplicityOfTwo = 0;
                    while (division.Remainder.Sign == 0)
                    {
                        denominator = division.Quotient;
                        division = denominator.EuclideanDivideBy(two);
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
                        cosine = (cosine * half + half).Exponentiate(half);
                        if (four < three * denominator && denominator < four)
                            cosine = cosine.Negative();
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
                    Number output = ComplexNumber.Create(angleMultipleCosine,
                        (One - angleMultipleCosine * angleMultipleCosine).Exponentiate(half));
                    Integer t = four * numerator;
                    if (denominator < t && t < two * denominator ||
                        three * denominator < t && t < four * denominator)
                        return output.Negative();
                    return output;
                }
                return new ComplexExponential(exponent);
            }
            protected override Number Add(Number number)
            {
                if (number is Integer integer)
                {
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Fraction fraction)
                    return new Sum(new List<Term> { this, fraction });
                if (number is Exponentiation exponentiation)
                    return new Sum(new List<Term> { this, exponentiation });
                if (number is ComplexExponential complexExponential)
                    if (Exponent.CompareTo(complexExponential.Exponent) == 0)
                        return Product.Create(new Integer(2), new List<Factor> { this });                
                return number + this;
            }
            public override Number Negative()
            {
                return Product.Create(new Integer(-1), new List<Factor> { this });
            }
            public override Number Times(Number number)
            {
                if (number is Integer integer)
                    return Product.Create(integer, new List<Factor> { this });
                if (number is Fraction fraction)
                    return Product.Create(fraction, new List<Factor> { this, });
                if (number is Exponentiation exponentiation)
                    return Product.Create(One, new List<Factor> { this, exponentiation });
                if (number is ComplexExponential complexExponential)
                    return Create(Exponent + complexExponential.Exponent);
                return number * this;
            }
            public override Number Reciprocal()
            {
                return Create(Exponent.Negative());
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                if (Exponent is Fraction fraction)
                {
                    List<Integer> coefficients = new List<Integer> { new Integer(-1), One };
                    List<Polynomial> dividends = new List<Polynomial>();
                    for (Integer i = One; i < fraction.Denominator; ++i)
                    {
                        dividends.Add(new Polynomial(new List<Integer>(coefficients)));
                        coefficients.Insert(1, Zero);
                    }
                    Polynomial annullingPolynomial = dividends[dividends.Count - 1];
                    dividends.RemoveAt(dividends.Count - 1);
                    List<Polynomial> factors = annullingPolynomial.GetFactors();
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
                }
                return null;
            }
            public override List<Number> GetConjugates()
            {
                if (Exponent is Fraction fraction)
                {
                    List<Number> conjugates = new List<Number>();
                    for (Integer i = Zero; i < fraction.Denominator; ++i)
                        conjugates.Add(Create(Fraction.Create(i, fraction.Denominator)));
                    return conjugates;
                }
                return null;
            }
            public override Number Exponentiate(Number exponent)
            {
                return Create(exponent * exponent);
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
            public override string InsertString(string str)
            {
                return ToString() + '*' + str;
            }
            public override string ToString()
            {
                Number exponentWith_i = Exponent * ComplexNumber.Create(Zero, One);
                return "e^(" + exponentWith_i.InsertString("tau") + ')';
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
            public static Term Create(Rational coefficient, List<Factor> factors)
            {
                if (coefficient is Integer coefficientInteger)
                {
                    if (coefficientInteger.Sign == 0)
                        return Zero;
                    if (factors.Count == 1 && coefficientInteger == One)
                        return factors[0];
                }
                if (factors.Count == 0)
                    return coefficient;
                return new Product(coefficient, factors);
            }
            protected override Number Add(Number number)
            {
                if (number is Integer integer)
                {
                    if (integer.Sign == 0)
                        return this;
                    return new Sum(new List<Term> { this, integer });
                }
                if (number is Factor factor)
                {
                    if (Factors.Count == 1 && factor == Factors[0])
                        return Create((Rational)(Coefficient + One), Factors);
                    return new Sum(new List<Term> { this, factor });
                }
                if (number is Product product)
                {
                    if (Factors.Count == product.Factors.Count)
                    {
                        for (int i = 0; i < Factors.Count; ++i)
                            if (Factors[i].CompareTo(product.Factors[i]) != 0)
                                return new Sum(new List<Term> { this, product });
                        return Create((Rational)(Coefficient + product.Coefficient), Factors);
                    }
                    return new Sum(new List<Term> { this, product });
                }
                return number + this;
            }
            public override Number Negative()
            {
                return Create((Rational)Coefficient.Negative(), Factors);
            }
            public override Number Times(Number number)
            {
                if (number is Rational)
                    return Create((Rational)(Coefficient * number), Factors);
                if (number is Factor factor)
                {
                    List<Factor> factors = new List<Factor>(Factors);
                    for (int i = 0; i < factors.Count; ++i)
                    {
                        Number p = number * factors[i];
                        if (p != new Product(Coefficient,
                            new List<Factor> { factor, factors[i] }))
                        {
                            factors.RemoveAt(i);
                            return Create(Coefficient, factors) * p;
                        }
                    }
                    factors.Add(factor);
                    return Create(Coefficient, factors);
                }
                if (number is Product product)
                {
                    Number output =
                        new Product((Rational)(Coefficient * product.Coefficient), Factors);
                    for (int i = 0; i < product.Factors.Count; ++i)
                        output = output * product.Factors[i];
                    return output;
                }
                return number * this;
            }
            public override Number Reciprocal()
            {
                List<Factor> factors = new List<Factor>();
                foreach (Factor factor in Factors)
                    factors.Add((Factor)factor.Reciprocal());
                return Create((Rational)Coefficient.Reciprocal(), factors);
            }
            protected override Polynomial CalculateMinimalPolynomial()
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
                variableForm.SetCoefficient(indices, Coefficient);
                return variableForm.GetMinimalPolynomial();
            }
            public override List<Number> GetConjugates()
            {
                List<List<Number>> factorConjugates = new List<List<Number>>();
                foreach (Factor factor in Factors)
                {
                    List<Number> conjugates = factor.GetConjugates();
                    if (conjugates == null)
                        return null;
                    factorConjugates.Add(conjugates);
                }
                List<List<Number>> conjugateCombinations =
                    GenerateCartesianProduct(factorConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = One;
                    foreach (Number number in combination)
                        conjugate *= number;
                    conjugateCandidates.Add(conjugate);
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            public override Integer GetDenominatorLCM()
            {
                return Coefficient.GetDenominatorLCM();
            }
            public override Integer GetGreatestIntegerFactor()
            {
                return Coefficient.GetGreatestIntegerFactor();
            }
            public override Number Exponentiate(Number exponent)
            {
                Number output = Coefficient.Exponentiate(exponent);
                foreach (Factor factor in Factors)
                    output *= factor.Exponentiate(exponent);
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
            StringBuilder FormatFactors()
            {
                StringBuilder output = new StringBuilder(Factors[0].ToString());
                for (int i = 1; i < Factors.Count; ++i)
                    output.Append('*' + Factors[i].ToString());
                return output;
            }
            string AddCoefficient(StringBuilder factors)
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
            public override string InsertString(string str)
            {
                StringBuilder output = FormatFactors().Append('*' + str);
                return AddCoefficient(output);
            }
            public override string ToString()
            {
                StringBuilder output = FormatFactors();
                return AddCoefficient(output);
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
            protected override Number Add(Number number)
            {
                if (number is Term term)
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
                    terms.Add(term);
                    return new Sum(terms);
                }
                else if (number is Sum sum)
                {
                    Number output = this;
                    foreach (Term t in sum.Terms)
                        output += t;
                    return output;
                }
                return number + this;
            }
            public override Number Negative()
            {
                List<Term> terms = new List<Term>();
                foreach (Term term in Terms)
                    terms.Add((Term)term.Negative());
                return new Sum(terms);
            }
            public override Number Times(Number number)
            {
                List<Term> terms = new List<Term>();
                if (number is Term)
                {
                    Number output = Zero;
                    foreach (Term term in Terms)
                        output += term * number;
                    return output;
                }
                if (number is Sum multiplier)
                {
                    Number output = Zero;
                    foreach (Term multiplicandTerm in Terms)
                        foreach (Term multiplierTerm in multiplier.Terms)
                            output += multiplicandTerm * multiplierTerm;
                    return output;
                }
                return number * this;
            }
            public override Number Reciprocal()
            {
                List<Number> conjugates = GetConjugates();
                if (conjugates == null)
                    return new Transcendental(this, new Integer(-1));
                conjugates.Remove(this);
                Number numerator = One;
                foreach (Number conjugate in conjugates)
                    numerator *= conjugate;
                return numerator / (numerator * this);
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                Polynomial[] minimalPolynomials = new Polynomial[Terms.Count];
                Rational constant;
                int startIndex;
                if (Terms[0] is Rational rational)
                {
                    constant = rational;
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
                if (constant != Zero) 
                    variableForm.SetCoefficient(new int[minimalPolynomials.Length], constant);
                for (int i = 0; i < minimalPolynomials.Length; ++i)
                {
                    int[] indices = new int[minimalPolynomials.Length];
                    indices[i] = 1;
                    variableForm.SetCoefficient(indices, One);
                }
                return variableForm.GetMinimalPolynomial();
            }
            public override List<Number> GetConjugates()
            {
                List<List<Number>> termConjugates = new List<List<Number>>();
                foreach (Term term in Terms)
                {
                    List<Number> conjugates = term.GetConjugates();
                    if (conjugates == null)
                        return null;
                    termConjugates.Add(conjugates);
                }
                List<List<Number>> conjugateCombinations = GenerateCartesianProduct(termConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = Zero;
                    foreach (Number number in combination)
                        conjugate += number;
                    conjugateCandidates.Add(conjugate);
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            public override Integer GetDenominatorLCM()
            {
                Integer LCM = One;
                foreach (Term term in Terms)
                    LCM = Integer.GetLCM(LCM, term.GetDenominatorLCM());
                return LCM;
            }
            public override Integer GetGreatestIntegerFactor()
            {
                Integer GCD = Terms[0].GetGreatestIntegerFactor();
                for (int i = 1; i < Terms.Count; ++i)
                    GCD = Integer.GetGCD(GCD, Terms[i].GetGreatestIntegerFactor());
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
            public override string InsertString(string str)
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
            public static Number Create(Number real, Number imaginary)
            {
                if (imaginary is Integer imaginaryInteger)                
                    if (imaginaryInteger.Sign == 0)
                        return real;                
                return new ComplexNumber(real, imaginary);
            }
            protected override Number Add(Number number)
            {
                if (number is ComplexNumber complexNumber)
                    return Create(Real + complexNumber.Real, Imaginary + complexNumber.Imaginary);                
                return Create(Real + number, Imaginary);
            }
            public override Number Negative()
            {
                return Create(Real.Negative(), Imaginary.Negative());
            }
            public override Number Times(Number number)
            {
                if (number is ComplexNumber complexNumber)
                    return Create(Real * complexNumber.Real - Imaginary * complexNumber.Imaginary,
                        Real * complexNumber.Imaginary + Imaginary * complexNumber.Real);                
                return Create(Real * number, Imaginary * number);
            }
            public override Number Reciprocal()
            {
                Number denominator = Real * Real + Imaginary * Imaginary;
                return Create(Real / denominator, Imaginary.Negative() / denominator);
            }
            protected override Polynomial CalculateMinimalPolynomial()
            {
                MultivariatePolynomial variableForm = new MultivariatePolynomial(
                    new Polynomial[] { Real.MinimalPolynomial, Imaginary.MinimalPolynomial,
                        new Polynomial(new List<Integer> { One, Zero, One })});
                variableForm.SetCoefficient(new int[] { 1, 0, 0 }, One);
                variableForm.SetCoefficient(new int[] { 0, 1, 1 }, One);
                return variableForm.GetMinimalPolynomial();
            }
            public override List<Number> GetConjugates()
            {
                List<Number> realPartConjugates = Real.GetConjugates();
                List<Number> imaginaryPartConjugates = Imaginary.GetConjugates();
                List<Number> conjugateCandidates = new List<Number>();
                foreach (Number a in realPartConjugates)
                    foreach (Number b in imaginaryPartConjugates)
                        conjugateCandidates.Add(a + Create(Zero, One) * b);
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            public override Integer GetDenominatorLCM()
            {
                return Integer.GetLCM(Real.GetDenominatorLCM(), Imaginary.GetDenominatorLCM());
            }
            public override Integer GetGreatestIntegerFactor()
            {
                return Integer.GetGCD(Real.GetGreatestIntegerFactor(),
                    Imaginary.GetGreatestIntegerFactor());
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
            public override string InsertString(string str)
            {
                if (Real is Integer && ((Integer)Real).Sign == 0)
                    return Imaginary.InsertString("i*" + str);
                return '(' + ToString() + ')' + str;
            }
            public override string ToString()
            {
                StringBuilder output =
                    new StringBuilder(Imaginary.InsertString("i").ToString());
                if (!(Real is Integer integer && integer.Sign == 0))
                {
                    if (output[0] != '-')
                        output.Insert(0, '+');
                    output.Insert(0, Real.ToString());
                }
                return output.ToString();
            }
        }
        static List<Integer> Primes = new List<Integer> { new Integer(2), new Integer(3) };
        class Polynomial : IComparable, IDivisible<Polynomial>
        {
            //The index of the coefficient represents the degree of its term.
            public List<Integer> Coefficients { get; }
            Integer Characteristic = Zero;
            public Polynomial(List<Integer> coefficients)
            {
                while (coefficients.Count != 0 && coefficients[coefficients.Count - 1] == Zero) 
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
            }
            public Polynomial(List<Integer> coefficients, Integer characteristic)
            {                
                Characteristic = characteristic;
                if (Characteristic != Zero) 
                {
                    Coefficients = new List<Integer>();
                    foreach (Integer coefficient in coefficients)
                    {
                        Integer remainder = coefficient.EuclideanDivideBy(characteristic).Remainder;
                        if (remainder.Sign < 0)
                            Coefficients.Add(remainder + characteristic);
                        else
                            Coefficients.Add(remainder);
                    }
                }
                else
                    Coefficients = coefficients;
                while (Coefficients.Count != 0 && Coefficients[Coefficients.Count - 1] == Zero) 
                    Coefficients.RemoveAt(Coefficients.Count - 1);
            }
            public Polynomial GetAdditiveIdentity()
            {
                return new Polynomial(new List<Integer>(), Characteristic);
            }
            public static Polynomial operator +(Polynomial a, Polynomial b)
            {
                if (a.Coefficients.Count < b.Coefficients.Count)
                    return b + a;
                List<Integer> output = new List<Integer>(a.Coefficients);
                for (int i = 0; i < b.Coefficients.Count; ++i)
                    output[i] = output[i] + b.Coefficients[i];
                return new Polynomial(output, a.Characteristic);
            }
            public Polynomial Minus(Polynomial a)
            {
                return this - a;
            }
            public Polynomial Negative()
            {
                List<Integer> output = new List<Integer>();
                foreach (Integer coefficient in Coefficients)
                    output.Add(-coefficient);
                return new Polynomial(output, Characteristic);
            }
            public static Polynomial operator -(Polynomial a, Polynomial b)
            {
                return a + b.Negative();
            }
            public Polynomial GetMultiplicativeIdentity()
            {
                return new Polynomial(new List<Integer> { One }, Characteristic);
            }
            public Polynomial Times(Polynomial a)
            {
                return this * a;
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
            delegate Number Divider(Integer a, Integer b);
            public Division<Polynomial> EuclideanDivideBy(Polynomial divisor)
            {
                Divider divide;
                if (Characteristic == Zero)
                    divide = delegate (Integer a, Integer b) { return a / b; };
                else
                    divide = delegate (Integer a, Integer b) {
                        return a * b.Exponentiate(Characteristic - new Integer(2)); };
                Division<Polynomial> division;
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
                return a.EuclideanDivideBy(b).Quotient;
            }
            public static Polynomial operator %(Polynomial a, Polynomial b)
            {
                return a.EuclideanDivideBy(b).Remainder;
            }
            public Polynomial GetDerivative()
            {
                List<Integer> derivativeCoefficients = new List<Integer>();
                for (int i = 1; i < Coefficients.Count; ++i)
                    derivativeCoefficients.Add(Coefficients[i] * new Integer(i));
                return new Polynomial(derivativeCoefficients, Characteristic);
            }
            Integer ReduceToPrimitivePart()
            {//Returns the content used in the reduction.
                Integer content = Integer.GetGCD(Coefficients);
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
            static Polynomial GetGCD(Polynomial a, Polynomial b)
            {
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    Polynomial t = a;
                    a = new Polynomial(new List<Integer>(b.Coefficients), b.Characteristic);
                    b = new Polynomial(new List<Integer>(t.Coefficients), t.Characteristic);
                }
                else
                {
                    a = new Polynomial(new List<Integer>(a.Coefficients), a.Characteristic);
                    b = new Polynomial(new List<Integer>(b.Coefficients), b.Characteristic);
                }
                if (b.Coefficients.Count == 0)
                    return a;
                Integer d = Integer.GetGCD(a.ReduceToPrimitivePart(), b.ReduceToPrimitivePart());
                Integer g = One;
                Number h = One;
                Integer degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                Polynomial remainder = ((Integer)b.Coefficients[b.Coefficients.Count - 1].
                    Exponentiate(degree + One) * a).EuclideanDivideBy(b).Remainder;
                while (remainder.Coefficients.Count > 1)
                {
                    a = b;
                    Number divisor = (g * h).Exponentiate(degree);
                    b = remainder;
                    for (int i = 0; i < b.Coefficients.Count; ++i)
                        b.Coefficients[i] = (Integer)(b.Coefficients[i] / divisor);
                    g = a.Coefficients[a.Coefficients.Count - 1];
                    h = h.Exponentiate(One - degree) * g.Exponentiate(degree);
                    degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    remainder = ((Integer)b.Coefficients[b.Coefficients.Count - 1].Exponentiate(
                        degree + One) * a).EuclideanDivideBy(b).Remainder;
                }
                if (remainder.Coefficients.Count == 1)
                    return new Polynomial(new List<Integer> { d }, a.Characteristic);
                b.ReduceToPrimitivePart();
                return d * b;
            }
            public List<Polynomial> GetFactors()
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
                        W = Exponentiate(W, Characteristic).EuclideanDivideBy(this).Remainder;
                        Polynomial degreeDFactorProduct =
                            GetGCD(W - new Polynomial(new List<Integer> { Zero, One }), V);
                        if (degreeDFactorProduct.Coefficients.Count > 1)
                        {
                            distinctDegreeFactors.Add(d, degreeDFactorProduct);
                            V = V.EuclideanDivideBy(degreeDFactorProduct).Quotient;
                            W = W.EuclideanDivideBy(V).Remainder;
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
                                C = T + (C * C).EuclideanDivideBy(factorProduct).Remainder;
                            B = GetGCD(factorProduct, C);
                            while (B.Coefficients.Count < 2 ||
                                B.Coefficients.Count == factorProduct.Coefficients.Count)
                            {
                                T *= new Polynomial(new List<Integer> { Zero, Zero, One },
                                    Characteristic);
                                C = new Polynomial(new List<Integer>(T.Coefficients),
                                    Characteristic);
                                for (int j = 1; j < degree; ++j)
                                    C = T + (C * C).EuclideanDivideBy(factorProduct).Remainder;
                                B = GetGCD(factorProduct, C);
                            }
                        }
                        else
                        {
                            List<Integer> coefficients = new List<Integer>();
                            Integer p = Characteristic - One;
                            for (int j = 1; j < 2 * degree; ++j)
                                coefficients.Add(p);
                            coefficients.Add(One);
                            Integer power = (Integer)((Characteristic.Exponentiate(
                                new Integer(degree)) - One) / new Integer(2));
                            B = GetGCD(factorProduct, Exponentiate(new Polynomial(
                                new List<Integer>(coefficients), Characteristic), power) -
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
                                B = GetGCD(factorProduct, Exponentiate(new Polynomial(
                                    new List<Integer>(coefficients), Characteristic), power) -
                                    new Polynomial(new List<Integer> { One }, Characteristic));
                            }
                        }
                        CZSplit(B, degree);
                        CZSplit(factorProduct.EuclideanDivideBy(B).Quotient, degree);
                    }
                    foreach (int degree in distinctDegreeFactors.Keys)
                        CZSplit(distinctDegreeFactors[degree], degree);
                    return irreducibleFactors;
                }
                List<Polynomial> squarefreeFactors = new List<Polynomial>();
                Polynomial derivative = GetDerivative();
                Polynomial a = GetGCD(this, derivative);
                Polynomial b = this / a;
                Polynomial c = derivative / a - b.GetDerivative();
                while (!(b.Coefficients.Count == 1 &&
                    (b.Coefficients[0] == One || b.Coefficients[0] == new Integer(-1))))
                {
                    a = GetGCD(b, c);
                    if (!squarefreeFactors.Contains(a))
                        squarefreeFactors.Add(a);
                    b = b / a;
                    c = c / a - b.GetDerivative();
                }
                List<Polynomial> splitFactor(Polynomial factor)
                {
                    Polynomial moddedFactor;
                    int i = 0;
                    while (true)
                    {
                        if (i == Primes.Count)
                        {
                            Integer primeCandidate = Primes[Primes.Count - 1] + new Integer(2);
                            while (true)
                            {
                                bool isDivisible = false;
                                for (int j = 0; Primes[j] <=
                                    primeCandidate.EuclideanDivideBy(new Integer(2)).Quotient; ++j)
                                    if (primeCandidate.EuclideanDivideBy(
                                        Primes[j]).Remainder == Zero)
                                    {
                                        isDivisible = true;
                                        break;
                                    }
                                if (!isDivisible)
                                    break;
                                primeCandidate += new Integer(2);
                            }
                            Primes.Add(primeCandidate);
                        }
                        moddedFactor = new Polynomial(factor.Coefficients, Primes[i]);
                        Polynomial GCD = GetGCD(moddedFactor, moddedFactor.GetDerivative());
                        if (GCD.Coefficients.Count == 1 && GCD.Coefficients[0] == One)
                        {
                            moddedFactor *=
                                moddedFactor.Coefficients[moddedFactor.Coefficients.Count - 1].
                                Reciprocal(moddedFactor.Characteristic);
                            break;
                        }
                        ++i;
                    }
                    List<Polynomial> irreducibleModdedFactors = moddedFactor.GetFactors();
                    Integer bound = Zero;
                    foreach (Integer coefficient in factor.Coefficients)
                        bound += coefficient * coefficient;
                    Integer squareRoot = One;
                    while (squareRoot * squareRoot < bound)
                        ++squareRoot;
                    Integer factorDegree = new Integer(factor.Coefficients.Count - 1);
                    Integer k = factorDegree.EuclideanDivideBy(new Integer(4)).Quotient;
                    bound = (squareRoot * factorDegree.ThisChooseK(k) +
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude() *
                        factorDegree.ThisChooseK(k - One)) * new Integer(2) *
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude();
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
                        while (characteristicExponent < e)
                        {                            
                            ++characteristicExponent;
                            characteristicToPower *= moddedFactor.Characteristic;
                            A.Characteristic = characteristicToPower;
                            B.Characteristic = characteristicToPower;
                            ExtendedGCDInfo<Polynomial> extendedGCD = ExtendedGCD(A, B);
                            Integer reciprocal = extendedGCD.GCD.Coefficients[0].Reciprocal(
                                moddedFactor.Characteristic);
                            extendedGCD.ACoefficient *= reciprocal;
                            extendedGCD.BCoefficient *= reciprocal;
                            A.Characteristic = Zero;
                            B.Characteristic = Zero;
                            List<Integer> fCoefficients = (factor - A * B).Coefficients;
                            for (int l = 0; l < fCoefficients.Count; ++l)
                                fCoefficients[l] =
                                    (Integer)(fCoefficients[l] / characteristicToPower);
                            Polynomial f =
                                new Polynomial(fCoefficients, moddedFactor.Characteristic);
                            Division<Polynomial> division =
                                (extendedGCD.BCoefficient * f).EuclideanDivideBy(A);
                            division.Quotient.Characteristic = Zero;
                            division.Remainder.Characteristic = Zero;
                            A += characteristicToPower * division.Remainder;
                            B += characteristicToPower *
                                (extendedGCD.ACoefficient * f + B * division.Quotient);                            
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
                                Integer remainder = V.Coefficients[j].EuclideanDivideBy(
                                    characteristicPower).Remainder;
                                if (remainder * two <= characteristicPower)
                                    V.Coefficients[j] = remainder;
                                else
                                    V.Coefficients[j] = remainder - characteristicPower;
                            }
                            Division<Polynomial> division = factor.EuclideanDivideBy(V);
                            if (division.Remainder.Coefficients.Count == 0)
                            {
                                while (division.Remainder.Coefficients.Count == 0) 
                                {
                                    factor = division.Quotient;
                                    division = factor.EuclideanDivideBy(V);
                                }
                                V.ReduceToPrimitivePart();
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
                        finalFactor.ReduceToPrimitivePart();
                        factorSplit.Add(finalFactor);
                    }
                    return factorSplit;
                }
                foreach (Polynomial factor in squarefreeFactors)
                {
                    Polynomial X = new Polynomial(new List<Integer> { Zero, One });
                    if (factor.Coefficients[0] == Zero)
                    {
                        irreducibleFactors.Add(X);
                        factor.Coefficients.RemoveAt(0);
                        while (factor.Coefficients[0] == Zero)
                            factor.Coefficients.RemoveAt(0);
                    }
                    if (factor.Coefficients[0].Magnitude() <
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude())
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
            public void SetCoefficient(int[] indices, Rational coefficient)
            {//The nth index is the degree with respect to the nth variable of the term being set.
                foreach (int[] keyIndices in Coefficients.Keys)
                    if (AreEqualByValue(indices, keyIndices))
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
                        if (AreEqualByValue(indicesA, indicesB))
                        {
                            Rational termSum =
                                (Rational)(sum.Coefficients[indicesA] + b.Coefficients[indicesB]);
                            if (termSum == Zero)
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
            static MultivariatePolynomial ReductionlessMultiply(MultivariatePolynomial a,
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
                            if (AreEqualByValue(productIndices, keyIndices))
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
                                    reducedDegreeFactor.SetCoefficient(termComponentIndices, One);
                                    productComponent = ReductionlessMultiply(productComponent,
                                        reducedDegreeFactor);
                                }
                                for (int j = 0; j < a.MinimalPolynomials[i].Coefficients.Count - 1;
                                    ++j)
                                {
                                    if (a.MinimalPolynomials[i].Coefficients[j] == Zero)
                                        continue;
                                    int[] termComponentIndices = new int[indicesA.Length];
                                    termComponentIndices[i] = j;
                                    termComponent.Coefficients[termComponentIndices] =
                                        (Rational)(a.MinimalPolynomials[i].Coefficients[j] /
                                        a.MinimalPolynomials[i].Coefficients[
                                        a.MinimalPolynomials[i].Coefficients.Count - 1].Negative());
                                }
                            }
                            productComponent =
                                ReductionlessMultiply(productComponent, termComponent);
                        }
                        product += productComponent;
                    }
                return product;
            }
            public static MultivariatePolynomial operator *(Rational a, MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(b.MinimalPolynomials);
                if (a == Zero)
                    return product;
                foreach (int[] indices in b.Coefficients.Keys)
                    product.Coefficients.Add(indices, (Rational)(a * b.Coefficients[indices]));
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a, Rational b)
            {
                return b * a;
            }
            public Polynomial GetMinimalPolynomial()
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
                            if (AreEqualByValue(indices, termsPresent[i]))
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
                    if (matrixRow[0] != Zero)
                        constantIsPresent = true;
                    augmentationRow.Insert(0, Zero);
                    augmentation.Add(new Polynomial(new List<Integer>(augmentationRow)));
                }
                for (int i = 1; i < matrix.Count; ++i)
                    for (int j = i; j < matrix.Count; ++j)
                    {
                        if (matrix[i - 1][matrix.Count - i] == Zero)
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
                                if (matrix[k][matrix.Count - i] != Zero)
                                {
                                    ExtendedGCDInfo<Integer> extendedGCD = ExtendedGCD(
                                        matrix[k][matrix.Count - i].GetDenominatorLCM() *
                                        matrix[i - 1][matrix.Count - i].GetGreatestIntegerFactor(),
                                        matrix[k][matrix.Count - i].GetGreatestIntegerFactor() *
                                        matrix[i - 1][matrix.Count - i].GetDenominatorLCM());
                                    for (int l = 0; l < matrix.Count; ++l)
                                        matrix[k][l] = (Rational)(matrix[k][l] * extendedGCD.
                                            AOverGCD - matrix[i - 1][l] * extendedGCD.BOverGCD);
                                    augmentation[k] = augmentation[k] * extendedGCD.AOverGCD -
                                        augmentation[i - 1] * extendedGCD.BOverGCD;
                                }                            
                            break;
                        }
                    }
                List<Integer> annullingPolynomialCoefficients =
                    augmentation[matrix.Count - 1].Coefficients;
                for (int i = 0; i < annullingPolynomialCoefficients.Count; ++i)
                    annullingPolynomialCoefficients[i] *=
                        matrix[matrix.Count - 1][0].GetDenominatorLCM();
                annullingPolynomialCoefficients[0] -=
                    matrix[matrix.Count - 1][0].GetGreatestIntegerFactor();
                Polynomial minimalPolynomial = new Polynomial(annullingPolynomialCoefficients);
                List<Polynomial> factors = minimalPolynomial.GetFactors();
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
            static bool AreEqualByValue(int[] a, int[] b)
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
        static Number EvaluateExpression(List<char> operations, List<Number> numbers)
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
                    numbers[i] = EvaluateExpression(operations.GetRange(i + 1, matchingParenIndex -
                        i - 1), numbers.GetRange(i + 1, matchingParenIndex - i - 1));
                    operations[i] = ' ';
                    numbers.RemoveRange(i + 1, matchingParenIndex - i);
                    operations.RemoveRange(i + 1, matchingParenIndex - i);
                }
                else
                    ++i;
            }
            try
            {
                for (int i = 0; i < operations.Count;)
                {
                    if (i == 0 || numbers[i - 1] == null)
                    {
                        if (operations[i] == '+')
                            if (numbers[i + 1] != null || operations[i + 1] == '+' ||
                                operations[i + 1] == '-')
                            {
                                numbers.RemoveAt(i);
                                operations.RemoveAt(i);
                            }
                            else
                                throw new InvalidUserInput("Operator missing operand.");
                        else if (operations[i] == '-')
                        {
                            if (numbers[i + 1] != null)
                            {
                                numbers[i + 1] = numbers[i + 1].Negative();
                                numbers.RemoveAt(i);
                                operations.RemoveAt(i);
                            }
                            else if (operations[i + 1] == '+') 
                            {
                                operations[i + 1] = '-';
                                numbers.RemoveAt(i);
                                operations.RemoveAt(i);
                            }
                            else if(operations[i + 1] == '-')
                            {
                                operations[i + 1] = '+';
                                numbers.RemoveAt(i);
                                operations.RemoveAt(i);
                            }
                            else
                                throw new InvalidUserInput("Operator missing operand.");
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
                        numbers[i - 1] = numbers[i - 1].Exponentiate(numbers[i + 1]);
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
                                    numbers.Add(ComplexNumber.Create(Zero,
                                        Integer.Parse(numberCollector.ToString())));
                                    operations.Add(' ');
                                }
                                else if (!"()+-*/^".Contains(c.ToString()))
                                    throw new InvalidUserInput(c + " is an invalid character.");
                                else
                                {
                                    numbers.Add(Integer.Parse(numberCollector.ToString()));
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
                            numbers.Add(ComplexNumber.Create(Zero, One));
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
                        numbers.Add(Integer.Parse(numberCollector.ToString()));
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
                    Number expression = EvaluateExpression(operations, numbers);
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
