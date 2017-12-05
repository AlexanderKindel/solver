using System;
using System.Collections.Generic;
using System.Text;

namespace Calculator
{ 
    abstract class Number : IComparable, IEquatable<Number>
    {
        protected abstract Number add(Number number);
        public abstract Number negative();
        protected abstract Number multiply(Number number);
        public abstract Number reciprocal();
        public virtual int getDenominatorLCM()
        {
            return 1;
        }
        public virtual Number exponentiate(Number exponent)
        {
            if (exponent is Integer)
            {
                int exponentInt = ((Integer)exponent).Value;
                Number output = new Integer(1);
                if (exponentInt < 0)
                {
                    for (int i = 0; i > exponentInt; --i)
                        output = this * output;
                    return output.reciprocal();
                }
                for (int i = 0; i < exponentInt; ++i)
                    output = this * output;
                return output;
            }
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                Integer denominatorLCM = new Integer(getDenominatorLCM());
                if (denominatorLCM.Value != 1)
                    return (this * denominatorLCM.exponentiate(new Integer(
                        exponentFraction.Denominator - 1))).exponentiate(exponentFraction) /
                        denominatorLCM.exponentiate(new Integer(exponentFraction.Numerator));
            }
            if (exponent is NumberList)
            {
                List<Number> output = new List<Number>();
                List<Number> list = ((NumberList)exponent).Numbers;
                foreach (Number number in list)
                    output.Add(exponentiate(number));
                return new NumberList(output);
            }
            return new Exponentiation(this, exponent);
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
        protected int getTypeIndex()
        {
            if (this is Integer)
                return 1;
            if (this is GaussianInteger)
                return 2;
            if (this is Exponentiation)
                return 3;
            if (this is ComplexExponential)
                return 4;
            if (this is Product)
                return 5;
            if (this is Fraction)
                return 6;
            if (this is Sum)
                return 7;
            if (this is ComplexNumber)
                return 8;
            return 9;
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

        //Treats str as the String representation of a numerical constant, and returns the String
        //representation of the reference object multiplied by that constant.
        public abstract String insertString(String str);
        public abstract override String ToString();
    }
    abstract class Term : Number
    { }
    abstract class Factor : Term
    { }
    abstract class Rational : Factor
    { }
    abstract class WholeNumber : Rational
    {
        public abstract int getNormSquared();
        public Dictionary<WholeNumber, int>getFactorization()
        {//The keys represent the factors, and the values, their multiplicities.
            WholeNumber number = this;
            Dictionary<WholeNumber, int> factors = new Dictionary<WholeNumber, int>();
            int sumOfRealAndImaginary = 2;
            while (true)
            {
                int imaginaryPart = sumOfRealAndImaginary / 2;
                int realPart = sumOfRealAndImaginary - imaginaryPart;
                for (int i = 0; imaginaryPart - i >= 0;)
                {
                    int normSquared = number.getNormSquared();
                    if ((realPart + i) * (realPart + i) + (imaginaryPart - i) *
                        (imaginaryPart - i) / 2 > normSquared) 
                    {
                        if (factors.ContainsKey(number))
                            factors[number] += 1;
                        else
                            factors.Add(number, 1);
                        return factors;
                    }
                    bool testFactor(WholeNumber factor)
                    {
                        Number quotient = number / factor;
                        if (quotient is WholeNumber)
                        {
                            if (factors.ContainsKey(factor))
                                factors[factor] += 1;
                            else
                                factors.Add(factor, 1);
                            number = (WholeNumber)quotient;
                            return true;
                        }
                        return false;
                    }
                    if (imaginaryPart - i == 0)
                    {
                        if (testFactor(new Integer(realPart + i)))
                            continue;
                    }
                    else
                    {
                        if (testFactor(GaussianInteger.create(realPart + i, imaginaryPart - i)))
                            continue;
                        if (testFactor(GaussianInteger.create(-realPart - i, imaginaryPart - i)))
                            continue;
                    }
                    ++i;
                }
                ++sumOfRealAndImaginary;
            }
        }
        static protected bool returnAllRoots = true;
        public Number exponentiateWholeNumber(Number exponent)
        {
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                Dictionary<WholeNumber, int> radicandFactors = getFactorization();
                List<int> exponentDivisors = new List<int>();
                int divisor = 1;
                while (divisor <= exponentFraction.Denominator / 2)
                {
                    if (exponentFraction.Denominator % divisor == 0)
                        exponentDivisors.Add(divisor);
                    ++divisor;
                }
                exponentDivisors.Add(exponentFraction.Denominator);
                Dictionary<int, WholeNumber> termComponents = new Dictionary<int, WholeNumber>();
                termComponents.Add(1, new Integer(1));
                termComponents.Add(exponentFraction.Denominator, new Integer(1));
                foreach (WholeNumber factor in radicandFactors.Keys)
                    for (int i = exponentDivisors.Count - 1; i >= 0; --i)
                        if (exponentDivisors[i] <= radicandFactors[factor])
                        {
                            int index = exponentFraction.Denominator / exponentDivisors[i];
                            if (!termComponents.ContainsKey(index))
                                termComponents.Add(index, new Integer(1));
                            for (int j = 0;
                                j < radicandFactors[factor] / exponentDivisors[i]; ++j)
                                termComponents[index] =
                                    (WholeNumber)(termComponents[index] * factor);
                            for (int j = 0;
                                j < radicandFactors[factor] % exponentDivisors[i]; ++j)
                                termComponents[exponentFraction.Denominator] = (WholeNumber)(
                                    termComponents[exponentFraction.Denominator] * factor);
                            break;
                        }
                Dictionary<int, WholeNumber> termComponentsToPower =
                    new Dictionary<int, WholeNumber>();
                foreach (int index in termComponents.Keys)
                {
                    termComponentsToPower.Add(index, new Integer(1));
                    for (int i = 0; i < exponentFraction.Numerator; ++i)
                        termComponentsToPower[index] =
                            (WholeNumber)(termComponentsToPower[index] * termComponents[index]);
                }
                Number coefficient = termComponentsToPower[1];
                termComponentsToPower.Remove(1);
                if (termComponentsToPower[exponentFraction.Denominator] is Integer &&
                    ((Integer)termComponentsToPower[exponentFraction.Denominator]).Value == 1)
                    termComponentsToPower.Remove(exponentFraction.Denominator);
                List<int> indices = new List<int>();
                foreach (int index in termComponentsToPower.Keys)
                    indices.Add(index);
                foreach (int index in indices)
                {
                    if (termComponentsToPower[index] is Integer)
                    {
                        Integer radicandInteger = (Integer)termComponentsToPower[index];
                        if (radicandInteger.Value == 1)
                            termComponentsToPower.Remove(index);
                        else if (radicandInteger.Value < 0)
                        {
                            if (index % 2 == 1)
                                coefficient = coefficient.negative();
                            else if (index % 4 == 0)
                            {
                                if (termComponentsToPower.ContainsKey(2))
                                    termComponentsToPower[2] = (WholeNumber)(
                                        termComponentsToPower[2] * new Integer(2));
                                else
                                    termComponentsToPower.Add(2, new Integer(2));
                                coefficient *= ComplexNumber.create(
                                    Fraction.create(1, 2), Fraction.create(1, 2));
                            }
                            else
                                coefficient *= new GaussianInteger(0, 1);
                            if (radicandInteger.Value == -1)
                                termComponentsToPower.Remove(index);
                        }
                    }
                    else
                    {
                        GaussianInteger gaussianRadicand =
                            (GaussianInteger)termComponentsToPower[index];
                        if (gaussianRadicand.Real == 0)
                        {
                            Integer sign;
                            if (gaussianRadicand.Imaginary < 0)
                                sign = new Integer(-1);
                            else
                                sign = new Integer(1);
                            if (gaussianRadicand.Imaginary == sign.Value)
                                termComponentsToPower.Remove(index);
                            if (index % 4 == 0)
                            {
                                Integer a = new Integer(2);
                                Rational b = Fraction.create(1, 2);
                                Exponentiation c = new Exponentiation(a, b);
                                coefficient *= ComplexNumber.create(new Exponentiation(a + c, b) /
                                    a, new Exponentiation(a - c, b) * sign / a);
                            }
                            else if (index % 4 == 1)
                                coefficient *= new GaussianInteger(0, 1) * sign;
                            else if (index % 4 == 2)
                            {
                                if (termComponentsToPower.ContainsKey(2))
                                    termComponentsToPower[2] = (WholeNumber)(
                                        termComponentsToPower[2] * new Integer(2));
                                else
                                    termComponentsToPower.Add(2, new Integer(2));
                                coefficient *= ComplexNumber.create(Fraction.create(1, 2),
                                    Fraction.create(1, 2) * sign);
                            }
                            else
                                coefficient *= new GaussianInteger(0, -1) * sign;
                        }
                    }
                }
                List<Factor> factors = new List<Factor>();
                int largestIndex = 1;
                foreach (int index in termComponentsToPower.Keys)
                {
                    factors.Add(new Exponentiation(termComponentsToPower[index],
                        Fraction.create(1, index)));
                    if (index > largestIndex)
                        largestIndex = index;
                }
                if (returnAllRoots && largestIndex < exponentFraction.Denominator)
                {
                    returnAllRoots = false;
                    Number rootOfUnity =
                        ComplexExponential.create(Fraction.create(1, exponentFraction.Denominator));
                    List<Number> roots =
                        new List<Number> { Product.create(new Integer(1), factors) * coefficient };
                    for (int i = 1; i < exponentFraction.Denominator / largestIndex; ++i)
                        roots.Add(roots[i - 1] * rootOfUnity);
                    returnAllRoots = true;
                    return new NumberList(roots);
                }
                return Product.create(new Integer(1), factors) * coefficient;
            }
            return base.exponentiate(exponent);
        }
    }
    class Integer : WholeNumber
    {
        public int Value { get; }
        public Integer(int value)
        {
            Value = value;
        }
        public override int getNormSquared()
        {
            return Value * Value;
        }
        protected override Number add(Number number)
        {
            if (number is Integer)
            {
                return new Integer(Value + ((Integer)number).Value);
            }
            return number + this;
        }
        public override Number negative()
        {
            return new Integer(-Value);
        }
        protected override Number multiply(Number number)
        {
            if (number is Integer)
                return new Integer(Value * ((Integer)number).Value);
            return number * this;
        }
        public override Number reciprocal()
        {
            return Fraction.create(1, Value);
        }
        public override Number exponentiate(Number exponent)
        {
            if (Value == 0)
                return new Integer(0);
            return exponentiateWholeNumber(exponent);
        }
        public override int CompareTo(object obj)
        {
            int comparison = base.CompareTo(obj);
            if (comparison != 0)
                return comparison;
            return Value - ((Integer)obj).Value;
        }
        public String format()
        {
            if (Value != 1)
                if (Value == -1)
                    return "-";
                else
                    return Value.ToString();
            return "";
        }
        public override int GetHashCode()
        {
            return Value;
        }
        public override String insertString(String str)
        {
            return format() + str;
        }
        public override String ToString()
        {
            return Value.ToString();
        }
        public static int getGCD(int a, int b)
        {
            while (b != 0)
            {
                int c = b;
                b = a % b;
                a = c;
            }
            if (a < 0)
                a *= -1;
            return a;
        }
        public static int getLCM(int x, int y)
        {
            return x / getGCD(x, y) * y;
        }
    }
    class GaussianInteger : WholeNumber
    {
        public int Real { get; }
        public int Imaginary { get; }
        public GaussianInteger(int real, int imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }
        public static WholeNumber create(int real, int imaginary)
        {
            if (imaginary == 0)
                return new Integer(real);
            return new GaussianInteger(real, imaginary);
        }
        public override int getNormSquared()
        {
            return Real * Real + Imaginary * Imaginary;
        }
        protected override Number add(Number number)
        {
            if (number is Integer)
                return new GaussianInteger(Real + ((Integer)number).Value, Imaginary);
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return create(Real + gaussianInteger.Real, Imaginary + gaussianInteger.Imaginary);
            }
            return number + this;
        }
        public override Number negative()
        {
            return new GaussianInteger(-Real, -Imaginary);
        }
        protected override Number multiply(Number number)
        {
            if (number is Integer)
            {
                Integer integer = (Integer)number;
                return create(Real * integer.Value, Imaginary * integer.Value);
            }
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return create(Real * gaussianInteger.Real - Imaginary * gaussianInteger.Imaginary,
                    Imaginary * gaussianInteger.Real + Real * gaussianInteger.Imaginary);
            }
            return number * this;
        }
        public override Number reciprocal()
        {
            int denominator = Real * Real + Imaginary * Imaginary;
            return ComplexNumber.create(Fraction.create(Real, denominator),
                Fraction.create(-Imaginary, denominator));
        }
        public override Number exponentiate(Number exponent)
        {
            return exponentiateWholeNumber(exponent);
        }
        public override int CompareTo(object obj)
        {
            int comparison = base.CompareTo(obj);
            if (comparison != 0)
                return comparison;
            GaussianInteger gaussianInteger = (GaussianInteger)obj;
            comparison = Real - gaussianInteger.Real;
            if (comparison != 0)
                return comparison;
            return Imaginary - gaussianInteger.Imaginary;
        }
        public override int GetHashCode()
        {
            return Real ^ Imaginary;
        }
        public override String insertString(String str)
        {
            if (Real == 0)
                return new Integer(Imaginary).insertString("i*" + str);
            return '(' + ToString() + ')' + str;
        }
        public override String ToString()
        {
            StringBuilder output = new StringBuilder();
            if (Real != 0)
            {
                output.Append(Real);
                if (Imaginary > 0)
                    output.Append('+');
            }
            return output.Append(new Integer(Imaginary).insertString("i")).ToString();
        }
    }
    class Fraction : Rational
    {
        public int Numerator { get; }
        public int Denominator { get; }
        Fraction(int numerator, int denominator)
        {
            Numerator = numerator;
            Denominator = denominator;
        }
        public static Rational create(int numerator, int denominator)
        {
            if (denominator == 0)
                throw new DivideByZeroException();
            int GCD = Integer.getGCD(numerator, denominator);
            numerator /= GCD;
            denominator /= GCD;
            if (denominator < 0)
            {
                numerator *= 1;
                denominator *= 1;
            }
            if (denominator == 1)
                return new Integer(numerator);
            return new Fraction(numerator, denominator);
        }
        protected override Number add(Number number)
        {
            if (number is Integer)
                return create(Numerator + ((Integer)number).Value * Denominator, Denominator);
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(create(Numerator + gaussianInteger.Real * Denominator,
                    Denominator), new Integer(gaussianInteger.Imaginary));
            }
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
                return create(Numerator * ((Integer)number).Value, Denominator);
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(create(Numerator * gaussianInteger.Real, Denominator),
                    create(Numerator * gaussianInteger.Imaginary, Denominator));
            }
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
        public override int getDenominatorLCM()
        {
            return Denominator;
        }
        public override int CompareTo(object obj)
        {
            int comparison = base.CompareTo(obj);
            if (comparison != 0)
                return comparison;
            Fraction fraction = (Fraction)obj;
            comparison = Numerator - fraction.Numerator;
            if (comparison != 0)
                return comparison;
            return Denominator - fraction.Denominator;
        }
        public override int GetHashCode()
        {
            return Numerator ^ Denominator;
        }
        public override String insertString(String str)
        {
            return new Integer(Numerator).insertString(str) + "/" + Denominator;
        }
        public override String ToString()
        {
            return Numerator + "/" + Denominator;
        }
    }
    class Exponentiation : Factor
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
            if (number is Integer)
            {
                Integer integer = (Integer)number;
                if (integer.Value == 0)
                    return this;
                return new Sum(new List<Term> { this, integer });
            }
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(this + new Integer(gaussianInteger.Real),
                    new Integer(gaussianInteger.Imaginary));
            }
            if (number is Fraction)
                return new Sum(new List<Term> { this, (Fraction)number });
            if (number is Exponentiation)
            {
                Exponentiation exponentiation = (Exponentiation)number;
                if (CompareTo(exponentiation) == 0)
                    return Product.create(new Integer(2), new List<Factor> { this });
                if (CompareTo(exponentiation.negative()) == 0)
                    return new Integer(0);
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
            if (number is Integer)
                return Product.create((Integer)number, new List<Factor> { this });
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(Product.create(new Integer(gaussianInteger.Real),
                    new List<Factor> { this }), Product.create(new Integer(
                    gaussianInteger.Imaginary), new List<Factor> { this }));
            }
            if (number is Fraction)
                return Product.create((Fraction)number, new List<Factor> { this });
            if (number is Exponentiation)
            {
                Exponentiation exponentiation = (Exponentiation)number;
                if (Base.CompareTo(exponentiation.Base) == 0) 
                    return Base.exponentiate(Exponent + exponentiation.Exponent);
                if (Exponent.CompareTo(exponentiation.Exponent) == 0)
                {
                    Number outputBase = Base * exponentiation.Base;
                    if (outputBase is Exponentiation)
                        return new Exponentiation(((Exponentiation)outputBase).Base,
                            Exponent * Exponent);
                    return (outputBase).exponentiate(Exponent);
                }
                return Product.create(new Integer(1), new List<Factor> { this, exponentiation });
            }
            return number * this;
        }
        public override Number reciprocal()
        {
            if (Exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)Exponent;
                return exponentiate(Fraction.create(exponentFraction.Denominator -
                    exponentFraction.Numerator % exponentFraction.Denominator,
                    exponentFraction.Denominator)) * exponentiate(new Integer(
                    exponentFraction.Numerator / exponentFraction.Denominator + 1)).reciprocal();
            }
            return new Exponentiation(Base, Exponent.negative());
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
                if (number is Integer)
                {
                    Integer integer = (Integer)number;
                    if (integer.Value >= 0)
                        return integer.ToString();
                }
                if (number is GaussianInteger)
                {
                    GaussianInteger gaussian = (GaussianInteger)number;
                    if (gaussian.Real == 0 && gaussian.Imaginary == 1)
                        return "i";
                }
                return '(' + number.ToString() + ')';
            }
            return encloseNumber(Base) + '^' + encloseNumber(Exponent);
        }
    }
    class ComplexExponential : Factor
    {//Respresents e^(tau*i*Exponent).
        public Number Exponent { get; }
        ComplexExponential(Number exponent)
        {
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                Exponent = exponentFraction -
                    new Integer(exponentFraction.Numerator / exponentFraction.Denominator);
                if (exponentFraction.Numerator < 0)
                    Exponent += new Integer(1);
            }
            else
                Exponent = exponent;
        }
        public static Number create(Number exponent)
        {
            if (exponent is Integer)
                return new Integer(1);
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                int denominator = exponentFraction.Denominator;
                int numerator = exponentFraction.Numerator % denominator;
                int multiplicityOfTwo = 0;
                while (denominator % 2 == 0)
                {
                    denominator /= 2;
                    ++multiplicityOfTwo;
                }
                Number cosine;
                Rational a = Fraction.create(1, 2);
                switch (denominator)
                {
                    case 1:
                        cosine = new Integer(1);
                        break;
                    case 3:
                        cosine = Fraction.create(-1, 2);
                        break;
                    case 5:
                        cosine = (new Exponentiation(new Integer(5), a) - new Integer(1)) /
                            new Integer(4);
                        break;
                    case 15:
                        {
                            Exponentiation b = new Exponentiation(new Integer(5), a);
                            cosine = (new Integer(1) + b + new Exponentiation(new Integer(30) -
                                new Integer(6) * b, a)) / new Integer(8);
                            break;
                        }
                    case 17:
                        {
                            Integer b = new Integer(17);
                            Exponentiation c = new Exponentiation(b, a);
                            Integer d = new Integer(2);
                            Exponentiation e = new Exponentiation(b * d - c * new Integer(2), a);
                            cosine = (new Integer(-1) + c + e + new Exponentiation(b + c *
                                new Integer(3) - e - new Exponentiation((b + c) * d, a) * d, a) *
                                d) * Fraction.create(1, 16);
                            break;
                        }
                    default:
                        return new ComplexExponential(exponent);
                }
                for (int i = 0; i < multiplicityOfTwo; ++i)
                {
                    denominator *= 2;
                    if (4 < 3 * denominator && denominator < 4)
                        cosine = (cosine * a + a).exponentiate(a).negative();
                    else
                        cosine = (cosine * a + a).exponentiate(a);
                }
                Number s = new Integer(1);
                Number t = cosine;
                Number angleMultipleCosine = cosine;
                for (int i = 1; i < numerator; ++i)
                {
                    angleMultipleCosine = new Integer(2) * cosine * t - s;
                    s = t;
                    t = angleMultipleCosine;
                }
                return ComplexNumber.create(angleMultipleCosine, (new Integer(1) -
                    angleMultipleCosine * angleMultipleCosine).exponentiate(a).negative());
            }
            return new ComplexExponential(exponent);
        }
        protected override Number add(Number number)
        {
            if (number is Integer)
            {
                Integer integer = (Integer)number;
                if (integer.Value == 0)
                    return this;
                return new Sum(new List<Term> { this, integer });
            }
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return new ComplexNumber(this + new Integer(gaussianInteger.Real),
                    new Integer(gaussianInteger.Imaginary));
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
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(Product.create(new Integer(gaussianInteger.Real),
                    new List<Factor> { this }), Product.create(new Integer(
                    gaussianInteger.Imaginary), new List<Factor> { this }));
            }
            if (number is Fraction)
                return Product.create((Fraction)number, new List<Factor> { this, });
            if (number is Exponentiation)
                return Product.create(new Integer(1), new List<Factor> { this, (Factor)number });
            if (number is ComplexExponential)
                return create(Exponent + ((ComplexExponential)number).Exponent);
            return number * this;
        }
        public override Number reciprocal()
        {
            return create(Exponent.negative());
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
            Number exponentWith_i = Exponent * new GaussianInteger(0, 1);
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
                int coefficientInt = ((Integer)coefficient).Value;
                if (coefficientInt == 0)
                    return new Integer(0);
                if (factors.Count == 1 && coefficientInt == 1)
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
                if (integer.Value == 0)
                    return this;
                return new Sum(new List<Term> { this, integer });
            }
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return new ComplexNumber(this + new Integer(gaussianInteger.Real),
                    new Integer(gaussianInteger.Imaginary));
            }
            if (number is Factor) 
                return new Sum(new List<Term> { this, (Term)number });
            if (number is Product)
            {
                Product product = (Product)number;
                if (Factors.Count == product.Factors.Count)
                {
                    for (int i = 0; i < Factors.Count; ++i)
                        if (Factors[i].CompareTo(product.Factors[i]) != 0)
                            return new Sum(new List<Term> { this, product });
                    return new Product((Rational)(Coefficient + product.Coefficient), Factors);
                }
                return new Sum(new List<Term> { this, product });
            }
            return number + this;
        }
        public override Number negative()
        {
            return new Product((Rational)Coefficient.negative(), Factors);
        }
        protected override Number multiply(Number number)
        {
            if (number is Rational)
                return create((Rational)(Coefficient * number), Factors);
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return ComplexNumber.create(new Product((Rational)(Coefficient *
                    new Integer(gaussianInteger.Real)), Factors), new Product((Rational)(
                    Coefficient * new Integer(gaussianInteger.Imaginary)), Factors));
            }
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
                Number output = new Product((Rational)(Coefficient * product.Coefficient), Factors);
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
        public override int getDenominatorLCM()
        {
            return Coefficient.getDenominatorLCM();
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
            {
                string factorString = Factors[i].ToString();
                if (output[output.Length - 1] != ')' && factorString[0] != '(')
                    output.Append('*');
                output.Append(factorString);
            }
            return output;
        }
        string addCoefficient(StringBuilder factors)
        {
            int numerator;
            if (Coefficient is Integer)
                numerator = ((Integer)Coefficient).Value;
            else
            {
                Fraction coefficientFraction = (Fraction)Coefficient;
                factors.Append("/" + coefficientFraction.Denominator);
                numerator = coefficientFraction.Numerator;
            }
            if (numerator != 1)
            {
                if (numerator == -1)
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
                Number output = new Integer(0);
                foreach (Term term in Terms)
                    output += term * number;
                return output;
            }
            if (number is Sum)
            {
                Sum multiplier = (Sum)number;
                Number output = new Integer(0);
                foreach (Term multiplicandTerm in Terms)
                    foreach (Term multiplierTerm in multiplier.Terms)
                        output += multiplicandTerm * multiplierTerm;
                return output;
            }
            return number * this;
        }
        public override Number reciprocal()
        {//Placeholder.
            return new Exponentiation(this, new Integer(-1));
        }
        public override int getDenominatorLCM()
        {
            int LCM = 1;
            foreach (Term term in Terms)
                LCM = Integer.getLCM(LCM, term.getDenominatorLCM());
            return LCM;
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
        public ComplexNumber(Number real, Number imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }
        public static Number create(Number real, Number imaginary)
        {
            if (imaginary is Integer)
            {
                int imaginaryInt = ((Integer)imaginary).Value;
                if (imaginaryInt == 0)
                    return real;
                if (real is Integer)
                    return new GaussianInteger(((Integer)real).Value, imaginaryInt);
            }
            return new ComplexNumber(real, imaginary);
        }
        protected override Number add(Number number)
        {
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                return create(Real + new Integer(gaussianInteger.Real),
                    Imaginary + new Integer(gaussianInteger.Imaginary));
            }
            if (number is ComplexNumber)
            {
                ComplexNumber complexNumber = (ComplexNumber)number;
                return create(Real + complexNumber.Real, Imaginary + complexNumber.Imaginary);
            }
            if (!(number is NumberList))
                return new ComplexNumber(Real + number, Imaginary);
            return number + this;
        }
        public override Number negative()
        {
            return new ComplexNumber(Real.negative(), Imaginary.negative());
        }
        protected override Number multiply(Number number)
        {
            if (number is GaussianInteger)
            {
                GaussianInteger gaussianInteger = (GaussianInteger)number;
                Integer real = new Integer(gaussianInteger.Real);
                Integer imaginary = new Integer(gaussianInteger.Imaginary);
                return create(Real * real - Imaginary * imaginary,
                    Real * imaginary + Imaginary * real);
            }
            if (number is ComplexNumber)
            {
                ComplexNumber complexNumber = (ComplexNumber)number;
                return create(Real * complexNumber.Real - Imaginary * complexNumber.Imaginary,
                    Real * complexNumber.Imaginary + Imaginary * complexNumber.Real);
            }
            if (!(number is NumberList))
                return create(Real * number, Imaginary * number);
            return number * this;
        }
        public override Number reciprocal()
        {
            Number denominator = Real * Real + Imaginary * Imaginary;
            return create(Real / denominator, Imaginary.negative() / denominator);
        }
        public override int getDenominatorLCM()
        {
            return Integer.getLCM(Real.getDenominatorLCM(),
                Imaginary.getDenominatorLCM());
        }
        public override Number exponentiate(Number exponent)
        {
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                Integer LCM = new Integer(getDenominatorLCM());
                if (LCM.Value != 1)
                    return (this * LCM.exponentiate(new Integer(
                        exponentFraction.Denominator - 1))).exponentiate(exponentFraction) /
                        LCM.exponentiate(new Integer(exponentFraction.Numerator));
            }
            return base.exponentiate(exponent);
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
            if (Real is Integer && ((Integer)Real).Value == 0) 
                return Imaginary.insertString("i*" + str);
            return '(' + ToString() + ')' + str;
        }
        public override string ToString()
        { 
            StringBuilder output =
                new StringBuilder(Imaginary.insertString("i").ToString());
            if (!(Real is Integer && ((Integer)Real).Value == 0)) 
            {
                if (output[0] != '-')
                    output.Insert(0, '+');
                output.Insert(0, Real.ToString());
            }
            return output.ToString();
        }
    }
    class NumberList : Number
    {
        public List<Number> Numbers { get; }
        public NumberList(List<Number> numbers)
        {
            Numbers = numbers;
            Numbers.Sort();
        }
        protected override Number add(Number number)
        {
            List<Number> outputList = new List<Number>();
            if (number is NumberList)
            {
                List<Number> list = ((NumberList)number).Numbers;
                foreach (Number n in Numbers)
                    foreach (Number m in list)
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
        public override Number negative()
        {
            List<Number> output = new List<Number>();
            foreach (Number n in Numbers)
                output.Add(n.negative());
            return new NumberList(output);
        }
        protected override Number multiply(Number number)
        {
            List<Number> outputList = new List<Number>();
            if (number is NumberList)
            {
                List<Number> list = ((NumberList)number).Numbers;
                foreach (Number n in Numbers)
                    foreach (Number m in list)
                    {
                        Number product = n * m;
                        if (!outputList.Contains(product))
                            outputList.Add(product);
                    }
            }
            else
            {
                if (number is Integer)
                {
                    if (((Integer)number).Value == 0)
                        return number;
                }
                foreach (Number n in Numbers)
                {
                    Number product = n * number;
                    if (!outputList.Contains(product))
                        outputList.Add(product);
                }
            }
            return new NumberList(outputList);
        }
        public override Number reciprocal()
        {
            List<Number> output = new List<Number>();
            foreach (Number n in Numbers)
                output.Add(n.reciprocal());
            return new NumberList(output);
        }
        public override Number exponentiate(Number exponent)
        {
            List<Number> output = new List<Number>();
            foreach (Number n in Numbers)
                output.Add(n.exponentiate(exponent));
            return new NumberList(output);
        }
        public override int CompareTo(object obj)
        {
            int comparison = base.CompareTo(obj);
            if (comparison != 0)
                return comparison;
            NumberList list = (NumberList)obj;
            comparison = Numbers.Count - list.Numbers.Count;
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
        public override string insertString(string str)
        {
            throw new NotImplementedException("NumberList.insertString should never be called.");
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
    class Program
    {
        static Number evaluateExpression(List<char> operations, List<Number> numbers)
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
                if (i - 1 < 0 || numbers[i - 1] == null)
                {
                    if (operations[i] == '+')
                    {
                        numbers.RemoveAt(i);
                        operations.RemoveAt(i);
                    }
                    else if (operations[i] == '-')
                    {
                        numbers[i] = new Integer(-1);
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
                                int number = int.Parse(numberCollector.ToString());
                                numbers.Add(GaussianInteger.create(0, number));
                                operations.Add(' ');
                            }
                            else
                            {
                                numbers.Add(new Integer(int.Parse(numberCollector.ToString())));
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
                            numbers.Add(new GaussianInteger(0, 1));
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
                    numbers.Add(new Integer(int.Parse(numberCollector.ToString())));
                }
                for (int i = 0; i < operations.Count;)
                {
                    if (operations[i] == '(' && i > 0 && numbers[i - 1] != null)
                    {
                        operations.Insert(i, '*');
                        numbers.Insert(i, null);
                        i += 2;
                    }
                    else if ((operations[i] == ')' || numbers[i] is GaussianInteger) &&
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
