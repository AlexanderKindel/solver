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
        public virtual int getGreatestIntegerFactor()
        {
            return 1;
        }
        public virtual Number exponentiate(Number exponent)
        {
            if (exponent is Integer)
            {
                int intExponent = ((Integer)exponent).Value;
                Number output = new Integer(1);
                if (intExponent < 0)
                {
                    for (int i = 0; i > intExponent; --i)
                        output = this * output;
                    return output.reciprocal();
                }
                for (int i = 0; i < intExponent; ++i)
                    output = this * output;
                return output;
            }
            Integer exponentIntegerFactor = new Integer(exponent.getGreatestIntegerFactor());
            if (exponentIntegerFactor.Value > 1)
                return exponentiate(exponentIntegerFactor).exponentiate(
                    exponent / exponentIntegerFactor);
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                Integer denominatorLCM = new Integer(getDenominatorLCM());
                if (denominatorLCM.Value != 1)
                    return (this * denominatorLCM).exponentiate(exponent) *
                        denominatorLCM.exponentiate(Fraction.create(exponentFraction.Denominator -
                        exponentFraction.Numerator % exponentFraction.Denominator,
                        exponentFraction.Denominator)) / denominatorLCM.exponentiate(new Integer(
                        exponentFraction.Numerator / exponentFraction.Denominator + 1));
                Integer numeratorGCD = new Integer(getGreatestIntegerFactor());
                if (numeratorGCD.Value != 1)
                {
                    Number factor = numeratorGCD.exponentiate(exponent);
                    if (!(factor is Exponentiation))
                        return (this / numeratorGCD).exponentiate(exponent) * factor;
                }
            }
            return Exponentiation.create(this, exponent);
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
            if (this is GaussianInteger)
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
    class Integer : Rational
    {
        public int Value { get; }
        public Integer(int value)
        {
            Value = value;
        }
        public int getNormSquared()
        {
            return Value * Value;
        }
        protected override Number add(Number number)
        {
            if (number is Integer)
                return new Integer(Value + ((Integer)number).Value);
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
            Integer exponentIntegerFactor = new Integer(exponent.getGreatestIntegerFactor());
            if (exponentIntegerFactor.Value > 1)
                return base.exponentiate(exponentIntegerFactor).exponentiate(
                    exponent / exponentIntegerFactor);
            if (exponent is Fraction)
            {
                Fraction exponentFraction = (Fraction)exponent;
                int radicand = Value;
                Dictionary<Integer, int> radicandFactors = new Dictionary<Integer, int> { };
                if (radicand < 0)
                {
                    radicand *= -1;
                    radicandFactors.Add(new Integer(-1), 1);
                }
                Integer factor;
                for (int i = 2; i <= radicand / 2;)                
                    if (radicand % i == 0) 
                    {
                        factor = new Integer(i);
                        if (radicandFactors.ContainsKey(factor))
                            radicandFactors[factor] += 1;
                        else
                            radicandFactors.Add(factor, 1);
                        radicand /= i;
                    }
                    else
                        ++i;
                factor = new Integer(radicand);
                if (radicandFactors.ContainsKey(factor))
                    radicandFactors[factor] += 1;
                else
                    radicandFactors.Add(factor, 1);
                List<int> exponentDivisors = new List<int>();
                int divisor = 1;
                while (divisor <= exponentFraction.Denominator / 2)
                {
                    if (exponentFraction.Denominator % divisor == 0)
                        exponentDivisors.Add(divisor);
                    ++divisor;
                }
                exponentDivisors.Add(exponentFraction.Denominator);
                Dictionary<int, Number> termComponents = new Dictionary<int, Number>();
                termComponents.Add(1, new Integer(1));
                termComponents.Add(exponentFraction.Denominator, new Integer(1));
                foreach (Integer n in radicandFactors.Keys)
                    for (int i = exponentDivisors.Count - 1; i >= 0; --i)
                        if (exponentDivisors[i] <= radicandFactors[n])
                        {
                            int index = exponentFraction.Denominator / exponentDivisors[i];
                            if (!termComponents.ContainsKey(index))
                                termComponents.Add(index, new Integer(1));
                            for (int j = 0;
                                j < radicandFactors[n] / exponentDivisors[i]; ++j)
                                termComponents[index] = termComponents[index] * n;
                            for (int j = 0;
                                j < radicandFactors[n] % exponentDivisors[i]; ++j)
                                termComponents[exponentFraction.Denominator] =
                                    termComponents[exponentFraction.Denominator] * n;
                            break;
                        }
                Number coefficient = termComponents[1];
                termComponents.Remove(1);
                int highestIndexComponent =
                    ((Integer)termComponents[exponentFraction.Denominator]).Value;
                if (highestIndexComponent == 1)
                    termComponents.Remove(exponentFraction.Denominator);
                else if (highestIndexComponent < 0)
                {
                    if (highestIndexComponent == -1)
                        termComponents.Remove(exponentFraction.Denominator);
                    else
                        termComponents[exponentFraction.Denominator] =
                            termComponents[exponentFraction.Denominator].negative();
                    if (exponentFraction.Denominator % 2 == 0)                                            
                        coefficient *= ComplexExponential.create(
                            Fraction.create(1, 2 * exponentFraction.Denominator));                    
                    else                    
                        coefficient = coefficient.negative();                    
                }                
                List<Factor> factors = new List<Factor>();
                int largestIndex = 1;
                foreach (int index in termComponents.Keys)
                {
                    factors.Add(new Surd(termComponents[index], index));
                    if (index > largestIndex)
                        largestIndex = index;
                }
                Number output = Product.create(new Integer(1), factors) * coefficient;
                return output;                
            }
            return base.exponentiate(exponent);
        }
        public override int getGreatestIntegerFactor()
        {
            return Value;
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
        public override int getGreatestIntegerFactor()
        {
            return Numerator;
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
                if (exponentFraction.Numerator == 1)
                    return new Surd(expBase, exponentFraction.Denominator);
            }
            return new Transcendental(expBase, exponent);
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
                return Product.create(new Integer(1), new List<Factor> { this, exponentiation });
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
                if (number is Integer)
                {
                    Integer integer = (Integer)number;
                    if (integer.Value >= 0)
                        return integer.ToString();
                }
                if (number is GaussianInteger)
                {
                    GaussianInteger gaussian = (GaussianInteger)number;
                    if (((Integer)gaussian.Real).Value == 0 &&
                        ((Integer)gaussian.Imaginary).Value == 1)
                        return "i";
                }
                return '(' + number.ToString() + ')';
            }
            return encloseNumber(Base) + '^' + encloseNumber(Exponent);
        }
    }
    class Surd : Exponentiation
    {
        public int Index { get; }
        public Surd(Number radicand, int index) : base(radicand, Fraction.create(1, index))
        {
            Index = index;
        }
        protected override Number multiply(Number number)
        {
            if (number is Surd)
            {
                Surd surd = (Surd)number;
                int indexLCM = Integer.getLCM(Index, surd.Index);
                return (Base.exponentiate(new Integer(indexLCM / Index)) * surd.Base.exponentiate(
                    new Integer(indexLCM / surd.Index))).exponentiate(Fraction.create(1, indexLCM));
            }
            return base.multiply(number);
        }
        public override Number reciprocal()
        {
            return exponentiate(Fraction.create(Index - 1, Index)) / this;
        }
    }
    class Transcendental : Exponentiation
    {
        public Transcendental(Number expBase, Number exponent) : base(expBase, exponent)
        { }
        public override Number reciprocal()
        {
            return Exponentiation.create(Base, Exponent.negative());
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
                switch (denominator)
                {
                    case 1:
                        cosine = new Integer(1);
                        break;
                    case 3:
                        cosine = Fraction.create(-1, 2);
                        break;
                    case 5:
                        cosine = (new Surd(new Integer(5), 2) - new Integer(1)) /
                            new Integer(4);
                        break;
                    case 15:
                        {
                            Surd a = new Surd(new Integer(5), 2);
                            cosine = (new Integer(1) + a + new Surd(new Integer(30) -
                                new Integer(6) * a, 2)) / new Integer(8);
                            break;
                        }
                    case 17:
                        {
                            Integer a = new Integer(17);
                            Surd b = new Surd(a, 2);
                            Integer c = new Integer(2);
                            Surd d = new Surd(a * c - b * new Integer(2), 2);
                            cosine = (new Integer(-1) + b + d + new Surd(a + b * new Integer(3) -
                                d - new Surd((a + b) * c, 2) * c, 2) * c) * Fraction.create(1, 16);
                            break;
                        }
                    default:
                        return new ComplexExponential(exponent);
                }
                Rational e = Fraction.create(1, 2);
                for (int i = 0; i < multiplicityOfTwo; ++i)
                {
                    denominator *= 2;
                    cosine = (cosine * e + e).exponentiate(e);
                    if (4 < 3 * denominator && denominator < 4)
                        cosine = cosine.negative();
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
                Number output = ComplexNumber.create(angleMultipleCosine, (new Integer(1) -
                    angleMultipleCosine * angleMultipleCosine).exponentiate(e));
                if (denominator < 4 * numerator && 4 * numerator < 2 * denominator ||
                    3 * denominator < 4 * numerator && 4 * numerator < 4 * denominator)
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
                if (integer.Value == 0)
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
            Number exponentWith_i = Exponent * new GaussianInteger(new Integer(0), new Integer(1));
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
            if (number is Factor)
            {
                if (Factors.Count == 1 && number.Equals(Factors[0]))
                    return create((Rational)(Coefficient + new Integer(1)), Factors);
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
        public override int getGreatestIntegerFactor()
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
            return new Transcendental(this, new Integer(-1));
        }
        public override int getDenominatorLCM()
        {
            int LCM = 1;
            foreach (Term term in Terms)
                LCM = Integer.getLCM(LCM, term.getDenominatorLCM());
            return LCM;
        }
        public override int getGreatestIntegerFactor()
        {
            int GCD = Terms[0].getGreatestIntegerFactor();
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
                if (imaginaryInteger.Value == 0)
                    return real;
                if (real is Integer)
                    return new GaussianInteger((Integer)real, imaginaryInteger);
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
        public override int getDenominatorLCM()
        {
            return Integer.getLCM(Real.getDenominatorLCM(), Imaginary.getDenominatorLCM());
        }
        public override int getGreatestIntegerFactor()
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
    class GaussianInteger : ComplexNumber
    {
        public GaussianInteger(Integer real, Integer imaginary) : base(real, imaginary)
        { }
        public int getNormSquared()
        {
            int real = ((Integer)Real).Value;
            int imaginary = ((Integer)Imaginary).Value;
            return real * real + imaginary * imaginary;
        }
    }
    class Polynomial
    {
        //The index of the coefficient represents the degree of its term.
        List<Number> Coefficients { get; }
        public Polynomial(List<Number> coefficients)
        {
            Integer zero = new Integer(0);
            while (coefficients.Count != 0 && coefficients[coefficients.Count - 1].Equals(zero))
                coefficients.RemoveAt(coefficients.Count - 1);
            Coefficients = coefficients;
        }
        public static Polynomial operator *(Polynomial a, Polynomial b)
        {
            List<Number> output = new List<Number>();
            Integer zero = new Integer(0);
            for (int i = 0; i < a.Coefficients.Count + b.Coefficients.Count - 1; ++i)
                output.Add(zero);
            for (int i = 0; i < a.Coefficients.Count; ++i)
                for (int j = 0; j < b.Coefficients.Count; ++j)
                    output[i + j] += a.Coefficients[i] * b.Coefficients[j];
            return new Polynomial(output);
        }
        public Polynomial modulo(Polynomial divisor)
        {
            List<Number> remainder = Coefficients;
            while (remainder.Count >= divisor.Coefficients.Count)
            {
                Number multiplier = remainder[remainder.Count - 1] /
                    divisor.Coefficients[divisor.Coefficients.Count - 1];
                remainder.RemoveAt(remainder.Count - 1);
                for (int i = 1; i < divisor.Coefficients.Count - 1; ++i)
                    remainder[remainder.Count - i] -=
                        multiplier * divisor.Coefficients[divisor.Coefficients.Count - 1 - i];
            }
            return new Polynomial(remainder);
        }
    }
    class Program
    {
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
                        numbers[i] = new Integer(1);
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
                                    numbers.Add(ComplexNumber.create(new Integer(0),
                                        new Integer(number)));
                                    operations.Add(' ');
                                }
                                else if (!"()+-*/^".Contains(c.ToString()))
                                    throw new InvalidUserInput(c + " is an invalid character.");
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
                        else if (c == 'i')
                        {
                            numbers.Add(new GaussianInteger(new Integer(0), new Integer(1)));
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
                    Number expression = evaluateExpression(operations, numbers);
                    if (expression == null)
                        Console.WriteLine("=\n");
                    else
                        Console.Write(
                            "=\n" + evaluateExpression(operations, numbers).ToString() + "\n\n");
                }
                catch (InvalidUserInput e)
                {
                    Console.WriteLine(e.Message);
                }
                catch (DivideByZeroException e)
                {
                    Console.WriteLine(e.Message);
                }
            }
        }
    }
}
