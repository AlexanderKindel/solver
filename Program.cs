using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Text;

namespace Solver
{
    static class IntervalConversion
    {
        public static Solver.Interval<Solver.Rational> ToRational(
            this Solver.Interval<Solver.Float> a)
        {
            return new Solver.Interval<Solver.Rational>((Solver.Rational)a.Min,
                (Solver.Rational)a.Max);
        }
    }
    static class Solver
    {
        public struct Interval<T>
        {
            public T Min;
            public T Max;
            public Interval(T min, T max)
            {
                Min = min;
                Max = max;
            }
        }
        public struct RectangularEstimate
        {
            public Interval<Float> RealPart;
            public Interval<Float> ImaginaryPart;
            public RectangularEstimate(Interval<Float> realPart, Interval<Float> imaginaryPart)
            {
                RealPart = realPart;
                ImaginaryPart = imaginaryPart;
            }
        }
        public struct Division<T>
        {
            public T Quotient;
            public T Remainder;
            public Division(T quotient, T remainder)
            {
                Quotient = quotient;
                Remainder = remainder;
            }
        }
        struct ExtendedGCDInfo<T>
        {
            public T GCD;
            public T ACoefficient;
            public T BCoefficient;
            public T AOverGCD;
            public T BOverGCD;
        }
        public delegate Interval<Float> EstimateGetter<T>(T number, Rational errorIntervalSize);
        delegate T Multiplier<T>(T a, T b);
        interface IMultipliable<S, T>
        {
            S Times(T n);
        }
        interface IRingElement<T> : IMultipliable<T, T>
        {
            T GetMultiplicativeIdentity();
        }
        interface IArithmetic<T> : IRingElement<T>, IEquatable<T>
        {
            T GetAdditiveIdentity();
            T GetNegative();
            T Plus(T n);
            T Minus(T n);
            Division<T> EuclideanDivideBy(T divisor);
        }
        public abstract class Primitive
        {
            RationalPolynomial MinimalPolynomial_ = null;
            public RationalPolynomial MinimalPolynomial
            {
                get
                {
                    if (MinimalPolynomial_ == null)
                    {                        
                        MinimalPolynomial_ = CalculateMinimalPolynomial();
                    }
                    return MinimalPolynomial_;
                }
            }
            protected abstract RationalPolynomial CalculateMinimalPolynomial();
            protected void RemoveNonConjugates(List<Number> conjugateCandidates)
            {
                for (int i = 1;
                    conjugateCandidates.Count > MinimalPolynomial.Coefficients.Count - 1;)
                {
                    if (!conjugateCandidates[i].MinimalPolynomial.Equals(MinimalPolynomial))
                    {
                        conjugateCandidates.RemoveAt(i);
                    }
                    else
                    {
                        ++i;
                    }
                }
            }
            public List<Number> GetSumConjugates(List<Term> terms)
            {
                List<List<Number>> termConjugates = new List<List<Number>>();
                foreach (Term term in terms)
                {
                    termConjugates.Add(term.GetConjugates());
                }
                List<List<Number>> conjugateCombinations = GetCartesianProduct(termConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = Number.Zero;
                    foreach (Number number in combination)
                    {
                        conjugate += number;
                    }
                    conjugateCandidates.Add(conjugate);
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            public abstract List<Number> GetConjugates();//The first conjugate is *this.
            public RationalPolynomial ArgumentInTermsOfThis(Primitive a)
            {
                if (a is Rational rational)
                {
                    return new RationalPolynomial(MinimalPolynomial, rational);
                }
                if ((MinimalPolynomial.Coefficients.Count - 1) %
                    (a.MinimalPolynomial.Coefficients.Count - 1) != 0)
                {
                    return null;
                }
                NestedPolynomial minimalPolynomial = new NestedPolynomial(null,
                    a.MinimalPolynomial).SwitchVariables(MinimalPolynomial);
                List<NestedPolynomial> factors = minimalPolynomial.GetFactors();
                List<RationalPolynomial> candidateFactors = new List<RationalPolynomial>();
                foreach (NestedPolynomial polynomial in factors)
                {
                    if (polynomial.Coefficients.Count == 2)
                    {
                        candidateFactors.Add(-polynomial.Coefficients[0]);
                    }
                }
                List<Number> conjugates = a.GetConjugates();
                conjugates.RemoveAt(0);
                Rational errorIntervalSize = Number.One;
                for (int i = 0; i < candidateFactors.Count; ++i)
                {
                    List<Number> candidateConjugates = new List<Number>(conjugates);
                    bool factorEqualsA()
                    {
                        while (candidateConjugates.Count > 0)
                        {
                            Interval<Float> aRealPartEstimate =
                                a.GetRealPartEstimate(errorIntervalSize);
                            Interval<Float> aImaginaryPartEstimate =
                                a.GetImaginaryPartEstimate(errorIntervalSize);
                            RectangularEstimate factorEvaluatedAtThis =
                                candidateFactors[i].EstimateEvaluation(this, errorIntervalSize);
                            if (factorEvaluatedAtThis.RealPart.Min > aRealPartEstimate.Max ||
                                factorEvaluatedAtThis.RealPart.Max < aRealPartEstimate.Min ||
                                factorEvaluatedAtThis.ImaginaryPart.Min >
                                aImaginaryPartEstimate.Max ||
                                factorEvaluatedAtThis.ImaginaryPart.Max <
                                aImaginaryPartEstimate.Min)
                            {
                                return false;
                            }
                            for (int j = 0; j < candidateConjugates.Count;)
                            {
                                Interval<Float> candidateRealPartEstimate =
                                    candidateConjugates[j].GetRealPartEstimate(errorIntervalSize);
                                Interval<Float> candidateImaginaryPartEstimate =
                                    candidateConjugates[j].GetImaginaryPartEstimate(
                                    errorIntervalSize);
                                if (factorEvaluatedAtThis.RealPart.Min >
                                    candidateRealPartEstimate.Max ||
                                    factorEvaluatedAtThis.RealPart.Max <
                                    candidateRealPartEstimate.Min ||
                                    factorEvaluatedAtThis.ImaginaryPart.Min >
                                    candidateImaginaryPartEstimate.Max ||
                                    factorEvaluatedAtThis.ImaginaryPart.Max <
                                    candidateImaginaryPartEstimate.Min)
                                {
                                    candidateConjugates.RemoveAt(j);
                                }
                                else
                                {
                                    ++j;
                                }
                            }
                            errorIntervalSize /= Number.Two;
                        }
                        return true;
                    }
                    if (factorEqualsA())
                    {
                        return new RationalPolynomial(candidateFactors[i].Coefficients,
                            MinimalPolynomial);
                    }
                }
                return null;
            }
            public abstract Interval<Float> GetRealPartEstimate(Rational errorIntervalSize);
            public static Interval<Float> GetRealPartEstimate(Primitive a,
                Rational errorIntervalSize)
            {
                return a.GetRealPartEstimate(errorIntervalSize);
            }
            public abstract Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize);
            public static Interval<Float> GetImaginaryPartEstimate(Primitive a,
                Rational errorIntervalSize)
            {
                return a.GetImaginaryPartEstimate(errorIntervalSize);
            }
            protected delegate Interval<Rational> RationalEstimateGetter(
                Rational errorIntervalSize);
            protected Interval<Float> GetRefinedErrorInterval(Rational errorIntervalSize,
                RationalEstimateGetter getEstimate)
            {
                Rational placeValue = Float.GetEstimatePlaceValue(errorIntervalSize);
                Interval<Rational> rationalEstimate = getEstimate(placeValue);
                return new Interval<Float>(rationalEstimate.Min.GetRealPartEstimate(placeValue).Min,
                    rationalEstimate.Max.GetRealPartEstimate(placeValue).Max);
            }
            protected Interval<Float> EstimatePartSum<T>(List<T> terms, Rational errorIntervalSize,
                EstimateGetter<T> getPartEstimate) where T : Primitive
            {
                Interval<Float> sumEstimate = new Interval<Float>(Float.Zero, Float.Zero);
                for (int i = 0; i < terms.Count; ++i)
                {
                    Interval<Float> partEstimate =
                        getPartEstimate(terms[i], errorIntervalSize / new Integer(terms.Count - i));
                    sumEstimate.Min += partEstimate.Min;
                    sumEstimate.Max += partEstimate.Max;
                    errorIntervalSize -= (Rational)(partEstimate.Max - partEstimate.Min);
                }
                return sumEstimate;
            }
            public Interval<Float> EstimateRealPartSum<T>(List<T> terms, Rational errorIntervalSize)
                where T : Number
            {
                return EstimatePartSum(terms, errorIntervalSize, GetRealPartEstimate);
            }
            public Interval<Float> EstimateImaginaryPartSum<T>(List<T> terms,
                Rational errorIntervalSize) where T : Number
            {
                return EstimatePartSum(terms, errorIntervalSize, GetImaginaryPartEstimate);
            }
        }
        public abstract class Number : Primitive, IArithmetic<Number>, IComparable
        {
            public static Integer Zero = new Integer(0);
            public static Integer One = new Integer(1);
            public static Integer Two = new Integer(2);            
            public Number GetAdditiveIdentity()
            {
                return Zero;
            }
            public Number GetMultiplicativeIdentity()
            {
                return One;
            }
            public Number GetNegative()
            {
                return this * new Integer(-1);
            }
            public abstract Number Reciprocal();
            public abstract Rational RationalFactor();
            public abstract Number Plus(Number a);
            public Number Minus(Number a)
            {
                return this + a * new Integer(-1);
            }
            public abstract Number Times(Number a);
            public Division<Number> EuclideanDivideBy(Number a)
            {
                return new Division<Number>(this / a, Zero);
            }
            public abstract Number Exponentiate(Rational exponent);
            protected abstract int GetTypeIndex();
            public virtual int CompareTo(object obj)
            {//Essentially a cheap way to check ToString() == obj.ToString(), which is guaranteed to
             //correspond to numerical equality for Rationals, but not for other derived types. To
             //test for numerical equality of arbitrary Numbers, check whether their difference
             //evaluates to Zero, and to numerically compare unequal Numbers, compare their Estimate
             //fields.
                Number number = obj as Number;
                if (ReferenceEquals(number, null))
                {
                    return 1;
                }
                return GetTypeIndex() - number.GetTypeIndex();
            }
            public abstract override int GetHashCode();
            public override bool Equals(object obj)
            {
                return CompareTo(obj) == 0;
            }
            public bool Equals(Number number)
            {
                return CompareTo(number) == 0;
            }
            public static Number operator +(Number a, Number b)
            {
                return a.Plus(b);
            }
            public static Number operator -(Number a, Number b)
            {
                return a + b * new Integer(-1);
            }
            public static Number operator -(Number a)
            {
                return a * new Integer(-1);
            }
            public static Number operator *(Number a, Number b)
            {
                return a.Times(b);
            }
            public static Number operator /(Number a, Number b)
            {
                return a * b.Reciprocal();
            }
            public static bool operator ==(Number a, Number b)
            {
                if (ReferenceEquals(a, null))
                {
                    if (ReferenceEquals(b, null))
                    {
                        return true;
                    }
                    return false;
                }
                return a.CompareTo(b) == 0;
            }
            public static bool operator !=(Number a, Number b)
            {
                return !(a == b);
            }
            public abstract Interval<Float> GetArgumentEstimate(Rational errorIntervalSize);
            public static Interval<Float> GetArgumentEstimate(Number number,
                Rational errorIntervalSize)
            {
                return number.GetArgumentEstimate(errorIntervalSize);
            }
            public abstract Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize);
            protected Interval<Rational> EstimateCosineOfArgument(Rational errorIntervalSize)
            {
                errorIntervalSize /= new Integer(3);
                Interval<Rational> argumentEstimate =
                    GetArgumentEstimate(errorIntervalSize).ToRational();                
                Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimate.Min);
                Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimate.Max);
                Interval<Rational> cosine;
                if (argumentEstimate.Min <= Pi.LowEstimate &&
                    Pi.HighEstimate <= argumentEstimate.Max)
                {
                    cosine.Min = -One;
                    cosine.Max = argumentEstimate.Min.EstimateCosine(errorIntervalSize).Max;
                    Rational possibleCosineMax =
                        argumentEstimate.Max.EstimateCosine(errorIntervalSize).Max;
                    if (possibleCosineMax > cosine.Max)
                    {
                        cosine.Max = possibleCosineMax;
                    }
                }
                else
                {
                    cosine = argumentEstimate.Min.EstimateCosine(errorIntervalSize);
                    Interval<Rational> otherEstimate =
                        argumentEstimate.Max.EstimateCosine(errorIntervalSize);
                    if (otherEstimate.Min < cosine.Min)
                    {
                        cosine.Min = otherEstimate.Min;
                    }
                    if (otherEstimate.Max > cosine.Max)
                    {
                        cosine.Max = otherEstimate.Max;
                    }
                }
                return cosine;
            }
            protected Interval<Rational> EstimateSineOfArgument(Rational errorIntervalSize)
            {
                errorIntervalSize /= new Integer(3);
                Interval<Rational> argumentEstimate =
                    GetArgumentEstimate(errorIntervalSize).ToRational();                
                Interval<Rational> argumentEstimateMultiple =
                    new Interval<Rational>(Two * argumentEstimate.Min, Two * argumentEstimate.Max);
                Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimateMultiple.Min);
                Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimateMultiple.Max);
                Interval<Rational> sine;
                if (argumentEstimateMultiple.Min <= Pi.LowEstimate &&
                    Pi.HighEstimate <= argumentEstimateMultiple.Max)
                {
                    sine.Max = One;
                    sine.Min = argumentEstimate.Min.EstimateSine(errorIntervalSize).Min;
                    Rational possibleSineMin =
                        argumentEstimate.Max.EstimateSine(errorIntervalSize).Min;
                    if (possibleSineMin < sine.Min)
                    {
                        sine.Min = possibleSineMin;
                    }
                }
                else
                {
                    Integer three = new Integer(3);
                    argumentEstimateMultiple.Min /= three;
                    argumentEstimateMultiple.Max /= three;
                    Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimateMultiple.Min);
                    Pi.ShrinkErrorIntervalToOneSideOfValue(argumentEstimateMultiple.Max);
                    if (argumentEstimateMultiple.Min <= Pi.LowEstimate &&
                        Pi.LowEstimate <= argumentEstimateMultiple.Max)
                    {
                        sine.Min = -One;
                        sine.Max = argumentEstimate.Min.EstimateSine(errorIntervalSize).Max;
                        Rational possibleSineMax =
                            argumentEstimate.Max.EstimateSine(errorIntervalSize).Max;
                        if (possibleSineMax > sine.Max)
                        {
                            sine.Max = possibleSineMax;
                        }
                    }
                    else
                    {
                        sine = argumentEstimate.Min.EstimateSine(errorIntervalSize);
                        Interval<Rational> otherEstimate =
                            argumentEstimate.Min.EstimateSine(errorIntervalSize);
                        if (otherEstimate.Min < sine.Min)
                        {
                            sine.Min = otherEstimate.Min;
                        }
                        if (otherEstimate.Max > sine.Max)
                        {
                            sine.Max = otherEstimate.Max;
                        }
                    }
                }
                return sine;
            }
            protected delegate Interval<Rational> TrigFunction(Rational errorIntervalSize);
            protected Interval<Rational> EstimateRectangularPartFromPolarForm(
                Rational errorIntervalSize, TrigFunction sineOrCosine)
            {
                errorIntervalSize /= (Rational)GetMagnitudeEstimate(One).Max.Magnitude() + Two;
                Interval<Rational> magnitudeEstimate =
                    GetMagnitudeEstimate(errorIntervalSize).ToRational();
                Interval<Rational> trigValue = sineOrCosine(errorIntervalSize);
                if (trigValue.Min >= Zero)
                {
                    return new Interval<Rational>(trigValue.Min * magnitudeEstimate.Min,
                        trigValue.Max * magnitudeEstimate.Max);
                }
                else if (trigValue.Max <= Zero)
                {
                    return new Interval<Rational>(trigValue.Min * magnitudeEstimate.Max,
                        trigValue.Max * magnitudeEstimate.Min);
                }
                return new Interval<Rational>(trigValue.Min * magnitudeEstimate.Max,
                    trigValue.Max * magnitudeEstimate.Max);
            }
            public bool HasZeroImaginaryPart()
            {
                Rational errorIntervalSize = One;
                Interval<Float> imaginaryPartEstimate = GetImaginaryPartEstimate(errorIntervalSize);
                if (imaginaryPartEstimate.Max < Float.Zero ||
                    Float.Zero < imaginaryPartEstimate.Min)
                {
                    return false;
                }
                List<Number> conjugates = GetConjugates();
                conjugates.RemoveAt(0);
                while (conjugates.Count > 1)
                {
                    Interval<Float> realPartEstimate = GetRealPartEstimate(errorIntervalSize);
                    for (int i = 0; i < conjugates.Count; ++i)
                    {
                        Interval<Float> conjugateRealPartEstimate =
                            conjugates[i].GetRealPartEstimate(errorIntervalSize);
                        if (conjugateRealPartEstimate.Max < realPartEstimate.Min ||
                            realPartEstimate.Max < conjugateRealPartEstimate.Min)
                        {
                            conjugates.RemoveAt(i);
                        }
                        else
                        {
                            imaginaryPartEstimate = GetImaginaryPartEstimate(errorIntervalSize);
                            Interval<Float> conjugateImaginaryPartEstimate =
                                conjugates[i].GetImaginaryPartEstimate(errorIntervalSize);
                            if (-conjugateImaginaryPartEstimate.Min < imaginaryPartEstimate.Min ||
                                imaginaryPartEstimate.Max < -conjugateImaginaryPartEstimate.Max)
                            {
                                conjugates.RemoveAt(i);
                            }
                            else
                            {
                                ++i;
                            }
                        }
                    }
                    errorIntervalSize /= Two;
                }
                while (true)
                {
                    Interval<Float> realPartEstimate = GetRealPartEstimate(errorIntervalSize);
                    Interval<Float> conjugateRealPartEstimate =
                        conjugates[0].GetRealPartEstimate(errorIntervalSize);
                    if (conjugateRealPartEstimate.Max < realPartEstimate.Min ||
                        realPartEstimate.Max < conjugateRealPartEstimate.Min)
                    {
                        return true;
                    }
                    imaginaryPartEstimate = GetImaginaryPartEstimate(errorIntervalSize);
                    if (imaginaryPartEstimate.Min < Float.Zero ||
                        Float.Zero < imaginaryPartEstimate.Max)
                    {
                        return false;
                    }
                    errorIntervalSize /= Two;
                }
            }
            public abstract override string ToString();
        }
        public abstract class Term : Number
        {
            public override Number Plus(Number a)
            {
                if (a == Zero)
                {
                    return this;
                }
                RationalPolynomial x = new RationalPolynomial(null, Zero, One);
                if (a is Rational rational)
                {
                    return new LinearlyIndependentSum(new List<Term> { this, rational },
                        new List<RationalPolynomial> { x, new RationalPolynomial(null, rational) },
                        this);
                }
                if (a is Term term)
                {
                    RationalPolynomial termInTermsOfThis = ArgumentInTermsOfThis(term);
                    if (termInTermsOfThis != null)
                    {
                        if (termInTermsOfThis.Coefficients.Count == 2 &&
                            termInTermsOfThis.Coefficients[0] == Zero)
                        {
                            return (termInTermsOfThis.Coefficients[1] + One) * this;
                        }
                        return new LinearlyIndependentSum(new List<Term> { this, term },
                            new List<RationalPolynomial> { x, termInTermsOfThis }, this);
                    }
                    return new LinearlyIndependentSum(this, term);
                }
                return a + this;
            }
            public abstract Term CheapPlus(Term a);
            protected abstract Term TimesRational(Rational a);
            public static Term operator *(Term a, Rational b)
            {
                return a.TimesRational(b);
            }
            public static Term operator *(Rational a, Term b)
            {
                return b.TimesRational(a);
            }
        }
        public abstract class Rational : Term, IArithmetic<Rational>
        {
            public abstract Integer Numerator { get; }
            public abstract Integer Denominator { get; }
            protected static Rational Create(Integer numerator, Integer denominator)
            {
                if (denominator == Zero)
                {
                    throw new DivideByZeroException();
                }
                if (numerator == Zero)
                {
                    return Zero;
                }
                ExtendedGCDInfo<Integer> extendedGCD = ExtendedGCD(numerator, denominator);
                numerator = extendedGCD.AOverGCD;
                denominator = extendedGCD.BOverGCD;
                if (denominator < Zero)
                {
                    numerator = -numerator;
                    denominator = -denominator;
                }
                if (denominator == One)
                {
                    return numerator;
                }
                return new Fraction(numerator, denominator);
            }
            new public Rational GetAdditiveIdentity()
            {
                return Zero;
            }
            new public Rational GetMultiplicativeIdentity()
            {
                return One;
            }
            new public Rational GetNegative()
            {
                return this * new Integer(-1);
            }
            public Rational Magnitude()
            {
                if (this >= Zero)
                {
                    return this;
                }
                return -this;
            }
            public override Rational RationalFactor()
            {
                return this;
            }
            public Rational Plus(Rational a)
            {
                return Create(Numerator * a.Denominator + a.Numerator * Denominator,
                    Denominator * a.Denominator);
            }
            public override Term CheapPlus(Term a)
            {
                if (a is Rational rational)
                {
                    return Plus(rational);
                }
                return a.CheapPlus(this);
            }
            public Rational Minus(Rational a)
            {
                return this - a;
            }
            protected override Term TimesRational(Rational a)
            {
                return Create(Numerator * a.Numerator, Denominator * a.Denominator);
            }
            public Rational Times(Rational a)
            {
                return Create(Numerator * a.Numerator, Denominator * a.Denominator);
            }
            public Division<Rational> EuclideanDivideBy(Rational divisor)
            {
                return new Division<Rational>(this / divisor, Zero);
            }
            public bool Equals(Rational a)
            {
                return CompareTo(a) == 0;
            }
            public static Rational operator +(Rational a, Rational b)
            {
                return a.Plus(b);
            }
            public static Rational operator -(Rational a, Rational b)
            {
                return a + -b;
            }
            public static Rational operator -(Rational a)
            {
                return a.GetNegative();
            }
            public static Rational operator *(Rational a, Rational b)
            {
                return a.Times(b);
            }
            public static Rational operator /(Rational a, Rational b)
            {
                return Create(a.Numerator * b.Denominator, a.Denominator * b.Numerator);
            }
            public static bool operator <(Rational a, Rational b)
            {
                return a.Numerator * b.Denominator < b.Numerator * a.Denominator;
            }
            public static bool operator >(Rational a, Rational b)
            {
                return a.Numerator * b.Denominator > b.Numerator * a.Denominator;
            }
            public static bool operator <=(Rational a, Rational b)
            {
                return !(a > b);
            }
            public static bool operator >=(Rational a, Rational b)
            {
                return !(a < b);
            }
            public static explicit operator Rational(Float a)
            {
                return a.Significand / Solver.Exponentiate(Two, a.RadixPointPosition);
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                return new RationalPolynomial(null, -this, One);
            }
            public override List<Number> GetConjugates()
            {
                return new List<Number> { this };
            }
            public Interval<Rational> EstimateCosine(Rational errorIntervalSize)
            {
                Rational cosineValue = One;
                Rational thisSquared = this * this;
                Integer factorialComponent = Two;
                Rational delta = -thisSquared / factorialComponent;
                while (delta.Magnitude() > errorIntervalSize)
                {
                    cosineValue += delta;
                    ++factorialComponent;
                    delta *= thisSquared / factorialComponent;
                    ++factorialComponent;
                    delta /= -factorialComponent;
                }
                if (delta > Zero)
                {
                    return new Interval<Rational>(cosineValue, cosineValue + delta);
                }
                return new Interval<Rational>(cosineValue + delta, cosineValue);
            }
            public Interval<Rational> EstimateSine(Rational errorIntervalSize)
            {
                Rational sineValue = this;
                Rational thisSquared = this * this;
                Integer factorialComponent = new Integer(3);
                Rational delta = -thisSquared * this / new Integer(6);
                while (delta.Magnitude() > errorIntervalSize)
                {
                    sineValue += delta;
                    ++factorialComponent;
                    delta *= thisSquared / factorialComponent;
                    ++factorialComponent;
                    delta /= -factorialComponent;
                }
                if (delta > Zero)
                {
                    return new Interval<Rational>(sineValue, sineValue + delta);
                }
                return new Interval<Rational>(sineValue + delta, sineValue);
            }
            public Interval<Rational> EstimateArctangent(Rational errorIntervalSize)
            {
                Debug.Assert(Magnitude() <= One, "EstimateArctangent was called on *this whose " +
                    "magnitude was greater than One. Use an identity to evaluate an equivalent " +
                    "expression involving its reciprocal instead.");
                Rational arctangentValue = this;
                Rational thisSquared = -this * this;
                Rational deltaNumerator = this * thisSquared;
                Integer deltaDenominator = new Integer(3);
                Rational delta = deltaNumerator / deltaDenominator;
                while (delta.Magnitude() > errorIntervalSize)
                {
                    arctangentValue += delta;
                    deltaNumerator *= thisSquared;
                    deltaDenominator += Two;
                    delta = deltaNumerator / deltaDenominator;
                }
                if (delta > Zero)
                {
                    return new Interval<Rational>(arctangentValue, arctangentValue + delta);
                }
                return new Interval<Rational>(arctangentValue + delta, arctangentValue);
            }
            Interval<Rational> GetRationalArgumentEstimate(Rational errorIntervalSize)
            {
                if (Numerator > Zero)
                {
                    return new Interval<Rational>(Zero, Zero);
                }
                Pi.RefineErrorInterval(errorIntervalSize);
                return new Interval<Rational>(Pi.LowEstimate, Pi.HighEstimate);
            }
            public override Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize)
            {
                return new Interval<Float>(Float.Zero, Float.Zero);
            }
            public override Interval<Float> GetArgumentEstimate(Rational errorIntervalSize)
            {
                return GetRefinedErrorInterval(errorIntervalSize, GetRationalArgumentEstimate);
            }

            //Treats str as the string representation of a numerical constant, and returns the
            //string representation of the reference object multiplied by that constant.
            public abstract string InsertString(string str);
        }
        public class Integer : Rational, IArithmetic<Integer>
        {
            uint[] Values;//little-endian
            sbyte Sign;
            public override Integer Numerator { get => this; }
            public override Integer Denominator { get => One; }
            Integer(uint[] values, sbyte sign)
            {
                int lastNonzeroValueIndex = -1;
                for (int i = values.Length - 1; i >= 0; --i)
                {
                    if (values[i] != 0)
                    {
                        lastNonzeroValueIndex = i;
                        break;
                    }
                }
                Values = new uint[lastNonzeroValueIndex + 1];
                for (int i = 0; i <= lastNonzeroValueIndex; ++i)
                {
                    Values[i] = values[i];
                }
                if (lastNonzeroValueIndex < 0)
                {
                    Sign = 0;
                }
                else
                {
                    Sign = sign;
                }
            }
            public Integer(int value)
            {
                if (value == 0)
                {
                    Sign = 0;
                    Values = new uint[] { };
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
                    output += multiplier *
                        new Integer(int.Parse(decimalString.Substring(startIndex, 9)));
                    multiplier *= multiplierStepSize;
                    startIndex -= 9;
                }
                if (startIndex > -9)
                {
                    output += multiplier *
                        new Integer(int.Parse(decimalString.Substring(0, startIndex + 9)));
                }
                return output;
            }
            new public Integer GetAdditiveIdentity()
            {
                return Zero;
            }
            new public Integer GetMultiplicativeIdentity()
            {
                return One;
            }
            new public Integer GetNegative()
            {
                return new Integer(Values, (sbyte)-Sign);
            }
            new public Integer Magnitude()
            {
                if (Sign < 0)
                {
                    return new Integer(Values, 1);
                }
                return this;
            }
            public override Number Reciprocal()
            {
                if (Sign == 0)
                {
                    throw new DivideByZeroException();
                }
                Integer magnitude = Magnitude();
                if (magnitude == One)
                {
                    return this;
                }
                return new Fraction(new Integer(new uint[] { 1 }, Sign), magnitude);
            }
            public Integer Reciprocal(Integer characteristic)
            {
                Integer inverse = ExtendedGCD(this, characteristic).ACoefficient;
                if (inverse < Zero)
                {
                    inverse += characteristic;
                }
                return inverse;
            }
            public Integer Plus(Integer a)
            {
                uint[] aValues;
                uint[] bValues;
                uint[] sumValues;
                void calculateValues()
                {
                    for (int i = 0; i < aValues.Length; ++i)
                    {
                        if ((aValues[i] & 0x80000000) != 0)
                        {
                            aValues[i] -= 0x80000000;
                            if ((bValues[i] & 0x80000000) != 0)
                            {
                                bValues[i] -= 0x80000000;
                                sumValues[i] += aValues[i] + bValues[i];
                                sumValues[i + 1] += 1;
                            }
                            else
                            {
                                sumValues[i] += aValues[i] + bValues[i];
                                if ((sumValues[i] & 0x80000000) != 0)
                                {
                                    sumValues[i] -= 0x80000000;
                                    sumValues[i + 1] += 1;
                                }
                                else
                                {
                                    sumValues[i] += 0x80000000;
                                }
                            }
                        }
                        else
                        {
                            if ((bValues[i] & 0x80000000) != 0)
                            {
                                bValues[i] -= 0x80000000;
                                sumValues[i] += aValues[i] + bValues[i];
                                if ((sumValues[i] & 0x80000000) != 0)
                                {
                                    sumValues[i] -= 0x80000000;
                                    sumValues[i + 1] += 1;
                                }
                                else
                                {
                                    sumValues[i] += 0x80000000;
                                }
                            }
                            else
                            {
                                sumValues[i] += aValues[i] + bValues[i];
                            }
                        }
                    }
                }
                uint[] takeTwosComplement(uint[] values)
                {
                    uint[] outputValues = new uint[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        outputValues[i] = ~values[i];
                    }
                    for (int i = 0; i < values.Length; ++i)
                    {
                        uint power = 1;
                        for (int j = 0; j < 32; ++j)
                        {
                            outputValues[i] ^= power;
                            if ((values[i] & power) != 0)
                            {
                                return outputValues;
                            }
                            power = power << 1;
                        }
                    }
                    return outputValues;
                }
                if (Values.Length > a.Values.Length)
                {
                    aValues = new uint[Values.Length];
                }
                else
                {
                    aValues = new uint[a.Values.Length];
                }
                Values.CopyTo(aValues, 0);
                bValues = new uint[aValues.Length];
                a.Values.CopyTo(bValues, 0);
                sumValues = new uint[aValues.Length + 1];
                Integer handleSingleNegativeCase()
                {
                    sumValues[aValues.Length] = 0xffffffff;
                    calculateValues();
                    if (sumValues[aValues.Length] == 0)
                    {
                        return new Integer(sumValues, 1);
                    }
                    return new Integer(takeTwosComplement(sumValues), -1);
                }
                if (Sign < 0)
                {
                    aValues = takeTwosComplement(aValues);
                    if (a.Sign < 0)
                    {
                        bValues = takeTwosComplement(bValues);
                        sumValues[aValues.Length] = 0xfffffffe;
                        calculateValues();
                        return new Integer(takeTwosComplement(sumValues), -1);
                    }
                    return handleSingleNegativeCase();
                }
                if (a.Sign < 0)
                {
                    bValues = takeTwosComplement(bValues);
                    return handleSingleNegativeCase();
                }
                calculateValues();
                return new Integer(sumValues, 1);
            }
            public override Number Plus(Number a)
            {
                if (a is Integer integer)
                {
                    return integer + this;
                }
                return a + this;
            }
            public Integer Minus(Integer a)
            {
                return this - a;
            }
            public Integer Times(Integer a)
            {
                Integer product = Zero;
                for (int i = 0; i < Values.Length; ++i)
                {
                    for (int j = 0; j < a.Values.Length; ++j)
                    {
                        ulong productComponent = (ulong)(Values[i]) * a.Values[j];
                        product = product + new Integer(ShiftValuesLeft(new uint[] {
                            (uint)(productComponent & 0x00000000ffffffff),
                            (uint)((productComponent & 0xffffffff00000000) >> 32) }, i + j, 0), 1);
                    }
                }
                product.Sign = (sbyte)(Sign * a.Sign);
                return product;
            }
            public override Number Times(Number a)
            {
                if (a is Integer integer)
                {
                    return integer * this;
                }
                return a * this;
            }            
            public Division<Integer> EuclideanDivideBy(Integer divisor)
            {                
                if (divisor.Sign == 0)
                {
                    throw new DivideByZeroException();
                }
                Division<Integer> division = new Division<Integer>();
                Integer positiveDivisor = divisor.Magnitude();
                if (positiveDivisor > Magnitude())
                {
                    division.Quotient = Zero;
                    division.Remainder = this;
                    return division;
                }
                division.Remainder = new Integer(Values, 1);
                int divisorLeadingDigitPlace = 0;
                uint power = 0x80000000;
                for (int i = 32; i > 0; --i)
                {
                    if ((divisor.Values[divisor.Values.Length - 1] & power) != 0)
                    {
                        divisorLeadingDigitPlace = i;
                        break;
                    }
                    power = power >> 1;
                }
                uint[] quotient = new uint[Values.Length];
                void calculateValue(int valuePlace, int stoppingDigitPlace)
                {
                    power = 0x80000000;
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
                {
                    calculateValue(i, 1);
                }
                calculateValue(positiveDivisor.Values.Length - 1, divisorLeadingDigitPlace);
                division.Quotient = new Integer(quotient, (sbyte)(Sign * divisor.Sign)).ShiftLeft(
                    1 - positiveDivisor.Values.Length, 1 - divisorLeadingDigitPlace);
                if (division.Remainder.Sign != 0)
                {
                    division.Remainder.Sign = Sign;
                }
                return division;
            }
            public override Number Exponentiate(Rational exponent)
            {
                if (this == Zero)
                {
                    return Zero;
                }
                if (exponent.Numerator < Zero)
                {
                    return Reciprocal().Exponentiate(-exponent);
                }
                Integer power = Solver.Exponentiate(this, exponent.Numerator);
                if (exponent.Denominator == One)
                {
                    return power;
                }
                Integer coefficient = One;
                Integer takeRootOfPositiveRadicand(Integer radicand)
                {
                    for (Integer i = radicand.EuclideanDivideBy(Two).Quotient; i >= One; --i)
                    {
                        Division<Integer> division = radicand.EuclideanDivideBy(
                            Solver.Exponentiate(i, exponent.Denominator));
                        if (division.Remainder == Zero)
                        {
                            radicand = division.Quotient;
                            coefficient = i;
                            break;
                        }
                    }
                    return radicand;
                }
                Surd surd;
                if (power < Zero)
                {
                    power = -takeRootOfPositiveRadicand(-power);
                    surd = new Surd(power, exponent.Denominator);
                }
                else
                {
                    power = takeRootOfPositiveRadicand(power);
                    surd = new Surd(power, exponent.Denominator);
                }
                if (power == One)
                {
                    return coefficient;
                }
                if (coefficient == One)
                {
                    return surd;
                }
                return new Product(coefficient, new List<Surd> { surd });
            }
            protected override int GetTypeIndex()
            {
                return 1;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                {
                    return comparison;
                }
                return (this - (Integer)obj).Sign;
            }
            public override int GetHashCode()
            {
                uint output = 0;
                foreach (uint value in Values)
                {
                    output ^= value;
                }
                return (int)(output * Sign);
            }
            public bool Equals(Integer a)
            {
                return CompareTo(a) == 0;
            }
            public static Integer operator +(Integer a, Integer b)
            {
                return a.Plus(b);
            }
            public static Integer operator ++(Integer a)
            {
                return a + One;
            }
            public static Integer operator -(Integer a, Integer b)
            {
                return a + -b;
            }
            public static Integer operator -(Integer a)
            {
                return a.GetNegative();
            }
            public static Integer operator --(Integer a)
            {
                return a + new Integer(-1);
            }
            public static Integer operator *(Integer a, Integer b)
            {
                return a.Times(b);
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
            public static explicit operator int(Integer a)
            {
                return (int)a.Values[0] * a.Sign;
            }
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                Float floatThis = new Float(this, Zero);
                return new Interval<Float>(floatThis, floatThis);
            }
            public override Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize)
            {
                Float magnitude = new Float(Magnitude(), Zero);
                return new Interval<Float>(magnitude, magnitude);
            }
            static uint[] ShiftValuesLeft(uint[] values, int valuePlaces, int digitPlaces)
            {//Negative valuePlaces and digitPlaces values shift to the right. Left shifts preserve
             //all digits by adding extra uints, while right shifts truncate the rightmost digits.
                uint[] shiftedValues;
                int smallestValuePlace;
                if (valuePlaces < 0)
                {
                    smallestValuePlace = -valuePlaces;
                }
                else
                {
                    smallestValuePlace = 0;
                }
                if (digitPlaces == 0)
                {
                    shiftedValues = new uint[values.Length + valuePlaces];
                    for (int i = smallestValuePlace; i < values.Length; ++i)
                    {
                        shiftedValues[valuePlaces + i] = values[i];
                    }
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
                        {
                            shiftedValues[valuePlaces + i - 1] += values[i] << 32 + digitPlaces;
                        }
                        shiftedValues[valuePlaces + i] += values[i] >> -digitPlaces;
                    }
                }
                return shiftedValues;
            }
            public Integer ShiftLeft(int valuePlaces, int digitPlaces)
            {
                return new Integer(ShiftValuesLeft(Values, valuePlaces, digitPlaces), Sign);
            }
            public Integer ThisChooseK(Integer k)
            {
                Integer numerator = One;
                for (Integer i = this - k + One; i <= this; ++i)
                {
                    numerator *= i;
                }
                Integer denominator = One;
                for (Integer i = Two; i <= k; ++i)
                {
                    denominator *= i;
                }
                return numerator.EuclideanDivideBy(denominator).Quotient;
            }
            public static Integer GetGCD(Integer a, Integer b)
            {
                return Solver.GetGCD(a, b).Magnitude();
            }
            public static Integer GetLCM(Integer a, Integer b)
            {
                return a.EuclideanDivideBy(GetGCD(a, b)).Quotient * b;
            }
            public override string InsertString(string str)
            {
                if (this != One)
                {
                    if (this == -One)
                    {
                        return '-' + str;
                    }
                    else if (str.Length == 0 || str[0] == '(')
                    {
                        return ToString() + str;
                    }
                    else
                    {
                        return ToString() + '*' + str;
                    }
                }
                return str;
            }
            public override string ToString()
            {
                if (this == Zero)
                {
                    return "0";
                }
                StringBuilder output = new StringBuilder();
                Integer quotient = this;
                Integer power = new Integer(10);
                while (quotient.Sign != 0)
                {
                    Division<Integer> division = quotient.EuclideanDivideBy(power);
                    quotient = division.Quotient;
                    if (division.Remainder.Sign != 0)
                    {
                        output.Insert(0, division.Remainder.Values[0]);
                    }
                    else
                    {
                        output.Insert(0, '0');
                    }
                }
                if (this < Zero)
                {
                    output.Insert(0, '-');
                }
                return output.ToString();
            }
        }
        public class Fraction : Rational
        {
            Integer Numerator_;
            Integer Denominator_;
            public override Integer Numerator { get => Numerator_; }
            public override Integer Denominator { get => Denominator_; }
            Integer EstimateNumerator = null;
            Integer EstimateDenominator = One;
            Integer EstimateRemainder = null;
            Integer RadixPointPosition = Zero;
            public Fraction(Integer numerator, Integer denominator)
            {//To be called only when the arguments are known a priori to form a canonical-form
             //fraction. Otherwise, use numerator / denominator instead.
                Numerator_ = numerator;
                Denominator_ = denominator;
            }
            new public Fraction GetNegative()
            {
                return new Fraction(-Numerator, Denominator);
            }
            public override Number Reciprocal()
            {
                return Create(Denominator, Numerator);
            }
            public override Number Plus(Number a)
            {
                if (a is Integer integer)
                {
                    return Create(integer * Denominator + Numerator, Denominator);
                }
                if (a is Fraction fraction)
                {
                    return Create(Numerator * fraction.Denominator +
                        fraction.Numerator * Denominator, Denominator * fraction.Denominator);
                }
                return a + this;
            }
            public override Number Times(Number a)
            {
                if (a is Integer integer)
                {
                    return Create(Numerator * integer, Denominator);
                }
                if (a is Fraction fraction)
                {
                    return Create(Numerator * fraction.Numerator,
                        Denominator * fraction.Denominator);
                }
                return a * this;
            }
            public override Number Exponentiate(Rational exponent)
            {
                return (Numerator * Denominator.Exponentiate(exponent.Denominator -
                    One)).Exponentiate(exponent) / Denominator.Exponentiate(exponent.Numerator);
            }
            protected override int GetTypeIndex()
            {
                return 2;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                {
                    return comparison;
                }
                Fraction fraction = (Fraction)obj;
                comparison = Numerator.CompareTo(fraction.Numerator);
                if (comparison != 0)
                {
                    return comparison;
                }
                return Denominator.CompareTo(fraction.Denominator);
            }
            public override int GetHashCode()
            {
                return Numerator.GetHashCode() ^ Denominator.GetHashCode();
            }
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                Interval<Float> magnitudeEstimate = GetMagnitudeEstimate(errorIntervalSize);
                if (Numerator < Zero)
                {
                    return new Interval<Float>(-magnitudeEstimate.Max, -magnitudeEstimate.Min);
                }
                return magnitudeEstimate;
            }
            public override Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize)
            {
                Division<Integer> division;
                Integer numeratorMagnitude = Numerator.Magnitude();
                if (EstimateNumerator == null)
                {
                    division = numeratorMagnitude.EuclideanDivideBy(Denominator);
                    EstimateNumerator = division.Quotient;
                    EstimateRemainder = division.Remainder;
                }
                Float estimate;
                while (errorIntervalSize.Denominator >
                    EstimateDenominator * errorIntervalSize.Numerator)
                {
                    division = (Two * EstimateRemainder).EuclideanDivideBy(Denominator);
                    EstimateNumerator = Two * EstimateNumerator + division.Quotient;
                    EstimateDenominator *= Two;
                    EstimateRemainder = division.Remainder;
                    ++RadixPointPosition;
                    if (EstimateRemainder == Zero)
                    {
                        estimate = new Float(EstimateNumerator, RadixPointPosition);
                        return new Interval<Float>(estimate, estimate);
                    }
                }
                estimate = Float.GetReducedForm(EstimateNumerator, RadixPointPosition);
                Integer crossMultiplicationA = EstimateNumerator * Denominator;
                Integer crossMultiplicationB = EstimateDenominator * numeratorMagnitude;
                if (crossMultiplicationA > crossMultiplicationB)
                {
                    return new Interval<Float>(Float.GetReducedForm(
                        EstimateNumerator - One, RadixPointPosition), estimate);
                }
                return new Interval<Float>(estimate,
                    Float.GetReducedForm(EstimateNumerator + One, RadixPointPosition));
            }
            public override string InsertString(string str)
            {
                return Numerator.InsertString(str) + "/" + Denominator;
            }
            public override string ToString()
            {
                return Numerator + "/" + Denominator;
            }
        }
        public class Surd : Term
        {
            public Number Radicand { get; }
            public Integer Index { get; }
            Interval<Float> RealPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ImaginaryPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ArgumentEstimate = new Interval<Float>(null, null);
            Interval<Float> MagnitudeEstimate = new Interval<Float>(null, null);
            public Surd(Number radicand, Integer index)
            {//To be called only when it's known a priori that Index > One, radicand is not an
             //Integer that is a perfect indexth power and, for odd index,
             //radicand.RationalFactor() == One, or even index,
             //radicand.RationalFactor().Magnitude() == One. Otherwise, use radicand.Exponentiate().
                Radicand = radicand;
                Index = index;
            }
            public override Number Reciprocal()
            {
                return Radicand.Exponentiate(new Fraction(Index - One, Index)) / Radicand;
            }
            public override Rational RationalFactor()
            {
                return One;
            }
            public override Term CheapPlus(Term a)
            {
                if (a is Rational)
                {
                    return null;
                }
                if (a is Surd surd && surd == this)
                {
                    return new Product(Two, new List<Surd> { this });
                }
                return a.CheapPlus(this);
            }
            protected override Term TimesRational(Rational a)
            {
                if (a == Zero)
                {
                    return Zero;
                }
                if (a == One)
                {
                    return this;
                }
                return new Product(a, new List<Surd> { this });
            }
            public override Number Times(Number a)
            {
                if (a is Rational rational)
                {
                    return TimesRational(rational);
                }
                if (a is Surd surd)
                {
                    ExtendedGCDInfo<Integer> GCDInfo = ExtendedGCD(Index, surd.Index);
                    Integer productIndex = ((Index * surd.Index) / GCDInfo.GCD).Numerator;
                    Number radicandFactor1 = Radicand.Exponentiate(GCDInfo.BOverGCD.Magnitude());
                    Number radicandFactor2 =
                        surd.Radicand.Exponentiate(GCDInfo.AOverGCD.Magnitude());
                    Interval<Float> radicandFactor1ArgumentEstimate =
                        radicandFactor1.GetArgumentEstimate(One);
                    Interval<Float> radicandFactor2ArgumentEstimate =
                        radicandFactor2.GetArgumentEstimate(One);
                    if (surd.Index != Index)
                    {
                        Rational radicandFactorArgumentErrorInterval =
                            Pi.LowEstimate / productIndex;
                        if ((Rational)radicandFactor1ArgumentEstimate.Max / productIndex <
                            (Rational)GetArgumentEstimate(radicandFactorArgumentErrorInterval).Min ||
                            (Rational)radicandFactor2ArgumentEstimate.Max / productIndex <
                            (Rational)surd.GetArgumentEstimate(
                            radicandFactorArgumentErrorInterval).Min)
                        {
                            return new Product(One, new List<Surd> { this, surd });
                        }
                    }
                    Number productRadicand = radicandFactor1 * radicandFactor2;
                    if (productRadicand.GetArgumentEstimate(Pi.LowEstimate).Max <
                        radicandFactor1ArgumentEstimate.Min + radicandFactor2ArgumentEstimate.Min)
                    {
                        return productRadicand.Exponentiate(new Fraction(One, productIndex)) *
                            RootsOfUnity.GetNthRoots(productIndex)[1];
                    }
                    return productRadicand.Exponentiate(new Fraction(One, productIndex));
                }
                return a * this;
            }
            public override Number Exponentiate(Rational exponent)
            {
                return Radicand.Exponentiate(exponent / Index);
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                {
                    return comparison;
                }
                Surd surd = (Surd)obj;
                comparison = Index.CompareTo(surd.Index);
                if (comparison != 0)
                {
                    return comparison;
                }
                return Radicand.CompareTo(surd.Radicand);
            }
            public override int GetHashCode()
            {
                return Radicand.GetHashCode() ^ Index.GetHashCode();
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                RationalPolynomial radicandMinimalPolynomial = Radicand.MinimalPolynomial;
                List<Rational> coefficients =
                    new List<Rational> { radicandMinimalPolynomial.Coefficients[0] };
                for (int i = 1; i < radicandMinimalPolynomial.Coefficients.Count; ++i)
                {
                    for (Integer j = One; j < Index; ++j)
                    {
                        coefficients.Add(Zero);
                    }
                    coefficients.Add(radicandMinimalPolynomial.Coefficients[i]);
                }
                List<IntegerPolynomial> candidates =
                    new RationalPolynomial(coefficients).GetPrimitivePart().GetFactors();
                if (candidates.Count == 1)
                {
                    return new RationalPolynomial(GetMonicCoefficients(candidates[0].Coefficients));
                }
                List<RationalPolynomial> rationalCandidates = new List<RationalPolynomial>();
                foreach (IntegerPolynomial candidate in candidates)
                {
                    rationalCandidates.Add(
                        new RationalPolynomial(GetMonicCoefficients(candidate.Coefficients)));
                }
                Rational precision = One;
                while (true)
                {
                    for (int i = 0; i < rationalCandidates.Count;)
                    {
                        RectangularEstimate candidateEstimate =
                            rationalCandidates[i].EstimateEvaluation(this, precision);
                        if (Float.Zero < candidateEstimate.RealPart.Min ||
                            candidateEstimate.RealPart.Max < Float.Zero ||
                            Float.Zero < candidateEstimate.ImaginaryPart.Min ||
                            candidateEstimate.ImaginaryPart.Max < Float.Zero)
                        {
                            rationalCandidates.RemoveAt(i);
                            if (rationalCandidates.Count == 1)
                            {
                                return rationalCandidates[0];
                            }
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    precision /= Two;
                }
            }
            public override List<Number> GetConjugates()
            {
                List<Number> radicandConjugates = Radicand.GetConjugates();
                List<Number> rootsOfUnity = RootsOfUnity.GetNthRoots(Index);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (Number conjugate in radicandConjugates)
                {
                    foreach (Number root in rootsOfUnity)
                    {
                        conjugateCandidates.Add(new Surd(conjugate, Index) * root);
                    }
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            protected override int GetTypeIndex()
            {
                return 3;
            }
            Interval<Rational> GetRationalRealEstimate(Rational errorIntervalSize)
            {
                return EstimateRectangularPartFromPolarForm(errorIntervalSize,
                    EstimateCosineOfArgument);
            }            
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                if (RealPartEstimate.Min == null || (Rational)(RealPartEstimate.Max -
                    RealPartEstimate.Min) > errorIntervalSize)
                {
                    RealPartEstimate =
                        GetRefinedErrorInterval(errorIntervalSize, GetRationalRealEstimate);
                }
                return RealPartEstimate;
            }
            Interval<Rational> GetRationalImaginaryEstimate(Rational errorIntervalSize)
            {
                return EstimateRectangularPartFromPolarForm(errorIntervalSize,
                    EstimateSineOfArgument);
            }            
            public override Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize)
            {
                if (ImaginaryPartEstimate.Min == null || (Rational)(ImaginaryPartEstimate.Max -
                    ImaginaryPartEstimate.Min) > errorIntervalSize)
                {
                    ImaginaryPartEstimate =
                        GetRefinedErrorInterval(errorIntervalSize, GetRationalImaginaryEstimate);
                }
                return ImaginaryPartEstimate;
            }
            Interval<Rational> GetRationalArgumentEstimate(Rational errorIntervalSize)
            {
                Interval<Rational> radicandArgumentEstimate =
                    Radicand.GetArgumentEstimate(Index * errorIntervalSize).ToRational();
                return new Interval<Rational>(radicandArgumentEstimate.Min / Index,
                    radicandArgumentEstimate.Max / Index);
            }            
            public override Interval<Float> GetArgumentEstimate(Rational errorIntervalSize)
            {
                if (ArgumentEstimate.Min == null || (Rational)(ArgumentEstimate.Max -
                    ArgumentEstimate.Min) > errorIntervalSize)
                {
                    ArgumentEstimate =
                        GetRefinedErrorInterval(errorIntervalSize, GetRationalArgumentEstimate);
                }
                return ArgumentEstimate;
            }            
            public override Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize)
            {
                if (MagnitudeEstimate.Min == null || (Rational)(MagnitudeEstimate.Max -
                    MagnitudeEstimate.Min) > errorIntervalSize)
                {
                    Integer three = new Integer(3);
                    Interval<Float> radicandMagnitudeEstimate = Radicand.GetMagnitudeEstimate(
                        Solver.Exponentiate(errorIntervalSize, Index) / three);
                    errorIntervalSize /= three;
                    MagnitudeEstimate = new Interval<Float>(
                        radicandMagnitudeEstimate.Min.EstimateRoot(errorIntervalSize, Index).Min,
                        radicandMagnitudeEstimate.Max.EstimateRoot(errorIntervalSize, Index).Max);
                }
                return MagnitudeEstimate;
            }
            public override string ToString()
            {
                if (Radicand is Integer integer && integer > Zero)
                {
                    return Radicand + "^(1/" + Index + ')';
                }
                return "(" + Radicand + ')' + "^(1/" + Index + ')';
            }
        }
        public class Product : Term
        {
            public Rational Coefficient { get; }
            public List<Surd> Factors { get; }
            Interval<Float> RealPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ImaginaryPartEstimate = new Interval<Float>(null, null);
            Interval<Float> MagnitudeEstimate = new Interval<Float>(null, null);
            Interval<Float> ArgumentEstimate = new Interval<Float>(null, null);
            public Product(Rational coefficient, List<Surd> factors)
            {
                Coefficient = coefficient;
                Factors = factors;
                Factors.Sort();
            }
            public override Number Reciprocal()
            {
                Number output = Coefficient.Reciprocal();
                foreach (Surd factor in Factors)
                {
                    output *= factor.Reciprocal();
                }
                return output;
            }
            public override Rational RationalFactor()
            {
                return Coefficient;
            }
            public override Term CheapPlus(Term a)
            {
                if (a is Rational)
                {
                    return null;
                }
                if (a is Surd surd)
                {
                    if (Factors.Count == 1 && Factors[0] == surd)
                    {
                        return new Product(Coefficient + One, Factors);
                    }
                    return null;
                }
                Product product = (Product)a;
                if (product.Factors.Count != Factors.Count)
                {
                    return null;
                }
                for (int i = 0; i < Factors.Count; ++i)
                {
                    if (product.Factors[i] != Factors[i])
                    {
                        return null;
                    }
                }
                Rational coefficient = Coefficient + product.Coefficient;
                if (coefficient == Zero)
                {
                    return Zero;
                }
                if (coefficient == One && Factors.Count == 1)
                {
                    return Factors[0];
                }
                return new Product(coefficient, Factors);
            }
            protected override Term TimesRational(Rational a)
            {
                if (a == Zero)
                {
                    return Zero;
                }
                if (a == One)
                {
                    return this;
                }
                Rational coefficient = Coefficient * a;
                if (Factors.Count == 1 && coefficient == One)
                {
                    return Factors[0];
                }
                return new Product(coefficient, Factors);
            }
            public override Number Times(Number a)
            {
                if (a is Rational rational)
                {
                    return TimesRational(rational);
                }
                if (a is Surd surd)
                {
                    List<Surd> factors = new List<Surd>(Factors);
                    for (int i = 0; i < Factors.Count; ++i)
                    {
                        Number factorProduct = surd * Factors[i];
                        if (factorProduct != new Product(One, new List<Surd> { surd, Factors[i] }))
                        {
                            factors.RemoveAt(i);
                            if (factors.Count == 0)
                            {
                                return Coefficient * factorProduct;
                            }
                            if (factors.Count == 1 && Coefficient == One)
                            {
                                return factors[0] * factorProduct;
                            }
                            return new Product(Coefficient, factors) * factorProduct;
                        }
                    }
                    factors.Add(surd);
                    return new Product(Coefficient, factors);
                }
                if (a is Product product)
                {
                    Number output = product.Coefficient * this;
                    foreach (Surd factor in product.Factors)
                    {
                        output *= factor;
                    }
                    return output;
                }
                return a * this;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                {
                    return comparison;
                }
                Product product = (Product)obj;
                comparison = Factors.Count - product.Factors.Count;
                if (comparison != 0)
                {
                    return comparison;
                }
                for (int i = 0; i < Factors.Count; ++i)
                {
                    comparison = Factors[i].CompareTo(product.Factors[i]);
                    if (comparison != 0)
                    {
                        return comparison;
                    }
                }
                return Coefficient.CompareTo(product.Coefficient);
            }
            public override Number Exponentiate(Rational exponent)
            {
                Number output = Coefficient.Exponentiate(exponent);
                foreach (Surd factor in Factors)
                {
                    output *= factor.Exponentiate(exponent);
                }
                return output;
            }
            public override int GetHashCode()
            {
                int output = Coefficient.GetHashCode();
                foreach (Surd surd in Factors)
                {
                    output ^= surd.GetHashCode();
                }
                return output;
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                RationalPolynomial[] minimalPolynomials = new RationalPolynomial[Factors.Count];
                int[] indicies = new int[Factors.Count];
                for (int i = 0; i < Factors.Count; ++i)
                {
                    minimalPolynomials[i] = Factors[i].MinimalPolynomial;
                    indicies[i] = 1;
                }
                MultivariatePolynomial variableForm = 
                    new MultivariatePolynomial(minimalPolynomials);
                variableForm.SetCoefficient(indicies, Coefficient);
                return variableForm.CalculateValueMinimalPolynomial(this);
            }
            public override List<Number> GetConjugates()
            {
                List<List<Number>> factorConjugates = new List<List<Number>>();
                foreach (Surd factor in Factors)
                {
                    factorConjugates.Add(factor.GetConjugates());
                }
                List<List<Number>> conjugateCombinations = GetCartesianProduct(factorConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = Coefficient;
                    foreach (Number number in combination)
                    {
                        conjugate *= number;
                    }
                    conjugateCandidates.Add(conjugate);
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            protected override int GetTypeIndex()
            {
                return 4;
            }            
            void RefineRectangularEstimateErrorInterval(Rational errorIntervalSize)
            {
                if (RealPartEstimate.Min != null &&
                    (Rational)(RealPartEstimate.Max - RealPartEstimate.Min) <= errorIntervalSize)
                {
                    return;
                }
                Interval<Float> coefficientRealPartEstimate = Coefficient.GetRealPartEstimate(One);
                Float errorScale =
                    GetMagnitudeInterval(coefficientRealPartEstimate).Max + Float.One;
                Float termToSubtract = Float.One;
                foreach (Surd factor in Factors)
                {
                    Float realBound = GetMagnitudeInterval(factor.GetRealPartEstimate(One)).Max;
                    Float imaginaryBound =
                        GetMagnitudeInterval(factor.GetImaginaryPartEstimate(One)).Max;
                    if (realBound > imaginaryBound)
                    {
                        errorScale *= realBound + Float.One;
                        termToSubtract *= realBound;
                    }
                    else
                    {
                        errorScale *= imaginaryBound + Float.One;
                        termToSubtract *= imaginaryBound;
                    }
                }
                errorIntervalSize /= (Rational)(errorScale - termToSubtract);
                coefficientRealPartEstimate = Coefficient.GetRealPartEstimate(errorIntervalSize);
                List<List<Interval<Float>>> terms =
                    new List<List<Interval<Float>>> { new List<Interval<Float>>() };
                List<int> termPowerOfIFactors = new List<int> { 0 };
                foreach (Surd factor in Factors)
                {
                    int currentTermCount = terms.Count;
                    for (int i = 0; i < currentTermCount; ++i)
                    {
                        List<Interval<Float>> term = new List<Interval<Float>>(terms[i]);
                        terms[i].Add(Factors[i].GetRealPartEstimate(errorIntervalSize));
                        term.Add(Factors[i].GetImaginaryPartEstimate(errorIntervalSize));
                        terms.Add(term);
                        termPowerOfIFactors.Add(termPowerOfIFactors[i] + 1);
                    }
                }
                Interval<Float> realEstimate = new Interval<Float>(Float.Zero, Float.Zero);
                Interval<Float> imaginaryEstimate = new Interval<Float>(Float.Zero, Float.Zero);
                for (int i = 0; i < terms.Count; ++i)
                {
                    List<Float> boundCandidates = new List<Float> {
                        coefficientRealPartEstimate.Min, coefficientRealPartEstimate.Max };
                    for (int j = 0; j < terms[i].Count; ++j)
                    {
                        int currentBoundCandidateCount = boundCandidates.Count;
                        for (int k = 0; k < currentBoundCandidateCount; ++k)
                        {
                            Float boundCandidate = boundCandidates[k];
                            boundCandidates[k] *= terms[i][j].Min;
                            boundCandidate *= terms[i][j].Max;
                            boundCandidates.Add(boundCandidate);
                        }
                    }
                    Interval<Float> termEstimate =
                        new Interval<Float>(boundCandidates[0], boundCandidates[0]);
                    for (int j = 1; j < boundCandidates.Count; ++j)
                    {
                        if (boundCandidates[j] < termEstimate.Min)
                        {
                            termEstimate.Min = boundCandidates[j];
                        }
                        else if (boundCandidates[j] > termEstimate.Max)
                        {
                            termEstimate.Max = boundCandidates[j];
                        }
                    }
                    switch (termPowerOfIFactors[i] % 4)
                    {
                        case 0:
                            realEstimate.Min += termEstimate.Min;
                            realEstimate.Max += termEstimate.Max;
                            break;
                        case 1:
                            imaginaryEstimate.Min += termEstimate.Min;
                            imaginaryEstimate.Max += termEstimate.Max;
                            break;
                        case 2:
                            realEstimate.Min -= termEstimate.Min;
                            realEstimate.Max -= termEstimate.Max;
                            break;
                        case 3:
                            imaginaryEstimate.Min -= termEstimate.Min;
                            imaginaryEstimate.Max -= termEstimate.Max;
                            break;
                    }
                }
                RealPartEstimate = realEstimate;
                ImaginaryPartEstimate = imaginaryEstimate;
            }
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                RefineRectangularEstimateErrorInterval(errorIntervalSize);
                return RealPartEstimate;
            }
            public override Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize)
            {
                RefineRectangularEstimateErrorInterval(errorIntervalSize);
                return ImaginaryPartEstimate;
            }            
            public override Interval<Float> GetArgumentEstimate(Rational errorIntervalSize)
            {
                if (ArgumentEstimate.Min == null ||
                    (Rational)(ArgumentEstimate.Max - ArgumentEstimate.Min) > errorIntervalSize)
                {
                    List<Term> terms = new List<Term>(Factors);
                    terms.Add(Coefficient);
                    ArgumentEstimate =
                        EstimatePartSum(terms, errorIntervalSize, GetArgumentEstimate);
                }
                return ArgumentEstimate;
            }            
            public override Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize)
            {
                if (MagnitudeEstimate.Min == null ||
                    (Rational)(MagnitudeEstimate.Max - MagnitudeEstimate.Min) > errorIntervalSize)
                {
                    Float errorScale = Coefficient.GetMagnitudeEstimate(One).Max + Float.One;
                    foreach (Surd factor in Factors)
                    {
                        errorScale *= factor.GetMagnitudeEstimate(One).Max + Float.One;
                    }
                    errorIntervalSize /= (Rational)errorScale;
                    Interval<Float> productMagnitude =
                        Coefficient.GetMagnitudeEstimate(errorIntervalSize);
                    foreach (Surd factor in Factors)
                    {
                        Interval<Float> factorMagnitudeEstimate =
                            factor.GetMagnitudeEstimate(errorIntervalSize);
                        productMagnitude.Min *= factorMagnitudeEstimate.Min;
                        productMagnitude.Max *= factorMagnitudeEstimate.Max;
                    }
                    MagnitudeEstimate = productMagnitude;
                }
                return MagnitudeEstimate;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder(Factors[0].ToString());
                for (int i = 1; i < Factors.Count; ++i)
                {
                    output.Append('*' + Factors[i].ToString());
                }
                return Coefficient.InsertString(output.ToString());
            }
        }
        public class PolynomialIndependentSum : Primitive
        {
            public List<Term> Terms { get; }
            Interval<Float> RealPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ImaginaryPartEstimate = new Interval<Float>(null, null);
            public PolynomialIndependentSum(List<Term> terms)
            {
                Terms = terms;
                Terms.Sort();
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                MultivariatePolynomial variableForm;
                RationalPolynomial[] minimalPolynomials;
                if (Terms[0] is Rational rational)
                {
                    minimalPolynomials = new RationalPolynomial[Terms.Count - 1];
                    for (int i = 1; i < Terms.Count; ++i)
                    {
                        RationalPolynomial polynomial = Terms[i].MinimalPolynomial;
                        minimalPolynomials[i - 1] = polynomial;
                    }
                    variableForm = new MultivariatePolynomial(minimalPolynomials);
                    variableForm.SetCoefficient(new int[minimalPolynomials.Length], rational);
                }
                else
                {
                    minimalPolynomials = new RationalPolynomial[Terms.Count];
                    for (int i = 0; i < Terms.Count; ++i)
                    {
                        RationalPolynomial polynomial = Terms[i].MinimalPolynomial;
                        minimalPolynomials[i] = polynomial;
                    }
                    variableForm = new MultivariatePolynomial(minimalPolynomials);
                }
                for (int i = 0; i < variableForm.MinimalPolynomials.Length; ++i)
                {
                    int[] indices = new int[variableForm.MinimalPolynomials.Length];
                    indices[i] = 1;
                    variableForm.SetCoefficient(indices, Number.One);
                }
                return variableForm.CalculateValueMinimalPolynomial(this);
            }
            public override List<Number> GetConjugates()
            {
                return GetSumConjugates(Terms);
            }
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                if (RealPartEstimate.Min == null ||
                    (Rational)(RealPartEstimate.Max - RealPartEstimate.Min) > errorIntervalSize)
                {
                    RealPartEstimate = EstimateRealPartSum(Terms, errorIntervalSize);
                }
                return RealPartEstimate;
            }
            public override Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize)
            {
                if (ImaginaryPartEstimate.Min == null || (Rational)(ImaginaryPartEstimate.Max -
                    ImaginaryPartEstimate.Min) > errorIntervalSize)
                {
                    ImaginaryPartEstimate = EstimateImaginaryPartSum(Terms, errorIntervalSize);
                }
                return ImaginaryPartEstimate;
            }
        }
        static Interval<Float> GetMagnitudeInterval(Interval<Float> interval)
        {
            Float magnitudeA = interval.Min.Magnitude();
            Float magnitudeB = interval.Max.Magnitude();
            if (magnitudeA < magnitudeB)
            {
                return new Interval<Float>(magnitudeA, magnitudeB);
            }
            return new Interval<Float>(magnitudeB, magnitudeA);
        }
        class LinearlyIndependentSum : Number
        {
            public List<Term> Terms { get; }
            public List<RationalPolynomial> TermsInTermsOfPrimitiveElement { get; }
            Interval<Float> RealPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ImaginaryPartEstimate = new Interval<Float>(null, null);
            Interval<Float> ArgumentEstimate = new Interval<Float>(null, null);
            Interval<Float> MagnitudeEstimate = new Interval<Float>(null, null);
            Term PrimitiveCandidate;
            Primitive PrimitiveElement_ = null;
            public Primitive PrimitiveElement
            {
                get
                {
                    if (PrimitiveElement_ != null)
                    {
                        return PrimitiveElement_;
                    }
                    bool candidateIsPrimitive(int candidateIndex, int otherTermIndex)
                    {
                        RationalPolynomial otherTermInTermsOfCandidate =
                            Terms[candidateIndex].ArgumentInTermsOfThis(Terms[otherTermIndex]);
                        if (otherTermInTermsOfCandidate != null)
                        {
                            TermsInTermsOfPrimitiveElement[candidateIndex] =
                                new RationalPolynomial(null, Zero, One);
                            TermsInTermsOfPrimitiveElement[otherTermIndex] =
                                otherTermInTermsOfCandidate;
                            PrimitiveElement_ = Terms[candidateIndex];
                            return true;
                        }
                        return false;
                    }
                    if (Terms[0] == PrimitiveCandidate && candidateIsPrimitive(0, 1))
                    {
                        return PrimitiveElement_;
                    }
                    else if (candidateIsPrimitive(1, 0))
                    {
                        return PrimitiveElement_;
                    }
                    RationalPolynomial x = new RationalPolynomial(null, Zero, One);
                    Integer k = One;
                    bool primitiveElementFound(int termIndexA, int termIndexB)
                    {
                        PrimitiveElement_ = new PolynomialIndependentSum(
                            new List<Term> { k * Terms[termIndexA], Terms[termIndexB] });
                        RationalPolynomial termAInTermsOfPrimitive =
                            PrimitiveElement_.ArgumentInTermsOfThis(Terms[termIndexA]);
                        if (termAInTermsOfPrimitive != null)
                        {
                            TermsInTermsOfPrimitiveElement[termIndexA] = termAInTermsOfPrimitive;
                            TermsInTermsOfPrimitiveElement[termIndexB] =
                                x - k * termAInTermsOfPrimitive;
                            return true;
                        }
                        return false;
                    }
                    while (true)
                    {
                        if (primitiveElementFound(0, 1))
                        {
                            return PrimitiveElement_;
                        }
                        if (primitiveElementFound(1, 0))
                        {
                            return PrimitiveElement_;
                        }
                        ++k;
                    }
                }
            }
            public LinearlyIndependentSum(List<Term> terms,
                List<RationalPolynomial> termsInTermsOfPrimitiveElement, Primitive primitiveElement)
            {
                List<Tuple<Term, RationalPolynomial>> termsWithPrimitiveElementRepresentations =
                    new List<Tuple<Term, RationalPolynomial>>();
                for (int i = 0; i < terms.Count; ++i)
                {
                    termsWithPrimitiveElementRepresentations.Add(new Tuple<Term,
                        RationalPolynomial>(terms[i], termsInTermsOfPrimitiveElement[i]));
                }
                termsWithPrimitiveElementRepresentations.Sort(delegate (
                    Tuple<Term, RationalPolynomial> a, Tuple<Term, RationalPolynomial> b)
                    { return a.Item1.CompareTo(b.Item1); });
                Terms = new List<Term>();
                TermsInTermsOfPrimitiveElement = new List<RationalPolynomial>();
                for (int i = 0; i < termsWithPrimitiveElementRepresentations.Count; ++i)
                {
                    Terms.Add(termsWithPrimitiveElementRepresentations[i].Item1);
                    TermsInTermsOfPrimitiveElement.Add(
                        termsWithPrimitiveElementRepresentations[i].Item2);
                }
                PrimitiveElement_ = primitiveElement;
            }
            public LinearlyIndependentSum(Term primitiveCandidateTerm, Term otherTerm)
            {
                if (primitiveCandidateTerm.CompareTo(otherTerm) < 0)
                {
                    Terms = new List<Term> { primitiveCandidateTerm, otherTerm };
                }
                else
                {
                    Terms = new List<Term> { otherTerm, primitiveCandidateTerm };
                }
                TermsInTermsOfPrimitiveElement = new List<RationalPolynomial> { null, null };
                PrimitiveCandidate = primitiveCandidateTerm;
            }
            public override Number Plus(Number a)
            {
                if (a is Term term)
                {
                    for (int i = 0; i < Terms.Count; ++i)
                    {
                        Term termSum = term.CheapPlus(Terms[i]);
                        if (!ReferenceEquals(termSum, null))
                        {
                            List<Term> outputTerms = new List<Term>(Terms);
                            if (termSum == Zero && Terms.Count == 2)
                            {                                
                                outputTerms.RemoveAt(i);
                                return outputTerms[0];
                            }
                            outputTerms[i] = termSum;
                            List<RationalPolynomial> outputTermsInTermsOfPrimitiveElement =
                                new List<RationalPolynomial>(TermsInTermsOfPrimitiveElement);
                            outputTermsInTermsOfPrimitiveElement[i] *=
                                termSum.RationalFactor() / Terms[i].RationalFactor();
                            return new LinearlyIndependentSum(outputTerms,
                                outputTermsInTermsOfPrimitiveElement, PrimitiveElement);
                        }
                    }
                    Number thisPlusTerm(Primitive newPrimitiveElement,
                        List<RationalPolynomial> oldTermsInTermsOfNewPrimitiveElement,
                        RationalPolynomial newTermInTermsOfNewPrimitiveElement)
                    {
                        Matrix matrix =
                            new Matrix(oldTermsInTermsOfNewPrimitiveElement).GetTranspose();
                        List<Term> outputTerms;
                        if (newTermInTermsOfNewPrimitiveElement.Coefficients.Count >
                            matrix.Rows.Count)
                        {
                            outputTerms = new List<Term>(Terms);
                            outputTerms.Add(term);
                            oldTermsInTermsOfNewPrimitiveElement.Add(
                                newTermInTermsOfNewPrimitiveElement);
                            return new LinearlyIndependentSum(outputTerms,
                                oldTermsInTermsOfNewPrimitiveElement, newPrimitiveElement);
                        }
                        List<Rational> augmentation =
                            new List<Rational>(newTermInTermsOfNewPrimitiveElement.Coefficients);
                        while (augmentation.Count < matrix.Rows.Count)
                        {
                            augmentation.Add(Zero);
                        }
                        matrix.GetRowEchelonForm(augmentation);
                        for (int i = augmentation.Count - 1; i >= Terms.Count; --i)
                        {
                            if (augmentation[i] != Zero)
                            {
                                outputTerms = new List<Term>(Terms);
                                outputTerms.Add(term);
                                oldTermsInTermsOfNewPrimitiveElement.Add(
                                    newTermInTermsOfNewPrimitiveElement);
                                return new LinearlyIndependentSum(outputTerms,
                                    oldTermsInTermsOfNewPrimitiveElement, newPrimitiveElement);
                            }
                            else
                            {
                                augmentation.RemoveAt(i);
                                matrix.Rows.RemoveAt(i);
                            }
                        }
                        matrix.Diagonalize(augmentation);
                        Integer minusOne = -One;
                        outputTerms = new List<Term>();
                        List<RationalPolynomial> outputTermsInTermsOfNewPrimitiveElement =
                            new List<RationalPolynomial>();
                        for (int i = 0; i < augmentation.Count; ++i)
                        {
                            if (augmentation[i] != minusOne)
                            {
                                Rational scale = One + augmentation[i];
                                Term outputTerm = scale * Terms[i];
                                outputTerms.Add(outputTerm);
                                outputTermsInTermsOfNewPrimitiveElement.Add(
                                    scale * oldTermsInTermsOfNewPrimitiveElement[i]);
                            }
                        }
                        if (outputTerms.Count == 0)
                        {
                            return Zero;
                        }
                        if (outputTerms.Count == 1)
                        {
                            return outputTerms[0];
                        }
                        return new LinearlyIndependentSum(outputTerms,
                            outputTermsInTermsOfNewPrimitiveElement, newPrimitiveElement);
                    }
                    RationalPolynomial termInTermsOfPrimitiveElement =
                        PrimitiveElement.ArgumentInTermsOfThis(term);
                    if (termInTermsOfPrimitiveElement != null)
                    {
                        return thisPlusTerm(PrimitiveElement, TermsInTermsOfPrimitiveElement,
                            termInTermsOfPrimitiveElement);
                    }
                    List<RationalPolynomial> convertTermsToNewPrimitiveElement(
                        RationalPolynomial newPrimitiveElementInTermsOfOld)
                    {
                        List<RationalPolynomial> termsInTermsOfNewPrimitiveElement =
                            new List<RationalPolynomial>();
                        List<RationalPolynomial> powers = new List<RationalPolynomial> {
                            newPrimitiveElementInTermsOfOld.GetMultiplicativeIdentity() };
                        for (int i = 0; i < Terms.Count; ++i)
                        {
                            RationalPolynomial termInTermsOfNewPrimitiveElement =
                                newPrimitiveElementInTermsOfOld.GetAdditiveIdentity();
                            for (int j = 0;
                                j < TermsInTermsOfPrimitiveElement[i].Coefficients.Count; ++j)
                            {
                                if (j >= powers.Count)
                                {
                                    powers.Add(powers[powers.Count - 1] *
                                        newPrimitiveElementInTermsOfOld);
                                }
                                termInTermsOfNewPrimitiveElement += powers[j] *
                                    TermsInTermsOfPrimitiveElement[i].Coefficients[j];
                            }
                            termsInTermsOfNewPrimitiveElement.Add(termInTermsOfNewPrimitiveElement);
                        }
                        return termsInTermsOfNewPrimitiveElement;
                    }
                    RationalPolynomial primitiveElementInTermsOfTerm =
                        term.ArgumentInTermsOfThis(PrimitiveElement);
                    if (primitiveElementInTermsOfTerm != null)
                    {
                        return thisPlusTerm(term,
                            convertTermsToNewPrimitiveElement(primitiveElementInTermsOfTerm),
                            new RationalPolynomial(term.MinimalPolynomial, Zero, One));
                    }
                    Integer k = One;
                    while (true)
                    {
                        List<Term> primitiveElementTerms = new List<Term>(Terms);
                        primitiveElementTerms.Add(k * term);
                        PrimitiveElement_ = new PolynomialIndependentSum(primitiveElementTerms);
                        RationalPolynomial termInTermsOfPrimitive =
                            PrimitiveElement_.ArgumentInTermsOfThis(term);
                        if (termInTermsOfPrimitive != null)
                        {
                            return thisPlusTerm(PrimitiveElement_,
                                convertTermsToNewPrimitiveElement(new RationalPolynomial(
                                PrimitiveElement_.MinimalPolynomial, Zero, One) -
                                k * termInTermsOfPrimitive), termInTermsOfPrimitive);
                        }
                        ++k;
                    }
                }
                LinearlyIndependentSum sum = (LinearlyIndependentSum)a;
                Number output = this;
                foreach (Term t in sum.Terms)
                {
                    output += t;
                }
                return output;
            }
            public override Number Times(Number a)
            {
                if (a == Zero)
                {
                    return Zero;
                }
                if (a == One)
                {
                    return this;
                }
                Number output = Zero;
                if (a is Term)
                {
                    foreach (Term term in Terms)
                    {
                        output += term * a;
                    }
                    return output;
                }
                LinearlyIndependentSum sum = (LinearlyIndependentSum)a;
                foreach (Term multiplicandTerm in Terms)
                {
                    foreach (Term multiplierTerm in sum.Terms)
                    {
                        output += multiplicandTerm * multiplierTerm;
                    }
                }
                return output;
            }
            public override Number Reciprocal()
            {
                List<Number> conjugates = GetConjugates();
                conjugates.Remove(this);
                Number output = One;
                foreach (Number conjugate in conjugates)
                {
                    output *= conjugate;
                }
                if (MinimalPolynomial.Coefficients.Count % 2 == 0)
                {
                    return -output / MinimalPolynomial.Coefficients[0];
                }
                return output / MinimalPolynomial.Coefficients[0];
            }
            public override Rational RationalFactor()
            {
                Rational termRationalFactor = Terms[0].RationalFactor();
                Integer GCD = termRationalFactor.Numerator;
                Integer LCM = termRationalFactor.Denominator;
                for (int i = 1; i < Terms.Count; ++i)
                {
                    termRationalFactor = Terms[i].RationalFactor();
                    GCD = Integer.GetGCD(GCD, termRationalFactor.Numerator);
                    LCM = Integer.GetLCM(LCM, termRationalFactor.Denominator);
                }
                return GCD / LCM;
            }
            public override Number Exponentiate(Rational exponent)
            {
                if (exponent == One)
                {
                    return this;
                }
                if (exponent.Numerator != One)
                {
                    return Exponentiate<Number>(this, exponent.Numerator).Exponentiate(
                        One / exponent.Denominator);
                }
                Rational rationalFactor = RationalFactor();
                Number coefficientExponentiation = rationalFactor.Exponentiate(exponent);
                if (coefficientExponentiation is Surd surd)
                {
                    return new Surd(this, exponent.Denominator);
                }
                return coefficientExponentiation *
                    new Surd(this / rationalFactor, exponent.Denominator);
            }
            protected override int GetTypeIndex()
            {
                return 5;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                {
                    return comparison;
                }
                LinearlyIndependentSum expression = (LinearlyIndependentSum)obj;
                comparison = Terms.Count - expression.Terms.Count;
                if (comparison != 0)
                {
                    return comparison;
                }
                for (int i = 0; i < Terms.Count; ++i)
                {
                    comparison = Terms[i].CompareTo(expression.Terms[i]);
                    if (comparison != 0)
                    {
                        return comparison;
                    }
                }
                return 0;
            }
            public override int GetHashCode()
            {
                int output = 0;
                foreach (Term term in Terms)
                {
                    output ^= term.GetHashCode();
                }
                return output;
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                List<Rational> variableFormCoefficients = new List<Rational>();
                Primitive ensurePrimitiveElementHasBeenFound = PrimitiveElement;
                foreach (RationalPolynomial polynomial in TermsInTermsOfPrimitiveElement)
                {
                    variableFormCoefficients =
                        CoefficientAdd(variableFormCoefficients, polynomial.Coefficients);
                }
                return new RationalPolynomial(variableFormCoefficients,
                    PrimitiveElement.MinimalPolynomial).GetValueMinimalPolynomial();
            }
            public override List<Number> GetConjugates()
            {
                return GetSumConjugates(Terms);
            }
            Interval<Rational> GetRationalArgumentEstimate(Rational errorIntervalSize)
            {
                Number thisTimesI = this * new Surd(-One, Two);
                if (HasZeroImaginaryPart())
                {
                    if (thisTimesI.HasZeroImaginaryPart())
                    {
                        return new Interval<Rational>();
                    }

                    //HasZeroImaginaryPart will already have refined thisTimesI's imaginary part
                    //estimate to a size of One or smaller, so passing One to this call of
                    //GetImaginaryEstimate will retrieve that estimate, whose bounds
                    //HasZeroImaginaryPart ensures will fall on the same side of Zero.
                    if (thisTimesI.GetImaginaryPartEstimate(One).Min > Float.Zero)
                    {
                        return new Interval<Rational>(Zero, Zero);
                    }
                    Pi.RefineErrorInterval(errorIntervalSize);
                    return new Interval<Rational>(Pi.LowEstimate, Pi.HighEstimate);
                }
                if (thisTimesI.HasZeroImaginaryPart())
                {
                    if (thisTimesI.GetRealPartEstimate(One).Max < Float.Zero)
                    {
                        Pi.RefineErrorInterval(Two * errorIntervalSize);
                        return new Interval<Rational>(Pi.LowEstimate / Two, Pi.HighEstimate / Two);
                    }
                    Rational threeOverTwo = new Fraction(new Integer(3), Two);
                    Pi.RefineErrorInterval(errorIntervalSize / threeOverTwo);
                    return new Interval<Rational>(threeOverTwo * Pi.LowEstimate,
                        threeOverTwo * Pi.HighEstimate);
                }
                Interval<Rational> realMagnitudeInterval =
                    GetMagnitudeInterval(GetRealPartEstimate(One)).ToRational();
                Interval<Rational> imaginaryMagnitudeInterval =
                    GetMagnitudeInterval(GetImaginaryPartEstimate(One)).ToRational();
                errorIntervalSize /= new Integer(3) *
                    (realMagnitudeInterval.Max + imaginaryMagnitudeInterval.Max + One);
                Interval<Rational> realPartEstimate =
                    GetRealPartEstimate(Min(errorIntervalSize, errorIntervalSize *
                    realMagnitudeInterval.Min * realMagnitudeInterval.Min)).ToRational();
                Interval<Rational> imaginaryPartEstimate =
                    GetImaginaryPartEstimate(Min(errorIntervalSize, errorIntervalSize * 
                    imaginaryMagnitudeInterval.Min * imaginaryMagnitudeInterval.Min)).ToRational();
                Interval<Rational> argumentEstimateA = EstimateAtan2(imaginaryPartEstimate.Max,
                    realPartEstimate.Min, errorIntervalSize);
                Interval<Rational> argumentEstimateB = EstimateAtan2(imaginaryPartEstimate.Min,
                    realPartEstimate.Max, errorIntervalSize);
                return new Interval<Rational>(Min(argumentEstimateA.Min, argumentEstimateB.Min),
                    Max(argumentEstimateA.Max, argumentEstimateB.Max));
            }
            public override Interval<Float> GetRealPartEstimate(Rational errorIntervalSize)
            {
                if (RealPartEstimate.Min == null ||
                    (Rational)(RealPartEstimate.Max - RealPartEstimate.Min) > errorIntervalSize)
                {
                    RealPartEstimate = EstimateRealPartSum(Terms, errorIntervalSize);
                }
                return RealPartEstimate;
            }
            public override Interval<Float> GetImaginaryPartEstimate(Rational errorIntervalSize)
            {
                if (ImaginaryPartEstimate.Min == null || (Rational)(ImaginaryPartEstimate.Max -
                    ImaginaryPartEstimate.Min) > errorIntervalSize)
                {
                    ImaginaryPartEstimate = EstimateImaginaryPartSum(Terms, errorIntervalSize);
                }
                return ImaginaryPartEstimate;
            }
            public override Interval<Float> GetMagnitudeEstimate(Rational errorIntervalSize)
            {
                if (MagnitudeEstimate.Min == null || (Rational)(MagnitudeEstimate.Max -
                    MagnitudeEstimate.Min) > errorIntervalSize)
                {
                    errorIntervalSize /= new Integer(3);
                    Interval<Float> getPartInterval(EstimateGetter<Number> getPartEstimate)
                    {
                        return GetMagnitudeInterval(getPartEstimate(this, errorIntervalSize /
                            (new Integer(4) * (Rational)GetMagnitudeInterval(
                            GetRealPartEstimate(One)).Max + Two)));
                    }
                    Interval<Float> realPartInterval = getPartInterval(GetRealPartEstimate);
                    Interval<Float> imaginaryPartInterval =
                        getPartInterval(GetImaginaryPartEstimate);
                    Interval<Float> getBoundInterval(Float realPart, Float imaginaryPart)
                    {
                        return new Surd(
                            (Rational)(realPart * realPart + imaginaryPart * imaginaryPart), Two).
                            GetRealPartEstimate(errorIntervalSize);
                    }
                    MagnitudeEstimate = new Interval<Float>(
                        getBoundInterval(realPartInterval.Min, imaginaryPartInterval.Min).Min,
                        getBoundInterval(realPartInterval.Max, imaginaryPartInterval.Max).Max);
                }
                return MagnitudeEstimate;
            }            
            public override Interval<Float> GetArgumentEstimate(Rational errorIntervalSize)
            {
                if (ArgumentEstimate.Min == null || (Rational)(ArgumentEstimate.Max -
                    ArgumentEstimate.Min) > errorIntervalSize)
                {
                    ArgumentEstimate =
                        GetRefinedErrorInterval(errorIntervalSize, GetRationalArgumentEstimate);
                }
                return ArgumentEstimate;
            }
            public override string ToString()
            {
                StringBuilder output = new StringBuilder(Terms[0].ToString());
                for (int i = 1; i < Terms.Count; ++i)
                {
                    string termString = Terms[i].ToString();
                    if (termString[0] != '-')
                    {
                        output.Append('+');
                    }
                    output.Append(termString);
                }
                return output.ToString();
            }
        }
        class Matrix
        {
            public List<List<Rational>> Rows { get; }
            public Matrix(List<List<Rational>> rows)
            {
                Rows = rows;
            }
            public Matrix(List<RationalPolynomial> rows)
            {
                Rows = new List<List<Rational>>();
                for (int i = 0; i < rows.Count; ++i)
                {
                    Rows.Add(new List<Rational>(rows[i].Coefficients));
                    for (int j = Rows[0].Count; j < Rows[i].Count; ++j)
                    {
                        for (int k = 0; k < i; ++k)
                        {
                            Rows[k].Add(Number.Zero);
                        }
                    }
                    while (Rows[i].Count < Rows[0].Count)
                    {
                        Rows[i].Add(Number.Zero);
                    }
                }
            }
            public Matrix GetTranspose()
            {
                List<List<Rational>> outputRows = new List<List<Rational>>();
                for (int i = 0; i < Rows[0].Count; ++i)
                {
                    List<Rational> row = new List<Rational>();
                    for (int j = 0; j < Rows.Count; ++j)
                    {
                        row.Add(Rows[j][i]);
                    }
                    outputRows.Add(row);
                }
                return new Matrix(outputRows);
            }
            public void GetRowEchelonForm<T>(List<T> augmentation)
                where T : IArithmetic<T>, IMultipliable<T, Rational>
            {
                int smallDimension;
                if (Rows.Count <= Rows[0].Count)
                {
                    smallDimension = Rows.Count;
                }
                else
                {
                    smallDimension = Rows[0].Count;
                }
                for (int i = 0; i < smallDimension; ++i)
                {
                    for (int j = i + 1; j < Rows.Count; ++j)
                    {
                        if (Rows[i][i] == Number.Zero)
                        {
                            List<Rational> tempRow = Rows[j];
                            Rows[j] = Rows[i];
                            Rows[i] = tempRow;
                            T tempAugmentationRow = augmentation[j];
                            augmentation[j] = augmentation[i];
                            augmentation[i] = tempAugmentationRow;
                        }
                        else
                        {
                            for (int k = i + 1; k < Rows.Count; ++k)
                            {
                                Rational scalar = Rows[k][i] / Rows[i][i];
                                for (int l = i; l < Rows[0].Count; ++l)
                                {
                                    Rows[k][l] -= Rows[i][l] * scalar;
                                }
                                augmentation[k] =
                                    augmentation[k].Minus(augmentation[i].Times(scalar));
                            }
                            break;
                        }
                    }
                }
            }
            public void Diagonalize(List<Rational> augmentation)
            {
                GetRowEchelonForm(augmentation);
                for (int i = Rows.Count - 1; i >= 0; --i)
                {
                    augmentation[i] /= Rows[i][i];
                    for (int j = 0; j <= i; ++j)
                    {
                        Rows[i][j] /= Rows[i][i];
                    }
                    for (int j = 0; j < i; ++j)
                    {
                        augmentation[j] -= augmentation[i] * Rows[j][i];
                        Rows[j][i] = Number.Zero;
                    }
                }
            }
        }
        static List<T> CoefficientAdd<T>(List<T> coefficientsA, List<T> coefficientsB)
            where T : IArithmetic<T>
        {
            if (coefficientsA.Count < coefficientsB.Count)
            {
                return CoefficientAdd(coefficientsB, coefficientsA);
            }
            List<T> output = new List<T>(coefficientsA);
            for (int i = 0; i < coefficientsB.Count; ++i)
            {
                output[i] = output[i].Plus(coefficientsB[i]);
            }
            return output;
        }
        static List<T> CoefficientMultiply<T>(List<T> coefficientsA,
            List<T> coefficientsB) where T : IArithmetic<T>
        {
            List<T> output = new List<T>();
            T instance;
            if (coefficientsA.Count > 0)
            {
                instance = coefficientsA[0];
            }
            else if (coefficientsB.Count > 0)
            {
                instance = coefficientsB[0];
            }
            else
            {
                return output;
            }
            for (int i = 0; i < coefficientsA.Count + coefficientsB.Count - 1; ++i)
            {
                output.Add(instance.GetAdditiveIdentity());
            }
            for (int i = 0; i < coefficientsA.Count; ++i)
            {
                for (int j = 0; j < coefficientsB.Count; ++j)
                {
                    output[i + j] = output[i + j].Plus(coefficientsA[i].Times(coefficientsB[j]));
                }
            }
            return output;
        }
        static List<Rational> GetMonicCoefficients<T>(List<T> coefficients) where T : Rational
        {
            List<Rational> monicCoefficients = new List<Rational>();
            foreach (Rational coefficient in coefficients)
            {
                monicCoefficients.Add(coefficient / coefficients[coefficients.Count - 1]);
            }
            return monicCoefficients;
        }
#if DEBUG
        static string CoefficientsToString<T>(List<T> coefficients) where T : Rational
        {
            if (coefficients.Count == 0)
            {
                return "0";
            }
            StringBuilder output = new StringBuilder();
            if (coefficients[0] != Number.Zero)
            {
                output.Append(coefficients[0]);
            }
            if (coefficients.Count > 1)
            {
                if (coefficients[1] < Number.Zero)
                {
                    output.Append(coefficients[1].InsertString("x"));
                }
                else if (coefficients[1] > Number.Zero)
                {
                    if (output.Length != 0)
                    {
                        output.Append("+");
                    }
                    output.Append(coefficients[1].InsertString("x"));
                }
                for (int i = 2; i < coefficients.Count; ++i)
                {
                    if (coefficients[i] < Number.Zero)
                    {
                        output.Append(coefficients[i].InsertString("x^" + i));
                    }
                    else if (coefficients[i] > Number.Zero)
                    {
                        if (output.Length != 0)
                        {
                            output.Append("+");
                        }
                        output.Append(coefficients[i].InsertString("x^" + i));
                    }
                }
            }
            return output.ToString();
        }
#endif
        public class IntegerPolynomial : IArithmetic<IntegerPolynomial>
        {
            //The index of each element corrosponds to the degree of the term of which it is the
            //coefficient.
            public List<Integer> Coefficients { get; }

            //See RationalPolynomial.MinimalPolynomial. The results of operations in the non-null
            //case satisfy the conditions described there only up to a scalar multiple, because
            //otherwise they would not in general have integer coefficients. The non-null case is
            //currently only used in RationalPolynomial.GetGCD() (via GetPrimitivePart()) where the
            //constant multiple isn't an issue.
            RationalPolynomial MinimalPolynomial;
            public IntegerPolynomial(RationalPolynomial minimalPolynomial,
                List<Integer> coefficients)
            {//Can mutate the coefficients parameter.
                while (coefficients.Count != 0 &&
                    coefficients[coefficients.Count - 1] == Number.Zero)
                {
                    coefficients.RemoveAt(coefficients.Count - 1);
                }
                Coefficients = coefficients;
                MinimalPolynomial = minimalPolynomial;
            }
            public IntegerPolynomial(RationalPolynomial minimalPolynomial,
                params Integer[] coefficients)
            {
                Coefficients = new List<Integer>(coefficients);
                while (Coefficients.Count != 0 &&
                    coefficients[Coefficients.Count - 1] == Number.Zero)
                {
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                }
                MinimalPolynomial = minimalPolynomial;
            }
            public IntegerPolynomial GetAdditiveIdentity()
            {
                return new IntegerPolynomial(MinimalPolynomial);
            }
            public IntegerPolynomial GetMultiplicativeIdentity()
            {
                return new IntegerPolynomial(MinimalPolynomial, Number.One);
            }
            public IntegerPolynomial GetNegative()
            {
                List<Integer> negativeCoefficients = new List<Integer>();
                foreach (Integer coefficient in Coefficients)
                {
                    negativeCoefficients.Add(-coefficient);
                }
                return new IntegerPolynomial(MinimalPolynomial, negativeCoefficients);
            }
            public IntegerPolynomial Plus(IntegerPolynomial a)
            {
                return new IntegerPolynomial(MinimalPolynomial,
                    CoefficientAdd(Coefficients, a.Coefficients));
            }
            public IntegerPolynomial Minus(IntegerPolynomial a)
            {
                return Plus(a.GetNegative());
            }
            public IntegerPolynomial Times(IntegerPolynomial a)
            {
                if (MinimalPolynomial == null)
                {
                    return new IntegerPolynomial(MinimalPolynomial,
                        CoefficientMultiply(Coefficients, a.Coefficients));
                }
                List<Rational> coefficients = CoefficientMultiply(
                    ((RationalPolynomial)this).Coefficients, ((RationalPolynomial)a).Coefficients);
                for (int i = coefficients.Count; i >= MinimalPolynomial.Coefficients.Count; --i)
                {
                    for (int j = 0; j < MinimalPolynomial.Coefficients.Count - 1; ++j)
                    {
                        coefficients[i - MinimalPolynomial.Coefficients.Count + j] -=
                            MinimalPolynomial.Coefficients[j];
                    }
                    coefficients.RemoveAt(coefficients.Count - 1);
                }
                return new RationalPolynomial(coefficients, MinimalPolynomial).GetPrimitivePart();
            }
            public Division<IntegerPolynomial> EuclideanDivideBy(IntegerPolynomial divisor)
            {
                Division<IntegerPolynomial> division;
                List<Integer> quotient = new List<Integer>();
                List<Integer> remainder = new List<Integer>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                {
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        Rational quotientCoefficient = remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1];
                        if (quotientCoefficient.Denominator == Number.One)
                        {
                            quotient.Insert(0, quotientCoefficient.Numerator);
                            remainder.RemoveAt(remainder.Count - 1);
                            for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            {
                                remainder[remainder.Count - j] -= quotient[0] *
                                    divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                            }
                        }
                        else
                        {
                            division.Quotient = null;
                            division.Remainder = null;
                            return division;
                        }
                    }
                }
                division.Quotient = new IntegerPolynomial(MinimalPolynomial, quotient);
                division.Remainder = new IntegerPolynomial(MinimalPolynomial, remainder);
                return division;
            }
            public bool Equals(IntegerPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                {
                    return false;
                }
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    if (Coefficients[i] != a.Coefficients[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            public static IntegerPolynomial operator +(IntegerPolynomial a, IntegerPolynomial b)
            {
                return a.Plus(b);
            }
            public static IntegerPolynomial operator -(IntegerPolynomial a, IntegerPolynomial b)
            {
                return a.Minus(b);
            }
            public static IntegerPolynomial operator *(IntegerPolynomial a, IntegerPolynomial b)
            {
                return a.Times(b);
            }
            public static IntegerPolynomial operator *(IntegerPolynomial a, Integer b)
            {
                return a.Times(new IntegerPolynomial(a.MinimalPolynomial, b));
            }
            public static IntegerPolynomial operator *(Integer a, IntegerPolynomial b)
            {
                return b * a;
            }
            public static IntegerPolynomial operator /(IntegerPolynomial a, IntegerPolynomial b)
            {
                return a.EuclideanDivideBy(b).Quotient;
            }
            public static IntegerPolynomial operator /(IntegerPolynomial a, Integer b)
            {
                List<Integer> quotientCoefficients = new List<Integer>();
                foreach (Integer coefficient in a.Coefficients)
                {
                    Division<Integer> division = coefficient.EuclideanDivideBy(b);
                    Debug.Assert(division.Remainder == Number.Zero,
                        "An element of a.Coefficients was not divisible by b.");
                    quotientCoefficients.Add(division.Quotient);
                }
                return new IntegerPolynomial(a.MinimalPolynomial, quotientCoefficients);
            }
            public static IntegerPolynomial operator %(IntegerPolynomial a, IntegerPolynomial b)
            {
                return a.EuclideanDivideBy(b).Remainder;
            }
            public IntegerPolynomial GetDerivative()
            {
                List<Integer> coefficients = new List<Integer>();
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    coefficients.Add(Coefficients[i] * new Integer(i));
                }
                return new IntegerPolynomial(null, coefficients);
            }
            public IntegerPolynomial GetPrimitivePart()
            {
                if (Coefficients.Count == 0)
                {
                    return this;
                }
                return this / Solver.GetGCD(Coefficients);
            }
            public static IntegerPolynomial GetGCD(IntegerPolynomial a, IntegerPolynomial b)
            {
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    if (a.Coefficients.Count == 0)
                    {
                        return b;
                    }
                    IntegerPolynomial t = a;
                    a = b;
                    b = t;
                }
                else if (b.Coefficients.Count == 0)
                {
                    return a;
                }
                Integer aContent = Solver.GetGCD(a.Coefficients);
                Integer bContent = Solver.GetGCD(b.Coefficients);
                a /= aContent;
                b /= bContent;
                Integer d = Solver.GetGCD(aContent, bContent);
                Integer g = Number.One;
                Integer h = Number.One;
                Integer degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                IntegerPolynomial remainder = (a * Exponentiate(
                    b.Coefficients[b.Coefficients.Count - 1], degree + Number.One)) % b;
                while (remainder.Coefficients.Count > 1)
                {
                    a = b;
                    b = remainder / (g * Exponentiate(h, degree));
                    g = a.Coefficients[a.Coefficients.Count - 1];
                    Integer hExponent = Number.One - degree;
                    if (hExponent >= Number.Zero)
                    {
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    }
                    else
                    {
                        h = (Exponentiate(g, degree) / Exponentiate(h, -hExponent)).Numerator;
                    }
                    degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    remainder = (a * Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                        degree + Number.One)) % b;
                }
                if (remainder.Coefficients.Count == 1)
                {
                    return new IntegerPolynomial(a.MinimalPolynomial, d);
                }
                return b.GetPrimitivePart() * d;
            }
            public List<IntegerPolynomial> GetFactors()
            {
                List<IntegerPolynomial> squarefreeFactors = new List<IntegerPolynomial>();
                IntegerPolynomial derivative = GetDerivative();
                IntegerPolynomial a = GetGCD(this, derivative);
                IntegerPolynomial b = this / a;
                IntegerPolynomial c = derivative / a - b.GetDerivative();
                while (!(b.Coefficients.Count == 1 && b.Coefficients[0].Magnitude() == Number.One))
                {
                    a = GetGCD(b, c);
                    if (!squarefreeFactors.Contains(a))
                    {
                        squarefreeFactors.Add(a);
                    }
                    b = b / a;
                    c = c / a - b.GetDerivative();
                }
                List<IntegerPolynomial> splitFactor(IntegerPolynomial factor)
                {
                    ModdedPolynomial moddedFactor = null;
                    ModdedPolynomial GCD;
                    foreach (Integer prime in Primes())
                    {
                        moddedFactor = new ModdedPolynomial(factor.Coefficients, prime);
                        if (moddedFactor.Coefficients.Count == factor.Coefficients.Count)
                        {
                            GCD = Solver.GetGCD(moddedFactor, moddedFactor.GetDerivative());
                            if (GCD.Coefficients.Count == 1)
                            {
                                break;
                            }
                        }
                    }
                    moddedFactor = moddedFactor.GetMonic();
                    List<ModdedPolynomial> irreducibleModdedFactors =
                        moddedFactor.GetFactorsOfSquarefree();
                    Integer bound = Number.Zero;
                    foreach (Integer coefficient in factor.Coefficients)
                        bound += coefficient * coefficient;
                    Float floatSquareRoot = new Float(bound, Number.Zero).EstimateRoot(
                        Number.One, Number.Two).Max;
                    Rational rationalRoot = (Rational)floatSquareRoot;
                    Integer squareRoot = rationalRoot.Numerator.EuclideanDivideBy(
                        rationalRoot.Denominator).Quotient + Number.One;
                    Integer factorDegreeMinusOne = new Integer(factor.Coefficients.Count - 2);
                    Integer k = factorDegreeMinusOne.EuclideanDivideBy(Number.Two).Quotient;
                    bound = (squareRoot * factorDegreeMinusOne.ThisChooseK(k) +
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude() *
                        factorDegreeMinusOne.ThisChooseK(k - Number.One)) * Number.Two *
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude();
                    Integer e = Number.One;
                    Integer characteristicPower = moddedFactor.Characteristic;
                    while (characteristicPower < bound)
                    {
                        ++e;
                        characteristicPower *= moddedFactor.Characteristic;
                    }
                    List<IntegerPolynomial> liftedFactors = new List<IntegerPolynomial>();
                    ModdedPolynomial productOfUnliftedFactors = new ModdedPolynomial(
                        factor.Coefficients, moddedFactor.Characteristic).GetMonic();
                    IntegerPolynomial unmoddedATimesB =
                        new IntegerPolynomial(null, new ModdedPolynomial(factor.Coefficients,
                        characteristicPower).GetMonic().Coefficients);
                    irreducibleModdedFactors.RemoveAt(irreducibleModdedFactors.Count - 1);
                    foreach (ModdedPolynomial irreducibleModdedFactor in irreducibleModdedFactors)
                    {
                        productOfUnliftedFactors /= irreducibleModdedFactor;
                        Integer liftCharacteristic = moddedFactor.Characteristic;
                        IntegerPolynomial A =
                            new IntegerPolynomial(null, irreducibleModdedFactor.Coefficients);
                        IntegerPolynomial B =
                            new IntegerPolynomial(null, productOfUnliftedFactors.Coefficients);
                        while (liftCharacteristic < characteristicPower)
                        {
                            ModdedPolynomial AModCharacteristic =
                                new ModdedPolynomial(A.Coefficients, moddedFactor.Characteristic);
                            ModdedPolynomial BModCharacteristic =
                                new ModdedPolynomial(B.Coefficients, moddedFactor.Characteristic);
                            ExtendedGCDInfo<ModdedPolynomial> extendedGCD =
                                ExtendedGCD(AModCharacteristic, BModCharacteristic);
                            Integer reciprocal = extendedGCD.GCD.Coefficients[0].Reciprocal(
                                moddedFactor.Characteristic);
                            extendedGCD.ACoefficient *= reciprocal;
                            extendedGCD.BCoefficient *= reciprocal;
                            ModdedPolynomial f = new ModdedPolynomial(((unmoddedATimesB - A * B) /
                                liftCharacteristic).Coefficients, moddedFactor.Characteristic);
                            Division<ModdedPolynomial> division = (extendedGCD.BCoefficient *
                                f).EuclideanDivideBy(AModCharacteristic);
                            A += liftCharacteristic *
                                new IntegerPolynomial(null, division.Remainder.Coefficients);
                            B += liftCharacteristic *
                                new IntegerPolynomial(null, (extendedGCD.ACoefficient * f +
                                BModCharacteristic * division.Quotient).Coefficients);
                            liftCharacteristic *= moddedFactor.Characteristic;
                        }
                        liftedFactors.Add(A);
                        unmoddedATimesB = B;
                    }
                    liftedFactors.Add(unmoddedATimesB);
                    IntegerPolynomial factorTimesLeadingCoefficient =
                        factor * factor.Coefficients[factor.Coefficients.Count - 1];
                    List<IntegerPolynomial> factorSplit = new List<IntegerPolynomial>();
                    void BoundCoefficients(IntegerPolynomial polynomial)
                    {
                        for (int i = 0; i < polynomial.Coefficients.Count; ++i)
                        {
                            Integer remainder = polynomial.Coefficients[i].EuclideanDivideBy(
                                characteristicPower).Remainder;
                            if (remainder < Number.Zero)
                            {
                                remainder += characteristicPower;
                            }
                            if (remainder * Number.Two <= characteristicPower)
                            {
                                polynomial.Coefficients[i] = remainder;
                            }
                            else
                            {
                                polynomial.Coefficients[i] = remainder - characteristicPower;
                            }
                        }
                    }
                    bool testFactorCombination(int elementsToAdd, int indexOfPreviousElement,
                        List<int> combinationIndices)
                    {
                        if (elementsToAdd == 0)
                        {
                            IntegerPolynomial v = new IntegerPolynomial(null,
                                factor.Coefficients[factor.Coefficients.Count - 1]);
                            foreach (int index in combinationIndices)
                            {
                                v *= liftedFactors[index];
                            }
                            BoundCoefficients(v);
                            Division<IntegerPolynomial> division =
                                factorTimesLeadingCoefficient.EuclideanDivideBy(v);
                            if (division.Remainder != null &&
                                division.Remainder.Coefficients.Count == 0)
                            {
                                while (division.Remainder != null &&
                                    division.Remainder.Coefficients.Count == 0)
                                {
                                    factor = division.Quotient;
                                    division = factor.EuclideanDivideBy(v);
                                }
                                factorSplit.Add(v.GetPrimitivePart());
                                return true;
                            }
                            return false;
                        }
                        for (int i = indexOfPreviousElement + 1;
                            i <= liftedFactors.Count - elementsToAdd; ++i)
                        {
                            List<int> enlargedCombinationIndices =
                                new List<int>(combinationIndices);
                            enlargedCombinationIndices.Add(i);
                            if (testFactorCombination(elementsToAdd - 1, i,
                                enlargedCombinationIndices))
                            {
                                liftedFactors.RemoveAt(i);
                                return true;
                            }
                        }
                        return false;
                    }
                    int combinationSize = 1;
                    while (2 * combinationSize < liftedFactors.Count)
                    {
                        while (testFactorCombination(combinationSize, -1, new List<int>()))
                        { }
                        ++combinationSize;
                    }
                    if (2 * combinationSize == liftedFactors.Count)
                    {
                        if (testFactorCombination(combinationSize - 1, 0, new List<int> { 0 }))
                        {
                            liftedFactors.RemoveAt(0);
                        }
                    }
                    IntegerPolynomial finalFactor = new IntegerPolynomial(null,
                        factor.Coefficients[factor.Coefficients.Count - 1]);
                    foreach (IntegerPolynomial liftedFactor in liftedFactors)
                    {
                        finalFactor *= liftedFactor;
                    }
                    if (finalFactor.Coefficients.Count > 1)
                    {
                        BoundCoefficients(finalFactor);
                        factorSplit.Add(finalFactor.GetPrimitivePart());
                    }
                    return factorSplit;
                }
                List<IntegerPolynomial> irreducibleFactors = new List<IntegerPolynomial>();
                foreach (IntegerPolynomial factor in squarefreeFactors)
                {
                    IntegerPolynomial x = new IntegerPolynomial(null, Number.Zero, Number.One);
                    if (factor.Coefficients[0] == Number.Zero)
                    {
                        irreducibleFactors.Add(x);
                        factor.Coefficients.RemoveAt(0);
                        while (factor.Coefficients[0] == Number.Zero)
                        {
                            factor.Coefficients.RemoveAt(0);
                        }
                    }
                    if (factor.Coefficients[0].Magnitude() <
                        factor.Coefficients[factor.Coefficients.Count - 1].Magnitude())
                    {
                        factor.Coefficients.Reverse();
                        List<IntegerPolynomial> splitFactors = splitFactor(factor);
                        for (int i = 0; i < splitFactors.Count; ++i)
                        {
                            splitFactors[i].Coefficients.Reverse();
                            irreducibleFactors.Add(splitFactors[i]);
                        }
                    }
                    else
                    {
                        irreducibleFactors.AddRange(splitFactor(factor));
                    }
                }
                return irreducibleFactors;
            }
#if DEBUG
            public override string ToString()
            {
                return CoefficientsToString(Coefficients);
            }
#endif
        }
        class ModdedPolynomial : IArithmetic<ModdedPolynomial>
        {
            public List<Integer> Coefficients { get; }
            public Integer Characteristic { get; }
            public ModdedPolynomial(List<Integer> coefficients, Integer characteristic)
            {
                Coefficients = new List<Integer>();
                foreach (Integer coefficient in coefficients)
                {
                    Integer remainder = coefficient.EuclideanDivideBy(characteristic).Remainder;
                    if (remainder < Number.Zero)
                    {
                        Coefficients.Add(remainder + characteristic);
                    }
                    else
                    {
                        Coefficients.Add(remainder);
                    }
                }
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1] == Number.Zero)
                {
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                }
                Characteristic = characteristic;
            }
            public ModdedPolynomial(Integer characteristic, params Integer[] coefficients)
            {
                Coefficients = new List<Integer>();
                foreach (Integer coefficient in coefficients)
                {
                    Integer remainder = coefficient.EuclideanDivideBy(characteristic).Remainder;
                    if (remainder < Number.Zero)
                    {
                        Coefficients.Add(remainder + characteristic);
                    }
                    else
                    {
                        Coefficients.Add(remainder);
                    }
                }
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1] == Number.Zero)
                {
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                }
                Characteristic = characteristic;
            }
            public ModdedPolynomial GetAdditiveIdentity()
            {
                return new ModdedPolynomial(Characteristic);
            }
            public ModdedPolynomial GetMultiplicativeIdentity()
            {
                return new ModdedPolynomial(Characteristic, Number.One);
            }
            public ModdedPolynomial GetNegative()
            {
                List<Integer> output = new List<Integer>();
                foreach (Integer coefficient in Coefficients)
                {
                    output.Add(-coefficient);
                }
                return new ModdedPolynomial(output, Characteristic);
            }
            public ModdedPolynomial Plus(ModdedPolynomial a)
            {
                Debug.Assert(Characteristic == a.Characteristic,
                    "Operation between ModdedPolynomials of unlike Characteristics.");
                return new ModdedPolynomial(CoefficientAdd(Coefficients, a.Coefficients),
                    Characteristic);
            }
            public ModdedPolynomial Minus(ModdedPolynomial a)
            {
                return Plus(a.GetNegative());
            }
            public ModdedPolynomial Times(ModdedPolynomial a)
            {
                Debug.Assert(Characteristic == a.Characteristic,
                    "Operation between ModdedPolynomials of unlike Characteristics.");
                return new ModdedPolynomial(CoefficientMultiply(Coefficients, a.Coefficients),
                    Characteristic);
            }
            public Division<ModdedPolynomial> EuclideanDivideBy(ModdedPolynomial divisor)
            {
                Debug.Assert(Characteristic == divisor.Characteristic,
                    "Operation between ModdedPolynomials of unlike Characteristics.");
                Division<ModdedPolynomial> division;
                List<Integer> quotient = new List<Integer>();
                List<Integer> remainder = new List<Integer>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                {
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] * divisor.Coefficients[
                            divisor.Coefficients.Count - 1].Reciprocal(Characteristic));
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                        {
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                        }
                    }
                }
                division.Quotient = new ModdedPolynomial(quotient, Characteristic);
                division.Remainder = new ModdedPolynomial(remainder, Characteristic);
                return division;
            }
            public bool Equals(ModdedPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                {
                    return false;
                }
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    if (Coefficients[i] != a.Coefficients[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            public static ModdedPolynomial operator +(ModdedPolynomial a, ModdedPolynomial b)
            {
                return a.Plus(b);
            }
            public static ModdedPolynomial operator -(ModdedPolynomial a, ModdedPolynomial b)
            {
                return a.Minus(b);
            }
            public static ModdedPolynomial operator *(ModdedPolynomial a, ModdedPolynomial b)
            {
                return a.Times(b);
            }
            public static ModdedPolynomial operator *(ModdedPolynomial a, Integer b)
            {
                return a.Times(new ModdedPolynomial(a.Characteristic, b));
            }
            public static ModdedPolynomial operator *(Integer a, ModdedPolynomial b)
            {
                return b * a;
            }
            public static ModdedPolynomial operator /(ModdedPolynomial a, ModdedPolynomial b)
            {
                return a.EuclideanDivideBy(b).Quotient;
            }
            public static ModdedPolynomial operator %(ModdedPolynomial a, ModdedPolynomial b)
            {
                return a.EuclideanDivideBy(b).Remainder;
            }
            public ModdedPolynomial GetDerivative()
            {
                List<Integer> coefficients = new List<Integer>();
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    coefficients.Add(Coefficients[i] * new Integer(i));
                }
                return new ModdedPolynomial(coefficients, Characteristic);
            }
            public ModdedPolynomial GetMonic()
            {
                return this * Coefficients[Coefficients.Count - 1].Reciprocal(Characteristic);
            }
            public List<ModdedPolynomial> GetFactorsOfSquarefree()
            {
                Dictionary<int, ModdedPolynomial> distinctDegreeFactors =
                    new Dictionary<int, ModdedPolynomial>();
                ModdedPolynomial V = this;
                ModdedPolynomial W = new ModdedPolynomial(Characteristic, Number.Zero, Number.One);
                int d = 0;
                while (2 * d <= Coefficients.Count - 3)
                {
                    ++d;
                    W = Exponentiate(W, Characteristic) % this;
                    ModdedPolynomial degreeDFactorProduct = GetGCD(W -
                        new ModdedPolynomial(Characteristic, Number.Zero, Number.One), V);
                    if (degreeDFactorProduct.Coefficients.Count > 1)
                    {
                        distinctDegreeFactors.Add(d, degreeDFactorProduct);
                        V /= degreeDFactorProduct;
                        W %= V;
                    }
                }
                if (V.Coefficients.Count > 1)
                {
                    distinctDegreeFactors.Add(V.Coefficients.Count - 1, V);
                }
                List<ModdedPolynomial> irreducibleFactors = new List<ModdedPolynomial>();
                void cantorZassenhausSplit(ModdedPolynomial factorProduct, int degree)
                {
                    if ((factorProduct.Coefficients.Count - 1) / degree == 1)
                    {
                        irreducibleFactors.Add(factorProduct * factorProduct.Coefficients[
                            factorProduct.Coefficients.Count - 1].Reciprocal(Characteristic));
                        return;
                    }
                    ModdedPolynomial B;
                    if (Characteristic == Number.Two)
                    {
                        ModdedPolynomial T =
                            new ModdedPolynomial(Characteristic, Number.Zero, Number.One);
                        ModdedPolynomial C =
                            new ModdedPolynomial(new List<Integer>(T.Coefficients), Characteristic);
                        for (int j = 1; j < degree; ++j)
                        {
                            C = T + (C * C) % factorProduct;
                        }
                        B = GetGCD(factorProduct, C);
                        while (B.Coefficients.Count < 2 ||
                            B.Coefficients.Count == factorProduct.Coefficients.Count)
                        {
                            T *= new ModdedPolynomial(Characteristic, Number.Zero, Number.Zero,
                                Number.One);
                            C = new ModdedPolynomial(new List<Integer>(T.Coefficients),
                                Characteristic);
                            for (int j = 1; j < degree; ++j)
                            {
                                C = T + (C * C) % factorProduct;
                            }
                            B = GetGCD(factorProduct, C);
                        }
                    }
                    else
                    {
                        List<Integer> coefficients = new List<Integer>();
                        Integer p = Characteristic - Number.One;
                        for (int j = 1; j < 2 * degree; ++j)
                        {
                            coefficients.Add(p);
                        }
                        coefficients.Add(Number.One);
                        Integer power = ((Exponentiate(Characteristic,
                            new Integer(degree)) - Number.One) / Number.Two).Numerator;
                        ModdedPolynomial one = factorProduct.GetMultiplicativeIdentity();
                        ModdedPolynomial exponentiate()
                        {
                            return Exponentiate(new ModdedPolynomial(
                                new List<Integer>(coefficients), Characteristic), power,
                                delegate (ModdedPolynomial a, ModdedPolynomial b)
                                { return (a * b) % factorProduct; }, one);
                        }
                        B = GetGCD(factorProduct, exponentiate() - one);
                        while (B.Coefficients.Count < 2 ||
                            B.Coefficients.Count == factorProduct.Coefficients.Count)
                        {
                            int j = coefficients.Count - 2;
                            while (coefficients[j] == Number.Zero)
                            {
                                coefficients[j] = p;
                                --j;
                            }
                            --coefficients[j];
                            B = GetGCD(factorProduct, exponentiate() - one);
                        }
                    }
                    cantorZassenhausSplit(B, degree);
                    cantorZassenhausSplit(factorProduct / B, degree);
                }
                foreach (int degree in distinctDegreeFactors.Keys)
                {
                    cantorZassenhausSplit(distinctDegreeFactors[degree], degree);
                }
                return irreducibleFactors;
            }
#if DEBUG
            public override string ToString()
            {
                return CoefficientsToString(Coefficients) + " mod " + Characteristic.ToString();
            }
#endif
        }
        public class RationalPolynomial :
            IArithmetic<RationalPolynomial>, IMultipliable<RationalPolynomial, Rational>
        {
            public List<Rational> Coefficients { get; }

            //When not null, indicates that the polynomial variable represents a root of
            //MinimalPolynomial. The result of an operation between instances with like
            //MinimalPolynomial values is the unique polynomial that has the same numeric value as
            //the polynomial that results from the operation in the ordinary polynomial arithmetic
            //sense, but whose degree is less than that of MinimalPolynomial. Operations between
            //instances with unlike MinimalPolynomial values are invalid.
            public RationalPolynomial MinimalPolynomial { get; }
            public RationalPolynomial(List<Rational> coefficients,
                RationalPolynomial minimalPolynomial = null)
            {//Can mutate the coefficients parameter.                
                while (coefficients.Count != 0 &&
                    coefficients[coefficients.Count - 1] == Number.Zero)
                {
                    coefficients.RemoveAt(coefficients.Count - 1);
                }
                Coefficients = coefficients;
                MinimalPolynomial = minimalPolynomial;
            }
            public RationalPolynomial(RationalPolynomial minimalPolynomial,
                params Rational[] coefficients)
            {
                Coefficients = new List<Rational>(coefficients);
                while (Coefficients.Count != 0 &&
                    coefficients[Coefficients.Count - 1] == Number.Zero)
                {
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                }
                MinimalPolynomial = minimalPolynomial;
            }
            public RationalPolynomial GetAdditiveIdentity()
            {
                return new RationalPolynomial(MinimalPolynomial);
            }
            public RationalPolynomial GetMultiplicativeIdentity()
            {
                return new RationalPolynomial(MinimalPolynomial, Number.One);
            }
            public RationalPolynomial GetNegative()
            {
                List<Rational> output = new List<Rational>();
                foreach (Rational coefficient in Coefficients)
                {
                    output.Add(-coefficient);
                }
                return new RationalPolynomial(output, MinimalPolynomial);
            }
            public RationalPolynomial Plus(RationalPolynomial a)
            {
                return new RationalPolynomial(CoefficientAdd(Coefficients, a.Coefficients),
                    MinimalPolynomial);
            }
            public RationalPolynomial Minus(RationalPolynomial a)
            {
                return Plus(a.GetNegative());
            }
            public RationalPolynomial Times(RationalPolynomial a)
            {
                List<Rational> coefficients = CoefficientMultiply(Coefficients, a.Coefficients);
                if (MinimalPolynomial != null)
                {
                    for (int i = coefficients.Count; i >= MinimalPolynomial.Coefficients.Count; --i)
                    {
                        for (int j = 0; j < MinimalPolynomial.Coefficients.Count - 1; ++j)
                        {
                            coefficients[i - MinimalPolynomial.Coefficients.Count + j] -=
                                coefficients[coefficients.Count - 1] *
                                MinimalPolynomial.Coefficients[j];
                        }
                        coefficients.RemoveAt(coefficients.Count - 1);
                    }
                }
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public RationalPolynomial Times(Rational a)
            {
                List<Rational> coefficients = new List<Rational>();
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    coefficients.Add(a * Coefficients[i]);
                }
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public Division<RationalPolynomial> EuclideanDivideBy(RationalPolynomial divisor)
            {
                Division<RationalPolynomial> division;
                List<Rational> quotient = new List<Rational>();
                List<Rational> remainder = new List<Rational>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                {
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1]);
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                        {
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                        }
                    }
                }
                division.Quotient = new RationalPolynomial(quotient, MinimalPolynomial);
                division.Remainder = new RationalPolynomial(remainder, MinimalPolynomial);
                return division;
            }
            public bool Equals(RationalPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                {
                    return false;
                }
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    if (Coefficients[i] != a.Coefficients[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            public static RationalPolynomial operator +(RationalPolynomial a, RationalPolynomial b)
            {
                return a.Plus(b);
            }
            public static RationalPolynomial operator -(RationalPolynomial a, RationalPolynomial b)
            {
                return a.Plus(b.GetNegative());
            }
            public static RationalPolynomial operator -(RationalPolynomial a)
            {
                return a.GetNegative();
            }
            public static RationalPolynomial operator *(RationalPolynomial a, RationalPolynomial b)
            {
                return a.Times(b);
            }
            public static RationalPolynomial operator *(Rational a, RationalPolynomial b)
            {
                return b.Times(a);
            }
            public static RationalPolynomial operator *(RationalPolynomial a, Rational b)
            {
                return a.Times(b);
            }
            public static RationalPolynomial operator /(RationalPolynomial a, RationalPolynomial b)
            {
                if (b.MinimalPolynomial == null)
                {
                    return a.EuclideanDivideBy(b).Quotient;
                }
                if (b.Coefficients.Count == 1)
                {
                    return a * (Number.One / b.Coefficients[0]);
                }
                RationalPolynomial[] minimalPolynomials =
                    new RationalPolynomial[b.MinimalPolynomial.Coefficients.Count];
                minimalPolynomials[0] = b.MinimalPolynomial;
                for (int i = 1; i < b.MinimalPolynomial.Coefficients.Count; ++i)
                {
                    minimalPolynomials[i] = null;
                }
                MultivariatePolynomial multivariateB =
                    new MultivariatePolynomial(minimalPolynomials);
                for (int i = 0; i < b.Coefficients.Count; ++i)
                {
                    int[] indices = new int[b.MinimalPolynomial.Coefficients.Count];
                    indices[0] = i;
                    multivariateB.SetCoefficient(indices, b.Coefficients[i]);
                }
                MultivariatePolynomial reciprocal = new MultivariatePolynomial(minimalPolynomials);
                for (int i = 1; i < b.MinimalPolynomial.Coefficients.Count; ++i)
                {
                    int[] indices = new int[b.MinimalPolynomial.Coefficients.Count];
                    indices[0] = i - 1;
                    indices[i] = 1;
                    reciprocal.SetCoefficient(indices, Number.One);
                }
                MultivariatePolynomial product = multivariateB * reciprocal;
                List<List<Rational>> matrixRows = new List<List<Rational>>();
                List<Rational> augmentation = new List<Rational>();
                for (int i = 1; i < a.MinimalPolynomial.Coefficients.Count; ++i)
                {
                    List<Rational> row = new List<Rational>();
                    for (int j = 1; j < b.MinimalPolynomial.Coefficients.Count; ++j)
                    {
                        row.Add(Number.Zero);
                    }
                    matrixRows.Add(row);
                    augmentation.Add(Number.Zero);
                }
                augmentation[0] = Number.One;
                foreach (int[] indices in product.Coefficients.Keys)
                {
                    int columnIndex = 0;
                    for (int i = 1; i < indices.Length; ++i)
                    {
                        if (indices[i] == 1)
                        {
                            columnIndex = i - 1;
                            break;
                        }
                    }
                    matrixRows[indices[0]][columnIndex] = product.Coefficients[indices];
                }
                Matrix matrix = new Matrix(matrixRows);
                matrix.Diagonalize(augmentation);
                return a * new RationalPolynomial(augmentation, b.MinimalPolynomial);
            }
            public static RationalPolynomial operator %(RationalPolynomial a, RationalPolynomial b)
            {
                return a.EuclideanDivideBy(b).Remainder;
            }
            public static explicit operator RationalPolynomial(IntegerPolynomial a)
            {
                List<Rational> coefficients = new List<Rational>();
                foreach (Rational coefficient in a.Coefficients)
                {
                    coefficients.Add(coefficient);
                }
                return new RationalPolynomial(coefficients);
            }
            public RectangularEstimate EstimateEvaluation(Primitive argument,
                Rational errorIntervalSize)
            {
                Debug.Assert(Coefficients.Count > 1,
                    "RationalPolynomial.EstimateEvaluation was called from *this of degree less " +
                    "than 1, which the method doesn't account for.");
                Integer numberOfRealTermsInEvaluation;
                Integer numberOfImaginaryTermsInEvaluation;
                Division<Integer> division =
                    new Integer(Coefficients.Count).EuclideanDivideBy(Number.Two);
                if (division.Remainder == Number.Zero)
                {
                    numberOfRealTermsInEvaluation =
                        division.Quotient * (division.Quotient + Number.One) - Number.One;
                    numberOfImaginaryTermsInEvaluation = division.Quotient * division.Quotient;
                }
                else
                {
                    Integer quotientPlusOne = division.Quotient + Number.One;
                    numberOfRealTermsInEvaluation = quotientPlusOne * quotientPlusOne - Number.One;
                    numberOfImaginaryTermsInEvaluation =
                        division.Quotient * (division.Quotient + Number.One);
                }
                Rational realErrorScale = Coefficients[0].Magnitude();
                Rational imaginaryErrorScale = Number.One;
                Interval<Float> argumentRealPartEstimate = argument.GetRealPartEstimate(Number.One);
                Interval<Float> argumentImaginaryPartEstimate =
                    argument.GetImaginaryPartEstimate(Number.One);
                Rational realMagnitudeBound =
                    (Rational)GetMagnitudeInterval(argumentRealPartEstimate).Max;
                Rational imaginaryMagnitudeBound =
                    (Rational)GetMagnitudeInterval(argumentImaginaryPartEstimate).Max;
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    Rational coefficientMagnitude = Coefficients[i].Magnitude();
                    Rational realMagnitudePower =
                        Exponentiate(realMagnitudeBound, new Integer(i - 1));
                    Rational realScaleCandidate =
                        Number.Two * coefficientMagnitude * realMagnitudePower;
                    if (realScaleCandidate > realErrorScale)
                    {
                        realErrorScale = realScaleCandidate;
                    }
                    Rational imaginaryMagnitudePower = Number.One;
                    Rational imaginaryScaleCandidate;
                    Integer activeTermCount = numberOfImaginaryTermsInEvaluation;
                    Integer inactiveTermCount = numberOfRealTermsInEvaluation;
                    for (int j = 1; j < i; ++j)
                    {
                        Rational nextImaginaryMagnitudePower =
                            imaginaryMagnitudePower * imaginaryMagnitudeBound;
                        Rational sharedScaleComponent = Number.Two * activeTermCount *
                            new Integer(i).ThisChooseK(new Integer(j)) * coefficientMagnitude *
                            (Number.One + realMagnitudePower + nextImaginaryMagnitudePower);
                        realMagnitudePower /= realMagnitudeBound;
                        realScaleCandidate = sharedScaleComponent * realMagnitudePower;
                        if (realScaleCandidate > realErrorScale)
                        {
                            realErrorScale = realScaleCandidate;
                        }
                        imaginaryScaleCandidate = sharedScaleComponent * imaginaryMagnitudePower;
                        imaginaryMagnitudePower = nextImaginaryMagnitudePower;
                        if (imaginaryScaleCandidate > imaginaryErrorScale)
                        {
                            imaginaryErrorScale = imaginaryScaleCandidate;
                        }
                        Integer temp = activeTermCount;
                        activeTermCount = inactiveTermCount;
                        inactiveTermCount = temp;
                    }
                    imaginaryScaleCandidate = Number.Two * coefficientMagnitude *
                        imaginaryMagnitudePower * activeTermCount;
                    if (imaginaryScaleCandidate > imaginaryErrorScale)
                    {
                        imaginaryErrorScale = imaginaryScaleCandidate;
                    }
                }
                errorIntervalSize = Float.GetEstimatePlaceValue(errorIntervalSize);
                argumentRealPartEstimate =
                    argument.GetRealPartEstimate(errorIntervalSize / realErrorScale);
                argumentImaginaryPartEstimate =
                    argument.GetImaginaryPartEstimate(errorIntervalSize / imaginaryErrorScale);
                Interval<Rational> evaluationRealPart =
                    new Interval<Rational>(Coefficients[0], Coefficients[0]);
                Interval<Rational> evaluationImaginaryPart =
                    new Interval<Rational>(Number.Zero, Number.Zero);
                Rational realPartMin = (Rational)argumentRealPartEstimate.Min;
                Rational realPartMax = (Rational)argumentRealPartEstimate.Max;
                Rational imaginaryPartMin = (Rational)argumentImaginaryPartEstimate.Min;
                Rational imaginaryPartMax = (Rational)argumentImaginaryPartEstimate.Max;
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    Integer termDegree = new Integer(i);
                    Rational realPartMinPower = Exponentiate(realPartMin, termDegree);
                    Rational realPartMaxPower = Exponentiate(realPartMax, termDegree);
                    Rational imaginaryPartMinPower = Number.One;
                    Rational imaginaryPartMaxPower = Number.One;
                    for (int j = 0; j <= i; ++j)
                    {
                        Rational coefficient =
                            Coefficients[i] * termDegree.ThisChooseK(new Integer(j));
                        Rational termEstimate =
                            coefficient * realPartMinPower * imaginaryPartMinPower;
                        Interval<Rational> termInterval =
                            new Interval<Rational>(termEstimate, termEstimate);
                        List<Rational> termBoundCandidates = new List<Rational> {
                            coefficient * realPartMinPower * imaginaryPartMaxPower,
                            coefficient * realPartMaxPower * imaginaryPartMinPower,
                            coefficient * realPartMaxPower * imaginaryPartMaxPower };
                        foreach (Rational candidate in termBoundCandidates)
                        {
                            if (candidate > termInterval.Max)
                            {
                                termInterval.Max = candidate;
                            }
                            else if (candidate < termInterval.Min)
                            {
                                termInterval.Min = candidate;
                            }
                        }
                        int remainder = j % 4;
                        switch (remainder)
                        {
                            case 0:
                                evaluationRealPart.Min += termInterval.Min;
                                evaluationRealPart.Max += termInterval.Max;
                                break;
                            case 1:
                                evaluationImaginaryPart.Min += termInterval.Min;
                                evaluationImaginaryPart.Max += termInterval.Max;
                                break;
                            case 2:
                                evaluationRealPart.Min -= termInterval.Max;
                                evaluationRealPart.Max -= termInterval.Min;
                                break;
                            case 3:
                                evaluationImaginaryPart.Min -= termInterval.Max;
                                evaluationImaginaryPart.Max -= termInterval.Min;
                                break;
                        }
                        realPartMinPower /= realPartMin;
                        realPartMaxPower /= realPartMax;
                        imaginaryPartMinPower *= imaginaryPartMin;
                        imaginaryPartMaxPower *= imaginaryPartMax;
                    }
                }
                return new RectangularEstimate(new Interval<Float>(
                    evaluationRealPart.Min.GetRealPartEstimate(errorIntervalSize).Min,
                    evaluationRealPart.Max.GetRealPartEstimate(errorIntervalSize).Max),
                    new Interval<Float>(
                    evaluationImaginaryPart.Min.GetRealPartEstimate(errorIntervalSize).Min,
                    evaluationImaginaryPart.Max.GetRealPartEstimate(errorIntervalSize).Max));
            }
            public RationalPolynomial GetDerivative()
            {
                List<Rational> coefficients = new List<Rational>();
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    coefficients.Add(Coefficients[i] * new Integer(i));
                }
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public IntegerPolynomial GetPrimitivePart()
            {
                Integer denomiatorProduct = Number.One;
                foreach (Rational coefficient in Coefficients)
                {
                    denomiatorProduct *= coefficient.Denominator;
                }
                List<Integer> coefficients = new List<Integer>();
                foreach (Rational coefficient in Coefficients)
                {
                    coefficients.Add((denomiatorProduct * coefficient).Numerator);
                }
                return new IntegerPolynomial(MinimalPolynomial, coefficients).GetPrimitivePart();
            }
            public static RationalPolynomial GetGCD(RationalPolynomial a, RationalPolynomial b)
            {
                return (RationalPolynomial)IntegerPolynomial.GetGCD(a.GetPrimitivePart(),
                    b.GetPrimitivePart());
            }
            public RationalPolynomial GetValueMinimalPolynomial()
            {
                List<RationalPolynomial> powers = new List<RationalPolynomial>();
                RationalPolynomial power = GetMultiplicativeIdentity();
                List<List<Rational>> matrixRows = new List<List<Rational>>();
                List<RationalPolynomial> augmentation = new List<RationalPolynomial>();
                List<Rational> augmentationRow = new List<Rational> { Number.One };
                for (int i = 1; i < MinimalPolynomial.Coefficients.Count; ++i)
                {
                    power *= this;
                    powers.Add(power);
                    List<Rational> matrixRow = new List<Rational>(power.Coefficients);
                    while (matrixRow.Count < MinimalPolynomial.Coefficients.Count - 1)
                    {
                        matrixRow.Add(Number.Zero);
                    }
                    matrixRow.Reverse();
                    matrixRows.Add(matrixRow);
                    augmentationRow.Insert(0, Number.Zero);
                    augmentation.Add(new RationalPolynomial(new List<Rational>(augmentationRow)));
                }
                Matrix matrix = new Matrix(matrixRows);
                matrix.GetRowEchelonForm(augmentation);
                List<Rational> annullingPolynomialCoefficients =
                    augmentation[augmentation.Count - 1].Coefficients;
                annullingPolynomialCoefficients[0] -=
                    matrix.Rows[matrix.Rows.Count - 1][matrix.Rows.Count - 1];
                IntegerPolynomial minimalPolynomial =
                    new RationalPolynomial(annullingPolynomialCoefficients).GetPrimitivePart();
                List<IntegerPolynomial> factors = minimalPolynomial.GetFactors();
                foreach (IntegerPolynomial factor in factors)
                {
                    RationalPolynomial sum =
                        new RationalPolynomial(MinimalPolynomial, factor.Coefficients[0]);
                    for (int i = 1; i < factor.Coefficients.Count; ++i)
                    {
                        sum += powers[i - 1] * factor.Coefficients[i];
                    }
                    if (sum.Coefficients.Count == 0)
                    {
                        minimalPolynomial = factor;
                        break;
                    }
                }
                return new RationalPolynomial(GetMonicCoefficients(minimalPolynomial.Coefficients));
            }
#if DEBUG
            public override string ToString()
            {
                return CoefficientsToString(Coefficients);
            }
#endif
        }
        class NestedPolynomial : IArithmetic<NestedPolynomial>
        {//I refer to the coefficients as the inner polynomials, and the whole object as the outer
         //polynomial.
            public List<RationalPolynomial> Coefficients { get; }
            public RationalPolynomial InnerMinimalPolynomial { get; set; }
            public NestedPolynomial(RationalPolynomial innerMinimalPolynomial,
                List<RationalPolynomial> coefficients)
            {//Can mutate the coefficients parameter.
                while (coefficients.Count != 0 &&
                    coefficients[coefficients.Count - 1].Coefficients.Count == 0)
                {
                    coefficients.RemoveAt(coefficients.Count - 1);
                }
                Coefficients = coefficients;
                InnerMinimalPolynomial = innerMinimalPolynomial;
            }
            public NestedPolynomial(RationalPolynomial innerMinimalPolynomial,
                params RationalPolynomial[] coefficients)
            {
                Coefficients = new List<RationalPolynomial>(coefficients);
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1].Coefficients.Count == 0)
                {
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                }
                InnerMinimalPolynomial = innerMinimalPolynomial;
            }
            public NestedPolynomial GetAdditiveIdentity()
            {
                return new NestedPolynomial(InnerMinimalPolynomial);
            }
            public NestedPolynomial GetMultiplicativeIdentity()
            {
                return new NestedPolynomial(InnerMinimalPolynomial,
                    new RationalPolynomial(InnerMinimalPolynomial, Number.One));
            }
            public NestedPolynomial GetNegative()
            {
                List<RationalPolynomial> negativeCoefficients = new List<RationalPolynomial>();
                foreach (RationalPolynomial coefficient in Coefficients)
                {
                    negativeCoefficients.Add(coefficient.GetNegative());
                }
                return new NestedPolynomial(InnerMinimalPolynomial, negativeCoefficients);
            }
            public NestedPolynomial Plus(NestedPolynomial a)
            {
                return new NestedPolynomial(InnerMinimalPolynomial,
                    CoefficientAdd(Coefficients, a.Coefficients));
            }
            public NestedPolynomial Minus(NestedPolynomial a)
            {
                return Plus(a.GetNegative());
            }
            public NestedPolynomial Times(NestedPolynomial a)
            {
                return new NestedPolynomial(InnerMinimalPolynomial,
                    CoefficientMultiply(Coefficients, a.Coefficients));
            }
            public NestedPolynomial Times(RationalPolynomial a)
            {
                List<RationalPolynomial> productCoefficients = new List<RationalPolynomial>();
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    productCoefficients.Add(a * Coefficients[i]);
                }
                return new NestedPolynomial(InnerMinimalPolynomial, productCoefficients);
            }
            public Division<NestedPolynomial> EuclideanDivideBy(NestedPolynomial divisor)
            {
                Division<NestedPolynomial> division;
                List<RationalPolynomial> quotient = new List<RationalPolynomial>();
                List<RationalPolynomial> remainder = new List<RationalPolynomial>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                {
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1]);
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                        {
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                        }
                    }
                }
                division.Quotient = new NestedPolynomial(InnerMinimalPolynomial, quotient);
                division.Remainder = new NestedPolynomial(InnerMinimalPolynomial, remainder);
                return division;
            }
            public bool Equals(NestedPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                {
                    return false;
                }
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    if (!Coefficients[i].Equals(a.Coefficients[i]))
                    {
                        return false;
                    }
                }
                return true;
            }
            public static NestedPolynomial operator +(NestedPolynomial a, NestedPolynomial b)
            {
                return a.Plus(b);
            }
            public static NestedPolynomial operator -(NestedPolynomial a, NestedPolynomial b)
            {
                return a.Minus(b);
            }
            public static NestedPolynomial operator *(NestedPolynomial a, NestedPolynomial b)
            {
                return a.Times(b);
            }
            public static NestedPolynomial operator *(RationalPolynomial a, NestedPolynomial b)
            {
                return b.Times(a);
            }
            public static NestedPolynomial operator *(NestedPolynomial a, RationalPolynomial b)
            {
                return a.Times(b);
            }
            public static NestedPolynomial operator /(NestedPolynomial a, NestedPolynomial b)
            {
                return a.EuclideanDivideBy(b).Quotient;
            }
            public static NestedPolynomial operator /(NestedPolynomial a, RationalPolynomial b)
            {
                List<RationalPolynomial> quotientCoefficients = new List<RationalPolynomial>();
                foreach (RationalPolynomial coefficient in a.Coefficients)
                {
                    quotientCoefficients.Add(coefficient / b);
                }
                return new NestedPolynomial(a.InnerMinimalPolynomial, quotientCoefficients);
            }
            public static NestedPolynomial operator %(NestedPolynomial a, NestedPolynomial b)
            {
                return a.EuclideanDivideBy(b).Remainder;
            }
            NestedPolynomial GetDerivative()
            {
                List<RationalPolynomial> derivativeCoefficients = new List<RationalPolynomial>();
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    derivativeCoefficients.Add(new Integer(i) * Coefficients[i]);
                }
                return new NestedPolynomial(InnerMinimalPolynomial, derivativeCoefficients);
            }
            NestedPolynomial GetPrimitivePart()
            {
                return this / Solver.GetGCD(Coefficients);
            }
            static RationalPolynomial GetResultant(NestedPolynomial a, NestedPolynomial b)
            {
                if (a.Coefficients.Count == 0 || b.Coefficients.Count == 0)
                {
                    return new RationalPolynomial(a.InnerMinimalPolynomial);
                }
                RationalPolynomial aContent = Solver.GetGCD(a.Coefficients);
                RationalPolynomial bContent = Solver.GetGCD(b.Coefficients);
                a /= aContent;
                b /= bContent;
                RationalPolynomial g = a.Coefficients[0].GetMultiplicativeIdentity();
                RationalPolynomial h = a.Coefficients[0].GetMultiplicativeIdentity();
                Rational s = Number.One;
                if (a.Coefficients.Count % 2 == 0 && b.Coefficients.Count % 2 == 0)
                {
                    s = -s;
                }
                RationalPolynomial t =
                    Exponentiate(aContent, new Integer(b.Coefficients.Count - 1)) *
                    Exponentiate(bContent, new Integer(a.Coefficients.Count - 1));
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    NestedPolynomial c = a;
                    a = b;
                    b = c;
                }
                Integer hExponent;
                while (b.Coefficients.Count > 1)
                {
                    Integer degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    if (a.Coefficients.Count % 2 == 0 && b.Coefficients.Count % 2 == 0)
                    {
                        s = -s;
                    }
                    NestedPolynomial remainder = (a * Exponentiate(
                        b.Coefficients[b.Coefficients.Count - 1], degree + Number.One)) % b;
                    a = b;
                    b = remainder / (g * Exponentiate(h, degree));
                    g = a.Coefficients[a.Coefficients.Count - 1];
                    hExponent = Number.One - degree;
                    if (hExponent >= Number.Zero)
                    {
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    }
                    else
                    {
                        h = Exponentiate(g, degree) / Exponentiate(h, -hExponent);
                    }
                }
                hExponent = new Integer(2 - a.Coefficients.Count);
                if (hExponent >= Number.Zero)
                {
                    return s * t * Exponentiate(h, hExponent) * Exponentiate(b.Coefficients[
                        b.Coefficients.Count - 1], new Integer(a.Coefficients.Count - 1));
                }
                return s * t * Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                    new Integer(a.Coefficients.Count - 1)) / Exponentiate(h, -hExponent);
            }
            static NestedPolynomial GetGCD(NestedPolynomial a, NestedPolynomial b)
            {
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    if (a.Coefficients.Count == 0)
                    {
                        return b;
                    }
                    NestedPolynomial t = a;
                    a = b;
                    b = t;
                }
                else if (b.Coefficients.Count == 0)
                {
                    return a;
                }
                RationalPolynomial aContent = Solver.GetGCD(a.Coefficients);
                RationalPolynomial bContent = Solver.GetGCD(b.Coefficients);
                a /= aContent;
                b /= bContent;
                RationalPolynomial d = Solver.GetGCD(aContent, bContent);
                RationalPolynomial g = a.Coefficients[0].GetMultiplicativeIdentity();
                RationalPolynomial h = a.Coefficients[0].GetMultiplicativeIdentity();
                Integer degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                NestedPolynomial remainder =
                    a.Times(Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                    degree + Number.One)) % b;
                while (remainder.Coefficients.Count > 1)
                {
                    a = b;
                    b = remainder / (g * Exponentiate(h, degree));
                    g = a.Coefficients[a.Coefficients.Count - 1];
                    Integer hExponent = Number.One - degree;
                    if (hExponent >= Number.Zero)
                    {
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    }
                    else
                    {
                        h = Exponentiate(g, degree) / Exponentiate(h, -hExponent);
                    }
                    degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    remainder = a.Times(Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                        degree + Number.One)) % b;
                }
                if (remainder.Coefficients.Count == 1)
                {
                    return new NestedPolynomial(a.InnerMinimalPolynomial, d);
                }
                return b.GetPrimitivePart().Times(d);
            }
            public NestedPolynomial SwitchVariables(
                RationalPolynomial innerMinimalPolynomial = null)
            {
                if (Coefficients.Count == 0)
                {
                    return this;
                }
                List<List<Rational>> coefficientCoefficients = new List<List<Rational>>();
                foreach (Rational coefficient in Coefficients[0].Coefficients)
                {
                    coefficientCoefficients.Add(new List<Rational> { coefficient });
                }
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    for (int j = coefficientCoefficients.Count;
                        j < Coefficients[i].Coefficients.Count; ++j)
                    {
                        coefficientCoefficients.Add(new List<Rational>());
                        for (int k = 0; k < i; ++k)
                        {
                            coefficientCoefficients[j].Add(Number.Zero);
                        }
                    }
                    for (int j = 0; j < Coefficients[i].Coefficients.Count; ++j)
                    {
                        coefficientCoefficients[j].Add(Coefficients[i].Coefficients[j]);
                    }
                }
                List<RationalPolynomial> coefficients = new List<RationalPolynomial>();
                foreach (List<Rational> list in coefficientCoefficients)
                {
                    coefficients.Add(new RationalPolynomial(list, innerMinimalPolynomial));
                }
                return new NestedPolynomial(innerMinimalPolynomial, coefficients);
            }
            public List<NestedPolynomial> GetFactors()
            {
                List<NestedPolynomial> squarefreeFactors = new List<NestedPolynomial>();
                NestedPolynomial derivative = GetDerivative();
                NestedPolynomial a = GetGCD(this, derivative);
                NestedPolynomial b = this / a;
                NestedPolynomial c = derivative / a - b.GetDerivative();
                while (!(b.Coefficients.Count == 1 && (b.Coefficients[0].Equals(
                    b.Coefficients[0].GetMultiplicativeIdentity()) || b.Coefficients[0].Equals(
                    b.Coefficients[0].GetMultiplicativeIdentity().GetNegative()))))
                {
                    a = GetGCD(b, c);
                    if (!squarefreeFactors.Contains(a))
                    {
                        squarefreeFactors.Add(a);
                    }
                    b = b / a;
                    c = c / a - b.GetDerivative();
                }
                List<NestedPolynomial> irreducibleFactors = new List<NestedPolynomial>();
                NestedPolynomial raisedMinimalPolynomial =
                    new NestedPolynomial(null, InnerMinimalPolynomial).SwitchVariables();
                foreach (NestedPolynomial squarefreeFactor in squarefreeFactors)
                {
                    RationalPolynomial resultant;
                    Integer k = Number.Zero;
                    while (true)
                    {
                        NestedPolynomial power =
                            new NestedPolynomial(null, new RationalPolynomial(null, Number.One));
                        NestedPolynomial d = new NestedPolynomial(null, new RationalPolynomial(null,
                            Number.Zero, k), new RationalPolynomial(null, Number.One));
                        NestedPolynomial e = new NestedPolynomial(null);
                        for (int i = 0; i < squarefreeFactor.Coefficients.Count; ++i)
                        {
                            e += squarefreeFactor.Coefficients[i] * power;
                            power *= d;
                        }
                        resultant = GetResultant(raisedMinimalPolynomial, e.SwitchVariables());
                        if (RationalPolynomial.GetGCD(resultant,
                            resultant.GetDerivative()).Coefficients.Count < 2)
                        {
                            break;
                        }
                        --k;
                    }
                    List<IntegerPolynomial> integerFactors =
                        resultant.GetPrimitivePart().GetFactors();
                    foreach (IntegerPolynomial factor in integerFactors)
                    {
                        NestedPolynomial power = GetMultiplicativeIdentity();
                        NestedPolynomial d = new NestedPolynomial(null,
                            new RationalPolynomial(InnerMinimalPolynomial, Number.Zero, -k),
                            new RationalPolynomial(InnerMinimalPolynomial, Number.One),
                            InnerMinimalPolynomial);
                        NestedPolynomial e = new NestedPolynomial(InnerMinimalPolynomial);
                        for (int i = 0; i < factor.Coefficients.Count; ++i)
                        {
                            e += new RationalPolynomial(InnerMinimalPolynomial,
                                factor.Coefficients[i]) * power;
                            power *= d;
                        }
                        NestedPolynomial irreducibleFactor = Solver.GetGCD(squarefreeFactor, e);
                        irreducibleFactors.Add(irreducibleFactor / irreducibleFactor.Coefficients[
                            irreducibleFactor.Coefficients.Count - 1]);
                    }
                }
                return irreducibleFactors;
            }
#if DEBUG
            public override string ToString()
            {
                if (Coefficients.Count == 0)
                {
                    return "0";
                }
                StringBuilder output = new StringBuilder();
                if (Coefficients[0].Coefficients.Count != 0)
                {
                    output.Append("(" + Coefficients[0] + ")");
                }
                if (Coefficients.Count > 1)
                {
                    if (Coefficients[1].Coefficients.Count != 0)
                    {
                        if (output.Length != 0)
                        {
                            output.Append("+");
                        }
                        output.Append("(" + Coefficients[1] + ")y");
                    }
                    for (int i = 2; i < Coefficients.Count; ++i)
                    {
                        if (Coefficients[i].Coefficients.Count != 0)
                        {
                            if (output.Length != 0)
                            {
                                output.Append("+");
                            }
                            output.Append("(" + Coefficients[i] + ")y^" + i);
                        }
                    }
                }
                return output.ToString();
            }
#endif
        }
        class MultivariatePolynomial : IRingElement<MultivariatePolynomial>
        {//The nth variable represents a root of MinimalPolynomials[n] if MinimalPolynomials[n] is
         //non-null, and member methods accordingly keep the degrees of all terms with respect to it
         //less than the degree of MinimalPolynomials[n]. Operations between MultivariatePolynomials
         //whose Coefficients fields don't match are invalid.
            public RationalPolynomial[] MinimalPolynomials { get; }
            public Dictionary<int[], Rational> Coefficients { get; }
            public MultivariatePolynomial(RationalPolynomial[] minimalPolynomials)
            {
                MinimalPolynomials = minimalPolynomials;
                Coefficients = new Dictionary<int[], Rational>();
            }
            static bool AreEqualByValue(int[] a, int[] b)
            {
                for (int i = 0; i < a.Length; ++i)
                {
                    if (a[i] != b[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            public void SetCoefficient(int[] indices, Rational coefficient)
            {//The nth index is the degree with respect to the nth variable of the term being set.
                if (coefficient == Number.Zero)
                {
                    return;
                }
                foreach (int[] keyIndices in Coefficients.Keys)
                {
                    if (AreEqualByValue(indices, keyIndices))
                    {
                        Coefficients[indices] = coefficient;
                        return;
                    }
                }
                Coefficients.Add(indices, coefficient);
            }
            public MultivariatePolynomial GetMultiplicativeIdentity()
            {
                MultivariatePolynomial one = new MultivariatePolynomial(MinimalPolynomials);
                one.SetCoefficient(new int[MinimalPolynomials.Length], Number.One);
                return one;
            }
            public MultivariatePolynomial Times(MultivariatePolynomial a)
            {
                return a * this;
            }
            public static MultivariatePolynomial operator +(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial sum = new MultivariatePolynomial(a.MinimalPolynomials);
                foreach (int[] indices in a.Coefficients.Keys)
                {
                    sum.Coefficients.Add(indices, a.Coefficients[indices]);
                }
                foreach (int[] indicesB in b.Coefficients.Keys)
                {
                    bool indicesFound = false;
                    foreach (int[] indicesA in sum.Coefficients.Keys)
                    {
                        if (AreEqualByValue(indicesA, indicesB))
                        {
                            Rational termSum =
                                sum.Coefficients[indicesA] + b.Coefficients[indicesB];
                            if (termSum == Number.Zero)
                            {
                                sum.Coefficients.Remove(indicesA);
                            }
                            else
                            {
                                sum.Coefficients[indicesA] = termSum;
                            }
                            indicesFound = true;
                            break;
                        }
                    }
                    if (!indicesFound)
                    {
                        sum.Coefficients.Add(indicesB, b.Coefficients[indicesB]);
                    }
                }
                return sum;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a,
                MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(a.MinimalPolynomials);
                foreach (int[] aIndices in a.Coefficients.Keys)
                {
                    foreach (int[] bIndices in b.Coefficients.Keys)
                    {
                        MultivariatePolynomial term =
                            new MultivariatePolynomial(a.MinimalPolynomials);
                        term.SetCoefficient(new int[aIndices.Length],
                            a.Coefficients[aIndices] * b.Coefficients[bIndices]);
                        for (int i = 0; i < aIndices.Length; ++i)
                        {
                            int termFactorDegree = aIndices[i] + bIndices[i];
                            if (a.MinimalPolynomials[i] == null ||
                                termFactorDegree < a.MinimalPolynomials[i].Coefficients.Count - 1)
                            {
                                foreach (int[] indices in term.Coefficients.Keys)
                                {
                                    indices[i] = termFactorDegree;
                                }
                            }
                            else
                            {
                                MultivariatePolynomial termFactor =
                                    new MultivariatePolynomial(a.MinimalPolynomials);
                                int[] termFactorIndices = new int[aIndices.Length];
                                termFactorIndices[i] = termFactorDegree -
                                    a.MinimalPolynomials[i].Coefficients.Count + 1;
                                termFactor.SetCoefficient(termFactorIndices, Number.One);
                                MultivariatePolynomial reducedFactorComponent =
                                    new MultivariatePolynomial(a.MinimalPolynomials);
                                for (int j = 0; j < a.MinimalPolynomials[i].Coefficients.Count - 1;
                                    ++j)
                                {
                                    int[] termIndices = new int[aIndices.Length];
                                    termIndices[i] = j;
                                    reducedFactorComponent.SetCoefficient(termIndices,
                                        -a.MinimalPolynomials[i].Coefficients[j]);
                                }
                                term *= termFactor * reducedFactorComponent;
                            }
                        }
                        product += term;
                    }
                }
                return product;
            }
            public static MultivariatePolynomial operator *(Rational a, MultivariatePolynomial b)
            {
                MultivariatePolynomial product = new MultivariatePolynomial(b.MinimalPolynomials);
                if (a == Number.Zero)
                {
                    return product;
                }
                foreach (int[] indices in b.Coefficients.Keys)
                {
                    product.Coefficients.Add(indices, a * b.Coefficients[indices]);
                }
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a, Rational b)
            {
                return b * a;
            }
            public RationalPolynomial CalculateValueMinimalPolynomial(Primitive value)
            {
                List<List<Rational>> matrixCoefficients = new List<List<Rational>>();
                List<int[]> termsPresent = new List<int[]> { new int[MinimalPolynomials.Length] };
                List<MultivariatePolynomial> powers = new List<MultivariatePolynomial>();
                MultivariatePolynomial power =
                    new MultivariatePolynomial(MinimalPolynomials);
                power.Coefficients[new int[MinimalPolynomials.Length]] = Number.One;
                bool constantIsPresent = false;
                while (!constantIsPresent || powers.Count < termsPresent.Count)
                {
                    power *= this;
                    powers.Add(power);
                    List<Rational> matrixRow = new List<Rational>();
                    for (int i = 0; i < termsPresent.Count; ++i)
                    {
                        matrixRow.Add(Number.Zero);
                    }
                    foreach (int[] indices in power.Coefficients.Keys)
                    {
                        bool indicesFound = false;
                        for (int i = 0; i < termsPresent.Count; ++i)
                        {
                            if (AreEqualByValue(indices, termsPresent[i]))
                            {
                                matrixRow[i] = power.Coefficients[indices];
                                indicesFound = true;
                                break;
                            }
                        }
                        if (!indicesFound)
                        {
                            termsPresent.Add(indices);
                            for (int i = 0; i < matrixCoefficients.Count; ++i)
                            {
                                matrixCoefficients[i].Add(Number.Zero);
                            }
                            matrixRow.Add(power.Coefficients[indices]);
                        }
                    }
                    matrixCoefficients.Add(matrixRow);
                    if (matrixRow[0] != Number.Zero)
                    {
                        constantIsPresent = true;
                    }
                }
                foreach (List<Rational> row in matrixCoefficients)
                {
                    row.Reverse();
                }
                Matrix matrix = new Matrix(matrixCoefficients);
                List<RationalPolynomial> augmentation = new List<RationalPolynomial>();
                List<Integer> augmentationRow = new List<Integer> { Number.One };
                for (int i = 0; i < powers.Count; ++i)
                {
                    augmentationRow.Insert(0, Number.Zero);
                    augmentation.Add(new RationalPolynomial(new List<Rational>(augmentationRow)));
                }
                matrix.GetRowEchelonForm(augmentation);
                List<Rational> annullingPolynomialCoefficients =
                    augmentation[matrixCoefficients.Count - 1].Coefficients;
                annullingPolynomialCoefficients[0] -=
                    matrixCoefficients[matrixCoefficients.Count - 1][matrixCoefficients.Count - 1];
                IntegerPolynomial minimalPolynomial =
                    new RationalPolynomial(annullingPolynomialCoefficients).GetPrimitivePart();
                List<IntegerPolynomial> factors = minimalPolynomial.GetFactors();
                List<RationalPolynomial> rationalFactors = new List<RationalPolynomial>();
                foreach (IntegerPolynomial factor in factors)
                {
                    rationalFactors.Add((RationalPolynomial)factor);
                }
                Rational errorIntervalSize = Number.One;
                while (rationalFactors.Count > 1)
                {
                    for (int i = 0; i < rationalFactors.Count;)
                    {
                        RectangularEstimate factorEvaluatedAtThis =
                            rationalFactors[i].EstimateEvaluation(value, errorIntervalSize);
                        if (factorEvaluatedAtThis.RealPart.Max < Float.Zero ||
                            factorEvaluatedAtThis.RealPart.Min > Float.Zero ||
                            factorEvaluatedAtThis.ImaginaryPart.Max < Float.Zero ||
                            factorEvaluatedAtThis.ImaginaryPart.Min > Float.Zero)
                        {
                            rationalFactors.RemoveAt(i);
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    errorIntervalSize /= Number.Two;
                }
                return new RationalPolynomial(
                    GetMonicCoefficients(rationalFactors[0].Coefficients));
            }
#if DEBUG
            public override string ToString()
            {
                if (Coefficients.Count == 0)
                {
                    return "0";
                }
                StringBuilder output = new StringBuilder();
                foreach (int[] indices in Coefficients.Keys)
                {
                    StringBuilder term = new StringBuilder();
                    for (int i = 0; i < indices.Length; ++i)
                    {
                        if (indices[i] != 0)
                        {
                            term.Append("(" + i + ")");
                            if (indices[i] != 1)
                            {
                                term.Append("^" + indices[i]);
                            }
                        }
                    }
                    output.Append(Coefficients[indices].InsertString(term.ToString()) + '+');
                }
                output.Remove(output.Length - 1, 1);
                return output.ToString();
            }
#endif
        }
        public class Float : IRingElement<Float>
        {//Behaves like a Rational whose denominator is limited to a power of two. 
            public Integer Significand { get; }
            public Integer RadixPointPosition { get; }
            public static Float Zero = new Float(Number.Zero, Number.Zero);
            public static Float One = new Float(Number.One, Number.Zero);
            public Float(Integer significand, Integer radixPointPosition)
            {
                Significand = significand;
                RadixPointPosition = radixPointPosition;
            }
            public static Float GetReducedForm(Integer significand, Integer radixPointPosition)
            {
                Division<Integer> division = significand.EuclideanDivideBy(Number.Two);
                while (radixPointPosition > Number.Zero && division.Remainder == Number.Zero)
                {
                    significand = division.Quotient;
                    --radixPointPosition;
                    division = significand.EuclideanDivideBy(Number.Two);
                }
                return new Float(significand, radixPointPosition);
            }
            public Float GetMultiplicativeIdentity()
            {
                return One;
            }
            public Float Magnitude()
            {
                if (Significand < Number.Zero)
                {
                    return new Float(-Significand, RadixPointPosition);
                }
                return this;
            }
            public Float Times(Float a)
            {
                return this * a;
            }
            public static Float operator +(Float a, Float b)
            {
                if (a.RadixPointPosition == b.RadixPointPosition)
                {
                    return GetReducedForm(a.Significand + b.Significand, a.RadixPointPosition);
                }
                if (b.RadixPointPosition > a.RadixPointPosition)
                {
                    return new Float(b.Significand + a.Significand * Exponentiate(Number.Two,
                        b.RadixPointPosition - a.RadixPointPosition), b.RadixPointPosition);
                }
                return new Float(a.Significand + b.Significand * Exponentiate(Number.Two,
                    a.RadixPointPosition - b.RadixPointPosition), a.RadixPointPosition);
            }
            public static Float operator -(Float a, Float b)
            {
                return a + -b;
            }
            public static Float operator -(Float a)
            {
                return new Float(-a.Significand, a.RadixPointPosition);
            }
            public static Float operator *(Float a, Float b)
            {
                return GetReducedForm(a.Significand * b.Significand,
                    a.RadixPointPosition + b.RadixPointPosition);
            }
            public static bool operator <(Float a, Float b)
            {
                return (a - b).Significand < Number.Zero;
            }
            public static bool operator >(Float a, Float b)
            {
                return (a - b).Significand > Number.Zero;
            }
            public static bool operator <=(Float a, Float b)
            {
                return !(a > b);
            }
            public static bool operator >=(Float a, Float b)
            {
                return !(a < b);
            }
            public static Float Max(Float a, Float b)
            {
                if (a > b)
                {
                    return a;
                }
                return b;
            }
            public Float ShiftRadixPointLeft(int placesToShift)
            {
                Debug.Assert(placesToShift >= 0);
                Integer significand = Significand;
                if (RadixPointPosition == Number.Zero)
                {
                    Division<Integer> division = significand.EuclideanDivideBy(Number.Two);
                    while (placesToShift > 0 && division.Remainder == Number.Zero)
                    {
                        significand = division.Quotient;
                        --placesToShift;
                        division = significand.EuclideanDivideBy(Number.Two);
                    }
                }
                return new Float(significand, RadixPointPosition + new Integer(placesToShift));
            }
            public static Rational GetEstimatePlaceValue(Rational errorIntervalSize)
            {
                Rational placeValue = Number.One;
                while (placeValue > errorIntervalSize)
                {
                    placeValue /= Number.Two;
                }
                return placeValue / Number.Two;
            }
            public Interval<Float> EstimateRoot(Rational errorIntervalSize, Integer index)
            {
                Debug.Assert(Significand >= Number.Zero, "Float.EstimateRoot was called from a " +
                    "negative *this, which the method doesn't account for.");
                Float rootEstimate;
                if (this < One)
                {
                    rootEstimate = One;
                }
                else
                {
                    rootEstimate = this;
                }
                Rational rationalRadicand = (Rational)this;
                Integer indexMinusOne = index - Number.One;
                Rational delta = (rationalRadicand / (Rational)Exponentiate(rootEstimate,
                    indexMinusOne) - (Rational)rootEstimate) / index;
                Interval<Float> deltaRealPartEstimate =
                    delta.GetRealPartEstimate(-delta / Number.Two);
                while ((Rational)deltaRealPartEstimate.Max - Number.Two * delta > errorIntervalSize)
                {
                    rootEstimate += deltaRealPartEstimate.Max;
                    delta = (rationalRadicand / (Rational)Exponentiate(rootEstimate,
                        indexMinusOne) - (Rational)rootEstimate) / index;
                    deltaRealPartEstimate = delta.GetRealPartEstimate(-delta / Number.Two);
                }
                rootEstimate += deltaRealPartEstimate.Max;
                return new Interval<Float>(rootEstimate + deltaRealPartEstimate.Max, rootEstimate);
            }
#if DEBUG
            public override string ToString()
            {
                return ((Rational)this).ToString();
            }
#endif
        }
        static class Pi
        {
            //K is the index of the next approximation iteration.            
            static Integer SixteenToTheK = new Integer(16);
            static Integer EightK = new Integer(8);
            public static Rational ErrorIntervalSize { get; private set; } =
                new Fraction(new Integer(1696), new Integer(12285));
            public static Rational LowEstimate { get; private set; } =
                new Fraction(new Integer(47), new Integer(15));
            public static Rational HighEstimate { get => LowEstimate + ErrorIntervalSize; }
            public static void RefineErrorInterval()
            {
                Integer four = new Integer(4);
                LowEstimate += (four / (EightK + Number.One) - Number.Two / (EightK + four) -
                    Number.One / (EightK + new Integer(5)) -
                    Number.One / (EightK + new Integer(6))) / SixteenToTheK;
                Integer sixteen = new Integer(16);
                ErrorIntervalSize = ErrorIntervalSize / sixteen;
                SixteenToTheK *= sixteen;
                EightK += new Integer(8);
            }
            public static void RefineErrorInterval(Rational errorIntervalSize)
            {
                if (!ReferenceEquals(errorIntervalSize, null) &&
                    ErrorIntervalSize > errorIntervalSize)
                {
                    RefineErrorInterval();
                    while (ErrorIntervalSize > errorIntervalSize)
                    {
                        RefineErrorInterval();
                    }
                }
            }
            public static void ShrinkErrorIntervalToOneSideOfValue(Rational value)
            {
                while (LowEstimate <= value && value <= HighEstimate)
                {
                    RefineErrorInterval();
                }
            }
        }
        static Interval<Rational> EstimateAtan2(Rational y, Rational x, Rational errorIntervalSize)
        {
            Rational multipleOfPiToAdd = Number.Zero;
            if (x < Number.Zero) 
            {
                multipleOfPiToAdd = Number.One;
            }
            else if (y < Number.Zero)
            {
                multipleOfPiToAdd = Number.Two;
            }
            Rational ratio = y / x;
            Interval<Rational> atan2Value;
            if (ratio.Magnitude() < Number.One)
            {
                atan2Value = ratio.EstimateArctangent(errorIntervalSize);
            }
            else
            {
                if (ratio > Number.Zero)
                {
                    multipleOfPiToAdd += new Fraction(Number.One, Number.Two);
                }
                else
                {
                    multipleOfPiToAdd -= new Fraction(Number.One, Number.Two);
                }
                Interval<Rational> arctanOfReciprocal =
                    (x / y).EstimateArctangent(errorIntervalSize);
                atan2Value.Min = -arctanOfReciprocal.Max;
                atan2Value.Max = -arctanOfReciprocal.Min;
            }
            if (multipleOfPiToAdd == Number.Zero)
            {
                return atan2Value;
            }
            errorIntervalSize /= Number.Two;
            Pi.RefineErrorInterval(errorIntervalSize / multipleOfPiToAdd.Magnitude());
            if (multipleOfPiToAdd > Number.Zero) 
            {
                atan2Value.Min += multipleOfPiToAdd * Pi.LowEstimate;
                atan2Value.Max += multipleOfPiToAdd * Pi.HighEstimate;
            }
            else
            {
                atan2Value.Min += multipleOfPiToAdd * Pi.HighEstimate;
                atan2Value.Max += multipleOfPiToAdd * Pi.LowEstimate;
            }
            return atan2Value;
        }
        static class RootsOfUnity
        {
            static Dictionary<Integer, List<Number>> NthRoots = new Dictionary<Integer,
                List<Number>> { { Number.One, new List<Number> { Number.One } },
                { Number.Two, new List<Number> { Number.One, -Number.One } } };
            public static List<Number> GetNthRoots(Integer n)
            {
                if (NthRoots.ContainsKey(n))
                {
                    return NthRoots[n];
                }
                List<Number> nthRoots = new List<Number>();
                Division<Integer> division;
                Integer half = n.EuclideanDivideBy(Number.Two).Quotient;
                void sortRoots(List<Number> roots, Integer rootIndex)
                {
                    Rational errorIntervalSize = Pi.LowEstimate / rootIndex;
                    roots.Sort(delegate (Number a, Number b)
                    {
                        if (a.GetArgumentEstimate(errorIntervalSize).Max <
                            b.GetArgumentEstimate(errorIntervalSize).Min)
                        {
                            return -1;
                        }
                        return 1;
                    });
                }
                foreach (Integer prime in Primes())
                {
                    if (prime > half)
                    {
                        break;
                    }
                    division = n.EuclideanDivideBy(prime);
                    if (division.Remainder == Number.Zero)
                    {
                        List<List<Number>> rootCombinations = GetCartesianProduct(
                            new List<List<Number>> { GetNthRoots(prime),
                        GetNthRoots(division.Quotient) });
                        foreach (List<Number> combination in rootCombinations)
                        {
                            nthRoots.Add(combination[0] * combination[1].Exponentiate(
                                new Fraction(Number.One, prime)));
                        }
                        sortRoots(nthRoots, n);
                        NthRoots.Add(n, nthRoots);
                        return nthRoots;
                    }
                }
                Integer nMinusOne = n - Number.One;
                Integer unfactoredComponent = nMinusOne;
                half = unfactoredComponent.EuclideanDivideBy(Number.Two).Quotient;
                List<Integer> factors = new List<Integer>();
                foreach (Integer prime in Primes())
                {
                    if (prime > half)
                    {
                        break;
                    }
                    division = unfactoredComponent.EuclideanDivideBy(prime);
                    if (division.Remainder == Number.Zero)
                    {
                        factors.Add(prime);
                        unfactoredComponent = division.Quotient;
                        division = unfactoredComponent.EuclideanDivideBy(prime);
                        while (division.Remainder == Number.Zero)
                        {
                            unfactoredComponent = division.Quotient;
                            division = unfactoredComponent.EuclideanDivideBy(prime);
                        }
                        half = unfactoredComponent.EuclideanDivideBy(Number.Two).Quotient;
                    }
                }
                if (unfactoredComponent != Number.One)
                {
                    factors.Add(unfactoredComponent);
                }
                Integer generator = Number.One;
                bool isGenerator = false;
                while (!isGenerator)
                {
                    ++generator;
                    isGenerator = true;
                    foreach (Integer factor in factors)
                    {
                        if (Exponentiate(generator,
                            nMinusOne.EuclideanDivideBy(factor).Quotient) == Number.One)
                        {
                            isGenerator = false;
                            break;
                        }
                    }
                }
                List<Rational> nMinusFirstRootAnnullingPolynomialCoefficients =
                    new List<Rational> { -Number.One };
                for (Integer i = Number.One; i < nMinusOne; ++i)
                {
                    nMinusFirstRootAnnullingPolynomialCoefficients.Add(Number.Zero);
                }
                List<Rational> nthRootAnnullingPolynomialCoefficients =
                    new List<Rational>(nMinusFirstRootAnnullingPolynomialCoefficients);
                nMinusFirstRootAnnullingPolynomialCoefficients.Add(Number.One);
                nthRootAnnullingPolynomialCoefficients.Add(Number.Zero);
                nthRootAnnullingPolynomialCoefficients.Add(Number.One);
                RationalPolynomial nMinusFirstRootAnnullingPolynomial =
                    new RationalPolynomial(nMinusFirstRootAnnullingPolynomialCoefficients);
                RationalPolynomial nthRootAnnullingPolynomial =
                    new RationalPolynomial(nthRootAnnullingPolynomialCoefficients);
                List<MultivariatePolynomial> resolvents = new List<MultivariatePolynomial>();
                int intNMinusOne = (int)n - 1;
                for (int i = 1; i < intNMinusOne; ++i)
                {
                    MultivariatePolynomial resolvent =
                        new MultivariatePolynomial(new RationalPolynomial[] {
                        nMinusFirstRootAnnullingPolynomial, nthRootAnnullingPolynomial });
                    resolvent.SetCoefficient(new int[] { 0, 1 }, Number.One);
                    resolvents.Add(resolvent);
                }
                Integer generatorPower = Number.One;
                List<List<Rational>> resolventMultiplesInTermsOfNMinusFirstRoots =
                    new List<List<Rational>>();
                for (Integer i = Number.One; i < nMinusOne; ++i)
                {
                    List<Rational> resolventMultipleInTermsOfNMinusFirstRoots =
                        new List<Rational> { Number.Zero };
                    generatorPower = (generatorPower * generator).EuclideanDivideBy(n).Remainder;
                    Integer nMinusFirstRootExponent = Number.One;
                    for (int j = 1; j < intNMinusOne; ++j)
                    {
                        resolvents[j - 1].SetCoefficient(new int[] { (int)nMinusFirstRootExponent,
                            (int)generatorPower }, Number.One);
                        resolventMultipleInTermsOfNMinusFirstRoots.Add(Number.Zero);
                        nMinusFirstRootExponent = (nMinusFirstRootExponent + i).EuclideanDivideBy(
                            nMinusOne).Remainder;
                    }
                    resolventMultiplesInTermsOfNMinusFirstRoots.Add(
                        resolventMultipleInTermsOfNMinusFirstRoots);
                }
                MultivariatePolynomial resolventPower =
                    new MultivariatePolynomial(new RationalPolynomial[] {
                    nMinusFirstRootAnnullingPolynomial, nthRootAnnullingPolynomial });
                resolventPower.SetCoefficient(new int[] { 0, 0 }, Number.One);
                List<MultivariatePolynomial> resolventProducts = new List<MultivariatePolynomial>();
                for (int i = resolvents.Count - 1; i >= 0; --i)
                {
                    resolventPower *= resolvents[0];
                    MultivariatePolynomial resolventProduct = resolventPower * resolvents[i];
                    foreach (int[] key in resolventProduct.Coefficients.Keys)
                    {
                        if (key[1] == 0)
                        {
                            resolventMultiplesInTermsOfNMinusFirstRoots[i][key[0]] +=
                                resolventProduct.Coefficients[key];
                        }
                        else if (key[1] == 1)
                        {
                            resolventMultiplesInTermsOfNMinusFirstRoots[i][key[0]] -=
                                resolventProduct.Coefficients[key];
                        }
                    }
                }
                List<Number> nMinusFirstRoots = GetNthRoots(nMinusOne);
                List<Number> resolventProductValues = new List<Number>();
                foreach (List<Rational> coefficients in resolventMultiplesInTermsOfNMinusFirstRoots)
                {
                    Number resolventProductValue = coefficients[0];
                    for (int i = 1; i < intNMinusOne; ++i)
                    {
                        resolventProductValue += coefficients[i] * nMinusFirstRoots[i];
                    }
                    resolventProductValues.Add(resolventProductValue);
                }
                List<Number> resolventValues = new List<Number> { new Integer(-1),
                    resolventProductValues[0].Exponentiate(new Fraction(Number.One, nMinusOne)) };
                for (int i = 1; i < resolventProductValues.Count; ++i)
                {
                    resolventValues.Add(resolventValues[1] / resolventProductValues[i]);
                }
                nthRoots = new List<Number> { Number.One };
                for (int i = 0; i < intNMinusOne; ++i)
                {
                    Number nthRoot = Number.Zero;
                    int nMinusFirstRootExponent = 0;
                    for (int j = 0; j < resolventValues.Count; ++j)
                    {
                        nthRoot += nMinusFirstRoots[nMinusFirstRootExponent] * resolventValues[j];
                        nMinusFirstRootExponent = (nMinusFirstRootExponent + i) % intNMinusOne;
                    }
                    nthRoots.Add(nthRoot / nMinusOne);
                }
                sortRoots(nthRoots, n);
                NthRoots.Add(n, nthRoots);
                return nthRoots;
            }
        }
        public class InvalidUserInput : Exception
        {
            public InvalidUserInput(string message) : base(message)
            { }
        }
        static List<Integer> PrimeList = new List<Integer> { Number.Two, new Integer(3) };
        static IEnumerable<Integer> Primes()
        {
            int i = 0;
            while (true)
            {
                if (i == PrimeList.Count)
                {
                    Integer primeCandidate = PrimeList[PrimeList.Count - 1] + Number.Two;
                    while (true)
                    {
                        bool isDivisible = false;
                        Integer halfCandidate =
                            primeCandidate.EuclideanDivideBy(Number.Two).Quotient;
                        for (int j = 0; PrimeList[j] <= halfCandidate; ++j)
                        {
                            if (primeCandidate.EuclideanDivideBy(
                                PrimeList[j]).Remainder == Number.Zero)
                            {
                                isDivisible = true;
                                break;
                            }
                        }
                        if (!isDivisible)
                        {
                            break;
                        }
                        primeCandidate += Number.Two;
                    }
                    PrimeList.Add(primeCandidate);
                }
                yield return PrimeList[i];
                ++i;
            }
        }
        static Rational Max(Rational a, Rational b)
        {
            if (a > b)
            {
                return a;
            }
            return b;
        }
        static Rational Min(Rational a, Rational b)
        {
            if (a < b)
            {
                return a;
            }
            return b;
        }
        static List<List<T>> GetCartesianProduct<T>(List<List<T>> sets)
        {//The first element of the return value is the list of the first element of each element of
         //the input value. Number.GetConjugates() depends on this fact in order to ensure that the
         //first element of its return value is the instance it was called from.
            List<List<T>> elements = new List<List<T>>();
            void generateElements(int setIndex, List<T> element)
            {
                if (setIndex < sets.Count)
                {
                    foreach (T t in sets[setIndex])
                    {
                        List<T> enlargedElement = new List<T>(element);
                        enlargedElement.Add(t);
                        generateElements(setIndex + 1, enlargedElement);
                    }
                }
                else
                {
                    elements.Add(element);
                }
            }
            generateElements(0, new List<T>());
            return elements;
        }
        static T Exponentiate<T>(T expBase, Integer exponent, Multiplier<T> multiply,
            T multiplicativeIdentity)
        {
            Debug.Assert(exponent >= Number.Zero, "Solver.Exponentiate was called with a " +
                "negative exponent, which the method doesn't account for.");
            T output = multiplicativeIdentity;
            T baseToAPowerOfTwo = expBase;
            while (exponent > Number.Zero)
            {
                Division<Integer> division = exponent.EuclideanDivideBy(Number.Two);
                if (division.Remainder == Number.One)
                {
                    output = multiply(output, baseToAPowerOfTwo);
                }
                baseToAPowerOfTwo = multiply(baseToAPowerOfTwo, baseToAPowerOfTwo);
                exponent = division.Quotient;
            }
            return output;
        }
        static T Exponentiate<T>(T expBase, Integer exponent) where T : IRingElement<T>
        {
            return Exponentiate(expBase, exponent, delegate (T a, T b) { return a.Times(b); },
                expBase.GetMultiplicativeIdentity());
        }
        static T GetGCD<T>(T a, T b) where T : IArithmetic<T>
        {
            T additiveIdentity = a.GetAdditiveIdentity();
            if (b.Equals(additiveIdentity))
            {
                return a;
            }
            T c;
            while (!b.Equals(additiveIdentity))
            {
                c = b;
                b = a.EuclideanDivideBy(b).Remainder;
                a = c;
            }
            return a;
        }
        static T GetGCD<T>(List<T> list) where T : IArithmetic<T>
        {//This procedure has no way to generate an instance of T if list has zero elements; the
         //caller must account for that case instead.
            T GCD = list[0];
            for (int i = 1; i < list.Count; ++i)
            {
                GCD = GetGCD(GCD, list[i]);
            }
            return GCD;
        }
        static ExtendedGCDInfo<T> ExtendedGCD<T>(T a, T b) where T : IArithmetic<T>
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
        abstract class Token
        {
            public virtual bool Equals(char a)
            {
                return false;
            }
        }
        class NumericToken : Token
        {
            public Number Value;
            public NumericToken(Number value)
            {
                Value = value;
            }
        }
        class TextToken : Token
        {
            public char Value;
            public TextToken(char value)
            {
                Value = value;
            }
            public override bool Equals(char a)
            {
                return Value == a;
            }
        }
        static void Main(string[] args)
        {
            Number EvaluateExpression(List<Token> tokens)
            {
                for (int i = 0; i < tokens.Count;)
                {
                    if (tokens[i].Equals(')'))
                    {
                        throw new InvalidUserInput("Unmatched parentheses.");
                    }
                    if (tokens[i].Equals('('))
                    {
                        int numOfUnmatchedParens = 1;
                        int matchingParenIndex = i;
                        while (numOfUnmatchedParens > 0)
                        {
                            matchingParenIndex += 1;
                            if (matchingParenIndex == tokens.Count)
                            {
                                throw new InvalidUserInput("Unmatched parentheses.");
                            }
                            if (tokens[matchingParenIndex].Equals('('))
                            {
                                numOfUnmatchedParens += 1;
                            }
                            else if (tokens[matchingParenIndex].Equals(')'))
                            {
                                numOfUnmatchedParens -= 1;
                            }
                        }
                        tokens[i] = new NumericToken(
                            EvaluateExpression(tokens.GetRange(i + 1, matchingParenIndex - i - 1)));
                        tokens.RemoveRange(i + 1, matchingParenIndex - i);
                    }
                    else
                    {
                        ++i;
                    }
                }
                void executeOperations()
                {
                    int i = 0;
                    while (i < tokens.Count)
                    {
                        if (i == 0 || tokens[i - 1] is TextToken)
                        {
                            if (tokens[i].Equals('+'))
                            {
                                if (tokens.Count < i + 2 || (tokens[i + 1] is TextToken &&
                                    !tokens[i + 1].Equals('+') && !tokens[i + 1].Equals('-')))
                                {
                                    throw new InvalidUserInput("Operator missing operand.");
                                }
                                tokens.RemoveAt(i);
                            }
                            else
                            {
                                ++i;
                            }
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    i = 0;
                    while (i < tokens.Count)
                    {
                        if (tokens[i].Equals('^'))
                        {
                            if (tokens.Count < i + 2)
                            {
                                throw new InvalidUserInput("Operator missing operand.");
                            }
                            if (tokens[i + 1] is NumericToken n && n.Value is Rational exponent)
                            {
                                tokens[i - 1] = new NumericToken(
                                    ((NumericToken)tokens[i - 1]).Value.Exponentiate(exponent));
                            }
                            else
                            {
                                throw new InvalidUserInput(
                                    "The input expression contains an exponentiation\n" +
                                    "whose exponent is not both real and rational;\n" +
                                    "this program doesn't handle transcendental numbers.");
                            }
                            tokens.RemoveRange(i, 2);
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    i = 0;
                    while (i < tokens.Count)
                    {
                        if (i == 0 || tokens[i - 1] is TextToken)
                        {
                            if (tokens[i].Equals('-'))
                            {
                                if (tokens.Count < i + 2)
                                {
                                    throw new InvalidUserInput("Operator missing operand.");
                                }
                                if (tokens[i + 1] is NumericToken n)
                                {
                                    tokens[i + 1] = new NumericToken(-n.Value);
                                    tokens.RemoveAt(i);
                                }
                                else if (tokens[i + 1].Equals('+'))
                                {
                                    tokens[i + 1] = new TextToken('-');
                                    tokens.RemoveAt(i);
                                }
                                else if (tokens[i + 1].Equals('-'))
                                {
                                    tokens[i + 1] = new TextToken('+');
                                    tokens.RemoveAt(i);
                                }
                                else
                                {
                                    throw new InvalidUserInput("Operator missing operand.");
                                }
                            }
                            else
                            {
                                ++i;
                            }
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    i = 0;
                    while (i < tokens.Count)
                    {
                        if (tokens[i].Equals('*'))
                        {
                            if (tokens.Count < i + 2)
                            {
                                throw new InvalidUserInput("Operator missing operand.");
                            }
                            tokens[i - 1] = new NumericToken(((NumericToken)tokens[i - 1]).Value *
                                ((NumericToken)tokens[i + 1]).Value);
                            tokens.RemoveRange(i, 2);
                        }
                        else if (tokens[i].Equals('/'))
                        {
                            if (tokens.Count < i + 2)
                            {
                                throw new InvalidUserInput("Operator missing operand.");
                            }
                            tokens[i - 1] = new NumericToken(((NumericToken)tokens[i - 1]).Value /
                                ((NumericToken)tokens[i + 1]).Value);
                            tokens.RemoveRange(i, 2);
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    i = 0;
                    while (i < tokens.Count)
                    {                        
                        if (tokens[i].Equals('+'))
                        {
                            if (tokens.Count < i + 2)
                            {
                                throw new InvalidUserInput("Operator missing operand.");
                            }
                            tokens[i - 1] = new NumericToken(((NumericToken)tokens[i - 1]).Value +
                                ((NumericToken)tokens[i + 1]).Value);
                            tokens.RemoveRange(i, 2);
                        }
                        else if (tokens[i].Equals('-'))
                        {
                            if (tokens.Count < i + 2)
                            {
                                throw new InvalidUserInput("Operator missing operand.");
                            }
                            tokens[i - 1] = new NumericToken(((NumericToken)tokens[i - 1]).Value -
                                ((NumericToken)tokens[i + 1]).Value);
                            tokens.RemoveRange(i, 2);
                        }
                        else
                        {
                            ++i;
                        }
                    }
                }
                executeOperations();
                return ((NumericToken)tokens[0]).Value;
            }
            void EvaluateString(string input)
            {
                List<Token> tokens = new List<Token>();
                StringBuilder numberCollector = new StringBuilder();
                bool lastCharWasDigit = false;
                try
                {
                    foreach (char c in input)
                    {
                        if (lastCharWasDigit)
                        {
                            if (char.IsDigit(c))
                            {
                                numberCollector.Append(c);
                            }
                            else
                            {
                                if (!"()+-*/^".Contains(c.ToString()))
                                {
                                    throw new InvalidUserInput(c + " is an invalid character.");
                                }
                                else
                                {
                                    tokens.Add(new NumericToken(
                                        Integer.Parse(numberCollector.ToString())));
                                    tokens.Add(new TextToken(c));
                                }
                                lastCharWasDigit = false;
                            }
                        }
                        else if (char.IsDigit(c))
                        {
                            numberCollector = new StringBuilder(c.ToString());
                            lastCharWasDigit = true;
                        }
                        else if (!"()+-*/^".Contains(c.ToString()))
                        {
                            throw new InvalidUserInput(c + " is an invalid character.");
                        }
                        else
                        {
                            tokens.Add(new TextToken(c));
                        }
                    }
                    if (lastCharWasDigit)
                    {
                        tokens.Add(new NumericToken(Integer.Parse(numberCollector.ToString())));
                    }
                    for (int i = 0; i < tokens.Count;)
                    {
                        if (tokens[i].Equals('(') && i > 0 && tokens[i - 1] is NumericToken)
                        {
                            tokens.Insert(i, new TextToken('*'));
                            i += 2;
                        }
                        else if (tokens[i].Equals(')') && i + 1 < tokens.Count &&
                            (tokens[i + 1].Equals('(') || tokens[i + 1] is NumericToken))
                        {
                            tokens.Insert(i + 1, new TextToken('*'));
                            i += 2;
                        }
                        else
                        {
                            ++i;
                        }
                    }
                    Console.Write("=\n" + EvaluateExpression(tokens).ToString() + "\n\n");
                }
                catch (InvalidUserInput e)
                {
                    Console.WriteLine(e.Message);
                }
                catch (DivideByZeroException)
                {
                    Console.WriteLine("Tried to divide by 0.");
                }
            }
#if DEBUG
            EvaluateString("1/(1+2^(1/3))");
#else
            while (true)
            {
                string input = Console.ReadLine();
                if (input != "")
                {
                    if (input[0] == 'q')
                    {
                        return;
                    }
                    EvaluateString(input);
                }
            }
#endif
        }
    }
}
