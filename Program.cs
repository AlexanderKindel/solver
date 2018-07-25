using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Text;

namespace Solver
{
    static class Solver
    {
        struct Interval<T>
        {
            public T Min;
            public T Max;
            public Interval(T min, T max)
            {
                Min = min;
                Max = max;
            }
        }
        struct RectangularEstimate
        {
            public Interval<Float> RealPart;
            public Interval<Float> ImaginaryPart;
            public RectangularEstimate(Interval<Float> realPart, Interval<Float> imaginaryPart)
            {
                RealPart = realPart;
                ImaginaryPart = imaginaryPart;
            }
        }
        struct Division<T>
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
        delegate Interval<Float> PartGetter(Primitive a);
        delegate void PartRefiner(Primitive a, Rational errorIntervalSize);
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
        abstract class Primitive
        {
            public Interval<Float> RealPartEstimate { get; protected set; }
            public Interval<Float> ImaginaryPartEstimate { get; protected set; }
            RationalPolynomial MinimalPolynomial_ = null;
            public RationalPolynomial MinimalPolynomial
            {
                get
                {
                    if (MinimalPolynomial_ == null)
                        MinimalPolynomial_ = CalculateMinimalPolynomial();
                    return MinimalPolynomial_;
                }
            }
            public static Interval<Float> GetRealPartEstimate(Primitive a)
            {
                return a.RealPartEstimate;
            }
            public static Interval<Float> GetImaginaryPartEstimate(Primitive a)
            {
                return a.ImaginaryPartEstimate;
            }
            protected abstract RationalPolynomial CalculateMinimalPolynomial();
            protected void RemoveNonConjugates(List<Number> conjugateCandidates)
            {
                for (int i = 1;
                    conjugateCandidates.Count > MinimalPolynomial.Coefficients.Count - 1;)
                    if (!conjugateCandidates[i].MinimalPolynomial.Equals(MinimalPolynomial))
                        conjugateCandidates.RemoveAt(i);
                    else
                        ++i;
            }
            public List<Number> GetSumConjugates(List<Term> terms)
            {
                List<List<Number>> termConjugates = new List<List<Number>>();
                foreach (Term term in terms)
                    termConjugates.Add(term.GetConjugates());
                List<List<Number>> conjugateCombinations =
                    GetCartesianProduct(termConjugates);
                List<Number> conjugateCandidates = new List<Number>();
                foreach (List<Number> combination in conjugateCombinations)
                {
                    Number conjugate = Number.Zero;
                    foreach (Number number in combination)
                        conjugate += number;
                    conjugateCandidates.Add(conjugate);
                }
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            public abstract List<Number> GetConjugates();//The first conjugate is *this.
            public RationalPolynomial ArgumentInTermsOfThis(Primitive a)
            {
                if ((MinimalPolynomial.Coefficients.Count - 1) %
                    (a.MinimalPolynomial.Coefficients.Count - 1) != 0)
                    return null;
                NestedPolynomial minimalPolynomial = new NestedPolynomial(null,
                    a.MinimalPolynomial).SwitchVariables(MinimalPolynomial);
                List<NestedPolynomial> factors = minimalPolynomial.GetFactors();
                List<RationalPolynomial> candidateFactors = new List<RationalPolynomial>();
                foreach (NestedPolynomial polynomial in factors)
                    if (polynomial.Coefficients.Count == 2)
                        candidateFactors.Add(-polynomial.Coefficients[0]);
                List<Number> conjugates = a.GetConjugates();
                conjugates.RemoveAt(0);
                Rational precision = Number.One;
                for (int i = 0; i < candidateFactors.Count; ++i)
                {
                    List<Number> candidateConjugates = new List<Number>(conjugates);
                    bool factorEqualsA()
                    {
                        Interval<Float> getDistanceInterval(Interval<Float> b, Interval<Float> c)
                        {
                            Float magnitude = (b.Min - c.Min).Magnitude();
                            Interval<Float> output = new Interval<Float>(magnitude, magnitude);
                            void updateBounds(Float boundCandidate)
                            {
                                if (boundCandidate < output.Min)
                                    output.Min = boundCandidate;
                                else if (boundCandidate > output.Max)
                                    output.Max = boundCandidate;
                            }
                            updateBounds((b.Min - c.Max).Magnitude());
                            updateBounds((b.Max - c.Min).Magnitude());
                            updateBounds((b.Max - c.Max).Magnitude());
                            return output;
                        }
                        while (candidateConjugates.Count > 0)
                        {
                            a.RefineRealPartErrorInterval(precision);
                            a.RefineImaginaryPartErrorInterval(precision);
                            foreach (Number conjugate in candidateConjugates)
                            {
                                conjugate.RefineRealPartErrorInterval(precision);
                                conjugate.RefineImaginaryPartErrorInterval(precision);
                            }
                            RectangularEstimate aEstimate =
                                candidateFactors[i].EstimateEvaluation(this, precision);
                            Interval<Float> realDistanceToThis =
                                getDistanceInterval(RealPartEstimate, aEstimate.RealPart);
                            Interval<Float> imaginaryDistanceToThis =
                                getDistanceInterval(ImaginaryPartEstimate, aEstimate.ImaginaryPart);
                            for (int j = 0; j < candidateConjugates.Count;)
                            {
                                Interval<Float> realDistanceToConjugate = getDistanceInterval(
                                    candidateConjugates[j].RealPartEstimate, aEstimate.RealPart);
                                Interval<Float> imaginaryDistanceToConjugate = getDistanceInterval(
                                    candidateConjugates[j].ImaginaryPartEstimate,
                                    aEstimate.ImaginaryPart);
                                if (realDistanceToConjugate.Max < realDistanceToThis.Min ||
                                    imaginaryDistanceToConjugate.Max < imaginaryDistanceToThis.Min)
                                    return false;
                                else if (realDistanceToConjugate.Min > realDistanceToThis.Max ||
                                    imaginaryDistanceToConjugate.Min > imaginaryDistanceToThis.Max)
                                    candidateConjugates.RemoveAt(j);
                                else
                                    ++j;
                            }
                            precision /= Number.Two;
                        }
                        return true;
                    }
                    if (factorEqualsA())
                        return candidateFactors[i];
                }
                return null;
            }
            protected delegate Interval<Rational> RationalEstimateGetter(
                Rational errorIntervalSize);
            protected Interval<Float> GetRefinedErrorInterval(Rational errorIntervalSize,
                RationalEstimateGetter getEstimate)
            {
                Rational placeValue = Float.GetEstimatePlaceValue(errorIntervalSize);
                Interval<Rational> rationalEstimate = getEstimate(placeValue);
                rationalEstimate.Min.RefineRealPartErrorInterval(placeValue);
                rationalEstimate.Max.RefineRealPartErrorInterval(placeValue);
                return new Interval<Float>(rationalEstimate.Min.RealPartEstimate.Min,
                    rationalEstimate.Max.RealPartEstimate.Max);
            }
            protected virtual Interval<Rational> GetRationalRealEstimate(
                Rational errorIntervalSize)
            {
                throw new NotImplementedException("Any derived class that doesn't override " +
                    "GetRationalRealEstimate should override RefineRealPartErrorInterval.");
            }
            protected virtual Interval<Rational> GetRationalImaginaryEstimate(
                Rational errorIntervalSize)
            {
                throw new NotImplementedException("Any derived class that doesn't override " +
                    "GetRationalImaginaryEstimate should override " +
                    "RefineImaginaryPartErrorInterval.");
            }
            public virtual void RefineRealPartErrorInterval(Rational errorIntervalSize)
            {
                if (RealPartEstimate.Min == null || (Rational)(RealPartEstimate.Max -
                    RealPartEstimate.Min) > errorIntervalSize)
                    RealPartEstimate = GetRefinedErrorInterval(errorIntervalSize,
                        GetRationalRealEstimate);
            }
            public static void RefineRealPartErrorInterval(Primitive a, Rational errorIntervalSize)
            {
                a.RefineRealPartErrorInterval(errorIntervalSize);
            }
            public virtual void RefineImaginaryPartErrorInterval(Rational errorIntervalSize)
            {
                if (ImaginaryPartEstimate.Min == null || (Rational)(ImaginaryPartEstimate.Max -
                    ImaginaryPartEstimate.Min) > errorIntervalSize)
                    ImaginaryPartEstimate = GetRefinedErrorInterval(errorIntervalSize,
                        GetRationalImaginaryEstimate);
            }
            public static void RefineImaginaryPartErrorInterval(Primitive a,
              Rational errorIntervalSize)
            {
                a.RefineImaginaryPartErrorInterval(errorIntervalSize);
            }
            Interval<Float> EstimateSumPart<T>(List<T> terms, Rational errorIntervalSize,
                PartGetter getPart, PartRefiner refinePart) where T : Number
            {
                Interval<Float> sumEstimate = new Interval<Float>(Float.Zero, Float.Zero);
                for (int i = 0; i < terms.Count; ++i)
                {
                    refinePart(terms[i], errorIntervalSize / new Integer(terms.Count - i));
                    Interval<Float> part = getPart(terms[i]);
                    sumEstimate.Min += part.Min;
                    sumEstimate.Max += part.Max;
                    errorIntervalSize -= (Rational)(part.Max - part.Min);
                }
                return sumEstimate;
            }
            public Interval<Float> EstimateSumRealPart<T>(List<T> terms, Rational errorIntervalSize)
                where T : Number
            {
                return EstimateSumPart(terms, errorIntervalSize, GetRealPartEstimate,
                    RefineRealPartErrorInterval);
            }
            public Interval<Float> EstimateSumImaginaryPart<T>(List<T> terms,
                Rational errorIntervalSize) where T : Number
            {
                return EstimateSumPart(terms, errorIntervalSize, GetImaginaryPartEstimate,
                    RefineImaginaryPartErrorInterval);
            }
        }
        abstract class Number : Primitive, IArithmetic<Number>, IComparable
        {
            public static Integer Zero = new Integer(0);
            public static Integer One = new Integer(1);
            public static Integer Two = new Integer(2);
            public Interval<Float> ArgumentEstimate { get; protected set; }
            public Interval<Float> MagnitudeEstimate { get; protected set; }
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
            {//Uses a cheap lexicographic scheme that is guaranteed to detect numerical equality for
             //Rationals, but for other derived types is only good for enabling superficial sorting
             //operations. To test for numerical equality of arbitrary Numbers, check whether the
             //minimal polynomial of their difference is x, and to numerically compare unequal
             //Numbers, compare their Estimate fields.
                Number number = obj as Number;
                if (ReferenceEquals(number, null))
                    return 1;
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
                        return true;
                    return false;
                }
                return a.CompareTo(b) == 0;
            }
            public static bool operator !=(Number a, Number b)
            {
                return !(a == b);
            }
            protected abstract Interval<Rational> GetRationalArgumentEstimate(
                Rational errorIntervalSize);
            protected virtual Interval<Rational> GetRationalMagnitudeEstimate(
                Rational errorIntervalSize)
            {
                throw new NotImplementedException("Any derived class that doesn't override " +
                    "GetRationalMagnitudeEstimate should override RefineMagnitudeErrorInterval.");
            }
            public void RefineArgumentErrorInterval(Rational errorIntervalSize)
            {
                if (ArgumentEstimate.Min == null || (Rational)(ArgumentEstimate.Max -
                    ArgumentEstimate.Min) > errorIntervalSize)
                    ArgumentEstimate = GetRefinedErrorInterval(errorIntervalSize,
                        GetRationalArgumentEstimate);
            }
            public virtual void RefineMagnitudeErrorInterval(Rational errorIntervalSize)
            {
                if (MagnitudeEstimate.Min == null || (Rational)(MagnitudeEstimate.Max -
                    MagnitudeEstimate.Min) > errorIntervalSize)
                    MagnitudeEstimate = GetRefinedErrorInterval(errorIntervalSize,
                        GetRationalMagnitudeEstimate);
            }
            protected Interval<Rational> EstimateCosineOfArgument(Rational errorIntervalSize)
            {
                errorIntervalSize /= new Integer(3);
                RefineArgumentErrorInterval(errorIntervalSize);
                Rational delta;
                Rational estimateCosine(Float argument)
                {
                    Rational cosineValue = One;
                    Rational argumentSquared = (Rational)(argument * argument);
                    Integer factorialComponent = Two;
                    delta = -argumentSquared / factorialComponent;
                    while (delta.Magnitude() > errorIntervalSize)
                    {
                        cosineValue += delta;
                        ++factorialComponent;
                        delta *= argumentSquared / factorialComponent;
                        ++factorialComponent;
                        delta /= -factorialComponent;
                    }
                    return cosineValue;
                }
                Rational rationalArgumentEstimateMin = (Rational)ArgumentEstimate.Min;
                while (Pi.LowEstimate <= rationalArgumentEstimateMin &&
                    rationalArgumentEstimateMin <= Pi.HighEstimate)
                    Pi.RefineErrorInterval();
                Rational rationalArgumentEstimateMax = (Rational)ArgumentEstimate.Max;
                while (Pi.LowEstimate <= rationalArgumentEstimateMax &&
                    rationalArgumentEstimateMax <= Pi.HighEstimate)
                    Pi.RefineErrorInterval();
                Interval<Rational> cosine;
                if (rationalArgumentEstimateMin <= Pi.LowEstimate &&
                    Pi.HighEstimate <= rationalArgumentEstimateMax)
                {
                    cosine.Min = -One;
                    cosine.Max = estimateCosine(ArgumentEstimate.Min);
                    if (delta > Zero)
                        cosine.Max += delta;
                    Rational cosineValue = estimateCosine(ArgumentEstimate.Max);
                    if (delta > Zero)
                        cosineValue += delta;
                    if (cosineValue > cosine.Max)
                        cosine.Max = cosineValue;
                }
                else
                {
                    cosine.Min = estimateCosine(ArgumentEstimate.Min);
                    cosine.Max = cosine.Min;
                    if (delta > Zero)
                        cosine.Max += delta;
                    else
                        cosine.Min += delta;
                    Rational cosineValue = estimateCosine(ArgumentEstimate.Max);
                    if (delta > Zero)
                        if (cosineValue < cosine.Min)
                            cosine.Min = cosineValue;
                        else
                        {
                            cosineValue += delta;
                            if (cosineValue > cosine.Max)
                                cosine.Max = cosineValue;
                        }
                    else if (cosineValue > cosine.Max)
                        cosine.Max = cosineValue;
                    else
                    {
                        cosineValue += delta;
                        if (cosineValue < cosine.Min)
                            cosine.Min = cosineValue;
                    }
                }
                return cosine;
            }
            protected Interval<Rational> EstimateSineOfArgument(Rational errorIntervalSize)
            {
                errorIntervalSize /= new Integer(3);
                RefineArgumentErrorInterval(errorIntervalSize);
                Rational delta;
                Rational estimateSine(Rational argument)
                {
                    Rational sineValue = argument;
                    Rational argumentSquared = argument * argument;
                    Integer factorialComponent = new Integer(3);
                    delta = -argumentSquared * argument / new Integer(6);
                    while (delta.Magnitude() > errorIntervalSize)
                    {
                        sineValue += delta;
                        ++factorialComponent;
                        delta *= argumentSquared / factorialComponent;
                        ++factorialComponent;
                        delta /= -factorialComponent;
                    }
                    return sineValue;
                }
                Rational rationalArgumentEstimateMin = (Rational)ArgumentEstimate.Min;
                while (Pi.LowEstimate / Two <= rationalArgumentEstimateMin &&
                    rationalArgumentEstimateMin <= Pi.HighEstimate / Two)
                    Pi.RefineErrorInterval();
                Rational rationalArgumentEstimateMax = (Rational)ArgumentEstimate.Max;
                while (Pi.LowEstimate / Two <= rationalArgumentEstimateMax &&
                    rationalArgumentEstimateMax <= Pi.HighEstimate / Two)
                    Pi.RefineErrorInterval();
                Interval<Rational> sine;
                if (rationalArgumentEstimateMin <= Pi.LowEstimate / Two &&
                    Pi.HighEstimate / Two <= rationalArgumentEstimateMax)
                {
                    sine.Max = One;
                    sine.Min = estimateSine(rationalArgumentEstimateMin);
                    if (delta < Zero)
                        sine.Min += delta;
                    Rational sineValue = estimateSine(rationalArgumentEstimateMax);
                    if (delta < Zero)
                        sineValue += delta;
                    if (sineValue < sine.Min)
                        sine.Min = sineValue;
                }
                else
                {
                    Fraction threeOverTwo = new Fraction(new Integer(3), Two);
                    while (threeOverTwo * Pi.LowEstimate <= rationalArgumentEstimateMin &&
                        rationalArgumentEstimateMin <= threeOverTwo * Pi.HighEstimate)
                        Pi.RefineErrorInterval();
                    while (threeOverTwo * Pi.LowEstimate <= rationalArgumentEstimateMax &&
                        rationalArgumentEstimateMax <= threeOverTwo * Pi.LowEstimate)
                        Pi.RefineErrorInterval();
                    if (rationalArgumentEstimateMin <= threeOverTwo * Pi.LowEstimate &&
                        threeOverTwo * Pi.LowEstimate <= rationalArgumentEstimateMax)
                    {
                        sine.Min = -One;
                        sine.Max = estimateSine(rationalArgumentEstimateMin);
                        if (delta > Zero)
                            sine.Max += delta;
                        Rational sineValue = estimateSine(rationalArgumentEstimateMax);
                        if (delta > Zero)
                            sineValue += delta;
                        if (sineValue > sine.Max)
                            sine.Max = sineValue;
                    }
                    else
                    {
                        sine.Min = estimateSine(rationalArgumentEstimateMin);
                        sine.Max = sine.Min;
                        if (delta > Zero)
                            sine.Max += delta;
                        else
                            sine.Min += delta;
                        Rational sineValue = estimateSine(rationalArgumentEstimateMax);
                        if (delta > Zero)
                            if (sineValue < sine.Min)
                                sine.Min = sineValue;
                            else
                            {
                                sineValue += delta;
                                if (sineValue > sine.Max)
                                    sine.Max = sineValue;
                            }
                        else if (sineValue > sine.Max)
                            sine.Max = sineValue;
                        else
                        {
                            sineValue += delta;
                            if (sineValue < sine.Min)
                                sine.Min = sineValue;
                        }
                    }
                }
                return sine;
            }
            protected delegate Interval<Rational> TrigFunction(Rational errorIntervalSize);
            protected Interval<Rational> EstimateRectangularPartFromPolarForm(
                Rational errorIntervalSize, TrigFunction sineOrCosine)
            {
                RefineMagnitudeErrorInterval(One);
                errorIntervalSize /= (Rational)MagnitudeEstimate.Max.Magnitude() + Two;
                RefineMagnitudeErrorInterval(errorIntervalSize);
                Interval<Rational> trigValue = sineOrCosine(errorIntervalSize);
                if (trigValue.Min >= Zero)
                    return new Interval<Rational>(trigValue.Min * (Rational)MagnitudeEstimate.Min,
                        trigValue.Max * (Rational)MagnitudeEstimate.Max);
                else if (trigValue.Max <= Zero)
                    return new Interval<Rational>(trigValue.Min * (Rational)MagnitudeEstimate.Max,
                        trigValue.Max * (Rational)MagnitudeEstimate.Min);
                return new Interval<Rational>(trigValue.Min * (Rational)MagnitudeEstimate.Max,
                    trigValue.Max * (Rational)MagnitudeEstimate.Max);
            }
            public abstract override string ToString();
        }
        abstract class Term : Number
        {
            public RationalPolynomial ThisInTermsOfParentSumPrimitiveElement = null;
            public abstract Term Copy();
            public override Number Plus(Number a)
            {
                if (a == Zero)
                    return this;
                RationalPolynomial x = new RationalPolynomial(null, Zero, One);
                if (a is Rational rational)
                {
                    Term termA = Copy();
                    Term termB = rational.Copy();
                    termA.ThisInTermsOfParentSumPrimitiveElement = x;
                    termB.ThisInTermsOfParentSumPrimitiveElement =
                        new RationalPolynomial(null, rational);
                    return new LinearlyIndependentSum(new List<Term> { termA, termB }, this);
                }
                if (a is Term term)
                {
                    Term termA;
                    Term termB;
                    RationalPolynomial termInTermsOfThis = ArgumentInTermsOfThis(term);
                    if (termInTermsOfThis != null)
                    {
                        if (termInTermsOfThis.Coefficients.Count == 2 &&
                            termInTermsOfThis.Coefficients[0] == Zero)
                            return (termInTermsOfThis.Coefficients[1] + One) * this;
                        termA = Copy();
                        termB = term.Copy();
                        termA.ThisInTermsOfParentSumPrimitiveElement = x;
                        termB.ThisInTermsOfParentSumPrimitiveElement = termInTermsOfThis;
                        return new LinearlyIndependentSum(new List<Term> { termA, termB }, this);
                    }
                    termA = Copy();
                    termB = term.Copy();
                    termA.ThisInTermsOfParentSumPrimitiveElement = null;
                    termB.ThisInTermsOfParentSumPrimitiveElement = null;
                    return new LinearlyIndependentSum(term, this);
                }
                return a + this;
            }
            protected abstract Term TermTimes(Rational a);
            public static Term operator *(Term a, Rational b)
            {
                return a.TermTimes(b);
            }
            public static Term operator *(Rational a, Term b)
            {
                return b.TermTimes(a);
            }
        }
        abstract class Rational : Term, IArithmetic<Rational>
        {
            public abstract Integer Numerator { get; }
            public abstract Integer Denominator { get; }
            protected static Rational Create(Integer numerator, Integer denominator)
            {
                if (denominator == Zero)
                    throw new DivideByZeroException();
                if (numerator == Zero)
                    return Zero;
                ExtendedGCDInfo<Integer> extendedGCD = ExtendedGCD(numerator, denominator);
                numerator = extendedGCD.AOverGCD;
                denominator = extendedGCD.BOverGCD;
                if (denominator < Zero)
                {
                    numerator = -numerator;
                    denominator = -denominator;
                }
                if (denominator == One)
                    return numerator;
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
                    return this;
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
            public Rational Minus(Rational a)
            {
                return this - a;
            }
            protected override Term TermTimes(Rational a)
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
            protected override Interval<Rational> GetRationalImaginaryEstimate(
                Rational errorIntervalSize)
            {
                return new Interval<Rational>(Zero, Zero);
            }
            protected override Interval<Rational> GetRationalArgumentEstimate(
                Rational errorIntervalSize)
            {
                if (Numerator > Zero)
                    return new Interval<Rational>(Zero, Zero);
                Pi.RefineErrorInterval(errorIntervalSize);
                return new Interval<Rational>(Pi.LowEstimate, Pi.HighEstimate);
            }

            //Treats str as the string representation of a numerical constant, and returns the
            //string representation of the reference object multiplied by that constant.
            public abstract string InsertString(string str);
        }
        class Integer : Rational, IArithmetic<Integer>
        {
            uint[] Values;
            sbyte Sign;
            public override Integer Numerator { get => this; }
            public override Integer Denominator { get => One; }
            Integer(uint[] values, sbyte sign)
            {
                int lastNonzeroValueIndex = -1;
                for (int i = values.Length - 1; i >= 0; --i)
                    if (values[i] != 0)
                    {
                        lastNonzeroValueIndex = i;
                        break;
                    }
                Values = new uint[lastNonzeroValueIndex + 1];
                for (int i = 0; i <= lastNonzeroValueIndex; ++i)
                    Values[i] = values[i];
                if (lastNonzeroValueIndex < 0)
                    Sign = 0;
                else
                    Sign = sign;
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
            public override Term Copy()
            {
                Integer copy = new Integer(Values, Sign);
                copy.ArgumentEstimate = ArgumentEstimate;
                return new Integer(Values, Sign);
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
                    output += multiplier *
                        new Integer(int.Parse(decimalString.Substring(0, startIndex + 9)));
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
                    return new Integer(Values, 1);
                return this;
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
                Integer inverse = ExtendedGCD(this, characteristic).ACoefficient;
                if (inverse < Zero)
                    inverse += characteristic;
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
                                    sumValues[i] += 0x80000000;
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
                                    sumValues[i] += 0x80000000;
                            }
                            else
                                sumValues[i] += aValues[i] + bValues[i];
                        }
                    }
                }
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
                if (Values.Length > a.Values.Length)
                    aValues = new uint[Values.Length];
                else
                    aValues = new uint[a.Values.Length];
                Values.CopyTo(aValues, 0);
                bValues = new uint[aValues.Length];
                a.Values.CopyTo(bValues, 0);
                sumValues = new uint[aValues.Length + 1];
                Integer handleSingleNegativeCase()
                {
                    sumValues[aValues.Length] = 0xffffffff;
                    calculateValues();
                    if (sumValues[aValues.Length] == 0)
                        return new Integer(sumValues, 1);
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
                    return integer + this;
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
                    for (int j = 0; j < a.Values.Length; ++j)
                    {
                        ulong productComponent = (ulong)(Values[i]) * a.Values[j];
                        product = product + new Integer(ShiftValuesLeft(new uint[] {
                            (uint)(productComponent & 0x00000000ffffffff),
                            (uint)((productComponent & 0xffffffff00000000) >> 32) }, i + j, 0), 1);
                    }
                product.Sign = (sbyte)(Sign * a.Sign);
                return product;
            }
            public override Number Times(Number a)
            {
                if (a is Integer integer)
                    return integer * this;
                return a * this;
            }
            public Division<Integer> EuclideanDivideBy(Integer divisor)
            {
                if (divisor.Sign == 0)
                    throw new DivideByZeroException();
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
                    calculateValue(i, 1);
                calculateValue(positiveDivisor.Values.Length - 1, divisorLeadingDigitPlace);
                division.Quotient = new Integer(quotient, (sbyte)(Sign * divisor.Sign)).ShiftLeft(
                    1 - positiveDivisor.Values.Length, 1 - divisorLeadingDigitPlace);
                if (division.Remainder.Sign != 0)
                    division.Remainder.Sign = Sign;
                return division;
            }
            public override Number Exponentiate(Rational exponent)
            {
                if (this == Zero)
                    return Zero;
                if (exponent.Numerator < Zero)
                    return Reciprocal().Exponentiate(-exponent);
                Integer power = Solver.Exponentiate(this, exponent.Numerator);
                if (exponent.Denominator == One)
                    return power;
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
                if (power < Zero)
                {
                    power = takeRootOfPositiveRadicand(-power);
                    if (exponent.Denominator.EuclideanDivideBy(Two).Remainder == Zero)
                        return new Surd(coefficient, -power, exponent.Denominator);
                    coefficient = -coefficient;
                }
                else
                    power = takeRootOfPositiveRadicand(power);
                if (power == One)
                    return coefficient;
                return new Surd(coefficient, power, exponent.Denominator);
            }
            protected override int GetTypeIndex()
            {
                return 1;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                return (this - (Integer)obj).Sign;
            }
            public override int GetHashCode()
            {
                uint output = 0;
                foreach (uint value in Values)
                    output ^= value;
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
            public override void RefineRealPartErrorInterval(Rational errorIntervalSize)
            {
                Float floatThis = new Float(this, Zero);
                RealPartEstimate = new Interval<Float>(floatThis, floatThis);
            }
            public override void RefineMagnitudeErrorInterval(Rational errorIntervalSize)
            {
                Float magnitude = new Float(Magnitude(), Zero);
                MagnitudeEstimate = new Interval<Float>(magnitude, magnitude);
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
            public Integer ThisChooseK(Integer k)
            {
                Integer numerator = One;
                for (Integer i = this - k + One; i <= this; ++i)
                    numerator *= i;
                Integer denominator = One;
                for (Integer i = Two; i <= k; ++i)
                    denominator *= i;
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
            public override String InsertString(String str)
            {
                if (this != One)
                    if (this == -One)
                        return '-' + str;
                    else if (str[0] == '(')
                        return ToString() + str;
                    else
                        return ToString() + '*' + str;
                return str;
            }
            public override string ToString()
            {
                if (this == Zero)
                    return "0";
                StringBuilder output = new StringBuilder();
                Integer quotient = this;
                Integer power = new Integer(10);
                while (quotient.Sign != 0)
                {
                    Division<Integer> division = quotient.EuclideanDivideBy(power);
                    quotient = division.Quotient;
                    if (division.Remainder.Sign != 0)
                        output.Insert(0, division.Remainder.Values[0]);
                    else
                        output.Insert(0, '0');
                }
                if (this < Zero)
                    output.Insert(0, '-');
                return output.ToString();
            }
        }
        class Fraction : Rational
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
            public override Term Copy()
            {
                Fraction copy = new Fraction(Numerator, Denominator);
                copy.ArgumentEstimate = ArgumentEstimate;
                return copy;
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
                    return Create(integer * Denominator + Numerator, Denominator);
                if (a is Fraction fraction)
                    return Create(Numerator * fraction.Denominator +
                        fraction.Numerator * Denominator, Denominator * fraction.Denominator);
                return a + this;
            }
            public override Number Times(Number a)
            {
                if (a is Integer integer)
                    return Create(Numerator * integer, Denominator);
                if (a is Fraction fraction)
                    return Create(Numerator * fraction.Numerator,
                        Denominator * fraction.Denominator);
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
                    return comparison;
                Fraction fraction = (Fraction)obj;
                comparison = Numerator.CompareTo(fraction.Numerator);
                if (comparison != 0)
                    return comparison;
                return Denominator.CompareTo(fraction.Denominator);
            }
            public override int GetHashCode()
            {
                return Numerator.GetHashCode() ^ Denominator.GetHashCode();
            }
            public override void RefineRealPartErrorInterval(Rational errorIntervalSize)
            {
                RefineMagnitudeErrorInterval(errorIntervalSize);
                if (Numerator < Zero)
                    RealPartEstimate =
                        new Interval<Float>(-MagnitudeEstimate.Max, -MagnitudeEstimate.Min);
                else
                    RealPartEstimate = MagnitudeEstimate;
            }
            public override void RefineMagnitudeErrorInterval(Rational errorIntervalSize)
            {
                Division<Integer> division;
                Integer numeratorMagnitude = Numerator.Magnitude();
                if (EstimateNumerator == null)
                {
                    division = numeratorMagnitude.EuclideanDivideBy(Denominator);
                    EstimateNumerator = division.Quotient;
                    EstimateRemainder = division.Remainder;
                }
                else if ((Rational)(MagnitudeEstimate.Max - MagnitudeEstimate.Min) <
                    errorIntervalSize)
                    return;
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
                        MagnitudeEstimate = new Interval<Float>(estimate, estimate);
                        return;
                    }
                }
                estimate = Float.GetReducedForm(EstimateNumerator, RadixPointPosition);
                Integer crossMultiplicationA = EstimateNumerator * Denominator;
                Integer crossMultiplicationB = EstimateDenominator * numeratorMagnitude;
                if (crossMultiplicationA > crossMultiplicationB)
                    MagnitudeEstimate = new Interval<Float>(Float.GetReducedForm(
                        EstimateNumerator - One, RadixPointPosition), estimate);
                else
                    MagnitudeEstimate = new Interval<Float>(estimate,
                        Float.GetReducedForm(EstimateNumerator + One, RadixPointPosition));
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
        class Surd : Term
        {
            public Rational Coefficient { get; }
            public Number Radicand { get; }
            public Integer Index { get; }
            public Surd(Rational coefficient, Number radicand, Integer index)
            {//To be called only when it's known a priori that Index > One, radicand is not an
             //Integer that is a perfect indexth power and, for odd index,
             //radicand.RationalFactor() == One, or even index,
             //radicand.RationalFactor().Magnitude() == One. Otherwise, use radicand.Exponentiate().
                Coefficient = coefficient;
                Radicand = radicand;
                Index = index;
            }
            public override Term Copy()
            {
                Surd copy = new Surd(Coefficient, Radicand, Index);
                copy.RealPartEstimate = RealPartEstimate;
                copy.ImaginaryPartEstimate = ImaginaryPartEstimate;
                copy.ArgumentEstimate = ArgumentEstimate;
                copy.MagnitudeEstimate = MagnitudeEstimate;
                return copy;
            }
            public override Number Reciprocal()
            {
                return Radicand.Exponentiate(new Fraction(Index - One, Index)) /
                    (Coefficient * Radicand);
            }
            public override Rational RationalFactor()
            {
                return Coefficient;
            }
            protected override Term TermTimes(Rational a)
            {
                if (a == Zero)
                    return Zero;
                if (a == One)
                    return this;
                return new Surd(a * Coefficient, Radicand, Index);
            }
            public override Number Times(Number a)
            {
                if (a is Rational rational)
                    return TermTimes(rational);
                if (a is Surd surd)
                {
                    ExtendedGCDInfo<Integer> GCDInfo = ExtendedGCD(Index, surd.Index);
                    Number product = Coefficient * surd.Coefficient * (Radicand.Exponentiate(
                        GCDInfo.BOverGCD.Magnitude()) * surd.Radicand.Exponentiate(
                        GCDInfo.AOverGCD.Magnitude())).Exponentiate(One / ((Index * surd.Index) /
                        GCDInfo.GCD).Numerator);
                    Rational errorIntervalSize = Pi.LowEstimate / Index;
                    product.RefineArgumentErrorInterval(errorIntervalSize);
                    errorIntervalSize /= new Integer(4) * Pi.HighEstimate;
                    RefineArgumentErrorInterval(errorIntervalSize);
                    surd.RefineArgumentErrorInterval(errorIntervalSize);
                    if (product.ArgumentEstimate.Max <
                        ArgumentEstimate.Min * surd.ArgumentEstimate.Min)
                        product *= RootsOfUnity.GetNthRoots(Index)[1];
                    return product;
                }
                return a * this;
            }
            public override Number Exponentiate(Rational exponent)
            {
                return Coefficient.Exponentiate(exponent) * Radicand.Exponentiate(exponent / Index);
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                Surd surd = (Surd)obj;
                comparison = Index.CompareTo(surd.Index);
                if (comparison != 0)
                    return comparison;
                comparison = Coefficient.CompareTo(surd.Coefficient);
                if (comparison != 0)
                    return comparison;
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
                        coefficients.Add(Zero);
                    coefficients.Add(radicandMinimalPolynomial.Coefficients[i]);
                }
                List<IntegerPolynomial> candidates =
                    new RationalPolynomial(coefficients).GetPrimitivePart().GetFactors();
                if (candidates.Count == 1)
                    return new RationalPolynomial(GetMonicCoefficients(candidates[0].Coefficients)).
                        GetMinimalPolynomialOfKTimesRootOfThis(Coefficient);
                List<RationalPolynomial> rationalCandidates = new List<RationalPolynomial>();
                foreach (IntegerPolynomial candidate in candidates)
                    rationalCandidates.Add(candidate);
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
                                return rationalCandidates[0].GetMinimalPolynomialOfKTimesRootOfThis(
                                    Coefficient);
                        }
                        else
                            ++i;
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
                    foreach (Number root in rootsOfUnity)
                        conjugateCandidates.Add(new Surd(Coefficient, conjugate, Index) * root);
                RemoveNonConjugates(conjugateCandidates);
                return conjugateCandidates;
            }
            protected override int GetTypeIndex()
            {
                return 3;
            }
            protected override Interval<Rational> GetRationalRealEstimate(
                Rational errorIntervalSize)
            {
                return EstimateRectangularPartFromPolarForm(errorIntervalSize,
                    EstimateCosineOfArgument);
            }
            protected override Interval<Rational> GetRationalImaginaryEstimate(
                Rational errorIntervalSize)
            {
                return EstimateRectangularPartFromPolarForm(errorIntervalSize,
                    EstimateSineOfArgument);
            }
            protected override Interval<Rational> GetRationalArgumentEstimate(
                Rational errorIntervalSize)
            {
                if (Coefficient > Zero)
                {
                    errorIntervalSize *= Index;
                    Radicand.RefineArgumentErrorInterval(errorIntervalSize);
                    return new Interval<Rational>((Rational)Radicand.ArgumentEstimate.Min / Index,
                        (Rational)Radicand.ArgumentEstimate.Max / Index);
                }
                errorIntervalSize /= Two;
                Coefficient.RefineArgumentErrorInterval(errorIntervalSize);
                Radicand.RefineArgumentErrorInterval(Index * errorIntervalSize);
                return new Interval<Rational>((Rational)Coefficient.ArgumentEstimate.Min +
                    (Rational)Radicand.ArgumentEstimate.Min / Index,
                    (Rational)Coefficient.ArgumentEstimate.Max +
                    (Rational)Radicand.ArgumentEstimate.Max / Index);
            }
            public override void RefineMagnitudeErrorInterval(Rational errorIntervalSize)
            {
                Radicand.RefineMagnitudeErrorInterval(One);
                if (Radicand.MagnitudeEstimate.Max <= Float.One)
                    errorIntervalSize /= Coefficient.Magnitude() + Two;
                else
                    errorIntervalSize /=
                        Coefficient.Magnitude() + (Rational)Radicand.MagnitudeEstimate.Max + One;
                Coefficient.RefineMagnitudeErrorInterval(errorIntervalSize);
                Integer three = new Integer(3);
                Radicand.RefineMagnitudeErrorInterval(
                    Solver.Exponentiate(errorIntervalSize, Index) / three);
                errorIntervalSize /= three;                
                MagnitudeEstimate = new Interval<Float>(Coefficient.MagnitudeEstimate.Min *
                    Radicand.MagnitudeEstimate.Min.EstimateRoot(errorIntervalSize, Index).Min,
                    Coefficient.MagnitudeEstimate.Max *
                    Radicand.MagnitudeEstimate.Max.EstimateRoot(errorIntervalSize, Index).Max);
            }
            public override string ToString()
            {
                if (Radicand is Integer integer && integer > Zero)
                    return Coefficient.InsertString(Radicand + "^(1/" + Index + ')');
                return Coefficient.InsertString("(" + Radicand + ')' + "^(1/" + Index + ')');
            }
        }
        class PolynomialIndependentSum : Primitive
        {
            public List<Term> Terms { get; }
            public PolynomialIndependentSum(List<Term> terms)
            {
                Terms = terms;
                Terms.Sort();
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                MultivariatePolynomial variableForm;
                if (Terms[0] is Rational rational)
                {
                    RationalPolynomial[] minimalPolynomials =
                        new RationalPolynomial[Terms.Count - 1];
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
                    RationalPolynomial[] minimalPolynomials = new RationalPolynomial[Terms.Count];
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
                return variableForm.GetMinimalPolynomial();
            }
            public override List<Number> GetConjugates()
            {
                return GetSumConjugates(Terms);
            }
            public override void RefineRealPartErrorInterval(Rational errorIntervalSize)
            {
                RealPartEstimate = EstimateSumRealPart(Terms, errorIntervalSize);
            }
            public override void RefineImaginaryPartErrorInterval(Rational errorIntervalSize)
            {
                ImaginaryPartEstimate = EstimateSumImaginaryPart(Terms, errorIntervalSize);
            }
        }
        static Interval<Float> GetMagnitudeInterval(Interval<Float> interval)
        {
            Float magnitudeA = interval.Min.Magnitude();
            Float magnitudeB = interval.Max.Magnitude();
            if (magnitudeA < magnitudeB)
                return new Interval<Float>(magnitudeA, magnitudeB);
            return new Interval<Float>(magnitudeB, magnitudeA);
        }
        class LinearlyIndependentSum : Number
        {
            public List<Term> Terms { get; }
            Term PrimitiveCandidate;
            Primitive PrimitiveElement_ = null;
            public Primitive PrimitiveElement
            {
                get
                {
                    if (PrimitiveElement_ != null)
                        return PrimitiveElement_;
                    bool candidateIsPrimitive(Term candidate, Term otherTerm)
                    {
                        RationalPolynomial otherTermInTermsOfCandidate =
                            candidate.ArgumentInTermsOfThis(otherTerm);
                        if (otherTermInTermsOfCandidate != null)
                        {
                            candidate.ThisInTermsOfParentSumPrimitiveElement =
                                new RationalPolynomial(null, Zero, One);
                            otherTerm.ThisInTermsOfParentSumPrimitiveElement =
                                otherTermInTermsOfCandidate;
                            PrimitiveElement_ = candidate;
                            return true;
                        }
                        return false;
                    }
                    if (Terms[0] == PrimitiveCandidate &&
                        candidateIsPrimitive(Terms[0], Terms[1]))
                        return PrimitiveElement_;
                    else if (candidateIsPrimitive(Terms[1], Terms[0]))
                        return PrimitiveElement_;
                    RationalPolynomial x = new RationalPolynomial(null, Zero, One);
                    Integer k = One;
                    bool primitiveElementFound(Term termA, Term termB)
                    {
                        PrimitiveElement_ =
                            new PolynomialIndependentSum(new List<Term> { k * termA, termB });
                        RationalPolynomial termAInTermsOfPrimitive =
                            PrimitiveElement_.ArgumentInTermsOfThis(termA);
                        if (termAInTermsOfPrimitive != null)
                        {
                            termA.ThisInTermsOfParentSumPrimitiveElement = termAInTermsOfPrimitive;
                            termB.ThisInTermsOfParentSumPrimitiveElement =
                                x - k * termAInTermsOfPrimitive;
                            return true;
                        }
                        return false;
                    }
                    while (true)
                    {
                        if (primitiveElementFound(Terms[0], Terms[1]))
                            return PrimitiveElement_;
                        if (primitiveElementFound(Terms[1], Terms[0]))
                            return PrimitiveElement_;
                        ++k;
                    }
                }
            }
            RationalPolynomial VariableForm_ = null;
            RationalPolynomial VariableForm
            {
                get
                {
                    if (VariableForm_ == null)
                    {
                        List<Rational> variableFormCoefficients = new List<Rational>();
                        foreach (Term term in Terms)
                            variableFormCoefficients = CoefficientAdd(variableFormCoefficients,
                                term.ThisInTermsOfParentSumPrimitiveElement.Coefficients);
                        VariableForm_ = new RationalPolynomial(variableFormCoefficients,
                            PrimitiveElement.MinimalPolynomial);
                    }
                    return VariableForm_;
                }
            }
            public LinearlyIndependentSum(List<Term> terms, Primitive primitiveElement)
            {
                Terms = terms;
                Terms.Sort();
                PrimitiveElement_ = primitiveElement;
            }
            public LinearlyIndependentSum(Term primitiveCandidateTerm, Term otherTerm)
            {
                Terms = new List<Term> { primitiveCandidateTerm, otherTerm };
                Terms.Sort();
                PrimitiveCandidate = primitiveCandidateTerm;
            }
            public override Number Plus(Number a)
            {
                if (a is Term term)
                {
                    Number thisPlusTerm(Primitive newPrimitiveElement,
                        List<Term> termsInTermsOfNewPrimitiveElement)
                    {
                        List<RationalPolynomial> termsInTermsOfPrimitiveElement =
                            new List<RationalPolynomial>();
                        foreach (Term t in termsInTermsOfNewPrimitiveElement)
                            termsInTermsOfPrimitiveElement.Add(
                                t.ThisInTermsOfParentSumPrimitiveElement);
                        Matrix matrix = new Matrix(termsInTermsOfPrimitiveElement).GetTranspose();
                        List<Rational> augmentation = new List<Rational>(
                            term.ThisInTermsOfParentSumPrimitiveElement.Coefficients);
                        while (augmentation.Count < termsInTermsOfPrimitiveElement.Count)
                            augmentation.Add(Zero);
                        matrix.GetRowEchelonForm(augmentation);
                        List<Term> outputTerms;
                        for (int i = Terms.Count; i < augmentation.Count; ++i)
                            if (augmentation[i] != Zero)
                            {
                                outputTerms = new List<Term>(Terms);
                                outputTerms.Add(term);
                                return new LinearlyIndependentSum(outputTerms, PrimitiveElement);
                            }
                        List<Rational> linearCombinationCoefficients = new List<Rational>();
                        for (int i = Terms.Count - 1; i >= 0; --i)
                        {
                            Rational nextCoefficient = augmentation[i];
                            for (int j = 0; j < linearCombinationCoefficients.Count; ++j)
                                nextCoefficient -= matrix.Rows[i][Terms.Count - 1 - j] *
                                    linearCombinationCoefficients[j];
                            linearCombinationCoefficients.Add(nextCoefficient / matrix.Rows[i][i]);
                        }
                        outputTerms = new List<Term>();
                        Integer minusOne = -One;
                        for (int i = 0; i < Terms.Count; ++i)
                            if (linearCombinationCoefficients[Terms.Count - 1 - i] != minusOne)
                            {
                                Term outputTerm = Terms[i] *
                                    (One - linearCombinationCoefficients[Terms.Count - 1 - i]);
                                outputTerm.ThisInTermsOfParentSumPrimitiveElement =
                                    Terms[i].ThisInTermsOfParentSumPrimitiveElement;
                                outputTerms.Add(outputTerm);
                            }
                        if (outputTerms.Count == 0)
                            return Zero;
                        if (outputTerms.Count == 1)
                            return outputTerms[0];
                        return new LinearlyIndependentSum(outputTerms, PrimitiveElement);
                    }
                    RationalPolynomial termInTermsOfPrimitiveElement =
                        PrimitiveElement.ArgumentInTermsOfThis(term);
                    if (termInTermsOfPrimitiveElement != null)
                    {
                        term.ThisInTermsOfParentSumPrimitiveElement = termInTermsOfPrimitiveElement;
                        return thisPlusTerm(PrimitiveElement, Terms);
                    }
                    RationalPolynomial x = new RationalPolynomial(null, Zero, One);
                    List<Term> convertTermsToNewPrimitiveElement(
                        RationalPolynomial newPrimitiveElementInTermsOfOld)
                    {
                        List<Term> convertedTerms = new List<Term>();
                        List<RationalPolynomial> powers = new List<RationalPolynomial> {
                            newPrimitiveElementInTermsOfOld.GetMultiplicativeIdentity() };
                        for (int i = 0; i < Terms.Count; ++i)
                        {
                            Term copiedTerm = Terms[i].Copy();
                            copiedTerm.ThisInTermsOfParentSumPrimitiveElement =
                                newPrimitiveElementInTermsOfOld.GetAdditiveIdentity();
                            for (int j = 0; j < Terms[i].
                                ThisInTermsOfParentSumPrimitiveElement.Coefficients.Count; ++j)
                            {
                                if (j >= powers.Count)
                                    powers.Add(powers[powers.Count - 1] *
                                        newPrimitiveElementInTermsOfOld);
                                copiedTerm.ThisInTermsOfParentSumPrimitiveElement += powers[j] *
                                    Terms[i].ThisInTermsOfParentSumPrimitiveElement.Coefficients[j];
                            }
                            convertedTerms.Add(copiedTerm);
                        }
                        return convertedTerms;
                    }
                    RationalPolynomial primitiveElementInTermsOfTerm =
                        term.ArgumentInTermsOfThis(PrimitiveElement);
                    if (primitiveElementInTermsOfTerm != null)
                    {
                        term.ThisInTermsOfParentSumPrimitiveElement = x;
                        return thisPlusTerm(term,
                            convertTermsToNewPrimitiveElement(primitiveElementInTermsOfTerm));
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
                            Term copiedTerm = term.Copy();
                            copiedTerm.ThisInTermsOfParentSumPrimitiveElement =
                                termInTermsOfPrimitive;
                            RationalPolynomial oldPrimitiveInTermsOfNewPrimitive =
                                x - k * termInTermsOfPrimitive;
                            return thisPlusTerm(copiedTerm, convertTermsToNewPrimitiveElement(
                                oldPrimitiveInTermsOfNewPrimitive));
                        }
                        ++k;
                    }
                }
                LinearlyIndependentSum sum = (LinearlyIndependentSum)a;
                Number output = this;
                foreach (Term t in sum.Terms)
                    output = output + t;
                return output;
            }
            public override Number Times(Number a)
            {
                Number output = Zero;
                if (a is Term)
                {
                    foreach (Term term in Terms)
                        output += term * a;
                    return output;
                }
                LinearlyIndependentSum sum = (LinearlyIndependentSum)a;
                foreach (Term multiplicandTerm in Terms)
                    foreach (Term multiplierTerm in sum.Terms)
                        output += multiplicandTerm * multiplierTerm;
                return output;
            }
            public override Number Reciprocal()
            {
                List<Number> conjugates = GetConjugates();
                conjugates.Remove(this);
                Number numerator = One;
                foreach (Number conjugate in conjugates)
                    numerator *= conjugate;
                List<List<Rational>> matrixRows = new List<List<Rational>> { new List<Rational>() };
                for (int i = 0; i < MinimalPolynomial.Coefficients.Count - 2; ++i)
                    matrixRows[0].Add(Zero);
                matrixRows[0].Add(-MinimalPolynomial.Coefficients[0]);
                for (int i = 1; i < MinimalPolynomial.Coefficients.Count - 1; ++i)
                {
                    List<Rational> row = new List<Rational>();
                    for (int j = 0; j < MinimalPolynomial.Coefficients.Count - 2; ++j)
                        row.Add(Zero);
                    row[i - 1] = One;
                    row.Add(-MinimalPolynomial.Coefficients[i]);
                    matrixRows.Add(row);
                }
                Matrix transformMatrix = new Matrix(matrixRows);
                return numerator / transformMatrix.GetDeterminant();
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
                if (exponent.Numerator != One)
                    return Exponentiate<Number>(this, exponent.Numerator).Exponentiate(
                        One / exponent.Denominator);
                Rational rationalFactor = RationalFactor();
                return rationalFactor.Exponentiate(exponent) *
                    new Surd(One, this / rationalFactor, exponent.Denominator);
            }
            protected override int GetTypeIndex()
            {
                return 4;
            }
            public override int CompareTo(object obj)
            {
                int comparison = base.CompareTo(obj);
                if (comparison != 0)
                    return comparison;
                LinearlyIndependentSum expression = (LinearlyIndependentSum)obj;
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
                    output ^= term.GetHashCode();
                return output;
            }
            protected override RationalPolynomial CalculateMinimalPolynomial()
            {
                return VariableForm.MinimalPolynomial;
            }
            public override List<Number> GetConjugates()
            {
                return GetSumConjugates(Terms);
            }
            public override void RefineRealPartErrorInterval(Rational errorIntervalSize)
            {
                RealPartEstimate = EstimateSumRealPart(Terms, errorIntervalSize);
            }
            public override void RefineImaginaryPartErrorInterval(Rational errorIntervalSize)
            {
                ImaginaryPartEstimate = EstimateSumImaginaryPart(Terms, errorIntervalSize);
            }
            public override void RefineMagnitudeErrorInterval(Rational errorIntervalSize)
            {
                errorIntervalSize /= new Integer(3);
                Integer four = new Integer(4);
                RefineRealPartErrorInterval(One);
                RefineRealPartErrorInterval(errorIntervalSize /
                    (four * (Rational)GetMagnitudeInterval(RealPartEstimate).Max + Two));
                RefineImaginaryPartErrorInterval(One);
                RefineImaginaryPartErrorInterval(errorIntervalSize /
                    (four * (Rational)GetMagnitudeInterval(ImaginaryPartEstimate).Max + Two));
                Interval<Float> realPartInterval = GetMagnitudeInterval(RealPartEstimate);
                Interval<Float> imaginaryPartInterval = GetMagnitudeInterval(ImaginaryPartEstimate);
                Fraction half = new Fraction(One, Two);
                Number magnitudeUnderestimate =
                    ((Rational)(realPartInterval.Min * realPartInterval.Min +
                    imaginaryPartInterval.Min * imaginaryPartInterval.Min)).Exponentiate(half);
                magnitudeUnderestimate.RefineRealPartErrorInterval(errorIntervalSize);
                Number magnitudeOverestimate =
                    ((Rational)(realPartInterval.Max * realPartInterval.Max +
                    imaginaryPartInterval.Max * imaginaryPartInterval.Max)).Exponentiate(half);
                magnitudeOverestimate.RefineRealPartErrorInterval(errorIntervalSize);
                MagnitudeEstimate = new Interval<Float>(magnitudeUnderestimate.RealPartEstimate.Min,
                    magnitudeOverestimate.RealPartEstimate.Max);
            }
            protected override Interval<Rational> GetRationalArgumentEstimate(
                Rational errorIntervalSize)
            {
                bool hasZeroImaginaryPart(Number a)
                {
                    Rational localErrorIntervalSize = One;
                    a.RefineImaginaryPartErrorInterval(localErrorIntervalSize);
                    if (a.ImaginaryPartEstimate.Max < Float.Zero ||
                        Float.Zero < a.ImaginaryPartEstimate.Min)
                        return false;
                    List<Number> conjugates = a.GetConjugates();
                    conjugates.RemoveAt(0);
                    while (conjugates.Count > 1)
                    {
                        a.RefineRealPartErrorInterval(localErrorIntervalSize);
                        for (int i = 0; i < conjugates.Count; ++i)
                        {
                            conjugates[i].RefineRealPartErrorInterval(localErrorIntervalSize);
                            if (conjugates[i].RealPartEstimate.Max < a.RealPartEstimate.Min ||
                                a.RealPartEstimate.Max < conjugates[i].RealPartEstimate.Min)
                                conjugates.RemoveAt(i);
                            else
                                ++i;
                        }
                        localErrorIntervalSize /= Two;
                    }
                    while (true)
                    {
                        conjugates[0].RefineRealPartErrorInterval(localErrorIntervalSize);
                        a.RefineRealPartErrorInterval(localErrorIntervalSize);
                        if (conjugates[0].RealPartEstimate.Max < a.RealPartEstimate.Min ||
                            a.RealPartEstimate.Max < conjugates[0].RealPartEstimate.Min)
                            return true;
                        a.RefineImaginaryPartErrorInterval(localErrorIntervalSize);
                        if (a.ImaginaryPartEstimate.Min < Float.Zero ||
                            Float.Zero < a.ImaginaryPartEstimate.Max)
                            return false;
                        localErrorIntervalSize /= Two;
                    }
                }
                Number thisTimesI = this * new Surd(One, -One, Two);
                if (hasZeroImaginaryPart(this))
                {
                    if (hasZeroImaginaryPart(thisTimesI))
                        return new Interval<Rational>();
                    if (thisTimesI.ImaginaryPartEstimate.Min > Float.Zero)
                        return new Interval<Rational>(Zero, Zero);
                    Pi.RefineErrorInterval(errorIntervalSize);
                    return new Interval<Rational>(Pi.LowEstimate, Pi.HighEstimate);
                }
                if (hasZeroImaginaryPart(thisTimesI))
                {
                    if (thisTimesI.RealPartEstimate.Max < Float.Zero)
                    {
                        Pi.RefineErrorInterval(Two * errorIntervalSize);
                        return new Interval<Rational>(Pi.LowEstimate / Two, Pi.HighEstimate / Two);
                    }
                    Rational threeOverTwo = new Fraction(new Integer(3), Two);
                    Pi.RefineErrorInterval(errorIntervalSize / threeOverTwo);
                    return new Interval<Rational>(threeOverTwo * Pi.LowEstimate,
                        threeOverTwo * Pi.HighEstimate);
                }
                Interval<Float> realMagnitudeInterval =
                    GetMagnitudeInterval(RealPartEstimate);
                Interval<Float> imaginaryMagnitudeInterval =
                    GetMagnitudeInterval(ImaginaryPartEstimate);
                errorIntervalSize /= new Integer(3) * ((Rational)(realMagnitudeInterval.Max +
                    imaginaryMagnitudeInterval.Max) + One);
                RefineRealPartErrorInterval(errorIntervalSize *
                    (Rational)(realMagnitudeInterval.Min * realMagnitudeInterval.Min));
                RefineRealPartErrorInterval(errorIntervalSize);
                RefineImaginaryPartErrorInterval(errorIntervalSize * (Rational)(
                    imaginaryMagnitudeInterval.Min * imaginaryMagnitudeInterval.Min));
                RefineImaginaryPartErrorInterval(errorIntervalSize);
                Rational delta = null;
                Rational estimateArctangent(Rational a, Rational localErrorIntervalSize)
                {
                    Rational arctangentValue = a;
                    Rational aSquared = a * a;
                    Integer denominator = new Integer(3);
                    delta = -aSquared * a / denominator;
                    while (delta.Magnitude() > localErrorIntervalSize)
                    {
                        arctangentValue += delta;
                        denominator += Two;
                        delta *= aSquared / denominator;
                        delta = -delta;
                    }
                    return arctangentValue;
                }
                Interval<Rational> argument = new Interval<Rational>();
                if (RealPartEstimate.Min > Float.Zero)
                    if (ImaginaryPartEstimate.Min > Float.Zero)
                    {
                        if (ImaginaryPartEstimate.Min < RealPartEstimate.Max)
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Min = estimateArctangent((Rational)RealPartEstimate.Max /
                                (Rational)ImaginaryPartEstimate.Min, errorIntervalSize -
                                Pi.ErrorIntervalSize / Two) + Pi.LowEstimate / Two;
                        }
                        else
                            argument.Min = estimateArctangent((Rational)ImaginaryPartEstimate.Min /
                                (Rational)RealPartEstimate.Max, errorIntervalSize);
                        if (delta < Zero)
                            argument.Min += delta;
                        if (ImaginaryPartEstimate.Max < RealPartEstimate.Min)
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Max = estimateArctangent((Rational)RealPartEstimate.Min /
                                (Rational)ImaginaryPartEstimate.Max, errorIntervalSize -
                                Pi.ErrorIntervalSize / Two) + Pi.HighEstimate / Two;
                        }
                        else
                            argument.Max = estimateArctangent((Rational)ImaginaryPartEstimate.Max /
                                (Rational)RealPartEstimate.Min, errorIntervalSize);
                    }
                    else
                    {
                        if (-ImaginaryPartEstimate.Min < RealPartEstimate.Min)
                        {
                            Integer three = new Integer(3);
                            Pi.RefineErrorInterval(errorIntervalSize / three);
                            Rational threeOverTwo = new Fraction(three, Two);
                            argument.Min = estimateArctangent((Rational)RealPartEstimate.Min /
                                (Rational)ImaginaryPartEstimate.Min, errorIntervalSize - threeOverTwo *
                                Pi.ErrorIntervalSize) + threeOverTwo * Pi.LowEstimate;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize / new Integer(4));
                            argument.Min = estimateArctangent((Rational)ImaginaryPartEstimate.Min /
                                (Rational)RealPartEstimate.Min, errorIntervalSize -
                                Two * Pi.ErrorIntervalSize) + Two * Pi.LowEstimate;
                        }
                        if (delta < Zero)
                            argument.Min += delta;
                        if (-ImaginaryPartEstimate.Max < RealPartEstimate.Max)
                        {
                            Integer three = new Integer(3);
                            Pi.RefineErrorInterval(errorIntervalSize / three);
                            Rational threeOverTwo = new Fraction(three, Two);
                            argument.Max = estimateArctangent((Rational)RealPartEstimate.Max /
                                (Rational)ImaginaryPartEstimate.Max,
                                errorIntervalSize - threeOverTwo * Pi.ErrorIntervalSize) +
                                threeOverTwo * Pi.HighEstimate;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize / new Integer(4));
                            argument.Max = estimateArctangent((Rational)ImaginaryPartEstimate.Max /
                                (Rational)RealPartEstimate.Max, errorIntervalSize -
                                Two * Pi.ErrorIntervalSize) + Two * Pi.HighEstimate;
                        }
                    }
                else if (RealPartEstimate.Max < Float.Zero)
                    if (ImaginaryPartEstimate.Min > Float.Zero)
                    {
                        if (ImaginaryPartEstimate.Max < -RealPartEstimate.Min)
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Min = estimateArctangent((Rational)RealPartEstimate.Min /
                                (Rational)ImaginaryPartEstimate.Max, errorIntervalSize -
                                Pi.ErrorIntervalSize / Two) + Pi.LowEstimate / Two;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize / Two);
                            argument.Min = estimateArctangent((Rational)ImaginaryPartEstimate.Max /
                                (Rational)RealPartEstimate.Min, errorIntervalSize -
                                Pi.HighEstimate + Pi.LowEstimate) + Pi.LowEstimate;
                        }
                        if (delta < Zero)
                            argument.Min += delta;
                        if (ImaginaryPartEstimate.Min < -RealPartEstimate.Max)
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Max = estimateArctangent((Rational)RealPartEstimate.Max /
                                (Rational)ImaginaryPartEstimate.Min, errorIntervalSize -
                                Pi.ErrorIntervalSize / Two) + Pi.HighEstimate / Two;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize / Two);
                            argument.Max = estimateArctangent((Rational)ImaginaryPartEstimate.Min /
                                (Rational)RealPartEstimate.Max, errorIntervalSize -
                                Pi.ErrorIntervalSize) + Pi.HighEstimate;
                        }
                    }
                    else
                    {
                        if (ImaginaryPartEstimate.Min > RealPartEstimate.Max)
                        {
                            Integer three = new Integer(3);
                            Pi.RefineErrorInterval(errorIntervalSize / three);
                            Rational threeOverTwo = new Fraction(three, Two);
                            argument.Min = estimateArctangent((Rational)RealPartEstimate.Max /
                                (Rational)ImaginaryPartEstimate.Min, errorIntervalSize -
                                threeOverTwo * Pi.ErrorIntervalSize) +
                                threeOverTwo * Pi.LowEstimate;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Min = estimateArctangent((Rational)ImaginaryPartEstimate.Min /
                                (Rational)RealPartEstimate.Max, errorIntervalSize -
                                Pi.ErrorIntervalSize) + Pi.LowEstimate;
                        }
                        if (delta < Zero)
                            argument.Min += delta;
                        if (ImaginaryPartEstimate.Max > RealPartEstimate.Min)
                        {
                            Integer three = new Integer(3);
                            Pi.RefineErrorInterval(errorIntervalSize / three);
                            Rational threeOverTwo = new Fraction(three, Two);
                            argument.Max = estimateArctangent((Rational)RealPartEstimate.Min /
                                (Rational)ImaginaryPartEstimate.Max,
                                errorIntervalSize - threeOverTwo * Pi.ErrorIntervalSize) +
                                threeOverTwo * Pi.HighEstimate;
                        }
                        else
                        {
                            Pi.RefineErrorInterval(errorIntervalSize);
                            argument.Max = estimateArctangent((Rational)ImaginaryPartEstimate.Max /
                                (Rational)RealPartEstimate.Min, errorIntervalSize -
                                Pi.ErrorIntervalSize) + Pi.HighEstimate;
                        }
                    }
                if (delta > Zero)
                    argument.Max += delta;
                return argument;
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
                    for (int j = Rows[i].Count; j < rows.Count; ++j)
                        Rows[i].Add(Number.Zero);
                }
            }
            public Matrix GetTranspose()
            {
                List<List<Rational>> outputRows = new List<List<Rational>>();
                for (int i = 0; i < Rows.Count; ++i)
                {
                    List<Rational> row = new List<Rational>();
                    for (int j = 0; j < Rows.Count; ++j)
                        row.Add(Rows[j][i]);
                    outputRows.Add(row);
                }
                return new Matrix(outputRows);
            }
            public void GetRowEchelonForm<T>(List<T> augmentation)
                where T : IArithmetic<T>, IMultipliable<T, Rational>
            {
                int smallDimension;
                if (Rows.Count <= Rows[0].Count)
                    smallDimension = Rows.Count;
                else
                    smallDimension = Rows[0].Count;
                for (int i = 0; i < smallDimension; ++i)
                    for (int j = i + 1; j < Rows.Count; ++j)
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
                                    Rows[k][l] -= Rows[i][l] * scalar;
                                augmentation[k] =
                                    augmentation[k].Minus(augmentation[i].Times(scalar));
                            }
                            break;
                        }
            }
            public void Diagonalize(List<Rational> augmentation)
            {
                GetRowEchelonForm(augmentation);
                for (int i = Rows.Count - 1; i >= 0; --i)
                {
                    augmentation[i] /= Rows[i][i];
                    for (int j = 0; j <= i; ++j)
                        Rows[i][j] /= Rows[i][i];
                    for (int j = 0; j < i; ++j)
                    {
                        augmentation[j] -= augmentation[i] * Rows[j][i];
                        Rows[j][i] = Number.Zero;
                    }
                }
            }
            public Rational GetDeterminant()
            {//Mutates *this.
                Rational determinant = -Number.One;
                for (int i = 1; i < Rows.Count; ++i)
                    for (int j = i; j < Rows.Count; ++j)
                        if (Rows[i - 1][Rows.Count - i] == Number.Zero)
                        {
                            List<Rational> tempRow = Rows[j];
                            Rows[j] = Rows[i - 1];
                            Rows[i - 1] = tempRow;
                            determinant = -determinant;
                        }
                        else
                        {
                            determinant *= Rows[i - 1][Rows.Count - i];
                            for (int k = i; k < Rows.Count; ++k)
                                if (Rows[k][Rows.Count - i] != Number.Zero)
                                {
                                    Rational scalar = Rows[k][Rows.Count - i] /
                                        Rows[i - 1][Rows.Count - i];
                                    for (int l = 0; l < Rows.Count; ++l)
                                        Rows[k][l] -= Rows[i - 1][l] * scalar;
                                    determinant *= Rows[i - 1][Rows.Count - i];
                                }
                            break;
                        }
                return determinant;
            }
        }
        static List<T> CoefficientAdd<T>(List<T> coefficientsA, List<T> coefficientsB)
            where T : IArithmetic<T>
        {
            if (coefficientsA.Count < coefficientsB.Count)
                return CoefficientAdd(coefficientsB, coefficientsA);
            List<T> output = new List<T>(coefficientsA);
            for (int i = 0; i < coefficientsB.Count; ++i)
                output[i] = output[i].Plus(coefficientsB[i]);
            return output;
        }
        static List<T> CoefficientMultiply<T>(List<T> coefficientsA,
            List<T> coefficientsB) where T : IArithmetic<T>
        {
            List<T> output = new List<T>();
            T instance;
            if (coefficientsA.Count > 0)
                instance = coefficientsA[0];
            else if (coefficientsB.Count > 0)
                instance = coefficientsB[0];
            else
                return output;
            for (int i = 0; i < coefficientsA.Count + coefficientsB.Count - 1; ++i)
                output.Add(instance.GetAdditiveIdentity());
            for (int i = 0; i < coefficientsA.Count; ++i)
                for (int j = 0; j < coefficientsB.Count; ++j)
                    output[i + j] = output[i + j].Plus(coefficientsA[i].Times(coefficientsB[j]));
            return output;
        }
        static List<Rational> GetMonicCoefficients<T>(List<T> coefficients) where T : Rational
        {
            List<Rational> monicCoefficients = new List<Rational>();
            foreach (Rational coefficient in coefficients)
                monicCoefficients.Add(coefficient / coefficients[coefficients.Count - 1]);
            return monicCoefficients;
        }
#if DEBUG
        static string CoefficientsToString<T>(List<T> coefficients) where T : Rational
        {
            if (coefficients.Count == 0)
                return "0";
            StringBuilder output = new StringBuilder();
            if (coefficients[0] != Number.Zero)
                output.Append(coefficients[0]);
            if (coefficients.Count > 1)
            {
                if (coefficients[1] < Number.Zero)
                    output.Append(coefficients[1].InsertString("x"));
                else if (coefficients[1] > Number.Zero)
                {
                    if (output.Length != 0)
                        output.Append("+");
                    output.Append(coefficients[1].InsertString("x"));
                }
                for (int i = 2; i < coefficients.Count; ++i)
                    if (coefficients[i] < Number.Zero)
                        output.Append(coefficients[i].InsertString("x^" + i));
                    else if (coefficients[i] > Number.Zero)
                    {
                        if (output.Length != 0)
                            output.Append("+");
                        output.Append(coefficients[i].InsertString("x^" + i));
                    }
            }
            return output.ToString();
        }
#endif
        class IntegerPolynomial : IArithmetic<IntegerPolynomial>
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
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
                MinimalPolynomial = minimalPolynomial;
            }
            public IntegerPolynomial(RationalPolynomial minimalPolynomial,
                params Integer[] coefficients)
            {
                Coefficients = new List<Integer>(coefficients);
                while (Coefficients.Count != 0 &&
                    coefficients[Coefficients.Count - 1] == Number.Zero)
                    Coefficients.RemoveAt(Coefficients.Count - 1);
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
                    negativeCoefficients.Add(-coefficient);
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
                    return new IntegerPolynomial(MinimalPolynomial,
                        CoefficientMultiply(Coefficients, a.Coefficients));
                List<Rational> coefficients = CoefficientMultiply(
                    ((RationalPolynomial)this).Coefficients, ((RationalPolynomial)a).Coefficients);
                for (int i = coefficients.Count; i >= MinimalPolynomial.Coefficients.Count; --i)
                {
                    for (int j = 0; j < MinimalPolynomial.Coefficients.Count - 1; ++j)
                        coefficients[i - MinimalPolynomial.Coefficients.Count + j] -=
                            MinimalPolynomial.Coefficients[j];
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
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        Rational quotientCoefficient = remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1];
                        if (quotientCoefficient.Denominator == Number.One)
                        {
                            quotient.Insert(0, quotientCoefficient.Numerator);
                            remainder.RemoveAt(remainder.Count - 1);
                            for (int j = 1; j < divisor.Coefficients.Count; ++j)
                                remainder[remainder.Count - j] -= quotient[0] *
                                    divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                        }
                        else
                        {
                            division.Quotient = null;
                            division.Remainder = null;
                            return division;
                        }
                    }
                division.Quotient = new IntegerPolynomial(MinimalPolynomial, quotient);
                division.Remainder = new IntegerPolynomial(MinimalPolynomial, remainder);
                return division;
            }
            public bool Equals(IntegerPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                    return false;
                for (int i = 0; i < Coefficients.Count; ++i)
                    if (Coefficients[i] != a.Coefficients[i])
                        return false;
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
                    coefficients.Add(Coefficients[i] * new Integer(i));
                return new IntegerPolynomial(null, coefficients);
            }
            public IntegerPolynomial GetPrimitivePart()
            {
                if (Coefficients.Count == 0)
                    return this;
                return this / Solver.GetGCD(Coefficients);
            }
            public static IntegerPolynomial GetGCD(IntegerPolynomial a, IntegerPolynomial b)
            {
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    if (a.Coefficients.Count == 0)
                        return b;
                    IntegerPolynomial t = a;
                    a = b;
                    b = t;
                }
                else if (b.Coefficients.Count == 0)
                    return a;
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
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    else
                        h = (Exponentiate(g, degree) / Exponentiate(h, -hExponent)).Numerator;
                    degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    remainder = (a * Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                        degree + Number.One)) % b;
                }
                if (remainder.Coefficients.Count == 1)
                    return new IntegerPolynomial(a.MinimalPolynomial, d);
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
                        squarefreeFactors.Add(a);
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
                                break;
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
                            if (remainder * Number.Two <= characteristicPower)
                                polynomial.Coefficients[i] = remainder;
                            else
                                polynomial.Coefficients[i] = remainder - characteristicPower;
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
                                v *= liftedFactors[index];
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
                            ;
                        ++combinationSize;
                    }
                    if (2 * combinationSize == liftedFactors.Count)                    
                        if (testFactorCombination(combinationSize - 1, 0, new List<int> { 0 }))
                            liftedFactors.RemoveAt(0);                    
                    IntegerPolynomial finalFactor = new IntegerPolynomial(null,
                        factor.Coefficients[factor.Coefficients.Count - 1]);
                    foreach (IntegerPolynomial liftedFactor in liftedFactors)
                        finalFactor *= liftedFactor;
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
                            factor.Coefficients.RemoveAt(0);
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
                        irreducibleFactors.AddRange(splitFactor(factor));
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
                        Coefficients.Add(remainder + characteristic);
                    else
                        Coefficients.Add(remainder);
                }
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1] == Number.Zero)
                    Coefficients.RemoveAt(Coefficients.Count - 1);
                Characteristic = characteristic;
            }
            public ModdedPolynomial(Integer characteristic, params Integer[] coefficients)
            {
                Coefficients = new List<Integer>();
                foreach (Integer coefficient in coefficients)
                {
                    Integer remainder = coefficient.EuclideanDivideBy(characteristic).Remainder;
                    if (remainder < Number.Zero)
                        Coefficients.Add(remainder + characteristic);
                    else
                        Coefficients.Add(remainder);
                }
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1] == Number.Zero)
                    Coefficients.RemoveAt(Coefficients.Count - 1);
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
                    output.Add(-coefficient);
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
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] * divisor.Coefficients[
                            divisor.Coefficients.Count - 1].Reciprocal(Characteristic));
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                    }
                division.Quotient = new ModdedPolynomial(quotient, Characteristic);
                division.Remainder = new ModdedPolynomial(remainder, Characteristic);
                return division;
            }
            public bool Equals(ModdedPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                    return false;
                for (int i = 0; i < Coefficients.Count; ++i)
                    if (Coefficients[i] != a.Coefficients[i])
                        return false;
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
                    coefficients.Add(Coefficients[i] * new Integer(i));
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
                while (d < (Coefficients.Count + 1) / 2)
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
                            C = T + (C * C) % factorProduct;
                        B = GetGCD(factorProduct, C);
                        while (B.Coefficients.Count < 2 ||
                            B.Coefficients.Count == factorProduct.Coefficients.Count)
                        {
                            T *= new ModdedPolynomial(Characteristic, Number.Zero, Number.Zero,
                                Number.One);
                            C = new ModdedPolynomial(new List<Integer>(T.Coefficients),
                                Characteristic);
                            for (int j = 1; j < degree; ++j)
                                C = T + (C * C) % factorProduct;
                            B = GetGCD(factorProduct, C);
                        }
                    }
                    else
                    {
                        List<Integer> coefficients = new List<Integer>();
                        Integer p = Characteristic - Number.One;
                        for (int j = 1; j < 2 * degree; ++j)
                            coefficients.Add(p);
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
                    cantorZassenhausSplit(distinctDegreeFactors[degree], degree);
                return irreducibleFactors;
            }
#if DEBUG
            public override string ToString()
            {
                return CoefficientsToString(Coefficients) + " mod " + Characteristic.ToString();
            }
#endif
        }
        class RationalPolynomial :
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
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
                MinimalPolynomial = minimalPolynomial;
            }
            public RationalPolynomial(RationalPolynomial minimalPolynomial,
                params Rational[] coefficients)
            {
                Coefficients = new List<Rational>(coefficients);
                while (Coefficients.Count != 0 &&
                    coefficients[Coefficients.Count - 1] == Number.Zero)
                    Coefficients.RemoveAt(Coefficients.Count - 1);
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
                    output.Add(-coefficient);
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
                    for (int i = coefficients.Count; i >= MinimalPolynomial.Coefficients.Count; --i)
                    {
                        for (int j = 0; j < MinimalPolynomial.Coefficients.Count - 1; ++j)
                            coefficients[i - MinimalPolynomial.Coefficients.Count + j] -=
                                coefficients[coefficients.Count - 1] *
                                MinimalPolynomial.Coefficients[j];
                        coefficients.RemoveAt(coefficients.Count - 1);
                    }
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public RationalPolynomial Times(Rational a)
            {
                List<Rational> coefficients = new List<Rational>();
                for (int i = 0; i < Coefficients.Count; ++i)
                    coefficients.Add(a * Coefficients[i]);
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public Division<RationalPolynomial> EuclideanDivideBy(RationalPolynomial divisor)
            {
                Division<RationalPolynomial> division;
                List<Rational> quotient = new List<Rational>();
                List<Rational> remainder = new List<Rational>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1]);
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                    }
                division.Quotient = new RationalPolynomial(quotient, MinimalPolynomial);
                division.Remainder = new RationalPolynomial(remainder, MinimalPolynomial);
                return division;
            }
            public bool Equals(RationalPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                    return false;
                for (int i = 0; i < Coefficients.Count; ++i)
                    if (Coefficients[i] != a.Coefficients[i])
                        return false;
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
                    return a.EuclideanDivideBy(b).Quotient;
                if (b.Coefficients.Count == 1)
                    return a * (Number.One / b.Coefficients[0]);
                RationalPolynomial[] minimalPolynomials =
                    new RationalPolynomial[b.MinimalPolynomial.Coefficients.Count];
                for (int i = 0; i < b.MinimalPolynomial.Coefficients.Count; ++i)
                    minimalPolynomials[i] = b.MinimalPolynomial;
                //Only the 0th value of minimalPolynomials is meaningful. Reducing the degree of a
                //term with respect to any of the other variables isn't relevant, but no term of
                //degree higher than one with respect to those variables will ever be produced, so
                //arbitrarily using b.MinimalPolynomial for them ensures that the reduction code
                //won't be invoked. It is necessary because I assumed all MinimalPolynomials values
                //are non-null when writing MultivariatePolynomial.
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
                        row.Add(Number.Zero);
                    matrixRows.Add(row);
                    augmentation.Add(Number.Zero);
                }
                augmentation[0] = Number.One;
                foreach (int[] indices in product.Coefficients.Keys)
                {
                    int columnIndex = 0;
                    for (int i = 1; i < indices.Length; ++i)
                        if (indices[i] == 1)
                        {
                            columnIndex = i - 1;
                            break;
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
            public static implicit operator RationalPolynomial(IntegerPolynomial a)
            {
                List<Rational> coefficients = new List<Rational>();
                foreach (Rational coefficient in a.Coefficients)
                    coefficients.Add(coefficient);
                return new RationalPolynomial(coefficients);
            }
            public RectangularEstimate EstimateEvaluation(Primitive argument,
                Rational errorIntervalSize)
            {
                argument.RefineRealPartErrorInterval(Number.One);
                argument.RefineImaginaryPartErrorInterval(Number.One);
                Rational realErrorScale = Number.One;
                Rational imaginaryErrorScale = Number.One;
                Rational realMagnitudeBound =
                    (Rational)GetMagnitudeInterval(argument.RealPartEstimate).Max;
                Rational imaginaryMagnitudeBound =
                    (Rational)GetMagnitudeInterval(argument.ImaginaryPartEstimate).Max;
                Division<Integer> division =
                    new Integer(Coefficients.Count - 1).EuclideanDivideBy(Number.Two);
                Integer numberOfTermsInRealEvaluation =
                    Number.One + division.Quotient * (new Integer(3) + division.Quotient);
                Integer numberOfTermsInImaginaryEvaluation;
                if (division.Remainder == Number.Zero)
                {
                    numberOfTermsInRealEvaluation -= division.Quotient + Number.One;
                    numberOfTermsInImaginaryEvaluation =
                        division.Quotient * (division.Quotient + Number.One);
                }
                else
                {
                    Integer quotientCeiling = division.Quotient + Number.One;
                    numberOfTermsInImaginaryEvaluation = quotientCeiling * quotientCeiling;
                }
                for (int i = 1; i < Coefficients.Count; ++i)
                {                    
                    Integer termDegree = new Integer(i);
                    Rational coefficientMagnitude = Coefficients[i].Magnitude();
                    Rational realMagnitudePower =
                        Exponentiate(realMagnitudeBound, termDegree - Number.One);
                    Rational realScaleCandidate = Number.Two * coefficientMagnitude *
                        realMagnitudePower * numberOfTermsInRealEvaluation;
                    if (realScaleCandidate > realErrorScale)
                        realErrorScale = realScaleCandidate;
                    Rational imaginaryMagnitudePower = Number.One;
                    Rational imaginaryScaleCandidate;
                    Integer activeTermCount = numberOfTermsInImaginaryEvaluation;
                    Integer inactiveTermCount = numberOfTermsInRealEvaluation;
                    for (Integer j = Number.One; j < termDegree; ++j)
                    {
                        Rational nextImaginaryMagnitudePower =
                            imaginaryMagnitudePower * imaginaryMagnitudeBound;
                        Rational sharedScaleComponent = Number.Two * activeTermCount *
                            termDegree.ThisChooseK(j) * coefficientMagnitude *
                            (Number.One + realMagnitudePower + nextImaginaryMagnitudePower);
                        realMagnitudePower /= realMagnitudeBound;
                        realScaleCandidate = sharedScaleComponent * realMagnitudePower;
                        if (realScaleCandidate > realErrorScale)
                            realErrorScale = realScaleCandidate;                        
                        imaginaryScaleCandidate = sharedScaleComponent * imaginaryMagnitudePower;
                        imaginaryMagnitudePower = nextImaginaryMagnitudePower;
                        if (imaginaryScaleCandidate > imaginaryErrorScale)
                            imaginaryErrorScale = imaginaryScaleCandidate;
                        Integer temp = activeTermCount;
                        activeTermCount = inactiveTermCount;
                        inactiveTermCount = temp;
                    }
                    imaginaryScaleCandidate = Number.Two * coefficientMagnitude *
                        imaginaryMagnitudePower * activeTermCount;
                    if (imaginaryScaleCandidate > imaginaryErrorScale)
                        imaginaryErrorScale = imaginaryScaleCandidate;
                }
                argument.RefineRealPartErrorInterval(errorIntervalSize / realErrorScale);
                argument.RefineImaginaryPartErrorInterval(errorIntervalSize / imaginaryErrorScale);
                Rational realEvaluationBoundA = Number.Zero;
                Rational realEvaluationBoundB = Number.Zero;
                Rational imaginaryEvaluationBoundA = Number.Zero;
                Rational imaginaryEvaluationBoundB = Number.Zero;                
                Rational realPartMin = (Rational)argument.RealPartEstimate.Min;
                Rational realPartMax = (Rational)argument.RealPartEstimate.Max;
                Rational imaginaryPartMin = (Rational)argument.ImaginaryPartEstimate.Min;
                Rational imaginaryPartMax = (Rational)argument.ImaginaryPartEstimate.Max;                
                for (int i = 0; i < Coefficients.Count; ++i)
                {
                    Integer termDegree = new Integer(i);
                    Rational realPartMinPower = Exponentiate(realPartMin, termDegree);
                    Rational realPartMaxPower = Exponentiate(realPartMax, termDegree);
                    Rational imaginaryPartMinPower = Number.One;
                    Rational imaginaryPartMaxPower = Number.One;
                    for (int j = 0; j <= i; ++j)
                    {
                        Integer binomialCoefficient = termDegree.ThisChooseK(new Integer(j));
                        Rational termBoundA = Coefficients[i] * binomialCoefficient *
                            realPartMinPower * imaginaryPartMinPower;
                        Rational termBoundB = Coefficients[i] * binomialCoefficient *
                            realPartMaxPower * imaginaryPartMaxPower;
                        int remainder = j % 4;
                        switch (remainder)
                        {
                            case 0:
                                realEvaluationBoundA += termBoundA;
                                realEvaluationBoundB += termBoundB;
                                break;
                            case 1:
                                imaginaryEvaluationBoundA += termBoundA;
                                imaginaryEvaluationBoundB += termBoundB;
                                break;
                            case 2:
                                realEvaluationBoundA -= termBoundA;
                                realEvaluationBoundB -= termBoundB;
                                break;
                            case 3:
                                imaginaryEvaluationBoundA -= termBoundA;
                                imaginaryEvaluationBoundB -= termBoundB;
                                break;
                        }
                        realPartMinPower /= realPartMin;
                        realPartMaxPower /= realPartMax;
                        imaginaryPartMinPower *= imaginaryPartMin;
                        imaginaryPartMaxPower *= imaginaryPartMax;
                    }
                }
                Rational placeValue = Float.GetEstimatePlaceValue(errorIntervalSize);
                realEvaluationBoundA.RefineRealPartErrorInterval(placeValue);
                realEvaluationBoundB.RefineRealPartErrorInterval(placeValue);
                imaginaryEvaluationBoundA.RefineRealPartErrorInterval(placeValue);
                imaginaryEvaluationBoundB.RefineRealPartErrorInterval(placeValue);
                Interval<Float> evaluationRealPart;
                if (realEvaluationBoundA > realEvaluationBoundB)
                    evaluationRealPart =
                        new Interval<Float>(realEvaluationBoundB.RealPartEstimate.Min,
                        realEvaluationBoundA.RealPartEstimate.Max);
                else
                    evaluationRealPart =
                        new Interval<Float>(realEvaluationBoundA.RealPartEstimate.Min,
                        realEvaluationBoundB.RealPartEstimate.Max);
                Interval<Float> evaluationImaginaryPart;
                if (imaginaryEvaluationBoundA > imaginaryEvaluationBoundB)
                    evaluationImaginaryPart =
                        new Interval<Float>(imaginaryEvaluationBoundB.RealPartEstimate.Min,
                        imaginaryEvaluationBoundA.RealPartEstimate.Max);
                else
                    evaluationImaginaryPart =
                        new Interval<Float>(imaginaryEvaluationBoundA.RealPartEstimate.Min,
                        imaginaryEvaluationBoundB.RealPartEstimate.Max);
                return new RectangularEstimate(evaluationRealPart, evaluationImaginaryPart);
            }
            public RationalPolynomial GetDerivative()
            {
                List<Rational> coefficients = new List<Rational>();
                for (int i = 1; i < Coefficients.Count; ++i)
                    coefficients.Add(Coefficients[i] * new Integer(i));
                return new RationalPolynomial(coefficients, MinimalPolynomial);
            }
            public IntegerPolynomial GetPrimitivePart()
            {
                Integer denomiatorProduct = Number.One;
                foreach (Rational coefficient in Coefficients)
                    denomiatorProduct *= coefficient.Denominator;
                List<Integer> coefficients = new List<Integer>();
                foreach (Rational coefficient in Coefficients)
                    coefficients.Add((denomiatorProduct * coefficient).Numerator);
                return new IntegerPolynomial(MinimalPolynomial, coefficients).GetPrimitivePart();
            }
            public static RationalPolynomial GetGCD(RationalPolynomial a, RationalPolynomial b)
            {
                return IntegerPolynomial.GetGCD(a.GetPrimitivePart(), b.GetPrimitivePart());
            }
            public RationalPolynomial GetMinimalPolynomial()
            {
                List<RationalPolynomial> powers = new List<RationalPolynomial>();
                RationalPolynomial power = GetMultiplicativeIdentity();
                List<List<Rational>> matrixRows = new List<List<Rational>>();
                List<RationalPolynomial> augmentation = new List<RationalPolynomial>();
                List<Rational> augmentationRow = new List<Rational> { Number.One };
                for (int i = 0; i < MinimalPolynomial.Coefficients.Count; ++i)
                {
                    power *= this;
                    powers.Add(power);
                    List<Rational> matrixRow = new List<Rational>(power.Coefficients);
                    while (matrixRow.Count < MinimalPolynomial.Coefficients.Count - 1)
                        matrixRow.Add(Number.Zero);
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
                        sum += powers[i - 1] * factor.Coefficients[i];
                    if (sum.Coefficients.Count == 0)
                    {
                        minimalPolynomial = factor;
                        break;
                    }
                }
                return new RationalPolynomial(GetMonicCoefficients(minimalPolynomial.Coefficients));
            }
            public RationalPolynomial GetMinimalPolynomialOfKTimesRootOfThis(Rational k)
            {
                List<Rational> coefficients = new List<Rational>();
                Rational power = Number.One;
                foreach (Rational coefficient in Coefficients)
                {
                    coefficients.Add(power * coefficient);
                    power /= k;
                }
                return new RationalPolynomial(GetMonicCoefficients(coefficients));
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
                    coefficients.RemoveAt(coefficients.Count - 1);
                Coefficients = coefficients;
                InnerMinimalPolynomial = innerMinimalPolynomial;
            }
            public NestedPolynomial(RationalPolynomial innerMinimalPolynomial,
                params RationalPolynomial[] coefficients)
            {
                Coefficients = new List<RationalPolynomial>(coefficients);
                while (Coefficients.Count != 0 &&
                    Coefficients[Coefficients.Count - 1].Coefficients.Count == 0)
                    Coefficients.RemoveAt(Coefficients.Count - 1);
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
                    negativeCoefficients.Add(coefficient.GetNegative());
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
                    productCoefficients.Add(a * Coefficients[i]);
                return new NestedPolynomial(InnerMinimalPolynomial, productCoefficients);
            }
            public Division<NestedPolynomial> EuclideanDivideBy(NestedPolynomial divisor)
            {
                Division<NestedPolynomial> division;
                List<RationalPolynomial> quotient = new List<RationalPolynomial>();
                List<RationalPolynomial> remainder = new List<RationalPolynomial>(Coefficients);
                if (divisor.Coefficients.Count <= Coefficients.Count)
                    for (int i = 0; i <= Coefficients.Count - divisor.Coefficients.Count; ++i)
                    {
                        quotient.Insert(0, remainder[remainder.Count - 1] /
                            divisor.Coefficients[divisor.Coefficients.Count - 1]);
                        remainder.RemoveAt(remainder.Count - 1);
                        for (int j = 1; j < divisor.Coefficients.Count; ++j)
                            remainder[remainder.Count - j] -= quotient[0] *
                                divisor.Coefficients[divisor.Coefficients.Count - j - 1];
                    }
                division.Quotient = new NestedPolynomial(InnerMinimalPolynomial, quotient);
                division.Remainder = new NestedPolynomial(InnerMinimalPolynomial, remainder);
                return division;
            }
            public bool Equals(NestedPolynomial a)
            {
                if (Coefficients.Count != a.Coefficients.Count)
                    return false;
                for (int i = 0; i < Coefficients.Count; ++i)
                    if (!Coefficients[i].Equals(a.Coefficients[i]))
                        return false;
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
                    quotientCoefficients.Add(coefficient / b);
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
                    derivativeCoefficients.Add(new Integer(i) * Coefficients[i]);
                return new NestedPolynomial(InnerMinimalPolynomial, derivativeCoefficients);
            }
            NestedPolynomial GetPrimitivePart()
            {
                return this / Solver.GetGCD(Coefficients);
            }
            static RationalPolynomial GetResultant(NestedPolynomial a, NestedPolynomial b)
            {
                if (a.Coefficients.Count == 0 || b.Coefficients.Count == 0)
                    return new RationalPolynomial(a.InnerMinimalPolynomial);
                RationalPolynomial aContent = Solver.GetGCD(a.Coefficients);
                RationalPolynomial bContent = Solver.GetGCD(b.Coefficients);
                a /= aContent;
                b /= bContent;
                RationalPolynomial g = a.Coefficients[0].GetMultiplicativeIdentity();
                RationalPolynomial h = a.Coefficients[0].GetMultiplicativeIdentity();
                Rational s = Number.One;
                if (a.Coefficients.Count % 2 == 0 && b.Coefficients.Count % 2 == 0)
                    s = -s;
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
                        s = -s;
                    NestedPolynomial remainder = (a * Exponentiate(
                        b.Coefficients[b.Coefficients.Count - 1], degree + Number.One)) % b;
                    a = b;
                    b = remainder / (g * Exponentiate(h, degree));
                    g = a.Coefficients[a.Coefficients.Count - 1];
                    hExponent = Number.One - degree;
                    if (hExponent >= Number.Zero)
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    else
                        h = Exponentiate(g, degree) / Exponentiate(h, -hExponent);
                }
                hExponent = new Integer(2 - a.Coefficients.Count);
                if (hExponent >= Number.Zero)
                    return s * t * Exponentiate(h, hExponent) * Exponentiate(b.Coefficients[
                        b.Coefficients.Count - 1], new Integer(a.Coefficients.Count - 1));
                return s * t * Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                    new Integer(a.Coefficients.Count - 1)) / Exponentiate(h, -hExponent);
            }
            static NestedPolynomial GetGCD(NestedPolynomial a, NestedPolynomial b)
            {
                if (b.Coefficients.Count > a.Coefficients.Count)
                {
                    if (a.Coefficients.Count == 0)
                        return b;
                    NestedPolynomial t = a;
                    a = b;
                    b = t;
                }
                else if (b.Coefficients.Count == 0)
                    return a;
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
                        h = Exponentiate(h, hExponent) * Exponentiate(g, degree);
                    else
                        h = Exponentiate(g, degree) / Exponentiate(h, -hExponent);
                    degree = new Integer(a.Coefficients.Count - b.Coefficients.Count);
                    remainder = a.Times(Exponentiate(b.Coefficients[b.Coefficients.Count - 1],
                        degree + Number.One)) % b;
                }
                if (remainder.Coefficients.Count == 1)
                    return new NestedPolynomial(a.InnerMinimalPolynomial, d);
                return b.GetPrimitivePart().Times(d);
            }
            public NestedPolynomial SwitchVariables(RationalPolynomial innerMinimalPolynomial =
                null)
            {
                if (Coefficients.Count == 0)
                    return this;
                List<List<Rational>> coefficientCoefficients = new List<List<Rational>>();
                foreach (Rational coefficient in Coefficients[0].Coefficients)
                    coefficientCoefficients.Add(new List<Rational> { coefficient });
                for (int i = 1; i < Coefficients.Count; ++i)
                {
                    for (int j = coefficientCoefficients.Count;
                        j < Coefficients[i].Coefficients.Count; ++j)
                    {
                        coefficientCoefficients.Add(new List<Rational>());
                        for (int k = 0; k < i; ++k)
                            coefficientCoefficients[j].Add(Number.Zero);
                    }
                    for (int j = 0; j < Coefficients[i].Coefficients.Count; ++j)
                        coefficientCoefficients[j].Add(Coefficients[i].Coefficients[j]);
                }
                List<RationalPolynomial> coefficients = new List<RationalPolynomial>();
                foreach (List<Rational> list in coefficientCoefficients)
                    coefficients.Add(new RationalPolynomial(list, innerMinimalPolynomial));
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
                        squarefreeFactors.Add(a);
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
                            break;
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
                    return "0";
                StringBuilder output = new StringBuilder();
                if (Coefficients[0].Coefficients.Count != 0)
                    output.Append("(" + Coefficients[0] + ")");
                if (Coefficients.Count > 1)
                {
                    if (Coefficients[1].Coefficients.Count != 0)
                    {
                        if (output.Length != 0)
                            output.Append("+");
                        output.Append("(" + Coefficients[1] + ")y");
                    }
                    for (int i = 2; i < Coefficients.Count; ++i)
                        if (Coefficients[i].Coefficients.Count != 0)
                        {
                            if (output.Length != 0)
                                output.Append("+");
                            output.Append("(" + Coefficients[i] + ")y^" + i);
                        }
                }
                return output.ToString();
            }
#endif
        }
        class MultivariatePolynomial : IRingElement<MultivariatePolynomial>
        {//The nth variable represents a root of MinimalPolynomials[n].
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
                    if (a[i] != b[i])
                        return false;
                return true;
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
                                    product.Coefficients[productIndices] +
                                    a.Coefficients[indicesA] * b.Coefficients[indicesB];
                                indicesFound = true;
                                break;
                            }
                        if (!indicesFound)
                            product.Coefficients.Add(productIndices,
                                a.Coefficients[indicesA] * b.Coefficients[indicesB]);
                    }
                return product;
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
                    sum.Coefficients.Add(indices, a.Coefficients[indices]);
                foreach (int[] indicesB in b.Coefficients.Keys)
                {
                    bool indicesFound = false;
                    foreach (int[] indicesA in sum.Coefficients.Keys)
                        if (AreEqualByValue(indicesA, indicesB))
                        {
                            Rational termSum =
                                sum.Coefficients[indicesA] + b.Coefficients[indicesB];
                            if (termSum == Number.Zero)
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
                            a.Coefficients[indicesA] * b.Coefficients[indicesB]);
                        for (int i = 0; i < indicesA.Length; ++i)
                        {
                            MultivariatePolynomial termComponent =
                                new MultivariatePolynomial(a.MinimalPolynomials);
                            int index = indicesA[i] + indicesB[i];
                            if (index < a.MinimalPolynomials[i].Coefficients.Count - 1)
                            {
                                int[] termComponentIndices = new int[indicesA.Length];
                                termComponentIndices[i] = index;
                                termComponent.Coefficients.Add(termComponentIndices, Number.One);
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
                                    reducedDegreeFactor.SetCoefficient(termComponentIndices,
                                        Number.One);
                                    productComponent = ReductionlessMultiply(productComponent,
                                        reducedDegreeFactor);
                                }
                                for (int j = 0; j < a.MinimalPolynomials[i].Coefficients.Count - 1;
                                    ++j)
                                {
                                    if (a.MinimalPolynomials[i].Coefficients[j] == Number.Zero)
                                        continue;
                                    int[] termComponentIndices = new int[indicesA.Length];
                                    termComponentIndices[i] = j;
                                    termComponent.Coefficients[termComponentIndices] =
                                        -a.MinimalPolynomials[i].Coefficients[j] /
                                        a.MinimalPolynomials[i].Coefficients[
                                        a.MinimalPolynomials[i].Coefficients.Count - 1];
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
                if (a == Number.Zero)
                    return product;
                foreach (int[] indices in b.Coefficients.Keys)
                    product.Coefficients.Add(indices, (Rational)(a * b.Coefficients[indices]));
                return product;
            }
            public static MultivariatePolynomial operator *(MultivariatePolynomial a, Rational b)
            {
                return b * a;
            }
            public RationalPolynomial GetMinimalPolynomial()
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
                        matrixRow.Add(Number.Zero);
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
                            for (int i = 0; i < matrixCoefficients.Count; ++i)
                                matrixCoefficients[i].Add(Number.Zero);
                            matrixRow.Add(power.Coefficients[indices]);
                        }
                    }
                    matrixRow.Reverse();
                    matrixCoefficients.Add(matrixRow);
                    if (matrixRow[0] != Number.Zero)
                        constantIsPresent = true;
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
                foreach (IntegerPolynomial factor in factors)
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
                return new RationalPolynomial(GetMonicCoefficients(minimalPolynomial.Coefficients));
            }
        }
        class Float : IRingElement<Float>
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
                    return new Float(-Significand, RadixPointPosition);
                return this;
            }
            public Float Times(Float a)
            {
                return this * a;
            }
            public static Float operator +(Float a, Float b)
            {
                if (a.RadixPointPosition == b.RadixPointPosition)
                    return GetReducedForm(a.Significand + b.Significand, a.RadixPointPosition);
                if (b.RadixPointPosition > a.RadixPointPosition)
                    return new Float(b.Significand + a.Significand * Exponentiate(Number.Two,
                        b.RadixPointPosition - a.RadixPointPosition), b.RadixPointPosition);
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
                    return a;
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
                    placeValue /= Number.Two;
                return placeValue / Number.Two;
            }
            public Interval<Float> EstimateRoot(Rational errorIntervalSize, Integer index)
            {
                Debug.Assert(Significand >= Number.Zero,
                    "*this was negative, which this algorithm doesn't account for.");
                Float rootEstimate;
                if (this < One)
                    rootEstimate = One;
                else
                    rootEstimate = this;
                Rational rationalRadicand = (Rational)this;
                Integer indexMinusOne = index - Number.One;
                Rational delta = (rationalRadicand / (Rational)Exponentiate(rootEstimate,
                    indexMinusOne) - (Rational)rootEstimate) / index;
                delta.RefineRealPartErrorInterval(-delta / Number.Two);
                while ((Rational)delta.RealPartEstimate.Max - Number.Two * delta >
                    errorIntervalSize)
                {
                    rootEstimate += delta.RealPartEstimate.Max;
                    delta = (rationalRadicand / (Rational)Exponentiate(rootEstimate,
                        indexMinusOne) - (Rational)rootEstimate) / index;
                    delta.RefineRealPartErrorInterval(-delta / Number.Two);
                }
                rootEstimate += delta.RealPartEstimate.Max;
                return new Interval<Float>(rootEstimate + delta.RealPartEstimate.Max, rootEstimate);
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
                        RefineErrorInterval();
                }
            }
        }
        static class RootsOfUnity
        {
            static Dictionary<Integer, List<Number>> NthRoots = new Dictionary<Integer,
                List<Number>> { { Number.One, new List<Number> { Number.One } },
                { Number.Two, new List<Number> { Number.One, -Number.One } } };
            public static List<Number> GetNthRoots(Integer n)
            {
                if (NthRoots.ContainsKey(n))
                    return NthRoots[n];
                List<Number> nthRoots = new List<Number>();
                Division<Integer> division;
                Integer half = n.EuclideanDivideBy(Number.Two).Quotient;
                foreach (Integer prime in Primes())
                {
                    if (prime > half)
                        break;
                    division = n.EuclideanDivideBy(prime);
                    if (division.Remainder == Number.Zero)
                    {
                        List<List<Number>> rootCombinations = GetCartesianProduct(
                            new List<List<Number>> { GetNthRoots(prime),
                        GetNthRoots(division.Quotient) });
                        foreach (List<Number> combination in rootCombinations)
                            nthRoots.Add(combination[0] * combination[1].Exponentiate(
                                new Fraction(Number.One, prime)));
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
                        break;
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
                    factors.Add(unfactoredComponent);
                Integer generator = Number.One;
                bool isGenerator = false;
                while (!isGenerator)
                {
                    ++generator;
                    isGenerator = true;
                    foreach (Integer factor in factors)
                        if (Exponentiate(generator,
                            nMinusOne.EuclideanDivideBy(factor).Quotient) == Number.One)
                        {
                            isGenerator = false;
                            break;
                        }
                }
                List<Rational> nMinusFirstRootAnnullingPolynomialCoefficients =
                    new List<Rational> { -Number.One };
                for (Integer i = Number.One; i < nMinusOne; ++i)
                    nMinusFirstRootAnnullingPolynomialCoefficients.Add(Number.Zero);
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
                        if (key[1] == 0)
                            resolventMultiplesInTermsOfNMinusFirstRoots[i][key[0]] +=
                                resolventProduct.Coefficients[key];
                        else if (key[1] == 1)
                            resolventMultiplesInTermsOfNMinusFirstRoots[i][key[0]] -=
                                resolventProduct.Coefficients[key];
                }
                List<Number> nMinusFirstRoots = GetNthRoots(nMinusOne);
                Pi.RefineErrorInterval(Pi.LowEstimate / (Number.Two * nMinusOne));
                foreach (Number root in nMinusFirstRoots)
                    root.RefineArgumentErrorInterval(Pi.LowEstimate);
                nMinusFirstRoots.Sort(delegate (Number a, Number b)
                {
                    if (a.ArgumentEstimate.Max <= b.ArgumentEstimate.Min)
                        return -1;
                    return 1;
                });
                List<Number> resolventProductValues = new List<Number>();
                foreach (List<Rational> coefficients in resolventMultiplesInTermsOfNMinusFirstRoots)
                {
                    Number resolventProductValue = coefficients[0];
                    for (int i = 1; i < intNMinusOne; ++i)
                        resolventProductValue += coefficients[i] * nMinusFirstRoots[i];
                    resolventProductValues.Add(resolventProductValue);
                }
                List<Number> resolventValues = new List<Number> { new Integer(-1),
                    resolventProductValues[0].Exponentiate(new Fraction(Number.One, nMinusOne)) };
                for (int i = 1; i < resolventProductValues.Count; ++i)
                    resolventValues.Add(resolventValues[1] / resolventProductValues[i]);
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
                            if (primeCandidate.EuclideanDivideBy(
                                PrimeList[j]).Remainder == Number.Zero)
                            {
                                isDivisible = true;
                                break;
                            }
                        if (!isDivisible)
                            break;
                        primeCandidate += Number.Two;
                    }
                    PrimeList.Add(primeCandidate);
                }
                yield return PrimeList[i];
                ++i;
            }
        }
        static List<List<T>> GetCartesianProduct<T>(List<List<T>> sets)
        {//The first element of the return value is the list of the first element of each element of
         //the input value. Number.GetConjugates() depends on this fact in order to ensure that the
         //first element of its return value is the instance it was called from.
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
        static T Exponentiate<T>(T expBase, Integer exponent, Multiplier<T> multiply,
            T multiplicativeIdentity)
        {
            Debug.Assert(exponent >= Number.Zero,
                "exponent was negative, which this algorithm doesn't account for.");
            T output = multiplicativeIdentity;
            T baseToAPowerOfTwo = expBase;
            while (exponent > Number.Zero)
            {
                Division<Integer> division = exponent.EuclideanDivideBy(Number.Two);
                if (division.Remainder == Number.One)
                    output = multiply(output, baseToAPowerOfTwo);
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
                return a;
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
                GCD = GetGCD(GCD, list[i]);
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
        static void Main(string[] args)
        {
            Number EvaluateExpression(List<char> operations, List<Number> numbers)
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
                        numbers[i] = EvaluateExpression(
                            operations.GetRange(i + 1, matchingParenIndex - i - 1),
                            numbers.GetRange(i + 1, matchingParenIndex - i - 1));
                        operations[i] = ' ';
                        numbers.RemoveRange(i + 1, matchingParenIndex - i);
                        operations.RemoveRange(i + 1, matchingParenIndex - i);
                    }
                    else
                        ++i;
                }
                void executeOperations()
                {
                    for (int i = 0; i < operations.Count;)
                        if (i == 0 || numbers[i - 1] == null)
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
                                if (numbers[i + 1] != null)
                                {
                                    numbers[i + 1] = -numbers[i + 1];
                                    numbers.RemoveAt(i);
                                    operations.RemoveAt(i);
                                }
                                else if (operations[i + 1] == '+')
                                {
                                    operations[i + 1] = '-';
                                    numbers.RemoveAt(i);
                                    operations.RemoveAt(i);
                                }
                                else if (operations[i + 1] == '-')
                                {
                                    operations[i + 1] = '+';
                                    numbers.RemoveAt(i);
                                    operations.RemoveAt(i);
                                }
                                else
                                    throw new InvalidUserInput("Operator missing operand.");
                            else
                                ++i;
                        else
                            ++i;
                    for (int i = 0; i < operations.Count;)
                        if (operations[i] == '^')
                        {
                            if (numbers[i + 1] is Rational exponent)
                                numbers[i - 1] = numbers[i - 1].Exponentiate(exponent);
                            else
                                throw new InvalidUserInput("The input expression contains an " +
                                    "exponentiation\nwhose exponent is not both real and " +
                                    "rational;\nthis program doesn't handle transcendental " +
                                    "numbers.");
                            numbers.RemoveRange(i, 2);
                            operations.RemoveRange(i, 2);
                        }
                        else
                            ++i;
                    for (int i = 0; i < operations.Count;)
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
                    for (int i = 0; i < operations.Count;)
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
#if DEBUG
                executeOperations();
#else
                try
                {
                    executeOperations();
                }
                catch (ArgumentOutOfRangeException)
                {
                    throw new InvalidUserInput("Operator missing operand.");
                }
#endif
                if (numbers.Count == 0)
                    return null;
                return numbers[0];
            }
            void EvaluateString(string input)
            {
                List<char> operations = new List<char>();
                List<Number> numbers = new List<Number>();
                StringBuilder numberCollector = new StringBuilder();
                bool lastCharWasDigit = false;
                try
                {
                    foreach (char c in input)
                    {
                        if (lastCharWasDigit)
                        {
                            if (char.IsDigit(c))
                                numberCollector.Append(c);
                            else
                            {
                                if (!"()+-*/^".Contains(c.ToString()))
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
                        else if (operations[i] == ')' && i + 1 < operations.Count &&
                            (numbers[i + 1] != null || operations[i + 1] == '('))
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
            }
#if DEBUG
            EvaluateString("1/(1+2^(1/3))");
#else
            while (true)
            {
                string input = Console.ReadLine();
                if (input[0] == 'q')
                    return;
                EvaluateString(input);
            }
#endif
        }
    }
}
