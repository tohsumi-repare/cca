#ifndef MOLBIOLIB_APPLIEDSTATISTICSALGORITHMS_H
#define MOLBIOLIB_APPLIEDSTATISTICSALGORITHMS_H 1

//////////////////////////////////////////////////////////////////////////////
// MolBioLib: A C++11 framework for rapidly developing bioinformatics tasks //
// Copyright (C)  2012  Massachusetts General Hospital                      //
//                                                                          //
// This program is free software: you can redistribute it and/or modify     //
// it under the terms of the GNU General Public License as published by     //
// the Free Software Foundation, version 3 of the License.                  //
//                                                                          //
// This program is distributed AS IS.                                       //
//////////////////////////////////////////////////////////////////////////////


/** \file AppliedStatisticsAlgorithms.hpp
 * Contains functions from StatLib.
 *
 * This is code that can be freely distributed, as indicated on
 *
 * http://lib.stat.cmu.edu/apstat/
 *
 * linked from
 *
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/cpp_src.html
 *
 */


#include "src/include/PrimitiveTypes.hpp"



double alnfac ( int n );
// Below is the successor to alogam; - computes ln(Gamma(xvalue)).
double alngam ( double xvalue, int *ifault );
double alnorm ( double x, bool upper );
// int ifault; double fx = gammad(x, p, &ifault); is equivalent to
// Numerical Recipes' double fx = gammp(p, x); - the normalized
// incomplete Gamma function.
double gammad ( double x, double p, int *ifault );
double r8_abs ( double x );
double r8_min ( double x, double y );
double chyper ( bool point, int kk, int ll, int mm, int nn, int *ifault );
void hypergeometric_cdf_values ( int *n_data, int *sam, int *suc, int *pop, int *n, double *fx );
void hypergeometric_pdf_values ( int *n_data, int *sam, int *suc, int *pop, int *n, double *fx );
int i4_min ( int i1, int i2 );
double r8_max ( double x, double y );
// int ifault;
// Note the - before the third term below.
// double beta_log = alngam(a, &ifault) +
//                   alngam(b, &ifault) - alngam(a+b, &ifault);
// double fx = betain(x, a, b, beta_log, &ifault); is equivalent to
// Numerical Recipes' double fx = betai(a, b, x);  - regularlized
// incomplete Beta function.  This is often denoted I_x(a, b).
// I_x(a, b) is expressed in terms of the usual incomplete Beta function:
// I_x(a, b) = B(x; a, b)/B(a, b).  Note betain = I, not B.
double betain ( double x, double p, double q, double beta, int *ifault );

// Below is the evaluation of the 2F1 hypergeometric function based on
// a recipe by Gosper.  See (below is primary).
// http://www.axelvogt.de/axalom/hyp2F1/hypergeometric_2F1_using_a_recipe_of_Gosper.mws.pdf
// and
// http://www.math.upenn.edu/~wilf/AeqB.html
// Used by Macsyma and others.  Works only if z < 1.
double hypg2F1(double a, double b, double c, double z,
               double epsilon, size_t maxIter);
// Below is from
// http://en.wikipedia.org/wiki/Student%27s_t-distribution CDF formula
double studentTDistCDF(double dof, double t);
// Below is from
// http://www.ncl.ucar.edu/Document/Functions/Built-in/betainc.shtml
double studentTDistPvalue(double dof, double t, bool twoTailed);
// Below is from
// http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
double spearmanZscore(double n, double r);
double spearmanPvalue(double n, double r, bool twoTailed);


// Below are some simple functions to compute the Student t value for
// various cases.  All are from
// http://en.wikipedia.org/wiki/Student's_t-test
double studentTOneSampleDOF(double sampleSize)
{
    return (sampleSize - 1);
}
double studentTOneSample(double mean, double expectedMean,
                         double stddev, double sampleSize)
{
    return ((mean - expectedMean)/(stddev/sqrt(sampleSize)));
}
double studentTTwoSamplesEqualSizesEqualVarianceDOF(double sampleSize)
{
    return (2.0*sampleSize - 2.0);
}
double studentTTwoSamplesEqualSizesEqualVariance(double mean1, double mean2,
        double stddev1, double stddev2,
        double sampleSize)
{
    double sx1x2 = sqrt(0.5*(stddev1*stddev1 + stddev2*stddev2));
    return ((mean1 - mean2)/(sx1x2*sqrt(2.0/sampleSize)));
}
double studentTTwoSamplesUnequalSizesEqualVarianceDOF(double sampleSize1,
        double sampleSize2)
{
    return (sampleSize1 + sampleSize2 - 2.0);
}
double studentTTwoSamplesUnequalSizesEqualVariance(double mean1, double mean2,
        double stddev1,
        double stddev2,
        double sampleSize1,
        double sampleSize2)
{
    double sx1x2 = sqrt(((sampleSize1 - 1.0)*stddev1*stddev1 +
                         (sampleSize2 - 1.0)*stddev2*stddev2)/
                        (sampleSize1 + sampleSize2 - 2.0));
    return ((mean1 - mean2)/(sx1x2*sqrt(1.0/sampleSize1 + 1.0/sampleSize2)));
}
double studentTTwoSamplesUnequalSizesUnequalVarianceDOF(double stddev1,
        double stddev2,
        double sampleSize1,
        double sampleSize2)
{
    double sdev12 = stddev1*stddev1, sdev22 = stddev2*stddev2;
    double sqrtnumer = sdev12/sampleSize1 + sdev22/sampleSize2;
    double denom1 = (sdev12/sampleSize1)*(sdev12/sampleSize1)/
                    (sampleSize1 - 1.0);
    double denom2 = (sdev22/sampleSize2)*(sdev22/sampleSize2)/
                    (sampleSize2 - 1.0);
    return (sqrtnumer*sqrtnumer/(denom1 + denom2));
}
double studentTTwoSamplesUnequalSizesUnequalVariance(double mean1, double mean2,
        double stddev1,
        double stddev2,
        double sampleSize1,
        double sampleSize2)
{
    double sx1x2 = sqrt(stddev1*stddev1/sampleSize1 +
                        stddev2*stddev2/sampleSize2);
    return ((mean1 - mean2)/sx1x2);
}





//****************************************************************************80

double hypg2F1(double a, double b, double c, double z,
               double epsilon = 1.0e-14, size_t maxIter = 1000)

//****************************************************************************80
{

#ifdef DEBUG
    if (fabs(z) >= 1.0)
    {
        cerr << "Warning!  In hypg2F1, z = " << z << " is >= 1.0, so the recurrence relation may not converge." << endl;
    }
#endif

    double dk = 0.0, ek = 1.0, fk = 0.0,
           dnext = 0.0, enext = 1.0, fnext = 0.0;
    double maxIterations = static_cast<double>(maxIter);
    double k = 0.0;

    // The k == 0 is there since the first time, the > epsilon will fail.
    while ((k < maxIterations) &&
            ((fabs(fnext - fk) > epsilon) || (k == 0.0)))
    {
        dk = dnext;
        ek = enext;
        fk = fnext;

        dnext = (k+a)*(k+b)*z*(ek - (k+c-b-a)*dk*z/(1.0-z))/
                (4.0*(k+1.0)*(k+c/2.0)*(k+(c+1.0)/2.0));
        enext = (k+a)*(k+b)*z*(a*b*dk*z/(1.0-z) + (k+c)*ek)/
                (4.0*(k+1.0)*(k+c/2.0)*(k+(c+1.0)/2.0));
        fnext = fk -
                dk*(k*((c-b-a)*z + k*(z-2.0)-c) - a*b*z)/
                (2.0*(k+c/2.0)*(1.0-z))                    +
                ek;
        k += 1.0;
    }

#ifdef DEBUG
    if (k == maxIterations)
        cerr << "hypg2F1 stopped on maximum iterations = " << maxIterations
             << endl;
#endif

    return fnext;
}




//****************************************************************************80

double studentTDistCDF(double dof, double t)

//****************************************************************************80
{

    // First compute f1 = Gamma((dof+1)/2)/Gamma(dof/2)
    // Since we have the ln(Gamma(x)) function, compute exp(ln(f)), which gives
    // f1 = exp(lnGamma((dof+1)/2) - lnGamma(dof/2))
    int ifault = 0;
    double t1 = alngam((dof + 1.0)/2.0, &ifault),
           t2 = alngam(dof/2.0, &ifault);
    double f1 = exp(t1 - t2);

    // Next, compute 2F1(0.5, (dof+1)/2, 1.5, -t*t/dof)
    double f2 = hypg2F1(0.5, (dof+1.0)/2.0, 1.5, -t*t/dof);

    // Now compute sqrt(pi*dof)
    double f3 = sqrt(M_PI*dof);

    // Finally, return CDF
    return (0.5 + t*f1*f2/f3);

}



//****************************************************************************80

double studentTDistPvalue(double dof, double t, bool twoTailed = true)

//****************************************************************************80
{

    double x = dof/(dof + t*t);
    double a = dof/2.0;
    double b = 0.5;
    int ifault = 0;
    double beta_log = alngam(a, &ifault) +
                      alngam(b, &ifault) - alngam(a+b, &ifault);
    double divFactor = twoTailed ? 1.0 : 2.0;
    double p = betain(x, a, b, beta_log, &ifault);
    return (p/divFactor);
}



//****************************************************************************80

double spearmanZscore(double n, double r)

//****************************************************************************80
{

#ifdef DEBUG
    if (r == 1.0)
    {
        cerr << "Error!  In spearmanZscore, r = 1, so the Fisher transform value will be infinite." << endl;
    }
#endif

    double F = 0.5*log((1.0+r)/(1.0-r));
    double z = sqrt((n-3)/1.06) * F;

    return z;
}



//****************************************************************************80

double spearmanPvalue(double n, double r, bool twoTailed = true)

//****************************************************************************80
{

    double t = r*sqrt((n - 2.0)/(1.0 - r*r));
    return studentTDistPvalue(n-2.0, t, twoTailed);
}



//****************************************************************************80

double alnfac ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    ALNFAC computes the logarithm of the factorial of N.
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial.
//
//    Output, double ALNFAC, the logarithm of the factorial of N.
//
{
    int ier;
    double value;

    value = alngam ( ( double ) ( n + 1 ), &ier );

    return value;
}



//****************************************************************************80

double alngam ( double xvalue, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    ALNGAM computes the logarithm of the gamma function.
//
//  Modified:
//
//    13 January 2008
//
//  Author:
//
//    Allan Macleod
//    C++ version by John Burkardt
//
//  Reference:
//
//    Allan Macleod,
//    Algorithm AS 245,
//    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
//    Applied Statistics,
//    Volume 38, Number 2, 1989, pages 397-402.
//
//  Parameters:
//
//    Input, double XVALUE, the argument of the Gamma function.
//
//    Output, int IFAULT, error flag.
//    0, no error occurred.
//    1, XVALUE is less than or equal to 0.
//    2, XVALUE is too big.
//
//    Output, double ALNGAM, the logarithm of the gamma function of X.
//
{
    double alr2pi = 0.918938533204673;
    double r1[9] =
    {
        -2.66685511495,
        -24.4387534237,
        -21.9698958928,
        11.1667541262,
        3.13060547623,
        0.607771387771,
        11.9400905721,
        31.4690115749,
        15.2346874070
    };
    double r2[9] =
    {
        -78.3359299449,
        -142.046296688,
        137.519416416,
        78.6994924154,
        4.16438922228,
        47.0668766060,
        313.399215894,
        263.505074721,
        43.3400022514
    };
    double r3[9] =
    {
        -2.12159572323E+05,
        2.30661510616E+05,
        2.74647644705E+04,
        -4.02621119975E+04,
        -2.29660729780E+03,
        -1.16328495004E+05,
        -1.46025937511E+05,
        -2.42357409629E+04,
        -5.70691009324E+02
    };
    double r4[5] =
    {
        0.279195317918525,
        0.4917317610505968,
        0.0692910599291889,
        3.350343815022304,
        6.012459259764103
    };
    double value;
    double x;
    double x1;
    double x2;
    double xlge = 510000.0;
    double xlgst = 1.0E+30;
    double y;

    x = xvalue;
    value = 0.0;
    //
    //  Check the input.
    //
    if ( xlgst <= x )
    {
        *ifault = 2;
        return value;
    }

    if ( x <= 0.0 )
    {
        *ifault = 1;
        return value;
    }

    *ifault = 0;
    //
    //  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
    //
    if ( x < 1.5 )
    {
        if ( x < 0.5 )
        {
            value = - log ( x );
            y = x + 1.0;
            //
            //  Test whether X < machine epsilon.
            //
            if ( y == 1.0 )
            {
                return value;
            }
        }
        else
        {
            value = 0.0;
            y = x;
            x = ( x - 0.5 ) - 0.5;
        }

        value = value + x * ((((
                                   r1[4]   * y
                                   + r1[3] ) * y
                               + r1[2] ) * y
                              + r1[1] ) * y
                             + r1[0] ) / ((((
                                                y
                                                + r1[8] ) * y
                                            + r1[7] ) * y
                                           + r1[6] ) * y
                                          + r1[5] );

        return value;
    }
    //
    //  Calculation for 1.5 <= X < 4.0.
    //
    if ( x < 4.0 )
    {
        y = ( x - 1.0 ) - 1.0;

        value = y * ((((
                           r2[4]   * x
                           + r2[3] ) * x
                       + r2[2] ) * x
                      + r2[1] ) * x
                     + r2[0] ) / ((((
                                        x
                                        + r2[8] ) * x
                                    + r2[7] ) * x
                                   + r2[6] ) * x
                                  + r2[5] );
    }
    //
    //  Calculation for 4.0 <= X < 12.0.
    //
    else if ( x < 12.0 )
    {
        value = ((((
                       r3[4]   * x
                       + r3[3] ) * x
                   + r3[2] ) * x
                  + r3[1] ) * x
                 + r3[0] ) / ((((
                                    x
                                    + r3[8] ) * x
                                + r3[7] ) * x
                               + r3[6] ) * x
                              + r3[5] );
    }
    //
    //  Calculation for 12.0 <= X.
    //
    else
    {
        y = log ( x );
        value = x * ( y - 1.0 ) - 0.5 * y + alr2pi;

        if ( x <= xlge )
        {
            x1 = 1.0 / x;
            x2 = x1 * x1;

            value = value + x1 * ( (
                                       r4[2]   *
                                       x2 + r4[1] ) *
                                   x2 + r4[0] ) / ( (
                                           x2 + r4[4] ) *
                                                    x2 + r4[3] );
        }
    }

    return value;
}
//****************************************************************************80

double alnorm ( double x, bool upper )

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    David Hill
//    C++ version by John Burkardt
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
    double a1 = 5.75885480458;
    double a2 = 2.62433121679;
    double a3 = 5.92885724438;
    double b1 = -29.8213557807;
    double b2 = 48.6959930692;
    double c1 = -0.000000038052;
    double c2 = 0.000398064794;
    double c3 = -0.151679116635;
    double c4 = 4.8385912808;
    double c5 = 0.742380924027;
    double c6 = 3.99019417011;
    double con = 1.28;
    double d1 = 1.00000615302;
    double d2 = 1.98615381364;
    double d3 = 5.29330324926;
    double d4 = -15.1508972451;
    double d5 = 30.789933034;
    double ltone = 7.0;
    double p = 0.398942280444;
    double q = 0.39990348504;
    double r = 0.398942280385;
    bool up;
    double utzero = 18.66;
    double value;
    double y;
    double z;

    up = upper;
    z = x;

    if ( z < 0.0 )
    {
        up = !up;
        z = - z;
    }

    if ( ltone < z && ( ( !up ) || utzero < z ) )
    {
        if ( up )
        {
            value = 0.0;
        }
        else
        {
            value = 1.0;
        }
        return value;
    }

    y = 0.5 * z * z;

    if ( z <= con )
    {
        value = 0.5 - z * ( p - q * y
                            / ( y + a1 + b1
                                / ( y + a2 + b2
                                    / ( y + a3 ))));
    }
    else
    {
        value = r * exp ( - y )
                / ( z + c1 + d1
                    / ( z + c2 + d2
                        / ( z + c3 + d3
                            / ( z + c4 + d4
                                / ( z + c5 + d5
                                    / ( z + c6 ))))));
    }

    if ( !up )
    {
        value = 1.0 - value;
    }

    return value;
}
//****************************************************************************80

// Below function is omitted, but the comments remain.
// void gamma_inc_values ( int *n_data, double *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    In Mathematica, the function can be evaluated by:
//
//      1 - GammaRegularized[A,X]
//
//  Modified:
//
//    20 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//

//****************************************************************************80

double gammad ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMAD computes the Incomplete Gamma Integral
//
//  Auxiliary functions:
//
//    ALOGAM = logarithm of the gamma function,
//    ALNORM = algorithm AS66
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    B Shea
//    C++ version by John Burkardt
//
//  Reference:
//
//    B Shea,
//    Algorithm AS 239:
//    Chi-squared and Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int IFAULT, error flag.
//    0, no error.
//    1, X < 0 or P <= 0.
//
//    Output, double GAMMAD, the value of the incomplete
//    Gamma integral.
//
{
    double a;
    double an;
    double arg;
    double b;
    double c;
    double elimit = - 88.0;
    double oflo = 1.0E+37;
    double plimit = 1000.0;
    double pn1;
    double pn2;
    double pn3;
    double pn4;
    double pn5;
    double pn6;
    double rn;
    double tol = 1.0E-14;
    bool upper;
    double value;
    double xbig = 1.0E+08;

    value = 0.0;
    //
    //  Check the input.
    //
    if ( x < 0.0 )
    {
        *ifault = 1;
        return value;
    }

    if ( p <= 0.0 )
    {
        *ifault = 1;
        return value;
    }

    *ifault = 0;

    if ( x == 0.0 )
    {
        value = 0.0;
        return value;
    }
    //
    //  If P is large, use a normal approximation.
    //
    if ( plimit < p )
    {
        pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 )
                                   + 1.0 / ( 9.0 * p ) - 1.0 );

        upper = false;
        value = alnorm ( pn1, upper );
        return value;
    }
    //
    //  If X is large set value = 1.
    //
    if ( xbig < x )
    {
        value = 1.0;
        return value;
    }
    //
    //  Use Pearson's series expansion.
    //  (Note that P is not large enough to force overflow in ALOGAM).
    //  No need to test IFAULT on exit since P > 0.
    //
    if ( x <= 1.0 || x < p )
    {
        arg = p * log ( x ) - x - alngam ( p + 1.0, ifault );
        c = 1.0;
        value = 1.0;
        a = p;

        for ( ; ; )
        {
            a = a + 1.0;
            c = c * x / a;
            value = value + c;

            if ( c <= tol )
            {
                break;
            }
        }

        arg = arg + log ( value );

        if ( elimit <= arg )
        {
            value = exp ( arg );
        }
        else
        {
            value = 0.0;
        }
    }
    //
    //  Use a continued fraction expansion.
    //
    else
    {
        arg = p * log ( x ) - x - alngam ( p, ifault );
        a = 1.0 - p;
        b = a + x + 1.0;
        c = 0.0;
        pn1 = 1.0;
        pn2 = x;
        pn3 = x + 1.0;
        pn4 = x * b;
        value = pn3 / pn4;

        for ( ; ; )
        {
            a = a + 1.0;
            b = b + 2.0;
            c = c + 1.0;
            an = a * c;
            pn5 = b * pn3 - an * pn1;
            pn6 = b * pn4 - an * pn2;

            if ( pn6 != 0.0 )
            {
                rn = pn5 / pn6;

                if ( r8_abs ( value - rn ) <= r8_min ( tol, tol * rn ) )
                {
                    break;
                }
                value = rn;
            }

            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;
            //
            //  Re-scale terms in continued fraction if terms are large.
            //
            if ( oflo <= abs ( pn5 ) )
            {
                pn1 = pn1 / oflo;
                pn2 = pn2 / oflo;
                pn3 = pn3 / oflo;
                pn4 = pn4 / oflo;
            }
        }

        arg = arg + log ( value );

        if ( elimit <= arg )
        {
            value = 1.0 - exp ( arg );
        }
        else
        {
            value = 1.0;
        }
    }

    return value;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
    double value;

    if ( 0.0 <= x )
    {
        value = x;
    }
    else
    {
        value = -x;
    }
    return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
    double value;

    if ( y < x )
    {
        value = y;
    }
    else
    {
        value = x;
    }
    return value;
}


//****************************************************************************80

double chyper ( bool point, int kk, int ll, int mm, int nn, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    CHYPER computes point or cumulative hypergeometric probabilities.
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    Richard Lund
//    C++ version by John Burkardt
//
//  Reference:
//
//    PR Freeman,
//    Algorithm AS 59:
//    Hypergeometric Probabilities,
//    Applied Statistics,
//    Volume 22, Number 1, 1973, pages 130-133.
//
//    Richard Lund,
//    Algorithm AS 152:
//    Cumulative hypergeometric probabilities,
//    Applied Statistics,
//    Volume 29, Number 2, 1980, pages 221-223.
//
//    BL Shea,
//    Remark AS R77:
//    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
//    Applied Statistics,
//    Volume 38, Number 1, 1989, pages 199-204.
//
//  Parameters:
//
//    Input, bool POINT, is TRUE if the point probability is desired,
//    and FALSE if the cumulative probability is desired.
//
//    Input, int KK, the sample size.
//    0 <= KK <= MM.
//
//    Input, int LL, the number of successes in the sample.
//    0 <= LL <= KK.
//
//    Input, int MM, the population size that was sampled.
//    0 <= MM.
//
//    Input, int NN, the number of "successes" in the population.
//    0 <= NN <= MM.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double CHYPER, the PDF (point probability) of
//    exactly LL successes out of KK samples, or the CDF (cumulative
//    probability) of up to LL successes out of KK samples.
//
// Modified Tue Feb 22 15:15:27 EST 2011 by Toshiro K. Ohsumi
// Changed any multiplication of ints (then cast to double) into
// cast first into double, then multiply, to prevent overflow.
{
    double arg;
    bool dir;
    double elimit = - 88.0;
    int i;
    int j;
    int k;
    int kl;
    int l;
    int m;
    int mbig = 600;
    double mean;
    int mnkl;
    int mvbig = 1000;
    int n;
    int nl;
    double p;
    double pt;
    double rootpi = 2.506628274631001;
    double scale = 1.0E+35;
    double sig;
    double value;

    *ifault = 0;

    k = kk + 1;
    l = ll + 1;
    m = mm + 1;
    n = nn + 1;

    dir = true;
    //
    //  Check arguments are within permitted limits.
    //
    value = 0.0;

    if ( n < 1 || m < n || k < 1 || m < k )
    {
        *ifault = 1;
        return value;
    }

    if ( l < 1 || m - n < k - l )
    {
        *ifault = 2;
        return value;
    }

    if ( !point )
    {
        value = 1.0;
    }

    if ( n < l || k < l )
    {
        *ifault = 2;
        return value;
    }

    *ifault = 0;
    value = 1.0;

    if ( k == 1 || k == m || n == 1 || n == m )
    {
        return value;
    }

    if ( !point && ll == i4_min ( kk, nn ) )
    {
        return value;
    }

    p = ( double ) ( nn ) / ( double ) ( mm - nn );

    if ( 16.0 * r8_max ( p, 1.0 / p )
            < ( double ) ( i4_min ( kk, mm - kk ) ) &&
            mvbig < mm && - 100.0 < elimit )
    {
        //
        //  Use a normal approximation.
        //
        // mean = ( double ) ( kk * nn ) / ( double ) ( mm );
        mean = ( double ) ( kk )  * ( double ) ( nn ) / ( double ) ( mm );

        sig = sqrt ( mean * ( ( double ) ( mm - nn ) / ( double ) ( mm ) )
                     * ( ( double ) ( mm - kk ) / ( ( double ) ( mm - 1 ) ) ) );

        if ( point )
        {
            arg = - 0.5 * ( pow ( ( ( double ) ( ll ) - mean ) / sig, 2 ) );
            if ( elimit <= arg )
            {
                value = exp ( arg ) / ( sig * rootpi );
            }
            else
            {
                value = 0.0;
            }
        }
        else
        {
            value = alnorm ( ( ( double ) ( ll ) + 0.5 - mean ) / sig, false );
        }
    }
    else
    {
        //
        //  Calculate exact hypergeometric probabilities.
        //  Interchange K and N if this saves calculations.
        //
        if ( i4_min ( n - 1, m - n ) < i4_min ( k - 1, m - k ) )
        {
            i = k;
            k = n;
            n = i;
        }

        if ( m - k < k - 1 )
        {
            dir = !dir;
            l = n - l + 1;
            k = m - k + 1;
        }

        if ( mbig < mm )
        {
            //
            //  Take logarithms of factorials.
            //
            p = alnfac ( nn )
                - alnfac ( mm )
                + alnfac ( mm - kk )
                + alnfac ( kk )
                + alnfac ( mm - nn )
                - alnfac ( ll )
                - alnfac ( nn - ll )
                - alnfac ( kk - ll )
                - alnfac ( mm - nn - kk + ll );

            if ( elimit <= p )
            {
                value = exp ( p );
            }
            else
            {
                value = 0.0;
            }
        }
        else
        {
            //
            //  Use Freeman/Lund algorithm.
            //
            for ( i = 1; i <= l - 1; i++ )
            {
                // value = value * ( double ) ( ( k - i ) * ( n - i ) )
                // / ( double ) ( ( l - i ) * ( m - i ) );

                value = value * (static_cast<double>( k - i )) * (static_cast<double>( n - i ))
                        / ( (static_cast<double>( l - i )) * (static_cast<double>( m - i )) ) ;
            }

            if ( l != k )
            {
                j = m - n + l;
                for ( i = l; i <= k - 1; i++ )
                {
                    value = value * static_cast<double>( j - i ) / static_cast<double>( m - i );
                }
            }
        }

        if ( point )
        {
            return value;
        }

        if ( value == 0.0 )
        {
            //
            //  We must recompute the point probability since it has underflowed.
            //
            if ( mm <= mbig )
            {
                p = alnfac ( nn )
                    - alnfac ( mm )
                    + alnfac ( kk )
                    + alnfac ( mm - nn )
                    - alnfac ( ll )
                    - alnfac ( nn - ll )
                    - alnfac ( kk - ll )
                    - alnfac ( mm - nn - kk + ll )
                    + alnfac ( mm - kk );
            }

            p = p + log ( scale );

            if ( p < elimit )
            {
                *ifault = 3;
                // if ( ( double ) ( nn * kk + nn + kk + 1 )
                if ( ( (static_cast<double>(nn)) * (static_cast<double>(kk)) + (static_cast<double>(nn)) + (static_cast<double>(kk)) + 1.0 )
                        / static_cast<double>( mm + 2 ) < static_cast<double>( ll ) )
                {
                    value = 1.0;
                }
                return value;
            }
            else
            {
                p = exp ( p );
            }
        }
        else
            //
            //  Scale up at this point.
            //
        {
            p = value * scale;
        }

        pt = 0.0;
        nl = n - l;
        kl = k - l;
        mnkl = m - n - kl + 1;

        if ( l <= kl )
        {
            for ( i = 1; i <= l - 1; i++ )
            {
                // p = p * ( double ) ( ( l - i ) * ( mnkl - i ) ) /
                // ( double ) ( ( nl + i ) * ( kl + i ) );
                p = p * ( ( (static_cast<double>(l)) - (static_cast<double>(i)) ) * ( (static_cast<double>(mnkl)) - (static_cast<double>(i)) ) ) /
                    ( ( (static_cast<double>(nl)) + (static_cast<double>(i)) ) * ( (static_cast<double>(kl)) + (static_cast<double>(i)) ) );
                pt = pt + p;
            }
        }
        else
        {
            dir = !dir;
            for ( j = 0; j <= kl - 1; j++ )
            {
                // p = p * ( double ) ( ( nl - j ) * ( kl - j ) )
                // / ( double ) ( ( l + j ) * ( mnkl + j ) );
                p = p * ( ( (static_cast<double>(nl)) - (static_cast<double>(j)) ) * ( (static_cast<double>(kl)) - (static_cast<double>(j)) ) )
                    / ( ( (static_cast<double>(l)) + (static_cast<double>(j)) ) * ( (static_cast<double>(mnkl)) + (static_cast<double>(j)) ) );
                pt = pt + p;
            }
        }

        if ( p == 0.0 )
        {
            *ifault = 3;
        }

        if ( dir )
        {
            value = value + ( pt / scale );
        }
        else
        {
            value = 1.0 - ( pt / scale );
        }
    }

    return value;
}
//****************************************************************************80

void hypergeometric_cdf_values ( int *n_data, int *sam, int *suc, int *pop,
                                 int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      CDF [ dist, n ]
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *SAM, int *SUC, int *POP, the sample size,
//    success size, and population parameters of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 16

    double fx_vec[N_MAX] =
    {
        0.6001858177500578E-01,
        0.2615284665839845E+00,
        0.6695237889132748E+00,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.5332595856827856E+00,
        0.1819495964117640E+00,
        0.4448047017527730E-01,
        0.9999991751316731E+00,
        0.9926860896560750E+00,
        0.8410799901444538E+00,
        0.3459800113391901E+00,
        0.0000000000000000E+00,
        0.2088888139634505E-02,
        0.3876752992448843E+00,
        0.9135215248834896E+00
    };

    int n_vec[N_MAX] =
    {
        7,  8,  9, 10,
        6,  6,  6,  6,
        6,  6,  6,  6,
        0,  0,  0,  0
    };

    int pop_vec[N_MAX] =
    {
        100, 100, 100, 100,
        100, 100, 100, 100,
        100, 100, 100, 100,
        90,  200, 1000, 10000
    };

    int sam_vec[N_MAX] =
    {
        10, 10, 10, 10,
        6,  7,  8,  9,
        10, 10, 10, 10,
        10, 10, 10, 10
    };

    int suc_vec[N_MAX] =
    {
        90, 90, 90, 90,
        90, 90, 90, 90,
        10, 30, 50, 70,
        90, 90, 90, 90
    };

    if ( *n_data < 0 )
    {
        *n_data = 0;
    }

    *n_data = *n_data + 1;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *sam = 0;
        *suc = 0;
        *pop = 0;
        *n = 0;
        *fx = 0.0;
    }
    else
    {
        *sam = sam_vec[*n_data-1];
        *suc = suc_vec[*n_data-1];
        *pop = pop_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return;
# undef N_MAX
}
//****************************************************************************80

void hypergeometric_pdf_values ( int *n_data, int *sam, int *suc, int *pop,
                                 int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      PDF [ dist, n ]
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *SAM, int *SUC, int *POP, the sample size,
//    success size, and population parameters of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 16

    double fx_vec[N_MAX] =
    {
        0.05179370533242827E+00,
        0.2015098848089788E+00,
        0.4079953223292903E+00,
        0.3304762110867252E+00,
        0.5223047493549780E+00,
        0.3889503452643453E+00,
        0.1505614239732950E+00,
        0.03927689321042477E+00,
        0.00003099828465518108E+00,
        0.03145116093938197E+00,
        0.2114132170316862E+00,
        0.2075776621999210E+00,
        0.0000000000000000E+00,
        0.002088888139634505E+00,
        0.3876752992448843E+00,
        0.9135215248834896E+00
    };

    int n_vec[N_MAX] =
    {
        7,  8,  9, 10,
        6,  6,  6,  6,
        6,  6,  6,  6,
        0,  0,  0,  0
    };

    int pop_vec[N_MAX] =
    {
        100, 100, 100, 100,
        100, 100, 100, 100,
        100, 100, 100, 100,
        90,  200, 1000, 10000
    };

    int sam_vec[N_MAX] =
    {
        10, 10, 10, 10,
        6,  7,  8,  9,
        10, 10, 10, 10,
        10, 10, 10, 10
    };

    int suc_vec[N_MAX] =
    {
        90, 90, 90, 90,
        90, 90, 90, 90,
        10, 30, 50, 70,
        90, 90, 90, 90
    };

    if ( *n_data < 0 )
    {
        *n_data = 0;
    }

    *n_data = *n_data + 1;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *sam = 0;
        *suc = 0;
        *pop = 0;
        *n = 0;
        *fx = 0.0;
    }
    else
    {
        *sam = sam_vec[*n_data-1];
        *suc = suc_vec[*n_data-1];
        *pop = pop_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return;
# undef N_MAX
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
    int value;

    if ( i1 < i2 )
    {
        value = i1;
    }
    else
    {
        value = i2;
    }
    return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
    double value;

    if ( x < y )
    {
        value = y;
    }
    else
    {
        value = x;
    }
    return value;
}



//****************************************************************************80

// Below function is omitted, but the comments remain.
// void beta_inc_values ( int *n_data, double *a, double *b, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    The incomplete Beta function is also sometimes called the
//    "modified" Beta function, or the "normalized" Beta function
//    or the Beta CDF (cumulative density function.
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//    The function can also be evaluated by using the Statistics package:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
//
//  Modified:
//
//    04 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//


// Note that beta below is
// double a, b, x, beta_log;
// int ifault;
// beta_log = alngam(a, &ifault) + alngam(b, &ifault) - alngam(a+b, &ifault);
// fx = betain(x, a, b, beta_log, &ifault);
//****************************************************************************80

double betain ( double x, double p, double q, double beta, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    KL Majumder, GP Bhattacharjee
//    C++ version by John Burkardt
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
//
{
    double acu = 0.1E-14;
    double ai;
    // Below variable never used.
    // double betain;
    double cx;
    bool indx;
    int ns;
    double pp;
    double psq;
    double qq;
    double rx;
    double temp;
    double term;
    double value;
    double xx;

    value = x;
    *ifault = 0;
    //
    //  Check the input arguments.
    //
    if ( p <= 0.0 || q <= 0.0 )
    {
        *ifault = 1;
        return value;
    }

    if ( x < 0.0 || 1.0 < x )
    {
        *ifault = 2;
        return value;
    }
    //
    //  Special cases.
    //
    if ( x == 0.0 || x == 1.0 )
    {
        return value;
    }
    //
    //  Change tail if necessary and determine S.
    //
    psq = p + q;
    cx = 1.0 - x;

    if ( p < psq * x )
    {
        xx = cx;
        cx = x;
        pp = q;
        qq = p;
        indx = true;
    }
    else
    {
        xx = x;
        pp = p;
        qq = q;
        indx = false;
    }

    term = 1.0;
    ai = 1.0;
    value = 1.0;
    ns = ( int ) ( qq + cx * psq );
    //
    //  Use the Soper reduction formula.
    //
    rx = xx / cx;
    temp = qq - ai;
    if ( ns == 0 )
    {
        rx = xx;
    }

    for ( ; ; )
    {
        term = term * temp * rx / ( pp + ai );
        value = value + term;;
        temp = r8_abs ( term );

        if ( temp <= acu && temp <= acu * value )
        {
            value = value * exp ( pp * log ( xx )
                                  + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

            if ( indx )
            {
                value = 1.0 - value;
            }
            break;
        }

        ai = ai + 1.0;
        ns = ns - 1;

        if ( 0 <= ns )
        {
            temp = qq - ai;
            if ( ns == 0 )
            {
                rx = xx;
            }
        }
        else
        {
            temp = psq;
            psq = psq + 1.0;
        }
    }

    return value;
}







#endif

