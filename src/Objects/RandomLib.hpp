#ifndef MOLBIOLIB_RANDOMLIB_H
#define MOLBIOLIB_RANDOMLIB_H 1

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


#include "src/include/PrimitiveTypes.hpp"


/** \file RandomLib.hpp
 * Contains the #RandomLib and #Permutation classes.
 */


/** This is a random number generator class by Marsaglia and Zaman.
 * See http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/random
 * This was modified from the above to be in a self-contained class.
 * Original comments are below:
 *
 * This Random Number Generator is based on the algorithm in a FORTRAN
 * version published by George Marsaglia and Arif Zaman, Florida State
 * University; ref.: see original comments below.
 * At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
 * Science, we have written sources in further languages (C, Modula-2
 * Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
 * results compared with the original FORTRAN version.
 * April 1989
 * Karl-L. Noell <NOELL@DWIFH1.BITNET>
 *    and  Helmut  Weber <WEBER@DWIFH1.BITNET>
 *
 * This random number generator originally appeared in "Toward a Universal
 * Random Number Generator" by George Marsaglia and Arif Zaman.
 * Florida State University Report: FSU-SCRI-87-50 (1987)
 * It was later modified by F. James and published in "A Review of Pseudo-
 * random Number Generators"
 * THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
 * (However, a newly discovered technique can yield
 * a period of 10^600. But that is still in the development stage.)
 * It passes ALL of the tests for random number generators and has a period
 * of 2^144, is completely portable (gives bit identical results on all
 * machines with at least 24-bit mantissas in the floating point
 * representation).
 * The algorithm is a combination of a Fibonacci sequence (with lags of 97
 * and 33, and operation "subtraction plus one, modulo one") and an
 * "arithmetic sequence" (using subtraction).
 *
 * Use IJ = 1802 & KL = 9373 to test the random number generator. The
 * subroutine RANMAR should be used to generate 20000 random numbers.
 * Then display the next six random numbers generated multiplied by 4096*4096
 * If the random number generator is working properly, the random numbers
 * should be:
 *         6533892.0  14220222.0  7275067.0
 *         6172232.0  8354498.0   10633180.0
 *
 * \section Usage Usage:
 * \code
 * long seed1 = 1802, seed2 = 9373;
 * RandomLib rand1(seed1, seed2);  // Default are the seed values shown.
 * double unifX = rand1.randomUniform();   // In [0, 1].
 * double gaussX = rand1.randomGaussian(0.0, 1.0);  // mean and stddev.
 *                                                  // Defaults are shown.
 * long randX = rand1.randomLong(1, 10);  // Inclusive of 1 and 10.
 * size_t randX = rand1.randomSize_t(2, 5);  // Inclusive of 2 and 5.
 * double rangeX = rand1.randomDouble();  // In [0.0, 1.0].
 * double rangeX = rand1.randomDouble(2.5, 8.6);  // In [2.5, 8.6].
 * \endcode
 */
class RandomLib
{
public:
    RandomLib(long ijl = 1802, long kll = 9373)
    {
        int ij = static_cast<int>(ijl);
        int kl = static_cast<int>(kll);
        double s, t;
        int ii, i, j, k, l, jj, m;

        //   Handle the seed range errors
        //      First random number seed must be between 0 and 31328
        //      Second seed must have a value between 0 and 30081
        if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081)
        {
#ifdef DEBUG
            if ((ij < 0) || (ij > 31328))
            {
                cerr << "Warning!  The first seed passed to RandomLib, ij = " << ij << " must be between 0 and 31328.   Using the default seed values of ij = 1802 and kl = 9373 and continuing." << endl;
            }
            if ((kl < 0) || (kl > 30081))
            {
                cerr << "Warning!  The second seed passed to RandomLib, kl = " << kl << " must be between 0 and 30081.   Using the default seed values of ij = 1802 and kl = 9373 and continuing." << endl;
            }
#endif
            ij = 1802;
            kl = 9373;
        }

        i = (ij / 177) % 177 + 2;
        j = (ij % 177)       + 2;
        k = (kl / 169) % 178 + 1;
        l = (kl % 169);

        for (ii=0; ii<97; ii++)
        {
            s = 0.0;
            t = 0.5;
            for (jj=0; jj<24; jj++)
            {
                m = (((i * j) % 179) * k) % 179;
                i = j;
                j = k;
                k = m;
                l = (53 * l + 1) % 169;
                if (((l * m % 64)) >= 32)
                {
                    s += t;
                }
                t *= 0.5;
            }
            u[ii] = s;
        }

        c    = 362436.0 / 16777216.0;
        cd   = 7654321.0 / 16777216.0;
        cm   = 16777213.0 / 16777216.0;
        i97  = 97;
        j97  = 33;
    }


    //   This is the random number generator proposed by George Marsaglia in
    //   Florida State University Report: FSU-SCRI-87-50
    double randomUniform()
    {
        double uni;

        uni = u[i97-1] - u[j97-1];
        if (uni <= 0.0)
        {
            uni++;
        }
        u[i97-1] = uni;
        i97--;
        if (i97 == 0)
        {
            i97 = 97;
        }
        j97--;
        if (j97 == 0)
        {
            j97 = 97;
        }
        c -= cd;
        if (c < 0.0)
        {
            c += cm;
        }
        uni -= c;
        if (uni < 0.0)
        {
            uni++;
        }

        return uni;
    }


    // ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    // THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    // VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
    // The function returns a normally distributed pseudo-random number
    // with a given mean and standard devaiation.  Calls are made to a
    // function subprogram which must return independent random
    // numbers uniform in the interval (0,1).
    // The algorithm uses the ratio of uniforms method of A.J. Kinderman
    // and J.F. Monahan augmented with quadratic bounding curves.
    double randomGaussian(double mean = 0.0, double stddev = 1.0)
    {
        double q, u, v, x, y;

        // Generate P = (u,v) uniform in rect. enclosing acceptance region
        // Make sure that any random numbers <= 0 are rejected, since
        // gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
        do
        {
            u = randomUniform();
            v = randomUniform();
            if (u <= 0.0 || v <= 0.0)
            {
                u = 1.0;
                v = 1.0;
            }
            v = 1.7156 * (v - 0.5);

            // Evaluate the quadratic form
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = x * x + y * (0.19600 * y - 0.25472 * x);

            // Accept P if inside inner ellipse
            if (q < 0.27597)
            {
                break;
            }
            //  Reject P if outside outer ellipse, or outside acceptance region
        }
        while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

        //  Return ratio of P's coordinates as the normal deviate
        return (mean + stddev * v / u);
    }


    // Return random long within a range, lower -> upper INCLUSIVE
    long randomLong(long lower, long upper)
    {
        return (static_cast<long>(floor((randomUniform() * (upper - lower + 1.0)))) + lower);
    }

    size_t randomSize_t(size_t lower, size_t upper)
    {
        return (static_cast<size_t>(floor((randomUniform() * (upper - lower + 1.0)))) + lower);
    }

    // Return random float within a range, lower -> upper
    double randomDouble(double lower = 0.0, double upper = 1.0)
    {
        return ((upper - lower) * randomUniform() + lower);
    }

private:
    double u[97], c, cd, cm;
    int i97, j97;
};


/** Permutation vector
 * This computates a permutation vector (0, 1, ..., length-1) rearranged
 * in some random way.  Usage:
 * \code
 * Permutation thePerm(10, 1802, 9373);   // length 10, and the two seeds.
 * vector<size_t> newPerm = thePerm.nextPerm();
 * // ...
 * newPerm = thePerm.nextPerm(3);   // Do k swaps
 * \endcode
 * Other operators included are <code>size()</code>, <code>reset()</code>
 * (where the permutation vector is reset to 0, 1, ..., length-1),
 * <code>resize(newLength)</code> (where the permutation is resized and
 * then reset).
 *
 */
class Permutation : public RandomLib
{
public:
    size_t size()
    {
        return currPerm.size();
    }
    void reset()
    {
        for (size_t i = 0; i < size(); ++i)
        {
            currPerm[i] = i;
        }
    }
    void resize(size_t length)
    {
        currPerm.resize(length);
        reset();
    }
    Permutation(size_t length = 1, long ijl = 1802, long kll = 9373) : RandomLib(ijl, kll)
    {
        resize(length);
    }
    // Below - specify k-permutation as defined in
    // http://www.techuser.net/randpermgen.html
    vector<size_t> nextPerm(size_t k = numeric_limits<size_t>::max())
    {
        // If not specified, do entire vector complete.
        if (k == numeric_limits<size_t>::max())
        {
            k = size() - 1;
        }
        for (size_t i = 0; i < k; ++i)
        {
            swap(currPerm[i], currPerm[randomSize_t(0, size()-1)]);
        }
        return currPerm;
    }
private:
    vector<size_t> currPerm;
};


#endif

