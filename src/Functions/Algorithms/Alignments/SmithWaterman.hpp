#ifndef MOLBIOLIB_SMITHWATERMAN_H
#define MOLBIOLIB_SMITHWATERMAN_H 1

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


#include "src/Objects/Sequence.hpp"



/** \file SmithWaterman.hpp
 * Various implementations of the Smith-Waterman alignment algorithm.
 * \fn void smithWatermanAffineGapPenalty(Sequence& input1, Sequence& input2, Sequence& output1, Sequence& output2, double matchWeight = 2.0, double mismatchWeight = -1.0, double delWeight = -1.0, double delMult = -1.0/3.0, double insWeight = -1.0, double insMult = -1.0/3.0)
 * This version uses an affine gap penalty.  The insertion is relative to 
 * input2.  The output is two sequences with dashes in them to indicate the 
 * alignment.
 */
void smithWatermanAffineGapPenalty(Sequence& input1, Sequence& input2,
                                   Sequence& output1, Sequence& output2,
                                   double matchWeight = 2.0,
                                   double mismatchWeight = -1.0,
                                   double delWeight = -1.0,
                                   double delMult = -1.0/3.0,
                                   double insWeight = -1.0,
                                   double insMult = -1.0/3.0)
{
    size_t i1len = input1.length(), i2len = input2.length();

    // First, initialize  the matrix
    vector< vector<double> > H(i2len+1);
    // Move entries can be d = diagonal, l = left, u = up.
    vector< vector<char> > move(i2len+1);
    for (size_t i = 0; i < i2len+1; ++i)
    {
        H[i].resize(i1len+1, 0.0);
        move[i].resize(i1len+1, ' ');
    }

    // Now, fill the matrix.
    double tempM = 0.0, tempDel = 0.0, tempIns = 0.0;
    double tempVal = -1.0, bestVal = -1.0;
    for (size_t i = 1; i < i2len+1; ++i)
    {
        for (size_t j = 1; j < i1len+1; ++j)
        {
            if (input1[j-1] == input2[i-1])
            {
                tempM = matchWeight;
            }
            else
            {
                tempM = mismatchWeight;
            }
            tempM += H[i-1][j-1];
            bestVal = -1.0;
            tempVal = -1.0;
            for (size_t k = 1; k <= i; ++k)
            {
                tempVal = H[i-k][j] + delWeight + delMult*static_cast<double>(k);
                if (tempVal > bestVal)
                {
                    bestVal = tempVal;
                }
            }
            tempDel = bestVal;
            bestVal = -1.0;
            tempVal = -1.0;
            for (size_t k = 1; k <= j; ++k)
            {
                tempVal = H[i][j-k] + insWeight + insMult*static_cast<double>(k);
                if (tempVal > bestVal)
                {
                    bestVal = tempVal;
                }
            }
            tempIns = bestVal;
            if ((tempM < 0.0) && (tempDel < 0.0) && (tempIns < 0.0))
            {
                H[i][j] = 0.0;
                move[i][j] = ' ';
            }
            else if ((tempM >= tempDel) && (tempM >= tempIns))
            {
                H[i][j] = tempM;
                move[i][j] = 'd';
            }
            else if ((tempDel >= tempM) && (tempDel >= tempIns))
            {
                H[i][j] = tempDel;
                move[i][j] = 'u';
            }
            else
            {
                H[i][j] = tempIns;
                move[i][j] = 'l';
            }
        }
    }

    // Find biggest value in matrix.
    double maxVal = -1.0;   // No such value in matrix.
    size_t maxI = numeric_limits<size_t>::max(),
           maxJ = numeric_limits<size_t>::max();
    for (size_t i = 0; i < i2len+1; ++i)
    {
        for (size_t j = 0; j < i1len+1; ++j)
        {
            if (H[i][j] >= maxVal)
            {
                maxVal = H[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    output1 = "";
    output2 = "";
    // Start from biggest value and work back.
    size_t i = maxI, j = maxJ;
    while (move[i][j] != ' ')
    {
        if (move[i][j] == 'd')
        {
            output1 = input1[j-1] + output1;
            output2 = input2[i-1] + output2;
            --i;
            --j;
        }
        else if (move[i][j] == 'l')
        {
            output1 = input1[j-1] + output1;
            output2 = "-" + output2;
            --j;
        }
        else if (move[i][j] == 'u')
        {
            output1 = "-" + output1;
            output2 = input2[i-1] + output2;
            --i;
        }
    }
}


#endif

