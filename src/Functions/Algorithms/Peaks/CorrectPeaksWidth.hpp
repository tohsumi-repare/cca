#ifndef MOLBIOLIB_CORRECTPEAKSWIDTH_H
#define MOLBIOLIB_CORRECTPEAKSWIDTH_H 1

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


#include "src/Objects/Peaks.hpp"




/** \file CorrectPeaksWidth.hpp
 * Correct peaks' width based on some coverage table.
 * \fn template<typename P, size_t N, size_t S, typename... T> bool midPointWidthAdjustment(Peaks<P>& thePeaks, Table<T...>& theTable, string contig = "", string strand = "", bool hasPosIndex = true, bool hasPosValue = true)
 * Adjust adjacent overlapping widths by making the endpoint of the widths
 * at the midpoint between the two peaks.
 * N is the column of T in which the amplitudes lie.  Amplitudes really not
 * needed.  N needed for S == N comparison.
 * If S == N, use index as distance, else use column S in T.
 * If S != N, it is still assumed the widthIndex will still have to be
 * adjusted - this may be slow, as we have to find the closest index to the
 * midpoint between the two peaks.  It is also assumed posIndex has been set,
 * unless hasPosIndex is set false.  If hasPosValue is set to false, do not
 * set the width, only widthIndex.
 * If specify contig and/or strand, then only work on the peaks with said
 * contig and/or strand.
 * Returns true if any adjustment took place.  Else return false.
 */
template<typename P, size_t N, size_t S, typename... T>
bool midPointWidthAdjustment(Peaks<P>& thePeaks, Table<T...>& theTable,
                             string contig = "",
                             string strand = "",
                             bool hasPosIndex = true,
                             bool hasPosValue = true)
{

#ifdef DEBUG
    assert(hasPosIndex || hasPosValue);
#endif

    bool hasAdjusted = false;
    size_t numPeaks = thePeaks.size();
    if ((numPeaks == 0) || (numPeaks == 1))
    {
        return hasAdjusted;
    }
    for (size_t i = 0; i < (numPeaks - 1); ++i)
    {
        if ((contig != "") && (thePeaks[i].pos.contig() != contig))
        {
            continue;
        }
        if ((strand != "") && (thePeaks[i].pos.strand() != strand))
        {
            continue;
        }
        bool foundOverlap = false;
        size_t breakPointSize_t = numeric_limits<size_t>::max();
        long breakPointLong = numeric_limits<long>::max();
        if (S != N)
        {
#ifdef DEBUG
            assert(hasPosValue);
#endif
            long rightCurr = thePeaks[i].width.upper();
            long leftNext = thePeaks[i+1].width.lower();
            if (leftNext <= rightCurr)
            {
                foundOverlap = true;
                // Below is the tentative break-point.  Need to adjust based on
                // what position values are in theTable, if have widthIndex.
                long rightCurrPos = thePeaks[i].pos.position().lower();
                long leftNextPos = thePeaks[i].pos.position().lower();
                breakPointLong = (rightCurrPos + leftNextPos)/2;
                if (hasPosIndex)
                {
                    size_t rightCurrIndex = thePeaks[i].posIndex;
                    size_t leftNextIndex = thePeaks[i+1].posIndex;
#ifdef PROG_DEBUG
                    assert(leftNextIndex <= rightCurrIndex);
#endif
                    size_t bestIndex = numeric_limits<size_t>::max();
                    long bestBreak = numeric_limits<long>::max();
                    long bestBreakDist = numeric_limits<long>::max();
                    for (size_t j = leftNextIndex; j <= rightCurrIndex; ++j)
                    {
                        long breakDist = labs(breakPointLong -
                                              static_cast<long>((get<S>(theTable.theData))[j]));
                        if (breakDist < bestBreakDist)
                        {
                            bestBreakDist = breakDist;
                            bestIndex = j;
                            bestBreak = static_cast<long>((get<S>(theTable.theData))[j]);
                        }
                    }
                    breakPointLong = bestBreak;
                    breakPointSize_t = bestIndex;
                }
            }
        }
        else
        {
#ifdef DEBUG
            assert(hasPosIndex);
#endif
            long rightCurr = thePeaks[i].widthIndex.upper();
            long leftNext = thePeaks[i+1].widthIndex.lower();
            if (leftNext <= rightCurr)
            {
                foundOverlap = true;
                size_t rightPosCurr = thePeaks[i].posIndex;
                size_t leftPosNext = thePeaks[i+1].posIndex;
                breakPointSize_t = (rightPosCurr + leftPosNext)/2;
                breakPointLong = static_cast<long>(breakPointSize_t);
            }
        }
        if (foundOverlap)
        {
            hasAdjusted = true;
            if (hasPosValue)
            {
                thePeaks[i].width.assign(thePeaks[i].width.lower(),
                                         breakPointLong);
                thePeaks[i+1].width.assign(breakPointLong + 1,
                                           thePeaks[i+1].width.upper());
            }
            if (hasPosIndex)
            {
                thePeaks[i].widthIndex.assign(thePeaks[i].widthIndex.lower(),
                                              breakPointSize_t);
                thePeaks[i+1].widthIndex.assign(breakPointSize_t + 1,
                                                thePeaks[i+1].widthIndex.upper());
            }
        }
    }
    return hasAdjusted;
}



#endif

