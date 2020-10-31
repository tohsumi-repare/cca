#ifndef MOLBIOLIB_CORRECTPEAKPOSITION_H
#define MOLBIOLIB_CORRECTPEAKPOSITION_H 1

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




/** \file CorrectPeakPosition.hpp
 * Correct peak position based on some coverage table.
 * \fn bool bestLocalPeak(P& currPeak, Table<T...>& theTable, size_t bumpLength = 1, size_t start = 0, size_t end = numeric_limits<size_t>::max())
 * Find best local peak.
 */
template<typename P, size_t N, size_t S, typename... T>
bool bestLocalPeak(P& currPeak, Table<T...>& theTable,
                   size_t bumpLength = 1,
                   size_t start = 0,
                   size_t end = numeric_limits<size_t>::max())
{
    if (end == numeric_limits<size_t>::max())
    {
        end = theTable.size() - 1;
    }

    size_t& peakIndex = currPeak.posIndex;
    long peakPos = currPeak.pos.position().lower();
    bool gotBetter = false;
    double bestVal = currPeak.amp;
    size_t bestIndex = peakIndex;
    size_t currIndex = peakIndex;
    bool exitNow = false;   // Only needed so do not try to do a decrement
    // on a size_t variable set to 0.
    if (N == S)
    {
        // Use index for bumpLength
        size_t numBump = 1;
        // First look to left.  start is always >= 0, so can do decrement.
        if (currIndex > start)
        {
            --currIndex;
        }
        while (!exitNow && (currIndex >= start) && (numBump <= bumpLength))
        {
            if (static_cast<double>(get<N>(theTable.theData)[currIndex]) >
                    bestVal)
            {
                bestVal = static_cast<double>(get<N>(theTable.theData)[currIndex]);
                bestIndex = currIndex;
                gotBetter = true;
            }
            if (currIndex == 0)
            {
                exitNow = true;
            }
            else
            {
                --currIndex;
            }
            ++numBump;
        }
        // Next look to right.
        numBump = 1;
        currIndex = peakIndex;
        // Below: unlikely we need to check against size_t's max().
        if (currIndex < end)
        {
            ++currIndex;
        }
        while ((currIndex <= end) && (numBump <= bumpLength))
        {
            if (static_cast<double>(get<N>(theTable.theData)[currIndex]) > bestVal)
            {
                bestVal = static_cast<double>(get<N>(theTable.theData)[currIndex]);
                bestIndex = currIndex;
                gotBetter = true;
            }
            ++currIndex;   // Again, not likely we need not check against max().
            ++numBump;
        }
    }
    else
    {
        // Use position for bumpLength
        long bumpLengthLong = static_cast<long>(bumpLength);
        // Look left
        if (currIndex > start)
        {
            --currIndex;
        }
        long bumpDist = labs(peakPos -
                             static_cast<long>(get<S>(theTable.theData)[currIndex]));
        while (!exitNow && (currIndex >= start) && (bumpDist <= bumpLengthLong))
        {
            if (static_cast<double>(get<N>(theTable.theData)[currIndex]) >
                    bestVal)
            {
                bestVal = get<N>(theTable.theData)[currIndex];
                bestIndex = currIndex;
                gotBetter = true;
            }
            if (currIndex == 0)
            {
                exitNow = true;
            }
            else
            {
                --currIndex;
            }
            bumpDist = labs(peakPos -
                            static_cast<long>(get<S>(theTable.theData)[currIndex]));
        }
        // Look right
        currIndex = peakIndex;
        if (currIndex < end)
        {
            ++currIndex;
        }
        bumpDist = labs(peakPos -
                        static_cast<long>(get<S>(theTable.theData)[currIndex]));
        while ((currIndex <= end) && (bumpDist <= bumpLengthLong))
        {
            if (static_cast<double>(get<N>(theTable.theData)[currIndex]) >
                    bestVal)
            {
                bestVal = static_cast<double>(get<N>(theTable.theData)[currIndex]);
                bestIndex = currIndex;
                gotBetter = true;
            }
            ++currIndex;
            if (currIndex < theTable.size())
                bumpDist = labs(peakPos -
                                static_cast<long>(get<S>(theTable.theData)[currIndex]));
        }
    }

    if (gotBetter)
    {
        // Found a higher valued peak than the current one, so replace.
        currPeak.posIndex = bestIndex;
        currPeak.amp = bestVal;
        if (N == S)
        {
            currPeak.pos.position().assign(static_cast<long>(bestIndex));
        }
        else
        {
            currPeak.pos.position().assign(static_cast<long>(get<S>(theTable.theData)[bestIndex]));
        }
    }

    return gotBetter;
}



/** \file CorrectPeakPosition.hpp
 * Correct peak position based on some coverage table.
 * \fn template<typename P, size_t N, size_t S, typename... T> long correctPeakPosition(P& aPeak, Table<T...>& theTable, size_t start = 0, size_t end = numeric_limits<size_t>::max(), bool findMiddleOfPlateau = true, size_t bumpLength = 1)
 * It is assumed P has an amplitude.
 * start and end fixes region in which to correct the value.
 * If one is using the table's position values, first determine from that
 * what the index start and end are before calling this function.
 * However, if one specifies S, then the new location is based on column S in
 * the table.
 * If contig and/or strand is specified, only change width on those.
 * Return value is the amount of shift a peak has in the table.
 * N is the column of T for the values and S (if != N) is the column of
 * T for the position (affects the return value) else use index.
 * This only changes the peak's position.
 * findMiddleOfPlateau will shift the maximum to the middle of a maximal
 * region where it is flat.
 * bumpLength is the distance (depends on S - index or length) to see if
 * there is a value increasing.
 */
template<typename P, size_t N, size_t S, typename... T>
long correctPeakPosition(P& aPeak,
                         Table<T...>& theTable,
                         size_t start = 0,
                         size_t end = numeric_limits<size_t>::max(),
                         bool findMiddleOfPlateau = true,
                         size_t bumpLength = 1)
{

    long orgPeakPos = static_cast<long>(aPeak.posIndex);   // N == S.
    if (N != S)
    {
        orgPeakPos = aPeak.pos.position().lower();
    }

    // Keep calling bestLocalPeak, until can't find anything better.
    while (bestLocalPeak<P, N, S>(aPeak, theTable, bumpLength, start, end)) { }

    if (findMiddleOfPlateau)
    {
        size_t currIndex = aPeak.posIndex;
        size_t leftIndex = currIndex, rightIndex = currIndex;
        double currAmp = aPeak.amp;
        if (currIndex > start)
        {
            --currIndex;
        }
        bool exitNow = false;
        while (!exitNow & (currIndex >= start) &&
                (currAmp == static_cast<double>(get<N>(theTable.theData)[currIndex])))
        {
            leftIndex = currIndex;
            if (currIndex == 0)
            {
                exitNow = true;
            }
            else
            {
                --currIndex;
            }
        }
        currIndex = aPeak.posIndex;
        if (currIndex < end)
        {
            ++currIndex;
        }
        while ((currIndex <= end) &&
                (currAmp == static_cast<double>(get<N>(theTable.theData)[currIndex])))
        {
            rightIndex = currIndex;
            ++currIndex;
        }
        aPeak.posIndex = (leftIndex + rightIndex)/2;
        if (N == S)
            // Use index.
        {
            aPeak.pos.position().assign(static_cast<long>(aPeak.posIndex));
        }
        else
        {
            aPeak.pos.position().assign(static_cast<long>((get<S>(theTable.theData)[leftIndex] + get<S>(theTable.theData)[rightIndex])/2));
        }
    }

    return (orgPeakPos - aPeak.pos.position().lower());
}


/** \file CorrectPeakPosition.hpp
 * Correct peaks position based on some coverage table.
 * \fn template<typename P, size_t N, size_t S, typename... T> void correctPeakPosition(Peaks<P>& thePeaks, Table<T...>& theTable, size_t start = 0, size_t end = numeric_limits<size_t>::max(), bool findMiddleOfPlateau = true, size_t bumpLength = 1, string contig = "", string strand = "")
 * See documentation for correctPeakPosition.  This applies correctPeakPosition
 * for all peaks in thePeaks.  The return result of correctPeakPosition is not
 * retained.  It is assumed this will be called before any peak width calling
 * methods are used.  If contig and/or strand is specified, then only correct
 * those on that particular contig and/or strand.
 */
template<typename P, size_t N, size_t S, typename... T>
void correctPeakPosition(Peaks<P>& thePeaks,
                         Table<T...>& theTable,
                         size_t start = 0,
                         size_t end = numeric_limits<size_t>::max(),
                         bool findMiddleOfPlateau = true,
                         size_t bumpLength = 1,
                         string contig = "",
                         string strand = "")
{
    for (size_t i = 0; i < thePeaks.size(); ++i)
    {
        if ((contig != "") && (thePeaks[i].pos.contig() != contig))
        {
            continue;
        }
        if ((strand != "") && (thePeaks[i].pos.strand() != strand))
        {
            continue;
        }
        correctPeakPosition<P, N, S>(thePeaks[i], theTable,
                                     start, end,
                                     findMiddleOfPlateau, bumpLength);
    }
}


/** \file CorrectPeakPosition.hpp
 * Correct peaks position based on some coverage table.
 * \fn template<typename P, size_t N, size_t S> size_t removeDuplicatePeaks(Peaks<P>& thePeaks, string contig = "", string strand = "")
 * This removes in a list of peaks any duplicate peaks.  Duplicates can occur
 * after correcting peak positions.
 * If S == N, then the posIndex is compared.
 * If S != N, the position is compared.
 * If contig and/or strand is specified, then only look for peaks on said
 * contig and/or strand.
 * Returns the number of peaks removed.
 * \section ProgNotes Programming Notes
 * Thought about instead of two template integers, use bool template parameter,
 * but this set of arguments is more consistent with what's being passed around
 * for Peaks.  Keep this way.
 */
template<typename P, size_t S, size_t N>
size_t removeDuplicatePeaks(Peaks<P>& thePeaks,
                            string contig = "",
                            string strand = "")
{
    set<size_t> posSetSize_t;
    set<long> posSetLong;
    IndexTable delVector;
    for (size_t i = 0; i < thePeaks.size(); ++i)
    {
        if ((contig != "") && (thePeaks[i].pos.contig() != contig))
        {
            continue;
        }
        if ((strand != "") && (thePeaks[i].pos.strand() != strand))
        {
            continue;
        }
        if (S == N)
        {
            size_t currPos = thePeaks[i].posIndex;
            if (posSetSize_t.find(currPos) == posSetSize_t.end())
            {
                posSetSize_t.insert(currPos);
            }
            else
            {
                delVector.push_back(i);
            }
        }
        else
        {
            long currPos = thePeaks[i].pos.position().lower();
            if (posSetLong.find(currPos) == posSetLong.end())
            {
                posSetLong.insert(currPos);
            }
            else
            {
                delVector.push_back(i);
            }
        }
    }
    vectorDeleteTable<P>(thePeaks.peaks, delVector);
    return delVector.size();
}



#endif


