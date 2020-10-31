#ifndef MOLBIOLIB_PEAKPERCENTAMPWIDTH_H
#define MOLBIOLIB_PEAKPERCENTAMPWIDTH_H 1

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


/** \file PeakPercentAmpWidth.hpp
 * Sets the width of peaks based on some percentage of the maximum amplitude.
 * \fn void peakPercentAmpWidth(Table<T...>& theTable, P& thePeak, double percent, size_t start = 0, size_t end = numeric_limits<size_t>::max())
 * If S = N, then do not use S for position - use index.
 * percent goes from 0 to 1.
 * \todo template and template parameters not included in Doxygen function
 *       definition because Doxygen cannot handle complex parameters.
 */
template<typename P, size_t N, size_t S, typename... T>
void peakPercentAmpWidth(Table<T...>& theTable, P& thePeak,
                         double percent,
                         size_t start = 0,
                         size_t end = numeric_limits<size_t>::max())
{
#ifdef DEBUG
    assert (percent <= 1.0);
#endif

    size_t currIndex = thePeak.posIndex;
    double maxAmp = thePeak.amp;
    // Used to check here in DEBUG
    // assert(static_cast<double>(get<N>(theTable.theData)[currIndex]) == maxAmp);
    // but this is a wrong assumption because the widths could be set to the
    // smoothed coverage table, but the position is set to the original
    // unsmoothed coverage table.
    double targetAmp = percent*maxAmp;

    size_t leftIndex = start;
    long leftPos = numeric_limits<long>::min();
    if (currIndex < start)
    {
#ifdef DEBUG
        cerr << "Warning, no left side found in peakPercentAmpWidth.  Using start = " << start << "." << endl;
#endif
    }
    else
    {
        bool stop = false;
        leftIndex = currIndex - 1;
        while (!stop && (leftIndex >= start))
        {
            if (static_cast<double>(get<N>(theTable.theData)[leftIndex]) <=
                    targetAmp)
            {
                stop = true;
            }
            else
            {
                --leftIndex;
            }
        }
#ifdef DEBUG
        if (!stop)
        {
            cerr << "Warning, no left side found in peakPercentAmpWidth.  Using start = " << start << "." << endl;
        }
#endif
    }
    if (S == N)
    {
        leftPos = static_cast<long>(leftIndex);
    }
    else
    {
        leftPos = static_cast<long>(get<S>(theTable.theData)[leftIndex]);
    }


    size_t rightIndex = end;
    long rightPos = numeric_limits<long>::min();
    if (currIndex > end)
    {
#ifdef DEBUG
        cerr << "Warning, no right side found in peakPercentAmpWidth.  Using end = " << end << "." << endl;
#endif
    }
    else
    {
        bool stop = false;
        rightIndex = currIndex + 1;
        while (!stop && (rightIndex <= end))
        {
            if (static_cast<double>(get<N>(theTable.theData)[rightIndex]) <=
                    targetAmp)
            {
                stop = true;
            }
            else
            {
                ++rightIndex;
            }
        }
#ifdef DEBUG
        if (!stop)
        {
            cerr << "Warning, no right side found in peakPercentAmpWidth.  Using end = " << end << "." << endl;
        }
#endif
    }
    if (S == N)
    {
        rightPos = static_cast<long>(rightIndex);
    }
    else
    {
        rightPos = static_cast<long>(get<S>(theTable.theData)[rightIndex]);
    }

    SequenceInterval tempWidth(leftPos, rightPos);
    Interval<size_t> tempWidthIndex(leftIndex, rightIndex);
    thePeak.width = tempWidth;
    thePeak.widthIndex = tempWidthIndex;
}



/** \file PeakPercentAmpWidth.hpp
 * Sets the width of peaks based on some percentage of the maximum amplitude.
 * \fn void peakPercentAmpWidth(Table<T...>& theTable, Peaks<P>& thePeaks, double percent, size_t start = 0, size_t end = numeric_limits<size_t>::max(), string contig = "", string strand = "")
 * start and end limits window in which to look for width.
 * If one is using the table's position values, first determine from that
 * what the index start and end are before calling this function.
 * If contig and/or strand is specified, only change width on those.
 * If S = N, then do not use S for position - use index.
 * percent goes from 0 to 1.
 * \todo template and template parameters not included in Doxygen function
 *       definition because Doxygen cannot handle complex parameters.
 */
template<typename P, size_t N, size_t S, typename... T>
void peakPercentAmpWidth(Table<T...>& theTable, Peaks<P>& thePeaks,
                         double percent,
                         size_t start = 0,
                         size_t end = numeric_limits<size_t>::max(),
                         string contig = "",
                         string strand = "")
{
    if (end == numeric_limits<size_t>::max())
    {
        end = theTable.size() - 1;
    }

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
        peakPercentAmpWidth<P, N, S>(theTable, thePeaks[i], percent, start, end);
    }
}

#endif



