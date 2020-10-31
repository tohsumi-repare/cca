#ifndef MOLBIOLIB_LOCALMAXIMUMWINDOW_H
#define MOLBIOLIB_LOCALMAXIMUMWINDOW_H 1

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


/** \file LocalMaximumWindow.hpp
 * Find maximums separated by at least windowLength.
 * \fn template<typename P, size_t N, size_t S, typename... T> void localMaximumWindow(Table<T...>& theTable, Peaks<P>& thePeaks, size_t windowLength, size_t start = 0, size_t end = numeric_limits<size_t>::max(), double minAmp = static_cast<double>(numeric_limits<long>::min()), string contig = "", string strand = "")
 * Maximum is assumed to be >= numeric_limits<long>::min().
 * Set contig and strand by passing in as parameters.
 * N is the column of T in which the value lies.  S is the column of T in which
 * the position lies.  Use index if S == N.
 * \todo template and template parameters not included in Doxygen function
 *       definition because Doxygen cannot handle complex parameters.
 */
template<typename P, size_t N, size_t S, typename... T>
void localMaximumWindow(Table<T...>& theTable, Peaks<P>& thePeaks,
                        size_t windowLength,
                        size_t start = 0,
                        size_t end = numeric_limits<size_t>::max(),
                        double minAmp = static_cast<double>(numeric_limits<long>::min()),
                        string contig = "",
                        string strand = "")
{

    if (end == numeric_limits<size_t>::max())
    {
        end = theTable.size() - 1;
    }

#ifdef DEBUG
    assert(start < theTable.size());
    assert(start <= end);
    assert(end < theTable.size());
    assert((end - start) >= windowLength);
#endif

    size_t i = start;
    while (i < end)
    {
        long maxPos = numeric_limits<long>::min();
        size_t maxIndex = numeric_limits<size_t>::max();
        double maxAmp = static_cast<double>(maxPos);
        size_t stop = min(i + windowLength - 1, end);
        for (size_t j = i; j <= stop; ++j)
        {
            double currAmp = static_cast<double>(get<N>(theTable.theData)[j]);
            if (currAmp > maxAmp)
            {
                maxAmp = currAmp;
                maxIndex = j;
                if (S == N)
                {
                    maxPos = static_cast<long>(j);
                }
                else
                {
                    maxPos = static_cast<long>(get<S>(theTable.theData)[j]);
                }
            }
        }
        if (maxAmp >= minAmp)
        {
            P newPeak;
            SequenceLocation newPos;
            newPos.contig() = contig;
            newPos.strand() = strand;
            newPos.position().assign(maxPos);
            newPeak.pos = newPos;
            newPeak.posIndex = maxIndex;
            newPeak.amp = maxAmp;
            thePeaks.push_back(newPeak);
            i = maxIndex + windowLength;
        }
        else
            // Can skip by windowLength since entire region has to be < minAmp.
        {
            i += windowLength;
        }
    }
}

#endif

