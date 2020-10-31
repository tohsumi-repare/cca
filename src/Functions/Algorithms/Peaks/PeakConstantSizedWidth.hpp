#ifndef MOLBIOLIB_PEAKCONSTANTSIZEDWIDTH_H
#define MOLBIOLIB_PEAKCONSTANTSIZEDWIDTH_H 1

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


/** \file PeakConstantSizedWidth.hpp
 * Set peak widths to be some constant value.  Does *not* set widthIndex.
 * \fn template<typename P> void peakConstantSizedWidth(P& thePeaks, long halfSize = 0, string contig = "", string strand = "")
 * If contig and/or strand is specified, only change width on those.
 */
template<typename P>
void peakConstantSizedWidth(P& thePeak, long halfSize,
                            string contig = "",
                            string strand = "")
{
    long currPos = thePeak.pos.position().lower();
    SequenceInterval tempI(currPos - halfSize, currPos + halfSize);
    thePeak.width = tempI;
}


/** \file PeakConstantSizedWidth.hpp
 * Set peak widths to be some constant value.  Does *not* set widthIndex.
 * \fn template<typename P> void peakConstantSizedWidth(Peaks<P>& thePeaks, long halfSize = 0, string contig = "", string strand = "")
 * If contig and/or strand is specified, only change width on those.
 */
template<typename P>
void peakConstantSizedWidth(Peaks<P>& thePeaks, long halfSize,
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
        peakConstantSizedWidth<P>(thePeaks[i], halfSize, contig, strand);
    }
}


#endif

