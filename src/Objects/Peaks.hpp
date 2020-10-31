#ifndef MOLBIOLIB_PEAKS_H
#define MOLBIOLIB_PEAKS_H 1

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


/** \file Peaks.hpp
 * The object to hold Peaks and its derived classes are defined.
 */

// Below loads Interval.hpp and Table.hpp
#include "src/Objects/Location.hpp"
#include "src/Functions/Algorithms/Tables/SortTable.hpp"
#include "src/Functions/Transformers/TableOperations/DistinctTable.hpp"


class PeakPos
{
public:
    SequenceLocation pos;
    size_t posIndex;
    string header()
    {
        return "Contig\tStrand\tPos";
    }
    string print()
    {
        return (pos.contig() + "\t" + pos.strand() + "\t" +
                convertToString<long>(pos.position().lower()));
    }
    PeakPos()
    {
        posIndex = numeric_limits<size_t>::max();
    }
};


class PeakWidth : public PeakPos
{
public:
    SequenceInterval width;
    Interval<size_t> widthIndex;
    string header()
    {
        return (PeakPos::header() + "\tWidth");
    }
    string print()
    {
        return (PeakPos::print() +
                "\t[" +
                convertToString<long>(width.lower()) + "," +
                convertToString<long>(width.upper()) +"]");
    }
    PeakWidth()
    {
        width.assign(0, 0);
        widthIndex.assign(0, 0);
    }
};


class PeakAmp : public PeakWidth
{
public:
    double amp;
    string header()
    {
        return (PeakWidth::header() + "\tAmplitude");
    }
    string print()
    {
        return (PeakWidth::print() + "\t" +
                convertToString<double>(amp));
    }
    PeakAmp()
    {
        amp = 0.0;
    }
};


class PeakRidgeLines : public PeakAmp
{
public:
    Table<size_t> ridgeIndices;
    size_t numGaps;
    size_t maxLevel;
    size_t levelLength;
    string header()
    {
        return (PeakAmp::header() +
                "\tLevel length\tNumber of gaps");
    }
    string print()
    {
        return (PeakAmp::print() + "\t" +
                convertToString<double>(levelLength) + "\t" +
                convertToString<double>(numGaps));
    }
    PeakRidgeLines()
    {
        numGaps = 0;
        maxLevel = 0;
        levelLength = 0;
    }
};



template <typename P>
bool compareTablePeaksPos(typename Table<P>::row_type& r1,
                          typename Table<P>::row_type& r2)
{
    if (get<0>(r1).pos.contig() != get<0>(r2).pos.contig())
    {
        return (get<0>(r1).pos.contig() < get<0>(r2).pos.contig());
    }
    else
        return (get<0>(r1).pos.position().lower() <
                get<0>(r2).pos.position().lower());
}


template <typename P>
bool compareTablePeaksPosIndex(typename Table<P>::row_type& r1,
                               typename Table<P>::row_type& r2)
{
    return (get<0>(r1).posIndex < get<0>(r2).posIndex);
}


template <typename P>
class Peaks
{
public:
    Table<P> peaks;   // public since deleters of peaks can directly access this.
    void sortPos(bool comparePos = true)
    {
        if (comparePos)
        {
            sortTable(peaks, compareTablePeaksPos<P>);
        }
        else // compare using indices
        {
            sortTable(peaks, compareTablePeaksPosIndex<P>);
        }
    }
    P& operator[](size_t i)
    {
#ifdef PROG_DEBUG
        assert(i < peaks.size());
#endif
        return peaks[i];
    }
    size_t size()
    {
        return peaks.size();
    }
    string header()
    {
        P temp;
        return temp.header();
    }
    void push_back(P& newPeak)
    {
        peaks.push_back(newPeak);
    }
    void clear()
    {
        peaks.clear();
    }
};





template <typename P>
class StrandedPeaks
{
public:
    Peaks<P> peaksPlus, peaksMinus;
};




#endif

