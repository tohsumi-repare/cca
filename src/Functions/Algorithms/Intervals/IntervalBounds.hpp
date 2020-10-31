#ifndef MOLBIOLIB_INTERVALBOUNDS_H
#define MOLBIOLIB_INTERVALBOUNDS_H 1

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


#include "src/Objects/Interval.hpp"
#include "src/Functions/Algorithms/Tables/TableBounds.hpp"


#ifndef INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING
/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \def INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING
 * Define minimum number of elements for lower and upper bound computing.
 * Otherwise, just give first and last elements as bounds.
 */ 
#define INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING 15
#endif



/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn bool intervalsTableLowerBoundsCompare(Interval<I>& x, Interval<I>& y)
 * Comparator to pass to lower_boundIndex in lower_boundIntervalsIndex.
 */
template<typename I>
bool intervalsTableLowerBoundsCompare(Interval<I>& x, Interval<I>& y)
{
    if (x.upper() < y.lower())
    {
        return true;
    }
    else
    {
        return false;
    }
}




/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn size_t lower_boundIntervalsIndex(vector< Interval<I> >& theIntervals, I& maxSize, I& value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool isUnsignedI = false, I unity = static_cast<I>(1), I zero = static_cast<I>(0))
 * Finds the lower bound index on a table of intervals.  This is the index for 
 * the first element in the sorted range [first,last) which does not compare 
 * less than value.  
 * \section Usage Usage:
 * \code
 * vector< Interval<long> > intervals;
 * ...
 * long maxSizeI = maxWidth(intervals) + 1;
 * size_t low = lower_boundIntervalsIndex<long>(intervals, maxSizeI, 55);
 * \endcode
 * If the vector size is zero, numeric_limits<size_t>::max() is returned.  One
 * must check that the first element is indeed >= value, as 0 is the lowest
 * return value, if 0 is returned.
 * If the type I is unsigned, one must pass true in the argument or else an
 * underflow wraparound error may occur.
 * last is one plus the last element to be used.
 *
 */
template<typename I>
size_t lower_boundIntervalsIndex(vector< Interval<I> >& theIntervals,
                                 I& maxSize,
                                 I& value,
                                 size_t first = 0,
                                 size_t last = numeric_limits<size_t>::max(),
                                 bool isUnsignedI = false,
                                 I unity = static_cast<I>(1),
                                 I zero = static_cast<I>(0))
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theIntervals.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    if (vecLength < INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING)
    {
        return 0;
    }

    Interval<I> theValue(value, value + unity);
    result = lower_boundIndex(theIntervals, theValue, first, last, 
                              intervalsTableLowerBoundsCompare);
    if (result == 0)
    {
        return result;
    }

    // Search for lower using maxSize.
    I lowest = theIntervals[result].lower() - maxSize;
    if ((isUnsignedI) && (theIntervals[result].lower() < maxSize))
    {
        lowest = zero;
    }
    size_t i = result - 1;
    bool breakNow = false;
    while (!breakNow && (theIntervals[i].upper() >= lowest))
    {
        if (in(value, theIntervals[i]))
        {
            result = i;
        }
        if (i == first) 
        {
            breakNow = true;
        }
        else
        {
            --i;
        }
    }

    return result;
}


/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn size_t lower_boundIntervalsIndex(Table< Interval<I> >& theIntervals, I& maxSize, I& value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool isUnsignedI = false, I unity = static_cast<I>(1), I zero = static_cast<I>(0))
 * Finds the lower bound index on a table of intervals.  This is the index for 
 * the first element in the sorted range [first,last) which does not compare 
 * less than value.  
 * \section Usage Usage:
 * \code
 * Table< Interval<long> > intervals;
 * ...
 * long maxSizeI = maxWidth(intervals) + 1;
 * size_t low = lower_boundIntervalsIndex<long>(intervals, maxSizeI, 55);
 * \endcode
 * If the table size is zero, numeric_limits<size_t>::max() is returned.  One
 * must check that the first element is indeed >= value, as 0 is the lowest
 * return value, if 0 is returned.
 * If the type I is unsigned, one must pass true in the argument or else an
 * underflow wraparound error may occur.
 * last is one plus the last element to be used.
 *
 */
template<typename I>
size_t lower_boundIntervalsIndex(Table< Interval<I> >& theIntervals,
                                 I& maxSize,
                                 I& value,
                                 size_t first = 0,
                                 size_t last = numeric_limits<size_t>::max(),
                                 bool isUnsignedI = false,
                                 I unity = static_cast<I>(1),
                                 I zero = static_cast<I>(0))
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theIntervals.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    if (vecLength < INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING)
    {
        return 0;
    }

    Interval<I> theValue(value, value + unity);
    result = lower_boundIndex<0>(theIntervals, theValue, first, last, 
                                 intervalsTableLowerBoundsCompare);
    if (result == 0)
    {
        return result;
    }

    // Search for lower using maxSize.
    I lowest = theIntervals[result].lower() - maxSize;
    if ((isUnsignedI) && (theIntervals[result].lower() < maxSize))
    {
        lowest = zero;
    }
    size_t i = result - 1;
    bool breakNow = false;
    while (!breakNow && (theIntervals[i].upper() >= lowest))
    {
        if (in(value, theIntervals[i]))
        {
            result = i;
        }
        if (i == first) 
        {
            breakNow = true;
        }
        else
        {
            --i;
        }
    }

    return result;
}








/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn bool intervalsTableUpperBoundsCompare(Interval<I>& x, Interval<I>& y)
 * Comparator to pass to lower_boundIndex in lower_boundIntervalsIndex.
 */
template<typename I>
bool intervalsTableUpperBoundsCompare(Interval<I>& x, Interval<I>& y)
{
    if (x.upper() < y.lower())
    {
        return true;
    }
    else
    {
        return false;
    }
}

/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn size_t upper_boundIntervalsIndex(vector< Interval<I> >& theIntervals, I& maxSize, I& value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool isUnsignedI = false, I unity = static_cast<I>(1))
 * Finds the upper bound index on a table.  This is the index for the first
 * element in the sorted range [first,last) which compares greater than value.
 * One must check that the first element is indeed > value, as 0 is the lowest
 * return value, if 0 is returned.
 * \section Usage Usage:
 * \code
 * vector< Interval<long> > intervals;
 * ...
 * long maxSizeI = maxWidth(intervals) + 1;
 * size_t low = upper_boundIntervalsIndex<long>(intervals, maxSizeI, 55);
 * \endcode
 * If the table size is zero, numeric_limits<size_t>::max() is returned.
 * Optionally, one may pass in a comparison function which compares the values.
 * The comparison is equivalent to operator<.
 * last is one plus the last element to be used.
 */
template<typename I>
size_t upper_boundIntervalsIndex(vector< Interval<I> >& theIntervals,
                                 I& maxSize,
                                 I& value,
                                 size_t first = 0,
                                 size_t last = numeric_limits<size_t>::max(),
                                 bool isUnsignedI = false,
                                 I unity = static_cast<I>(1))
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theIntervals.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    if (vecLength < INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING)
    {
        return vecLength;
    }

    if (last == numeric_limits<size_t>::max())
    {
        last = vecLength;
    }

    Interval<I> theValue(value, value + unity);
    result = upper_boundIndex(theIntervals, theValue, first, last, 
                              intervalsTableUpperBoundsCompare);

    // Search for lower using maxSize.
    I largest = theIntervals[result].upper() + maxSize;
    size_t i = result;
    while ((i < last) && (theIntervals[i].lower() <= largest))
    {
        if (in(value, theIntervals[i]))
        {
            result = i + 1;
        }
        ++i;
    }


    return result;
}


/** \file IntervalBounds.hpp
 * C++ analogs of lower_bound and upper_bound for a vector and Table of Intervals.
 * \fn size_t upper_boundIntervalsIndex(Table< Interval<I> >& theIntervals, I& maxSize, I& value, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool isUnsignedI = false, I unity = static_cast<I>(1))
 * Finds the upper bound index on a table.  This is the index for the first
 * element in the sorted range [first,last) which compares greater than value.
 * One must check that the first element is indeed > value, as 0 is the lowest
 * return value, if 0 is returned.
 * \section Usage Usage:
 * \code
 * Table< Interval<long> > intervals;
 * ...
 * long maxSizeI = maxWidth(intervals) + 1;
 * size_t low = upper_boundIntervalsIndex<long>(intervals, maxSize, 55);
 * \endcode
 * If the table size is zero, numeric_limits<size_t>::max() is returned.
 * Optionally, one may pass in a comparison function which compares the values.
 * The comparison is equivalent to operator<.
 * last is one plus the last element to be used.
 */
template<typename I>
size_t upper_boundIntervalsIndex(Table< Interval<I> >& theIntervals,
                                 I& maxSize,
                                 I& value,
                                 size_t first = 0,
                                 size_t last = numeric_limits<size_t>::max(),
                                 bool isUnsignedI = false,
                                 I unity = static_cast<I>(1))
{
    size_t result = numeric_limits<size_t>::max();
    size_t vecLength = theIntervals.size();
    if (vecLength == 0)
    {
        return result;    // No entries in table.
    }
    if (vecLength < INTERVAL_BOUNDS_MIN_ELEMENTS_FOR_COMPUTING)
    {
        return vecLength;
    }

    if (last == numeric_limits<size_t>::max())
    {
        last = vecLength;
    }

    Interval<I> theValue(value, value + unity);
    result = upper_boundIndex<0>(theIntervals, theValue, first, last, 
                                 intervalsTableUpperBoundsCompare);

    // Search for lower using maxSize.
    I largest = theIntervals[result].upper() + maxSize;
    size_t i = result;
    while ((i < last) && (theIntervals[i].lower() <= largest))
    {
        if (in(value, theIntervals[i]))
        {
            result = i + 1;
        }
        ++i;
    }


    return result;
}


#endif


