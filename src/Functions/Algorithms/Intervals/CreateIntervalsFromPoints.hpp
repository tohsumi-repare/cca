#ifndef MOLBIOLIB_CREATEINTERVALSFROMPOINTS_H
#define MOLBIOLIB_CREATEINTERVALSFROMPOINTS_H 1

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
#include "src/Functions/Algorithms/Tables/SortTable.hpp"




/** \file CreateIntervalsFromPoints.hpp
 * Given a Table of points, create intervals which represent those points.
 * \fn void void createIntervalsFromPoints(Table<I> points, Table< Interval<I> >& outIntervals, bool keepEmpty = true, I unity = static_cast<I>(1))
 * Combines points into intervals, then merges the intervals.
 * keepEmpty is if we want to keep an interval that spans one point, i.e [x, x].
 * unity is typename I's definition of 1.  This is needed since we want to 
 * know if a point is next to another point by adding 1 to it.  Thus, the +
 * operation must be defined for I, as well as operator< for sorting.
 * This method is only for countable classes I's.
 *
 * Programming note: passing by value points, since we are going to sort them.
 */
template<typename I>
void createIntervalsFromPoints(Table<I> points, 
                               Table< Interval<I> >& outIntervals,
                               bool keepEmpty = true,
                               I unity = static_cast<I>(1))
{
    sortTable<0>(points);

    Table< Interval<I> > intervals;
    outIntervals.clear();

    size_t n = points.size();

    if ((n == 0) || ((n == 1) && !keepEmpty))
    {
        return;
    }
    
    Interval<I> currI(points[0], points[0]);
    if (n == 1) 
    {
        intervals.push_back(currI);
        return;
    }
    
    for (size_t i = 1; i < n; ++i)
    {
        if ((points[i] != currI.upper()) &&
            (points[i] != (currI.upper() + unity)))
        {
            // Need to create new interval
            if (keepEmpty || (currI.lower() != currI.upper()))
            {
                intervals.push_back(currI);
            }
            currI.assign(points[i], points[i]);
            if (keepEmpty && (i == (n-1)))
            {
                intervals.push_back(currI);
            }
        }
        else
        {
            currI.assignUpper(points[i]);
            if (i == (n-1))
            {
                intervals.push_back(currI);
            }
        }
    }


    // No need to sort intervals - already in sorted order since points sorted.
    merge<I>(intervals, outIntervals, false, keepEmpty);
}



/** \file CreateIntervalsFromPoints.hpp
 * Given a vector of points, create intervals which represent those points.
 * \fn void createIntervalsFromPoints(vector<I> points, vector< Interval<I> >& outIntervals, bool keepEmpty = true, I unity = static_cast<I>(1))
 * Combines points into intervals, then merges the intervals.
 * keepEmpty is if we want to keep an interval that spans one point, i.e [x, x].
 * unity is typename I's definition of 1.  This is needed since we want to 
 * know if a point is next to another point by adding 1 to it.  Thus, the +
 * operation must be defined for I, as well as operator< for sorting.
 * This method is only for countable classes I's.
 *
 * Programming note: passing by value points, since we are going to sort them.
 */
template<typename I>
void createIntervalsFromPoints(vector<I> points, 
                               vector< Interval<I> >& outIntervals,
                               bool keepEmpty = true,
                               I unity = static_cast<I>(1))
{
    sort(points.begin(), points.end());

    vector< Interval<I> > intervals;
    outIntervals.clear();

    size_t n = points.size();

    if ((n == 0) || ((n == 1) && !keepEmpty))
    {
        return;
    }
    
    Interval<I> currI(points[0], points[0]);
    if (n == 1) 
    {
        intervals.push_back(currI);
        return;
    }
    
    for (size_t i = 1; i < n; ++i)
    {
        if ((points[i] != currI.upper()) &&
            (points[i] != (currI.upper() + unity)))
        {
            // Need to create new interval
            if (keepEmpty || (currI.lower() != currI.upper()))
            {
                intervals.push_back(currI);
            }
            currI.assign(points[i], points[i]);
            if (keepEmpty && (i == (n-1)))
            {
                intervals.push_back(currI);
            }
        }
        else
        {
            currI.assignUpper(points[i]);
            if (i == (n-1))
            {
                intervals.push_back(currI);
            }
        }
    }


    // No need to sort intervals - already in sorted order since points sorted.
    merge<I>(intervals, outIntervals, false, keepEmpty);
}


#endif


