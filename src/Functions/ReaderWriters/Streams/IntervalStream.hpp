#ifndef MOLBIOLIB_INTERVALSTREAM_H
#define MOLBIOLIB_INTERVALSTREAM_H 1

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



/** \file IntervalStream.hpp
 * This adds streams to #Interval.
 * \fn template<typename I> istream& operator>>(istream& in, Interval<I>& theInterval)
 * This adds input streams to #Interval.
 * \section Usage Usage:
 * \code
 * Interval i1;
 * \endcode
 * If the input stream is of the form [3,5], [3], or 5, then one can do
 * \code
 * cin >> i1;
 * \endcode
 * \section ProgNotes Programming Notes:
 * It is similar version of what's in Boost's example.
 */
template<typename I>
istream& operator>>(istream& in, Interval<I>& theInterval)
{

    I lower, upper;
    char c = 0;
    bool isEmpty = false;
    in >> c;
    if (c == '[')
    {
        c = in.get();
        if (c  == ']')
        {
            // Empty interval
            theInterval.setEmpty();
            isEmpty = true;
        }
        else
        {
            in.putback(c);
            in >> lower >> c;
            if (c == ',')
            {
                in >> upper >> c;
            }
            else
            {
                upper = lower;
            }
        }
    }
    else
    {
        in.putback(c);
        in >> lower;
        upper = lower;
    }

#ifdef DEBUG
    assert(in.bad() == false);
#endif
    // Only assign a value if the interval is not empty.
    if (!isEmpty)
    {
        theInterval.assign(lower, upper);
    }

    return in;
}




/** \file IntervalStream.hpp
 * This adds streams to #Interval.
 * \fn template<typename I> ostream& operator<<(ostream& out, Interval<I>& theInterval)
 * This adds output streams to #Interval.
 * \section Usage Usage:
 * \code
 * Interval i1(3, 5);
 * cout << i1 << endl;
 * \endcode
 * \section ProgNotes Programming Notes:
 * It is similar version of what's in Boost's example.
 */
template<typename I>
ostream& operator<<(ostream& out, Interval<I>& theInterval)
{
    if (theInterval.empty())
    {
        return out << "[]";
    }
    else if (theInterval.lower() == theInterval.upper())
    {
        return out << '[' << theInterval.lower() << ']';
    }
    else
    {
        return out << '[' << theInterval.lower() << "," << theInterval.upper() << ']';
    }

}

#endif

