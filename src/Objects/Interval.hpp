#ifndef MOLBIOLIB_INTERVAL_H
#define MOLBIOLIB_INTERVAL_H 1

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


#include "src/Objects/Table.hpp"
#include "src/Functions/Algorithms/Tables/SortTable.hpp"

/** \file Interval.hpp
 * A closed interval class.
 */

/** A closed interval class.
 * The class #Interval is just a closed interval
 * where the endpoints are of the same type as the index type of a vector.
 * A lot of the interface is similar or the same to Boost's, but more
 * specialized.
 * \section Usage Usage:
 * Suppose one has the declaration and assignment:
 *
 * \code
 * Interval<long> i1(3, 5), i2, i3;  // Always assign low to high.
 *     // Since Interval<long> happens so often, it is typedef'd to
 *     // SequenceInterval, so that can be used in the below examples instead.
 * \endcode
 *
 * Some things one can do with the Interval class are:
 * <ol>
 *    <li> Get the type of the #Interval:
 *         \code
           Interval<long>::element_type i0;   // Obviously of type long.
 *         \endcode
 *         This is useful when an #Interval is passed within another templated
 *         class, so the element type of the interval may be hidden.
 *    </li>
 *    <li> Another form of assignment:
 *         \code
 *         i2.assign(4, 8);
 *         \endcode
 *    </li>
 *    <li> Assign just the lower or upper values (no error checking done!):
 *         \code
 *         i2.assignLower(7);
 *         i2.assignUpper(15);
 *         \endcode
 *    </li>
 *    <li> Get the lower and upper value so one can do:
 *         \code
 *         cout << '[' << i2.lower() << ',' << i2.upper() << ']' << endl;
 *         \endcode
 *    </li>
 *    <li> Get the <code>width = upper - lower</code>.
 *         (This is *not* the size, which is +1 of this.)
 *         \code
 *         cout<< i1.width();
 *         \endcode
 *    </li>
 *    <li> Check if empty (returns a bool):
 *         \code
 *         if (i1.empty()) { ... }
 *         \endcode
 *         Can similarly call <code>i1.setEmpty(true);</code> to set the
 *         empty flag true (default) or false.
 *    </li>
 *    <li> Can add, subtract, multiply, or divide by a constant of
 *         the same type to both the upper and lower values:
 *         \code
 *         i1 += 2;
 *         i2 -= 5;
 *         i1 *= i2.lower();
 *         i2 /= i3.upper();
 *         \endcode
 *         This is usually done when doing some sort of translation or scaling.
 *    </li>
 *    <li> Widen
 *         \code
 *         i1.widen(2);
 *         \endcode
 *         This subtracts 2 from the lower and adds 2 to the upper value.
 *         Use negative numbers to shrink the interval.
 *    <li> Intersection
 *         \code
 *         i3 = intersection(i2, i3, false);
 *         \endcode
 *         In all of the below's boolean checks, such as subset, overlap,
 *         <i>etc.</i>, there is an optional boolean argument at the end
 *         (default is false) that if true will check the arguments for that
 *         operator's function anyway even if the interval(s) are empty.  This
 *         is contrary to standard set definitions, but still may be useful.
 *         <b>Note that overloaded operators such as <code>==</code> cannot
 *         and does not have this optional argument.</b>
 *    </li>
 *    <li> Subset - i1 subset of i2?
 *         \code
 *         if (subset(i1, i2)) { ... }
 *         \endcode
 *    </li>
 *    <li> Superset - i1 superset of i2?
 *         \code
 *         if (superset(i1, i2)) { ... }
 *         \endcode
 *    </li>
 *    <li> Overlap
 *         \code
 *         if (overlap(i1, i2)) { ... }
 *         \endcode
 *    </li>
 *    <li> Hull - like union if overlap.  The hull is the smallest interval
 *         that contains the two intervals.  So, here there the intervals need
 *         not intersect.
 *         \code
 *         i3 = hull(i1, i2);
 *         \endcode
 *         There is an option to ignore the fact that i1 or i2 may be empty.
 *    </li>
 *    <li> Contains, <i>e.g.</i> if 4 in the interval i1?
 *         \code
 *         if (in(4, i1)) { ... }
 *         \endcode
 *    </li>
 *    <li> Check for inequalities and equalities, where for
 *         equality both upper and lower values must be the same and for
 *         inequality, first the lower is checked and if equal, the
 *         upper one is checked.
 *    </li>
 *    <li> Given a vector of intervals, is a specified number (in the below
 *         example, 4) in one of them?  If so return
 *         return the first such interval, else return numeric_limits<size_t>::max()
 *         \code
 *         vector< Interval<long> > theInts;
 *         ...
 *         size_t n = in(4, theInts);
 *         \endcode
 *         <b> Typically, one may use <code>unsigned long</code> in place of
 *         <code>size_t</code> above and below.</b>
 *         This method has the option to not check for emptyness of intervals.
 *    </li>
 *    <li> Similarly, if one has an interval, then one can search for an
 *         overlap:
 *         \code
 *         size_t n = overlap(i1, theInts, false);
 *         \endcode
 *         This method has the option to not check for emptyness of intervals.
 *         <ul>
 *            <li> The above two are fairly inefficienct, as it does a linear
 *                 search.  The proposed TODO is more efficient for a number
 *                 of points and intervals:
 *            </li>
 *         </ul>
 *    </li>
 *    <li> If one has a vector of overlapping intervals, one may merge them as
 *         follows:
 *         \code
 *         vector<long> theOutput;
 *         merge(theInts, theOutput, false);
 *         // or, theOutput = merge(theInts, false);
 *         \endcode
 *         where the <code>false</code> (default) is to check for
 *         overlap (and thus, subsets).  <code>true</code> would only merge
 *         intervals which is a subset of some other interval.  The last
 *         parameter is optional.
 *    </li>
 *    <li> If one has a vector of intervals, one may do the
 *         following operations:
 *         \code
 *         vector< Interval<long> > theIntervals;
 *         // Fill theIntervals.
 *         long maxSize = maxWidth(theIntervals) + 1;
 *         \endcode
 *    </li>
 *    <li> All of the above <code>vector</code> works with <code>Table</code>
 *         too.  Replace vector[Name] by table[Name] and replace vector<type>
 *         by Table<type>.
 *    </li>
 * </ol>
 * \section ProgNotes Programming Notes:
 * <ul>
 *    <li> The lower and upper nomenclature comes from adopting Boost's
 *         interval library's nomenclature.
 *    </li>
 *    <li> Normallyfunctions methods outside of objects are put in the Functions
 *         subdirectory, but since all of the below methods pertain to
 *         #Interval, we put all the methods in this file.
 *    </li>
 *    <li> <em>The following is removed from the code.</em>  It seems
 *         confusing to put it in and just confuses issues.  Besides, it seems
 *         not to be needed generally.
 *
 *         If one has a vector of intervals and a larger interval, one may
 *         get the complement of the intervals relative to the larger
 *         interval as follows:
 *         \code
 *         Interval<long> intervalUniverse(2, 100);
 *         vector< Interval<long> > theOutput;
 *         complement(intervalUniverse, theInts, theOutput);
 *         // or, theOutput = complement(intervalUniverse, theInts);
 *         \endcode
 *         Note that this is the mathematical complement, so for example if the
 *         universe is [1, 10] and we have intervals [2, 4] and [6, 9], then
 *         the complement is the set of intervals [1, 2], [5, 6], and [9, 10].
 *    </li>
 * </ul>
 * \todo Write a function to find the position for each point or
 *       #Interval in a vector of such.
 */
template<typename I>
class Interval
{
public:
    typedef I element_type;

    const I lower() const
    {
        return _lower;
    }
    const I upper() const
    {
        return _upper;
    }
    const bool empty() const
    {
        return _empty;
    }
    void setEmpty(bool newValue = true)
    {
        _empty = newValue;
    }

    const bool singleton() const
    {
        return (!_empty && (_lower == _upper));
    }

    void assign(I low)
    {
        _lower = low;
        _upper = low;
        _empty = false;
    }
    void assign(I low, I up)
    {
        _lower = low;
        _upper = up;
        if (_upper < _lower)
        {
            _empty = true;
#ifdef PROG_DEBUG
            cerr << "Error in Interval.assign(I low, I up)!  up < low.  empty is set true.  Continuing execution." << endl;
#elif defined(DEBUG)
            cerr << "Warning in Interval.assign(I low, I up).  up < low.  Normally, empty would be set true.  Instead, swapping and Continuing execution." << endl;
            swap(_lower, _upper);
            _empty = false;
#endif
        }
        else
        {
            _empty = false;
        }
    }

    void assignLower(I low)
    {
        _lower = low;
    }

    void assignUpper(I up)
    {
        _upper = up;
    }

    I width() const
    {
        return (_upper - _lower);
    }

    Interval()
    {
        _empty = true;
    }
    Interval(I low)
    {
        assign(low);
    }
    Interval(I low, I up)
    {
        assign(low, up);
    }

    Interval<I>& operator+=(I rhs)
    {
        _lower += rhs;
        _upper += rhs;
        return *this;
    }

    Interval<I>& operator-=(I rhs)
    {
        _lower -= rhs;
        _upper -= rhs;
        return *this;
    }

    Interval<I>& operator*=(I rhs)
    {
        _lower *= rhs;
        _upper *= rhs;
        return *this;
    }

    Interval<I>& operator/=(I rhs)
    {
        _lower /= rhs;
        _upper /= rhs;
        return *this;
    }

    void widen(I amt)
    {
        _lower -= amt;
        _upper += amt;
        if (_upper < _lower)
        {
            _empty = true;
#ifdef DEBUG
            cerr << "Warning in Interval.widen(I up).  up < low.  Normally, empty would be set true.  Instead, swapping and Continuing execution." << endl;
            swap(_lower, _upper);
            _empty = false;
#endif
        }
    }

private:
    I _lower, _upper;
    bool _empty;
};


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool subset(const Interval<I>& i1, const Interval<I>& i2, bool checkEmpty = false)
 * Subset.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool subset(const Interval<I>& i1, const Interval<I>& i2,
            bool checkEmpty = false)
{
    if (!checkEmpty && i1.empty())
    {
        return true;
    }
    return ((!checkEmpty || !i2.empty()) &&
            (i2.lower() <= i1.lower()) &&
            (i1.upper() <= i2.upper()));
}


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool superset(const Interval<I>& i1, const Interval<I>& i2, bool checkEmpty = false)
 * Subset.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool superset(const Interval<I>& i1, const Interval<I>& i2,
              bool checkEmpty = false)
{
    return subset(i2, i1, checkEmpty);
}


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool overlap(const Interval<I>& i1, const Interval<I>& i2, bool checkEmpty = false)
 * Overlap.  True on overlap, or i1 or i2 is empty (assuming checkEmpty =
 * false), or one is superset of other.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool overlap(const Interval<I>& i1, const Interval<I>& i2,
             bool checkEmpty = false)
{
    if (!checkEmpty && (i1.empty() || i2.empty()))
    {
        return true;
    }
    return ((i1.lower() <= i2.upper()) && (i2.lower() <= i1.upper()));
}


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> Interval<I> intersection(const Interval<I>& i1, const Interval<I>& i2, bool checkEmpty = false)
 * Intersection.
 * See #Interval class documentation for usage.
 */
template<typename I>
Interval<I> intersection(const Interval<I>& i1, const Interval<I>& i2,
                         bool checkEmpty = false)
{
    if ((!checkEmpty && (i1.empty() || i2.empty() || !overlap(i1, i2))) ||
            (checkEmpty && !overlap(i1, i2, true)))
    {
        return Interval<I>();
    }

    const I& low = max(i1.lower(), i2.lower());
    const I& up = min(i1.upper(), i2.upper());

    // if up < low, then the empty flag will be set by the assign method.
    return Interval<I>(low, up);
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> Interval<I> hull(Interval<I>& i1, Interval<I>& i2, bool checkEmpty = false)
 * Hull.
 * See #Interval class documentation for usage.
 */
template<typename I>
Interval<I> hull(Interval<I>& i1, Interval<I>& i2, bool checkEmpty = false)
{
    Interval<I> out;
    if (i1.empty() && i2.empty())
    {
        return out;
    }
    else if (!checkEmpty && i1.empty())
    {
        out = i2;
    }
    else if (!checkEmpty && i2.empty())
    {
        out = i1;
    }
    else
    {
        const I& low = min(i1.lower(), i2.lower());
        const I& up = max(i1.upper(), i2.upper());
        out.assign(low, up);
    }
    return out;
}


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator==(const Interval<I>& i1, const Interval<I>& i2)
 * ==.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator==(const Interval<I>& i1, const Interval<I>& i2)
{
    if (i1.empty() && i2.empty())
    {
        return true;
    }
    else if (i1.empty() || i2.empty())
    {
        return false;
    }
    if ((i1.lower() == i2.lower()) && (i1.upper() == i2.upper()))
    {
        return true;
    }
    return false;
}
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator!=(const Interval<I>& i1, const Interval<I>& i2)
 * !=.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator!=(const Interval<I>& i1, const Interval<I>& i2)
{
    if (i1.empty() && i2.empty())
    {
        return false;
    }
    else if (i1.empty() || i2.empty())
    {
        return true;
    }
    if (!((i1.lower() == i2.lower()) && (i1.upper() == i2.upper())))
    {
        return true;
    }
    return false;
}
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator<(const Interval<I>& i1, const Interval<I>& i2)
 * <.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator<(const Interval<I>& i1, const Interval<I>& i2)
{
#ifdef PROG_DEBUG
    if (i1.empty() || i2.empty())
    {
        cerr << "Error!  < for an empty interval is undefined.  Exiting."
             << endl;
        exit(1);
    }
#endif
    if (i1.lower() < i2.lower())
    {
        return true;
    }
    else if ((i1.lower() == i2.lower()) && (i1.upper() < i2.upper()))
    {
        return true;
    }
    else
    {
        return false;
    }
}
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator<=(const Interval<I>& i1, const Interval<I>& i2)
 * <=.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator<=(const Interval<I>& i1, const Interval<I>& i2)
{
#ifdef PROG_DEBUG
    if (i1.empty() || i2.empty())
    {
        cerr << "Error!  <= for an empty interval is undefined.  Exiting."
             << endl;
        exit(1);
    }
#endif
    if (i1.lower() <= i2.lower())
    {
        return true;
    }
    else if ((i1.lower() == i2.lower()) && (i1.upper() <= i2.upper()))
    {
        return true;
    }
    else
    {
        return false;
    }
}
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator>(const Interval<I>& i1, const Interval<I>& i2)
 * >.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator>(const Interval<I>& i1, const Interval<I>& i2)
{
#ifdef PROG_DEBUG
    if (i1.empty() || i2.empty())
    {
        cerr << "Error!  > for an empty interval is undefined.  Exiting."
             << endl;
        exit(1);
    }
#endif
    if (i1.upper() > i2.upper())
    {
        return true;
    }
    else if ((i1.upper() == i2.upper()) && (i1.lower() > i2.lower()))
    {
        return true;
    }
    else
    {
        return false;
    }
}
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool operator>=(const Interval<I>& i1, const Interval<I>& i2)
 * >=.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool operator>=(const Interval<I>& i1, const Interval<I>& i2)
{
#ifdef PROG_DEBUG
    if (i1.empty() || i2.empty())
    {
        cerr << "Error!  >= for an empty interval is undefined.  Exiting."
             << endl;
        exit(1);
    }
#endif
    if (i1.upper() >= i2.upper())
    {
        return true;
    }
    else if ((i1.upper() == i2.upper()) && (i1.lower() >= i2.lower()))
    {
        return true;
    }
    else
    {
        return false;
    }
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> bool in(I x, Interval<I>& i1, bool checkEmpty = false)
 * In.
 * See #Interval class documentation for usage.
 */
template<typename I>
bool in(I x, Interval<I>& i1, bool checkEmpty = false)
{
    if (!checkEmpty && i1.empty())
    {
        return false;
    }
    if ((x >= i1.lower()) && (x <= i1.upper()))
    {
        return true;
    }
    else
    {
        return false;
    }
}




/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> template<typename I> size_t in(I x, vector< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for a scalar.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t in(I x, vector< Interval<I> >& theIntervals, bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (in(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> size_t in(I x, Table< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for a scalar.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t in(I x, Table< Interval<I> >& theIntervals, bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (in(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> size_t overlap(Interval<I>& x, vector< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t overlap(Interval<I>& x, vector< Interval<I> >& theIntervals,
               bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (overlap(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> size_t overlap(Interval<I>& x, Table< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t overlap(Interval<I>& x, Table< Interval<I> >& theIntervals,
               bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (overlap(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}




/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> size_t subset(Interval<I>& x, vector< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t subset(Interval<I>& x, vector< Interval<I> >& theIntervals,
              bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (subset(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> size_t subset(Interval<I>& x, Table< Interval<I> >& theIntervals, bool checkEmpty = false)
 * In for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
size_t subset(Interval<I>& x, Table< Interval<I> >& theIntervals,
              bool checkEmpty = false)
{
    size_t numIntervals = theIntervals.size();
    for (size_t i = 0; i < numIntervals; ++i)
    {
        if (subset(x, theIntervals[i], checkEmpty))
        {
            return i;
        }
    }
    return numeric_limits<size_t>::max();
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> void merge(vector< Interval<I> >& in, vector< Interval<I> >& out, bool isSubset = false, bool checkEmpty = false)
 * vector argument merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
void merge(vector< Interval<I> >& in, vector< Interval<I> >& out,
           bool isSubset = false, bool checkEmpty = false)
{
    if (in.size() <= 1)
    {
        out = in;
        return;
    }

    vector< Interval<I> > temp;
    for (size_t i = 0; i < in.size(); ++i)
    {
        if (checkEmpty || !(in[i].empty()))
        {
            temp.push_back(in[i]);
        }
    }
    sort(temp.begin(), temp.end());
    bool didMerge = true;
    while (didMerge)
    {
        vector<bool> merged(temp.size(), false);
        out.clear();
        didMerge = false;
        for (size_t i = 0; i < (temp.size() - 1); ++i)
        {
            if (merged[i])
            {
                continue;
            }
            bool found = false;
            size_t goodJ = 0;
            for (size_t j = i+1; !found && (j < temp.size()) &&
                    (temp[j].lower() <= temp[i].upper()); ++j)
            {
                if (merged[j])
                {
                    continue;
                }
                if ((!isSubset && overlap(temp[i], temp[j], checkEmpty)) ||
                        (isSubset && subset(temp[i], temp[j], checkEmpty)) ||
                        (isSubset && subset(temp[j], temp[i], checkEmpty)))
                {
                    found = true;
                    goodJ = j;
                }
            }
            if (found)
            {
                didMerge = true;
                Interval<I> newI = hull(temp[i], temp[goodJ], checkEmpty);
                out.push_back(newI);
                merged[i] = true;
                merged[goodJ] = true;
            }
            else
            {
                out.push_back(temp[i]);
            }
        }
        if (!merged[temp.size() - 1])
        {
            out.push_back(temp[temp.size() - 1]);
        }
        if (didMerge)
        {
            temp = out;
        }
    }
}

/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> void merge(Table< Interval<I> >& in, Table< Interval<I> >& out, bool isSubset = false, bool checkEmpty = false)
 * table arguments merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
void merge(Table< Interval<I> >& in, Table< Interval<I> >& out,
           bool isSubset = false, bool checkEmpty = false)
{
    if (in.size() <= 1)
    {
        out = in;
        return;
    }

    Table< Interval<I> > temp;
    for (size_t i = 0; i < in.size(); ++i)
    {
        if (checkEmpty || !(in[i].empty()))
        {
            temp.push_back(in[i]);
        }
    }
    sortTable<0>(temp);
    bool didMerge = true;
    while (didMerge)
    {
        vector<bool> merged(temp.size(), false);
        out.clear();
        didMerge = false;
        for (size_t i = 0; i < (temp.size() - 1); ++i)
        {
            if (merged[i])
            {
                continue;
            }
            bool found = false;
            size_t goodJ = 0;
            for (size_t j = i+1; !found && (j < temp.size()) &&
                    (temp[j].lower() <= temp[i].upper()); ++j)
            {
                if (merged[j])
                {
                    continue;
                }
                if ((!isSubset && overlap(temp[i], temp[j])) ||
                        (isSubset && subset(temp[i], temp[j])) ||
                        (isSubset && subset(temp[j], temp[i])))
                {
                    found = true;
                    goodJ = j;
                }
            }
            if (found)
            {
                didMerge = true;
                Interval<I> newI = hull(temp[i], temp[goodJ]);
                out.push_back(newI);
                merged[i] = true;
                merged[goodJ] = true;
            }
            else
            {
                out.push_back(temp[i]);
            }
        }
        if (!merged[temp.size() - 1])
        {
            out.push_back(temp[temp.size() - 1]);
        }
        if (didMerge)
        {
            temp = out;
        }
    }
}




/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> vector< Interval<I> > merge(vector< Interval<I> >& in, bool isSubset = false, bool checkEmpty = false)
 * vector merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
vector< Interval<I> > merge(vector< Interval<I> >& in,
                            bool isSubset = false, bool checkEmpty = false)
{
    vector< Interval<I> > theResult;
    merge(in, theResult, isSubset, checkEmpty);
    return theResult;
}

/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> Table< Interval<I> > merge(Table< Interval<I> >& in, bool isSubset = false, bool checkEmpty = false)
 * table merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
Table< Interval<I> > merge(Table< Interval<I> >& in,
                           bool isSubset = false, bool checkEmpty = false)
{
    Table< Interval<I> > theResult;
    merge(in, theResult, isSubset, checkEmpty);
    return theResult;
}



/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> I maxWidth(vector< Interval<I> >& in)
 * table merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
I maxWidth(vector< Interval<I> >& in)
{
    if (in.size() == 0)
    {
        // Cannot compute this on an empty set of intervals.
        assert(true == false);
    }
    I result = in[0].width();
    size_t numIn = in.size();
    for (size_t i = 1; i < numIn; ++i)
    {
        if (in[i].width() > result)
        {
            result = in[i].width();
        }
    }
    return result;
}


/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> I maxWidth(Table< Interval<I> >& in)
 * table merge for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
I maxWidth(Table< Interval<I> >& in)
{
    if (in.size() == 0)
    {
        // Cannot compute this on an empty set of intervals.
        assert(true == false);
    }
    I result = in[0].width();
    size_t numIn = in.size();
    for (size_t i = 1; i < numIn; ++i)
    {
        if (in[i].width() > result)
        {
            result = in[i].width();
        }
    }
    return result;
}



// The below is not included.  I think it is not the clearest function and is
// typically not needed.
#if 0
/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> void complement(Interval<I> universe, vector< Interval<I> >& in, vector< Interval<I> >& out)
 * complement for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
void complement(Interval<I> universe,
                vector< Interval<I> >& in, vector< Interval<I> >& out)
{
    out.clear();
    if (universe.width() < static_cast<I>(1))
    {
        return;    // Null universe.
    }
    if (in.size() == 0)
    {
        out.push_back(universe);
        return;  // Null set's complement.
    }

    // Sort and merge before taking complement.
    vector< Interval<I> > temp = in, tempout;
    sort(temp.begin(), temp.end());
    merge(temp, tempout);

    Interval<I> newI;
    if (tempout[0].lower() > universe.lower())
    {
        newI.assign(universe.lower(), tempout[0].lower());
        out.push_back(newI);
    }
    size_t numTemp = tempout.size();
    for (size_t i = 0; i < (numTemp - 1); ++i)
    {
#ifdef DEBUG
        assert(subset(tempout[i], universe));
#endif
        newI.assign(tempout[i].upper(), tempout[i+1].lower());
        out.push_back(newI);
    }
    if (tempout[numTemp - 1].upper() < universe.upper())
    {
        newI.assign(tempout[numTemp - 1].upper(), universe.upper());
        out.push_back(newI);
    }
}

/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> void complement(Interval<I> universe, Table< Interval<I> >& in, Table< Interval<I> >& out)
 * complement for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
void complement(Interval<I> universe,
                Table< Interval<I> >& in, Table< Interval<I> >& out)
{
    out.clear();
    if (universe.width() < static_cast<I>(1))
    {
        return;    // Null universe.
    }
    if (in.size() == 0)
    {
        out.push_back(universe);
        return;  // Null set's complement.
    }

    // Sort and merge before taking complement.
    Table< Interval<I> > temp = in, tempout;
    sortTable<0>(temp);
    merge(temp, tempout);

    Interval<I> newI;
    if (tempout[0].lower() > universe.lower())
    {
        newI.assign(universe.lower(), tempout[0].lower());
        out.push_back(newI);
    }
    size_t numTemp = tempout.size();
    for (size_t i = 0; i < (numTemp - 1); ++i)
    {
#ifdef DEBUG
        assert(subset(tempout[i], universe));
#endif
        newI.assign(tempout[i].upper(), tempout[i+1].lower());
        out.push_back(newI);
    }
    if (tempout[numTemp - 1].upper() < universe.upper())
    {
        newI.assign(tempout[numTemp - 1].upper(), universe.upper());
        out.push_back(newI);
    }
}




/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> vector< Interval<I> > complement(Interval<I> universe, vector< Interval<I> >& in)
 * complement for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
vector< Interval<I> > complement(Interval<I> universe,
                                 vector< Interval<I> >& in)
{
    vector< Interval<I> > theResult;
    complement(universe, in, theResult);
    return theResult;
}

/** \file Interval.hpp
 * #Interval functions.
 * \fn template<typename I> Table< Interval<I> > complement(Interval<I> universe, Table< Interval<I> >& in)
 * complement for an #Interval.
 * See #Interval class documentation for usage.
 */
template<typename I>
Table< Interval<I> > complement(Interval<I> universe, Table< Interval<I> >& in)
{
    Table< Interval<I> > theResult;
    complement(universe, in, theResult);
    return theResult;
}

// Below is endif for if 0 - removing the vector and table IntervalComplement
#endif




/** \file Interval.hpp
 * \var typedef Interval<long> SequenceInterval
 * Typedef'd since this is used as the standard interval container in MolBioLib.
 */
typedef Interval<long> SequenceInterval;



#endif

