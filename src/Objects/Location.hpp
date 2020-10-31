#ifndef MOLBIOLIB_LOCATION_H
#define MOLBIOLIB_LOCATION_H 1

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
#include "src/Objects/Qualifiers.hpp"



/** \file Location.hpp
 * Contains only the #Location.
 */

/** Holds contig, position, and other info.
 * This class simply contains both some position element type
 * (usually an integral type) which is used in the #Interval class and
 * additional contig-specifying information.
 * The template paramters are contig type, position type, and optionally
 * information type and whether it is zero-based in position (default true).
 * A contig is usually specfied by how it is specified in some FASTA file,
 * either by name (inefficient) or by the index of the entry within that file.
 * The strand information (usually either '+'/'-' or '0'/'1') of type
 * <code>string</code> is also stored.  One need not use this variable is no
 * such information exists.
 * The #Location object has an information member which is by default
 * <code>bool</code>, since it uses little space.  However, one can
 * store more complex information (as may be needed in a class
 * derived from #SequenceAlignmentFragmentReader) by passing
 * in a class with many members.
 *
 * \section Usage Usage:
 * \code
 * Location<int, double, false> theLoc;
 * Interval<int> x = theLoc.position();
 * // ...
 * theLoc.position() = x;
 * string whichContig = theLoc.contig();
 * \endcode
 * One can set contig information:
 * \code
 * theLoc.contig() = whichContig;
 * \endcode
 * One can access strand information:
 * \code
 * theLoc.strand() = "0";   // Strand is of type string
 * cout << theLoc.strand() << endl;
 * \endcode
 * The strand is purely optional, as there are objects with no strand
 * information.  One can access additional information:
 * \code
 * theLoc.info() = 0.75;
 * cout << theLoc.info() << endl;
 * \endcode
 * Finally, one can determine the origin of the position.  For example, for
 * 0-based indices, the origin is 0, for 1-based indices, the origin is 1, etc.
 * \code
 * // The below is of the same type as the element type.
 * cout << theLoc.origin() << endl;
 * \endcode
 */
template <typename TheElementType = size_t, typename TheInfoType = bool, TheElementType Origin = static_cast<TheElementType>(0)>
class Location : public tuple<string, Interval<TheElementType>, string, TheInfoType>
{
public:
    typedef Interval<TheElementType> position_type;
    typedef TheElementType element_type;
    typedef TheInfoType info_type;
    // There is no typedef string contig_type since we want to enforce that
    // all contigs are to be of type string.

    Location()
    {
        theOrigin = Origin;
    }
    // In C++11, any class which has a constructor automatically deletes the
    // default copy constructor.  Thus, one needs to be supplied.
    Location& operator=(Location const &rhs)
    {
        if (this != &rhs)
        {
            get<0>(*this) = get<0>(rhs);
            get<1>(*this) = get<1>(rhs);
            get<2>(*this) = get<2>(rhs);
            get<3>(*this) = get<3>(rhs);
            theOrigin = Origin;
        }
        return *this;
    }

    element_type origin()
    {
        return theOrigin;
    }

    string& contig()
    {
        return get<0>(*this);
    }
    position_type& position()
    {
        return get<1>(*this);
    }
    string& strand()
    {
        return get<2>(*this);
    }
    info_type& info()
    {
        return get<3>(*this);
    }

private:
    TheElementType theOrigin;

};

/** \file Location.hpp
 * \var typedef Location<long, StringQualifiers> SequenceLocation
 * Typedef'd since this is used as the standard location container in MolBioLib.
 * Since version on Sept. 2, 2009, the below is 0-indexed, not 1-indexed.
 */
typedef Location<long, StringQualifiers> SequenceLocation;

#endif

