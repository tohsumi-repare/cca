#ifndef MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTWRITER_H
#define MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTWRITER_H 1

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


#include "src/Objects/FileStream.hpp"
#include "src/Objects/Location.hpp"
#include "src/Objects/Table.hpp"
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"


/** \file SequenceAlignmentFragmentWriter.hpp
 * Defines base class #SequenceAlignmentFragmentWriter.
 * Documentation on how to use #BroadQltoutAlignmentFragmentWriter here.
 * Should always be using something derived from
 * #SequenceAlignmentFragmentWriter, not the class itself.
 */

/** Base class to write alignment files - not meant to be used - see docs.
 * An aligment fragment is a single entry in an alignment file.  It is true
 * that sometimes, an alignment consists of just one fragment.
 * This class is not meant to be used.  Rather a class derived from this
 * should be used where writeAlignmenFragment
 * is defined/overridden.  The expected usage of a class derived from
 * #SequenceAlignmentFragmentWriter is
 * \code
 * SequenceAlignmentFragment theFragment;  // fill in later...
 * // ...
 * BroadQltoutAlignmentFragmentWriter myWriter("output.qltout");
 *    // BroadQltoutAlignmentFragmentWriter is derived from
 *      // SequenceAlignmentFragmentWriter.
 * myWriter.write(theFragment);
 * \endcode
 */
class SequenceAlignmentFragmentWriter
{
public:
    typedef SequenceLocation location_type;
    typedef location_type::position_type position_type;
    typedef location_type::element_type element_type;
    typedef location_type::info_type info_type;

    // SequenceAlignmentFragmentWriter() = default;

    virtual ~SequenceAlignmentFragmentWriter() = default;


    virtual void write(SequenceAlignmentFragment& theFragment)
    {
        assert(true == false);  // Must override this function in derived class.
    }

    virtual void close()
    {
        // Default, do nothing.
    }

};

#endif

