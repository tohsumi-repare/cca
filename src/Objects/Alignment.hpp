#ifndef MOLBIOLIB_ALIGNMENT_H
#define MOLBIOLIB_ALIGNMENT_H 1

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


// Below loads Location.hpp
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Objects/Qualifiers.hpp"



/** \file Alignment.hpp
 * Contains only the #Alignment class.
 */

/** A table of #AlignmentFragment instances and a map of parameters to strings.
 * The parameters variable (map) can include the file names of relevant
 * files.
 */
template<typename LocationType>
class Alignment : public Table< AlignmentFragment<LocationType> >
{
public:
    StringQualifiers parameters;
};


/** \file Alignment.hpp
 * \var typedef Alignment<SequenceLocation> SequenceAlignment
 * Typedef'd since this is used as the standard alignment container
 * in MolBioLib.
 */
typedef Alignment<SequenceLocation> SequenceAlignment;


#endif

