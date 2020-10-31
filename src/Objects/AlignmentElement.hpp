#ifndef MOLBIOLIB_ALIGNMENTELEMENT_H
#define MOLBIOLIB_ALIGNMENTELEMENT_H 1

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


#include "src/Objects/Location.hpp"



/** \file AlignmentElement.hpp
 * Contains the #AlignmentElement class as well as associated enum and DEFINE.
 * \enum AlignmentElementOperation
 * An alignment is indicated by NoChange.  (Alignment is already taken by
 * the Alignment class.)  An Indel is either an insertion or deletion.  An
 * IndelSNPs is a combination of an Indel and a SNPs with SNP(s) and either
 * an insertion or deletion included. A Substitution event can also include a
 * change of codon.  The Combination type is never to be used as an actual
 * table but rather as a return code to signify an indivisible
 * alignment element.
 * \todo Replace <code>enum AlignmentElementOperation</code> with
 * <code>enum class AlignmentElementOperation</code> when strongly typed
 * enums work.
 */
enum AlignmentElementOperation
{
    NoChange, Deletion, HardClip, Insertion, Indel, IndelSNPs, Inversion, MatchOrMismatch, Padding, SNPs, SoftClip, Substitution, Translocation, Combination
};



/** \file AlignmentElement.hpp
 * Contains the #AlignmentElement class as well as associated enum and DEFINE.
 * \def ALIGNMENT_ELEMENT_ROW_TYPE
 * The row_type in #Table has the type of polymorphism, the source location,
 * the target location, and any inserted bases (if pertinent).
 */
#ifndef ALIGNMENT_ELEMENT_ROW_TYPE
#define ALIGNMENT_ELEMENT_ROW_TYPE AlignmentElementOperation, TheLocationType, TheLocationType
#endif


/** An #AlignmentElement holds a single alignment type for a pair of contigs.
 * Just is a tuple that will go into a row of the #AlignmentFragment class.
 * The #AlignmentElement class stores one instance of a single type
 * of alignment, <i>i.e.</i> SNPs, an insertion, a deletion, an inversion,
 * and/or a translocation.  Additional information about an alignment may be
 * stored by deriving a class from #Location that contains
 * more information than just the location (see #SequenceLocation).
 * \section Usage Usage:
 * Can get the alignment type:
 * \code
 * AlignmentElement<MyLocationType> anElement;
 *   // Note that if the location type is SequenceLocation, then one may do
 *   // SequenceAlignmentElement anElement;
 * ...
 * anElement.type();   // Get the type.
 * \endcode
 * Can then get and set the data (can put accessors on either side of =):
 * \code
 * MyLocationType theLoc;
 * ...
 * theLoc = anElement.queryLocation();
 * anElement.targetLocation() = theLoc;
 * \endcode
 *
 * Programmer's note: the sequence is not included in the element.  Rather, it
 * may be included if needed as a feature in the location type, such as is
 * possible with #SequenceLocation.
 */
template<typename TheLocationType>
class AlignmentElement : public tuple<ALIGNMENT_ELEMENT_ROW_TYPE>
{
    // While it is possible to derive PolymorphismEvent from
    // Table<POLYMORPHISM_EVENT_ROW_TYPE>, it does not yield much since
    // C++ is not specified to be able to resolve types and methods that are
    // template-based of a templated base class where template arguments are
    // passed from the derived class' template arguments (e.g.
    // template<typename T> class PolyE : public Table<POLYM..._TYPE> { ... }).
    // We would have to had put Table<POLY...>:: in front of almost everything.

public:
    typedef TheLocationType location_type;
    AlignmentElement() { }

    AlignmentElementOperation& type()
    {
        return get<0>(*this);
    }

    TheLocationType& queryLocation()
    {
        return get<1>(*this);
    }

    TheLocationType& targetLocation()
    {
        return get<2>(*this);
    }

};



/** \file AlignmentElement.hpp
 * \var typedef AlignmentElement<SequenceLocation> SequenceAlignmentElement
 * Typedef'd since this is used as the standard alignment container
 * in MolBioLib.
 */
typedef AlignmentElement<SequenceLocation> SequenceAlignmentElement;




#endif

