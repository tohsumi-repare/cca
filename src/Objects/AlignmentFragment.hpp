#ifndef MOLBIOLIB_ALIGNMENTFRAGMENT_H
#define MOLBIOLIB_ALIGNMENTFRAGMENT_H 1

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


// Below loads Location.hpp and thus Table.hpp
#include "src/Objects/AlignmentElement.hpp"




/** \file AlignmentFragment.hpp
 * Contains the #AlignmentFragment class
 */


/** An alignment of one contig to another.
 * The #AlignmentFragment consists of an alignment of one contig (typically
 * a single query) to another contig (typically a contig on a reference).
 * \section Usage Usage:
 * Can get the alignment type:
 * \code
 * AlignmentFragment<MyLocationType> anAlignFrag;
 *    // If location type is SequenceLocation, can use SequenceAlignmentFragment
 * size_t n, num;
 * ...
 * num = anAlignFrag.size();  // Number of alignment elements.
 * anAlignFrag.getAlignmentType();   // Get a consensus alignType.  May be Combination.
 * anAlignFrag.getAlignmentType(n);   // Get a specific alignType.
 * \endcode
 * Can then get data:
 * \code
 * MyLocationType theLoc;
 * theLoc = anAlignFrag.getQueryLocation();   // Get global location
 * theLoc = anAlignFrag.getTargetLocation();
 * theLoc = anAlignFrag.getQueryLocation(n);   // Get nth location
 * theLoc = anAlignFrag.getTargetLocation(n);
 * // Note that the location may contain information, e.g.
 * // the Broad alignment has number of errors in the info of type
 * // StringQualifier with key ERRS.
 * \endcode
 * \todo Someday when g++ correctly handles templated functions of templated
 * base classes, replace get<0>(anAlignFrag.theData) by
 * anAlignFrag.getColRef<0>() and similarly everywhere else in this file.
 */
template<typename TheLocationType>
class AlignmentFragment
{
public:

    // While it is possible to derive AlignmentFragment from
    // Table<ALIGNMENT_ELEMENT_ROW_TYPE>, it does not yield much since
    // C++ is not specified to be able to resolve types and methods that are
    // template-based of a templated base class where template arguments are
    // passed from the derived class' template arguments (e.g.
    // template<typename T> class PolyE : public Table<ALIGN..._TYPE> { ... }).
    // We would have to had put Table<ALIGN...>:: in front of almost everything.
    Table<ALIGNMENT_ELEMENT_ROW_TYPE> alignFragTable;

    typedef TheLocationType location_type;
    // Below must be done since in C++, a type that is dependent on the base
    // class' template parameters is not looked up when instantiating the
    // dervied class with a template parameter that is passed to a template
    // parameter of the base class, as we have here.
    typedef typename Table<ALIGNMENT_ELEMENT_ROW_TYPE>::row_type row_type;
    typedef typename Table<ALIGNMENT_ELEMENT_ROW_TYPE>::cols_type cols_type;
    typedef AlignmentElement<TheLocationType> alignment_element_type;

    TheLocationType globalQueryLocation, globalTargetLocation;
    TheLocationType& getQueryLocation()
    {
        return globalQueryLocation;
    }
    TheLocationType& getTargetLocation()
    {
        return globalTargetLocation;
    }


    AlignmentElementOperation getAlignmentType()
    {
        size_t numElements = alignFragTable.numRows();
        AlignmentElementOperation result = NoChange;
        // g++ 4.3 does not correctly handle templated functions of
        // templated base classes, so we resort to the below and similarly for
        // the below functions.
        // vector<AlignmentElementOperation>& theAlignFragTypes = alignFragTable.getColRef<0>();
        vector<AlignmentElementOperation>& theAlignFragTypes = get<0>(alignFragTable.theData);
        for (size_t i = 0; i < numElements && result != Combination; ++i)
        {
            if (i == 0)
            {
                result = theAlignFragTypes[0];
            }
            else if (result != theAlignFragTypes[i])
            {
                result = Combination;
            }
        }
        return result;
    }

    void push_back(row_type& theRow)
    {
        alignFragTable.push_back(theRow);
    }

    void clear()
    {
        alignFragTable.clear();
    }

    size_t size()
    {
        return alignFragTable.size();
    }

    AlignmentElementOperation& getAlignmentType(size_t n)
    {
        return get<0>(alignFragTable.theData)[n];
    }

    TheLocationType& getQueryLocation(size_t n)
    {
        return get<1>(alignFragTable.theData)[n];
    }

    TheLocationType& getTargetLocation(size_t n)
    {
        return get<2>(alignFragTable.theData)[n];
    }


};


/** \file AlignmentFragment.hpp
 * Standard alignment container in MolBioLib.  clear() overridden since
 * it is known in this case that the info() type has a clear method.
 */
class SequenceAlignmentFragment : public AlignmentFragment<SequenceLocation>
{
public:
    void clear()
    {
        alignFragTable.clear();
        globalQueryLocation.info().clear();
        globalTargetLocation.info().clear();
    }
};


#endif

