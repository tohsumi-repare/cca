#ifndef MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTTOFEATURETRACKERWRITER_H
#define MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTTOFEATURETRACKERWRITER_H 1

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
// Below loads Table.hpp
#include "src/Objects/Location.hpp"
#include "src/Objects/SequenceAlignmentFragmentToFeaturesTracker.hpp"



/** \file SequenceAlignmentFragmentToFeatureTrackerWriter.hpp
 * This file defines two functions to write output of the #SequenceAlignmentFragmentToFeaturesTracker class to a file.
 * \fn template<typename CountType> void writeSequenceAlignmentFragmentToFeatureTrackerCounts(SequenceAlignmentFragmentToFeaturesTracker<CountType>& theTracker, string filename, bool printAllLines = false, SequenceLocation::element_type offset = 1, bool separatePos = false, bool printCombined = false)
 * This prints out only genes that got hit (unless the below boolean parameter
 * is <code>true</code> - default is <code>false</code>):
 * \section Usage Usage:
 * \code
 * SequenceAlignmentFragmentToFeaturesTracker<CountType> theTracker;
 * writeSequenceAlignmentFragmentToFeatureTrackerCounts(theTracker, "theTracker.tsv", false, 1, false, false, false);
 * \endcode
 * where the last <code>1</code> is an applied index to the internal 0-based
 * index coordinates.  The second <code>false</code> means not to separate the
 * position.  Otherwise, output start [tab] stop [tab].  The last
 * <code>false</code> is to not print the addition of both strands.  Note that
 * if <code>true</code>, will print the total column only if the strand is
 * nog ignored.
 */
template<typename CountType>
void writeSequenceAlignmentFragmentToFeatureTrackerCounts(SequenceAlignmentFragmentToFeaturesTracker<CountType>& theTracker, string filename, bool printAllLines = false, SequenceLocation::element_type offset = 1, bool separatePos = false, bool printCombined = false)
{

    Ofstream fout(filename);
    SequenceFeatures& theFeatures = theTracker.getFeatures();
    StringQualifiers theHeader = theFeatures.getQualifiersOfIndex(0);
    if (!separatePos)
    {
        fout << "Feature ID\tChr\tStrand\tPosition\t" << theHeader.getStringKeys();
    }
    else
    {
        fout << "Feature ID\tChr\tStrand\tStart\tStop\t" << theHeader.getStringKeys();
    }
    if (!theTracker.ignoringStrand())
    {
        fout << "\tSame strand count";
        fout << "\tOpposite strand count";
    }
    else
    {
        fout << "\tCount";
    }
    if (!theTracker.ignoringStrand() && printCombined)
    {
        fout << "\tTotal count";
    }
    fout << endl;
    for (size_t i = 0; i < theFeatures.size(); ++i)
    {
        string geneID = theFeatures.getFeatureIDOfIndex(i);
        StringQualifiers theQualifiers = theFeatures.getQualifiersOfIndex(i);
        string theValues = theQualifiers.getStringValues();
        SequenceLocation currLoc = theFeatures.getLocationOfIndex(i);
        CountType senseNorm = theTracker.getNumQueryContigsAssociatedToFeatureIndex(i, true);
        CountType antisenseNorm = theTracker.getNumQueryContigsAssociatedToFeatureIndex(i, false);
        if (printAllLines || ((senseNorm + antisenseNorm) > 0))
        {
            SequenceLocation::position_type adjustedPos = currLoc.position();
            adjustedPos += offset;
            fout << geneID << "\t"
                 << currLoc.contig() << "\t"
                 << currLoc.strand() << "\t";
            if (!separatePos)
            {
                fout << adjustedPos << "\t";
            }
            else
            {
                fout << adjustedPos.lower() << "\t" << adjustedPos.upper() << "\t";
            }
            fout << theValues;
            fout << "\t" << senseNorm;
            if (!theTracker.ignoringStrand())
            {
                fout << "\t" << antisenseNorm;
            }
            if (!theTracker.ignoringStrand() && printCombined)
            {
                CountType temp = static_cast<CountType>(0);
                temp += senseNorm;
                temp += antisenseNorm;
                fout << "\t" << temp;
            }
            fout << endl;
        }
    }
    fout.close();
}




/** \file SequenceAlignmentFragmentToFeatureTrackerWriter.hpp
 * This file defines two functions to write output of the #SequenceAlignmentFragmentToFeaturesTracker class to a file.
 * \fn template<typename CountTypeNorm, typename CountTypeUnnorm> void writeSequenceAlignmentFragmentToFeatureTrackerCounts(SequenceAlignmentFragmentToFeaturesTracker<CountTypeNorm>& theNormTracker, SequenceAlignmentFragmentToFeaturesTracker<CountTypeUnnorm>& theUnnormTracker, string filename, bool printAllLines = false, SequenceLocation::element_type offset = 1, bool separatePos = false, bool printCombined = false)
 * This prints out only genes that got hit (unless the below boolean parameter
 * is <code>true</code> - default is <code>false</code>).  Here, we print both
 * the normalized and unnormalized versions of the tracker for the same #Features.
 * \section Usage Usage:
 * \code
 * // ...
 * SequenceAlignmentFragmentToFeaturesTracker<CountTypeNorm> theNormTracker;
 * SequenceAlignmentFragmentToFeaturesTracker<CountTypeUnnorm> theUnnormTracker;
 * writeSequenceAlignmentFragmentToFeatureTrackerCounts(
 *    theTracker, theUnnormTracker, "theTracker.tsv", false, 1, false, false, false);
 * \endcode
 * where the last <code>1</code> is an applied index to the internal 0-based
 * index coordinates.  The second <code>false</code> means not to separate the
 * position.  Otherwise, output start [tab] stop [tab].  The last
 * <code>false</code> is to print the addition of both strands.  Note that
 * if <code>true</code>, will print the total column only if the strand is
 * not ignored.
 */
template<typename CountTypeNorm, typename CountTypeUnnorm>
void writeSequenceAlignmentFragmentToFeatureTrackerCounts(SequenceAlignmentFragmentToFeaturesTracker<CountTypeNorm>& theNormTracker, SequenceAlignmentFragmentToFeaturesTracker<CountTypeUnnorm>& theUnnormTracker, string filename, bool printAllLines = false, SequenceLocation::element_type offset = 1, bool separatePos = false, bool printCombined = false)
{

    Ofstream fout(filename);
    SequenceFeatures& theFeatures = theNormTracker.getFeatures();
#ifdef DEBUG
    if (printCombined)
    {
        // Can only print total if either both ignore strand or both uses it.
        assert(theNormTracker.ignoringStrand() ==
               theUnnormTracker.ignoringStrand());
    }

    SequenceFeatures& theUnnormFeatures = theUnnormTracker.getFeatures();
    // Features for normalized and unnormalized counts must be the same.
    assert(&theFeatures == &theUnnormFeatures);
#endif
    StringQualifiers theHeader = theFeatures.getQualifiersOfIndex(0);
    if (!separatePos)
    {
        fout << "Feature ID\tChr\tStrand\tPosition\t" << theHeader.getStringKeys();
    }
    else
    {
        fout << "Feature ID\tChr\tStrand\tStart\tStop\t" << theHeader.getStringKeys();
    }
    if (!theNormTracker.ignoringStrand())
    {
        fout << "\tSame strand normalized count\tSame strand unnormalized count";
        fout << "\tOpposite strand normalized count\tOpposite strand unnormalized count";
    }
    else
    {
        fout << "\tNormalized count\tUnnormalized count";
    }
    if ((!theNormTracker.ignoringStrand() || !theUnnormTracker.ignoringStrand())
            && printCombined)
    {
        fout << "\tTotal normalized count\tTotal unnormalized count";
    }
    fout << endl;
    for (size_t i = 0; i < theFeatures.size(); ++i)
    {
        string geneID = theFeatures.getFeatureIDOfIndex(i);
        StringQualifiers theQualifiers = theFeatures.getQualifiersOfIndex(i);
        string theValues = theQualifiers.getStringValues();
        SequenceLocation currLoc = theFeatures.getLocationOfIndex(i);
        CountTypeNorm senseNorm = theNormTracker.getNumQueryContigsAssociatedToFeatureIndex(i, true);
        CountTypeNorm antisenseNorm = theNormTracker.getNumQueryContigsAssociatedToFeatureIndex(i, false);
        CountTypeUnnorm senseUnnorm = theUnnormTracker.getNumQueryContigsAssociatedToFeatureIndex(i, true);
        CountTypeUnnorm antisenseUnnorm = theUnnormTracker.getNumQueryContigsAssociatedToFeatureIndex(i, false);

        if (printAllLines ||
                ((senseNorm + senseUnnorm + antisenseNorm + antisenseUnnorm) > 0))
        {
            SequenceLocation::position_type adjustedPos = currLoc.position();
            adjustedPos.assign(adjustedPos.lower() + offset,
                               adjustedPos.upper() + offset);
            fout << geneID << "\t"
                 << currLoc.contig() << "\t"
                 << currLoc.strand() << "\t";
            if (!separatePos)
            {
                fout << adjustedPos << "\t";
            }
            else
            {
                fout << adjustedPos.lower() << "\t" << adjustedPos.upper() << "\t";
            }
            fout << theValues;
            ;
            fout << "\t" << senseNorm << "\t" << senseUnnorm;
            if (!theNormTracker.ignoringStrand())
            {
                fout << "\t" << antisenseNorm << "\t" << antisenseUnnorm;
            }
            if ((!theNormTracker.ignoringStrand() || !theUnnormTracker.ignoringStrand()) && printCombined)
            {
                CountTypeNorm temp1 = static_cast<CountTypeNorm>(0);
                CountTypeUnnorm temp2 = static_cast<CountTypeUnnorm>(0);
                if (!theNormTracker.ignoringStrand())
                {
                    temp1 += senseNorm;
                    temp2 += senseUnnorm;
                }
                if (!theUnnormTracker.ignoringStrand())
                {
                    temp1 += antisenseNorm;
                    temp2 += antisenseUnnorm;
                }
                fout << "\t" << temp1 << "\t" << temp2;
            }
            fout << endl;
        }
    }
    fout.close();
}

#endif


