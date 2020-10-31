#ifndef MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTTOFEATURESTRACKER_H
#define MOLBIOLIB_SEQUENCEALIGNMENTFRAGMENTTOFEATURESTRACKER_H 1

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


#include "src/Objects/Features.hpp"
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/FileSequenceAlignmentFragmentReader.hpp"


/** \file SequenceAlignmentFragmentToFeaturesTracker.hpp
 * Contains the #SequenceAlignmentFragmentToFeaturesTracker and
 * #SequenceAlignmentFragmentReaderToFeaturesTracker.
 */

/** Tracks which features are hit by which queries.
 * Track features hit by queries.  May do only counting, if necessary.
 * \section Usage Usage:
 * \code
 * SequenceFeatures refFeature;
 *    // SequenceFeatures is Features<SequenceLocation> typedef'd.
 * // ...
 * readFeaturesTableRefFeature(refFeature, refFeatureFile, refLinkFile);
 * SequenceAlignmentFragmentToFeaturesTracker<double> theTracker(refFeature, false, true, false, false, false);
 *   // The template parameters are optional.  Defaults shown.
 *   // The double is optional (default).  Can change how counted by changing
 *   // the double data type.
 *   // Above optional Boolean parameters are shown with default values.
 *   // first false = do not ignore strand information -
 *   //               match only if strands agree.
 *   // first true = do not reverse strands.  If false, then reverse strands.
 *   // second false = require only overlap - true means require query to be
 *   //                completely enveloped in the feature, i.e. subset.
 *   // third false = track features associated with queries (uses more memory)
 *   // fourth false = track queries associated with features (uses more memory)
 * ...
 * SequenceLocation queryLocation, targetLocation;
 * // ...
 * size_t numHits;
 * numHits = theTracker.addAlignmentLocationData(queryLocation,
 *                                               targetLocation, 1);
 *   // Returns the number of features hit.
 *   // 1 is the amount one counts the read.  One may instead use (when doing
 *   // double as the data type for counting) frequency/# places it aligns, i.e.
 *   // if theAligns is derived from AlignmentFragmentReader, then
 *   // string theContig = (theAligns().getQueryLocation()).contig();
 *   // double amtCount = frequency[theContig]/
 *   //                   (static_cast<double>(theAligns.numberAlignmentFragmentsRelatedToQueryContig(theContig));
 *   // theTracker.addAlignmentLocationData(queryLocation, targetLocation, amtCount);
 * // Then, there are the following functions:
 * StringTable contigsHit = theTracker.getContigsHit(true)
 *    // Above is a list of contigs hit by the alignment file
 *    // on the sense strand (false for antisense strand - default = true).
 *    // Subsequent boolean value of true means sense strand.
 * IndexTable featuresHit = theTracker.getIndiciesOfFeaturesHit(true)
 *    // Above is a list of feature indices (since one can do refFeature.getRow(i)) hit.
 * // One can then iterate through the above (use size() and [], like a vector)
 * // and use the methods
 * CountType numFeatures = getNumFeaturesAssociatedToQueryContig(contigsHit[i], true);
 * CountType numQueries = getNumQueryContigsAssociatedToFeatureIndex(featuresHit[i], true);
 * // If one did not pass true to the latter two parameters of the constructor,
 * // one could also call:
 * set<size_t>& theIndexSet = theTracker.getSetFeatureIndicesAssociatedToQueryContig(contigsHit[i], true);
 * IndexTable theIndexTable = theTracker.getFeatureIndicesAssociatedToQueryContig(contigsHit[i], true);
 * // where the difference between the two is that the latter requires a copy
 * // of a #Table whereas the former gives back a reference to a set of indices.
 * // Similarly,
 * set<string> theContigSet = theTracker.getSetContigsAssociatedToFeatureIndex(featuresHit[i], true);
 * StringTable theContigTable = theTracker.getContigsAssociatedToFeatureIndex(featuresHit[i], true);
 * // One may access the #Features being tracked by doing:
 * SequenceFeatures& theFeaturesUsed = theTracker.getFeatures();
 * // One may also access the type of analysis being done by doing the
 * // following.  This is more for programs that are using this class.
 * bool ignoreStrand = theTracker.ignoringStrand();
 * bool matchSense = theTracker.matchingSense();
 * bool useSubset = theTracker.usingSubset();
 * \endcode
 *
 * There is also a parameter, maxLoadFactor that comes after all the Boolean
 * parameters, which may also be set via the function
 * <code>setMaxLoadFactor(float maxLoadFactor)</code>.)  See the documentation
 * for this in Fasta.hpp.  The default is 1.0.
 *
 */
template<typename CountType = double>
class SequenceAlignmentFragmentToFeaturesTracker
{
public:
    void setMaxLoadFactor(float maxLoadFactor)
    {
        featuresAssocToQueryContigSense.max_load_factor(maxLoadFactor);
        featuresAssocToQueryContigAntisense.max_load_factor(maxLoadFactor);
        numFeaturesAssocToQueryContigSense.max_load_factor(maxLoadFactor);
        numFeaturesAssocToQueryContigAntisense.max_load_factor(maxLoadFactor);
        queryContigsAssocToFeatureSense.max_load_factor(maxLoadFactor);
        queryContigsAssocToFeatureAntisense.max_load_factor(maxLoadFactor);
        numQueryContigsAssocToFeatureSense.max_load_factor(maxLoadFactor);
        numQueryContigsAssocToFeatureAntisense.max_load_factor(maxLoadFactor);
    }

    SequenceAlignmentFragmentToFeaturesTracker(SequenceFeatures& inputFeatures,
            bool inputIgnoreStrand = false,
            bool inputKeepSense = true,
            bool inputIsSubset = false,
            bool inputDoNotTrackFeaturesAssocToQueries = false,
            bool inputDoNotTrackQueriesAssocToFeatures = false,
            bool inputDivideByFeatureLength = false,
            float maxLoadFactor = 1.0) :
        theFeatures(inputFeatures),
        ignoreStrand(inputIgnoreStrand), keepSense(inputKeepSense),
        isSubset(inputIsSubset),
        doNotTrackFeaturesAssocToQueries(inputDoNotTrackFeaturesAssocToQueries),
        doNotTrackQueriesAssocToFeatures(inputDoNotTrackQueriesAssocToFeatures),
        divideByFeatureLength(inputDivideByFeatureLength)
    {
        setMaxLoadFactor(maxLoadFactor);
    }

    size_t addAlignmentLocationData(SequenceLocation& queryLocation,
                                    SequenceLocation& targetLocation,
                                    CountType amtCount = static_cast<CountType>(1))
    {
        string& queryContig = queryLocation.contig();
        // First do sense (or lack thereof if ignoringStrand),
        // then do antisense if not ignoringStrand.
        size_t numAssocSense = 0, numAssocAntisense = 0;
        IndexTable associatedFeaturesSense;
        // Below keepSense explanation:
        //    If keepSense if true, then the below keepSense is
        //    true, meaning we look on the same strand as the feature.  On
        //    the other hand, if keepSense is true, then we look on the
        //    opposite strand as the feature.
        theFeatures.featuresWithLocation(targetLocation,
                                         associatedFeaturesSense,
                                         ignoreStrand, keepSense,
                                         isSubset);
        numAssocSense += associatedFeaturesSense.size();
        for (size_t i = 0; i < numAssocSense; ++i)
        {
            size_t currentFeatureIndex = associatedFeaturesSense[i];
            CountType realCount = amtCount;
            if (divideByFeatureLength)
            {
                SequenceFeatures::position_type theFeaturePos = theFeatures.getLocationOfIndex(currentFeatureIndex).position();
                long featureLength = theFeaturePos.width() + 1;
                realCount = realCount/static_cast<CountType>(featureLength);
            }
            if (numFeaturesAssocToQueryContigSense.find(queryContig) ==
                    numFeaturesAssocToQueryContigSense.end())
            {
                numFeaturesAssocToQueryContigSense[queryContig] = realCount;
            }
            else
            {
                numFeaturesAssocToQueryContigSense[queryContig] += realCount;
            }
            if (numQueryContigsAssocToFeatureSense.find(currentFeatureIndex) ==
                    numQueryContigsAssocToFeatureSense.end())
            {
                numQueryContigsAssocToFeatureSense[currentFeatureIndex] = realCount;
            }
            else
            {
                numQueryContigsAssocToFeatureSense[currentFeatureIndex] += realCount;
            }

            if (!doNotTrackFeaturesAssocToQueries)
            {
                featuresAssocToQueryContigSense[queryContig].insert(currentFeatureIndex);
            }
            if (!doNotTrackQueriesAssocToFeatures)
            {
                queryContigsAssocToFeatureSense[currentFeatureIndex].insert(queryContig);
            }
        }

        if (!ignoreStrand)
        {
            IndexTable associatedFeaturesAntisense;
            // See above's keepSense explanation.
            //    This is just the complement of it.
            theFeatures.featuresWithLocation(targetLocation,
                                             associatedFeaturesAntisense,
                                             ignoreStrand, !keepSense,
                                             isSubset);
            numAssocAntisense += associatedFeaturesAntisense.size();
            for (size_t i = 0; i < numAssocAntisense; ++i)
            {
                size_t currentFeatureIndex = associatedFeaturesAntisense[i];
                CountType realCount = amtCount;
                if (divideByFeatureLength)
                {
                    SequenceFeatures::position_type theFeaturePos = theFeatures.getLocationOfIndex(currentFeatureIndex).position();
                    long featureLength = theFeaturePos.width() + 1;
                    realCount = realCount/static_cast<CountType>(featureLength);
                }
                if (numFeaturesAssocToQueryContigAntisense.find(queryContig) ==
                        numFeaturesAssocToQueryContigAntisense.end())
                {
                    numFeaturesAssocToQueryContigAntisense[queryContig] = realCount;
                }
                else
                {
                    numFeaturesAssocToQueryContigAntisense[queryContig] += realCount;
                }
                if (numQueryContigsAssocToFeatureAntisense.find(currentFeatureIndex) ==
                        numQueryContigsAssocToFeatureAntisense.end())
                {
                    numQueryContigsAssocToFeatureAntisense[currentFeatureIndex] = realCount;
                }
                else
                {
                    numQueryContigsAssocToFeatureAntisense[currentFeatureIndex] += realCount;
                }

                if (!doNotTrackFeaturesAssocToQueries)
                {
                    featuresAssocToQueryContigAntisense[queryContig].insert(currentFeatureIndex);
                }
                if (!doNotTrackQueriesAssocToFeatures)
                {
                    queryContigsAssocToFeatureAntisense[currentFeatureIndex].insert(queryContig);
                }
            }
        }

        return (numAssocSense + numAssocAntisense);
    }


    StringTable getContigsHit(bool sense = true)
    {
        StringTable resultTable;
        if (sense)
        {
            for (typename unordered_map<string, CountType>::iterator i = numFeaturesAssocToQueryContigSense.begin();
                    i != numFeaturesAssocToQueryContigSense.end(); ++i)
            {
                resultTable.push_back(i->first);
            }
        }
        else
        {
            for (typename unordered_map<string, CountType>::iterator i = numFeaturesAssocToQueryContigAntisense.begin();
                    i != numFeaturesAssocToQueryContigAntisense.end(); ++i)
            {
                resultTable.push_back(i->first);
            }
        }
        return resultTable;
    }


    IndexTable getIndiciesOfFeaturesHit(bool sense = true)
    {
        IndexTable resultTable;
        if (sense)
        {
            for (typename unordered_map<size_t, CountType>::iterator i = numQueryContigsAssocToFeatureSense.begin();
                    i != numQueryContigsAssocToFeatureSense.end(); ++i)
            {
                resultTable.push_back(i->first);
            }
        }
        else
        {
            for (typename unordered_map<size_t, CountType>::iterator i = numQueryContigsAssocToFeatureAntisense.begin();
                    i != numQueryContigsAssocToFeatureAntisense.end(); ++i)
            {
                resultTable.push_back(i->first);
            }
        }
        return resultTable;
    }


    CountType getNumFeaturesAssociatedToQueryContig(string theContig,
            bool sense = true)
    {
        if (sense)
        {
            return numFeaturesAssocToQueryContigSense[theContig];
        }
        else
        {
            return numFeaturesAssocToQueryContigAntisense[theContig];
        }
    }


    set<size_t>& getSetFeatureIndicesAssociatedToQueryContig(string theContig,
            bool sense = true)
    {
#ifdef DEBUG
        assert(!doNotTrackFeaturesAssocToQueries);
#endif
        if (sense)
        {
            return featuresAssocToQueryContigSense[theContig];
        }
        else
        {
            return featuresAssocToQueryContigAntisense[theContig];
        }
    }


    IndexTable getFeatureIndicesAssociatedToQueryContig(string theContig,
            bool sense = true)
    {
#ifdef DEBUG
        assert(!doNotTrackFeaturesAssocToQueries);
#endif
        IndexTable theResult;
        set<size_t>& theSet = getSetFeatureIndicesAssociatedToQueryContig(theContig, sense);
        for (set<size_t>::iterator i = theSet.begin(); i != theSet.end(); ++i)
        {
            theResult.push_back(*i);
        }
        return theResult;
    }


    CountType getNumQueryContigsAssociatedToFeatureIndex(size_t theFeature,
            bool sense = true)
    {
        if (sense)
        {
            return numQueryContigsAssocToFeatureSense[theFeature];
        }
        else
        {
            return numQueryContigsAssocToFeatureAntisense[theFeature];
        }
    }


    set<string>& getSetContigsAssociatedToFeatureIndex(size_t theFeature,
            bool sense = true)
    {
#ifdef DEBUG
        assert(!doNotTrackQueriesAssocToFeatures);
#endif
        if (sense)
        {
            return queryContigsAssocToFeatureSense[theFeature];
        }
        else
        {
            return queryContigsAssocToFeatureAntisense[theFeature];
        }
    }


    StringTable getContigsAssociatedToFeatureIndex(size_t theFeature,
            bool sense = true)
    {
#ifdef DEBUG
        assert(!doNotTrackQueriesAssocToFeatures);
#endif
        StringTable theResult;
        set<string>& theSet = getSetContigsAssociatedToFeatureIndex(theFeature, sense);
        for (typename set<string>::iterator i = theSet.begin(); i != theSet.end(); ++i)
        {
            theResult.push_back(*i);
        }
        return theResult;
    }


    SequenceFeatures& getFeatures()
    {
        return theFeatures;
    }


    bool ignoringStrand()
    {
        return ignoreStrand;
    }

    bool keepingSense()
    {
        return keepSense;
    }

    bool usingSubset()
    {
        return isSubset;
    }



private:
    SequenceFeatures& theFeatures;
    bool ignoreStrand, keepSense, isSubset, doNotTrackFeaturesAssocToQueries, doNotTrackQueriesAssocToFeatures, divideByFeatureLength;
    unordered_map < string, set<size_t> > featuresAssocToQueryContigSense, featuresAssocToQueryContigAntisense;
    unordered_map <string, CountType> numFeaturesAssocToQueryContigSense, numFeaturesAssocToQueryContigAntisense;
    unordered_map < size_t, set<string> > queryContigsAssocToFeatureSense, queryContigsAssocToFeatureAntisense;
    unordered_map <size_t, CountType> numQueryContigsAssocToFeatureSense, numQueryContigsAssocToFeatureAntisense;
};





/** Tracks which features are hit by which queries on a #SequenceAlignmentFragmentReader.
 * Like #SequenceAlignmentFragmentToFeaturesTracker, but also stores an
 * #SequenceAlignmentFragmentReader.  Given the above, a constructor
equenceAlignmentFragmentReaderToFeaturesTracker
 * \code
 * // Here, define something derived from #AlignmentFragmentReader, e.g.
 * // Define string filename to be something then...
 * NCBIBlastnTableAlignmentFragmentReader theReader(filename);
 * SequenceAlignmentFragmentReaderToFeaturesTracker<double> theTracker(refFeature, theReader, false, false, false, false, false);
 * // where the Boolean parameters are the same as before and the double is
 * // optional.
 * \endcode
 * Instead of calling addAlignmentLocationData, one may do
 * \code
 * numHits = theTracker.addCurrentFragmentData(amtCount);
 * \endcode
 */
template<typename CountType = double>
class SequenceAlignmentFragmentReaderToFeaturesTracker : public SequenceAlignmentFragmentToFeaturesTracker<CountType>
{
public:
    SequenceAlignmentFragmentReaderToFeaturesTracker(SequenceFeatures& theFeatures,
            SequenceAlignmentFragmentReader& inputReader,
            bool ignoreStrand = false,
            bool keepSense = true,
            bool isSubset = false,
            bool inputDoNotTrackFeaturesAssocToQueries = false,
            bool inputDoNotTrackQueriesAssocToFeatures = false,
            bool divideByFeatureLength = false,
            float maxLoadFactor = 1.0) :
        SequenceAlignmentFragmentToFeaturesTracker<CountType>(
            theFeatures, ignoreStrand, keepSense, isSubset,
            inputDoNotTrackFeaturesAssocToQueries,
            inputDoNotTrackQueriesAssocToFeatures,
            divideByFeatureLength,
            maxLoadFactor
        ), theReader(inputReader)
    {
    }

    size_t addCurrentFragmentData(CountType amtCount = static_cast<CountType>(1))
    {
        SequenceAlignmentFragment theFragment = theReader();
        SequenceLocation& queryLocation = theFragment.getQueryLocation();
        SequenceLocation& targetLocation = theFragment.getTargetLocation();
        return SequenceAlignmentFragmentToFeaturesTracker<CountType>::addAlignmentLocationData(queryLocation, targetLocation, amtCount);
    }

private:
    SequenceAlignmentFragmentReader& theReader;
};





/** Tracks which features are hit by which queries on a #FileSequenceAlignmentFragmentReader.
 * Like #SequenceAlignmentFragmentToFeaturesTracker, but instead stores an
 * #FileSequenceAlignmentFragmentReader.  Given the above, a constructor
 * for this would look like
 * \code
 * // Here, define something derived from #AlignmentFragmentReader, e.g.
 * // Define string filename to be something then...
 * FileSequenceAlignmentFragmentReader theReader(alignmentFile,
 *                                               headerSizesFilename,
 *                                               theReads);
 * // where alignmentFile (string) contains the path and filename of the
 * // alignment file, the headerSizesFilename contains the path and filename
 * // of the output of readsHeaderSizes with KEEP_ONLY_FIRST_WORD=true and 
 * // HEADER_SIZES, and theReads are the initialized #Fasta data structure of 
 * // the reads.
 * FileSequenceAlignmentFragmentReaderToFeaturesTracker<double> theTracker(refFeature, theReader, false, true, false, false, false);
 * // where the Boolean parameters are the same as before and the double is
 * // optional.
 * \endcode
 * Instead of calling addAlignmentLocationData, one may do
 * \code
 * numHits = theTracker.addCurrentFragmentData(amtCount);
 * \endcode
 */
template<typename CountType = double>
class FileSequenceAlignmentFragmentReaderToFeaturesTracker : public SequenceAlignmentFragmentToFeaturesTracker<CountType>
{
public:
    FileSequenceAlignmentFragmentReaderToFeaturesTracker(SequenceFeatures& theFeatures,
            FileSequenceAlignmentFragmentReader& inputReader,
            bool ignoreStrand = false,
            bool keepSense = true,
            bool isSubset = false,
            bool inputDoNotTrackFeaturesAssocToQueries = false,
            bool inputDoNotTrackQueriesAssocToFeatures = false,
            bool divideByFeatureLength = false,
            float maxLoadFactor = 1.0) :
        SequenceAlignmentFragmentToFeaturesTracker<CountType>(
            theFeatures, ignoreStrand, keepSense, isSubset,
            inputDoNotTrackFeaturesAssocToQueries,
            inputDoNotTrackQueriesAssocToFeatures,
            divideByFeatureLength, maxLoadFactor
        ), theReader(inputReader)
    {
    }

    size_t addCurrentFragmentData(CountType amtCount = static_cast<CountType>(1))
    {
        SequenceAlignmentFragment theFragment = theReader();
        SequenceLocation& queryLocation = theFragment.getQueryLocation();
        SequenceLocation& targetLocation = theFragment.getTargetLocation();
        return SequenceAlignmentFragmentToFeaturesTracker<CountType>::addAlignmentLocationData(queryLocation, targetLocation, amtCount);
    }

private:
    FileSequenceAlignmentFragmentReader& theReader;
};



#endif

