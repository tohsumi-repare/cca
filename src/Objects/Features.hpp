#ifndef MOLBIOLIB_FEATURES_H
#define MOLBIOLIB_FEATURES_H 1

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


// Below loads Qualifiers.hpp and Table.hpp
#include "src/Objects/Location.hpp"



#ifndef FEATURES_ROW_TYPE
/** \file Features.hpp
 * Contains only the #Features class and methods and associated DEFINE.
 * \def FEATURES_ROW_TYPE
 * The assumption is that the contigs are of type string.
 * Integers are, of course, a special case.  LocationType is the template
 * through which the #Location class type is passed.
 */
#define FEATURES_ROW_TYPE string, LocationType, StringQualifiers
#endif



#ifndef MIN_FEATURES_SEARCH_FAST
/** \file Features.hpp
 * Minimum number of features before a fast search algorithm is chosen.
 * \def MIN_FEATURES_SEARCH_FAST
 * This is the minimum number of features in an array needed before a fast
 * bisection algorithm for a search is used.  If it is smaller, it is
 * faster to use a more direct algorithm, much like STL's introsort.
 */
#define MIN_FEATURES_SEARCH_FAST 10
#endif


/** Holds features information in a specialized #Table.
 * The below is a #Features class, derived from a #Table.  Also includes a map
 * for fast access of particular types of features.
 *
 * \section Usage Usage:
 * \code
 * typedef Location<unsigned long, unsigned long> FeaturesLocType;
 * Features<FeaturesLocType> setOfFeatures;
 *   // If location is of type SequenceLocation, then can use SequenceFeatures
 * FeaturesLocType newLoc;
 * newLoc.contig() = "chrX";
 * FeaturesLocType::position_type newI(3, 5);   // i.e. Interval<unsigned long>
 * newLoc.position() = newI;
 * newLoc.strand() = "+";
 * StringQualifiers newQualifiers;
 * newQualifiers.push_back("Repeat class", "SINE");
 * setOfFeatures.addFeature("feature1", newLoc, newQualifiers);
 *    // ...
 * StringTable setOfChrs = setOfFeatures.getContigs();
 *    // getContigs() and getQualifiers() return only a
 *    // single instance of a name.  Alternatively, avoiding a copy do
 *    // StringTable setOfChrs;
 *    // setOfFeatures.getContigs(setOfChrs);
 * IndexTable& featuresOnChrX = setOfFeatures.featuresOnContig("chrX");
 * StringTable setOfQualifiers = setOfFeatures.getQualifiers();
 *    // Alternatively, setOfFeatures.getQualifiers(setOfQualifiers);
 * // One may also get specific ID, positions, and qualifiers:
 * cout << setOfFeatures.getFeatureIDOfIndex(0) << endl;
 * cout << setOfFeatures.getLocationOfIndex(0) << endl;
 * StringQualifiers firstQualifiers = setOfFeatures.getQualifiersOfIndex(0);
 * pair<string, string> checkQualifier;
 * checkQualifier.first = "Repeat class";
 * checkQualifier.second = "SINE";
 * // or do just:
 * // pair<string, string> checkQualifier("Repeat class", "SINE");
 * IndexTable& sineFeatures = setOfFeatures.featuresWithQualifier(checkQualifier);
 * IndexTable locFeatures;
 * setOfFeatures.featuresWithLocation(newLoc, locFeatures, false, true, false);
 *    // Where the first false is that we do not ignore the strand and the
 *    // second false is that we do not require the query to be a subset of the
 *    // feature, but rather overlaps it.  The true means we want the sense -
 *    // false for antisense.
 *    // Alternatively (at a cost of a copy),
 *    // IndexTable locFeatures = setOfFeatures.featuresWithLocation(newLoc, false, true, false);
 * // If the table is sorted, one may then do
 * setOfFeatures.setSorted(true);
 * // Then, the method of accessing which features are at a particular location is
 * // much faster.
 * \endcode
 */
template<typename LocationType>
class Features : public Table<FEATURES_ROW_TYPE>
{
public:
    typedef LocationType location_type;
    typedef typename location_type::position_type position_type;
    typedef typename location_type::position_type::element_type element_type;

    // See the attributes in the private section below for details on what
    // they contain.

    Features()
    {
        vector<string> newLabels(3);
        newLabels[0] = "FeatureID";
        newLabels[1] = "Location";
        newLabels[2] = "Qualifiers";
        this->setLabels(newLabels);
        sorted = false;
        disableContigMaps = false;
        disableQualifierMaps = false;
        maxWidth = 0;
    }

    string getFeatureIDOfIndex(size_t i)
    {
        return get<0>(this->getRow(i));
    }

    location_type getLocationOfIndex(size_t i)
    {
        return get<1>(this->getRow(i));
    }

    StringQualifiers getQualifiersOfIndex(size_t i)
    {
        return get<2>(this->getRow(i));
    }

    // The below is not the most efficient.  However, it is assumed that
    // this will be called a relatively few amount of times, so we are not
    // going to try to optimize this too much.
    void addFeature(string newFeatureID,
                    location_type& newLoc,
                    StringQualifiers& newQualifiers)
    {
        // Recall StringQualifiers basically have the
        // data structure map< string, Table<string> >
        tuple<FEATURES_ROW_TYPE> newRow;
        get<0>(newRow) = newFeatureID;
        get<1>(newRow) = newLoc;
        get<2>(newRow) = newQualifiers;
        if ((newLoc.position().width() + 1) > maxWidth)
        {
            maxWidth = newLoc.position().width() + 1;
        }
        this->push_back(newRow);
        // currentIndex is the last index of the Table<FEATURES_ROW_TYPE>
        size_t currentIndex = this->size() - 1;
        if (!disableContigMaps)
        {
            contigMapFeatures[newLoc.contig()].push_back(currentIndex);
        }
        if (!disableQualifierMaps)
        {
            StringTable theKeys = newQualifiers.getKeys();
            size_t numKeys = theKeys.size();
            for (size_t i = 0; i < numKeys; ++i)
            {
                StringTable& theValues = newQualifiers[theKeys[i]];
                size_t numValues = theValues.size();
                for (size_t j = 0; j < numValues; ++j)
                {
                    pair<string, string> newPair(theKeys[i], theValues[j]);
                    qualifierMapFeatures[newPair].push_back(currentIndex);
                }
            }
        }
    }

    IndexTable& featuresOnContig(string theContig)
    {
#ifdef DEBUG
        assert(!disableContigMaps);
#endif
        return contigMapFeatures[theContig];
    }


    void getContigs(StringTable& theContigs)
    {
#ifdef DEBUG
        assert(!disableContigMaps);
#endif
        set<string> contigNames;
        for (unordered_map<string, IndexTable>::iterator i = contigMapFeatures.begin(); i != contigMapFeatures.end(); ++i)
        {
            contigNames.insert(i->first);
        }
        // Below is so we can make distinct the names.
        for (set<string>::iterator i = contigNames.begin(); i != contigNames.end(); ++i)
        {
            theContigs.push_back(*i);
        }
    }

    StringTable getContigs()
    {
#ifdef DEBUG
        assert(!disableContigMaps);
#endif
        StringTable theContigs;
        getContigs(theContigs);
        return theContigs;
    }



    IndexTable& featuresWithQualifier(pair<string, string>& theQualifier)
    {
#ifdef DEBUG
        assert(!disableQualifierMaps);
#endif
        return qualifierMapFeatures[theQualifier];
    }


    void getQualifiers(StringTable& theQualifiers)
    {
#ifdef DEBUG
        assert(!disableQualifierMaps);
#endif
        set<string> qualifierNames;
        for (map<pair<string, string>, IndexTable>::iterator i = qualifierMapFeatures.begin(); i != qualifierMapFeatures.end(); ++i)
        {
            qualifierNames.insert((i->first).first);
        }
        // Below is so we can make distinct the names.
        for (set<string>::iterator i = qualifierNames.begin(); i != qualifierNames.end(); ++i)
        {
            theQualifiers.push_back(*i);
        }
    }

    StringTable getQualifiers()
    {
#ifdef DEBUG
        assert(!disableQualifierMaps);
#endif
        StringTable theQualifiers;
        getQualifiers(theQualifiers);
        return theQualifiers;
    }


    void featuresWithLocation(location_type& theLoc,
                              IndexTable& theFeaturesIndices,
                              bool ignoreStrand = false,
                              bool sense = true,
                              bool isSubset = false,
                              bool isOverlap = true,
                              bool isSuperset = false)
    {
#ifdef DEBUG
        assert(!disableContigMaps);
#endif
        string& queryContig = theLoc.contig();
        IndexTable& contigFeatures = contigMapFeatures[queryContig];
        string& queryStrand = theLoc.strand();
        position_type& queryInterval = theLoc.position();
        element_type queryLow = queryInterval.lower();
        element_type queryUp = queryInterval.upper();
        size_t startIndex, endIndex;
        if (sorted && (contigFeatures.size() >= MIN_FEATURES_SEARCH_FAST))
        {
            // Find the first feature such that the upper end is smaller than
            // the low - maxWidth of the query.

            // Below finds the first loc such that the lower position' upper +
            // maxWidth is >= queryInterval.lower()
            size_t loc, first = 0, count = contigFeatures.size(), step;
            while (count > 0)
            {
                loc = first;
                step = count/2;
                loc += step;
                if (get<1>(this->theData)[contigFeatures[loc]].position().upper() +
                        maxWidth < queryLow)
                {
                    ++loc;
                    first = loc;
                    count -= step + 1;
                }
                else
                {
                    count = step;
                }
            }
            startIndex = first;
            if (startIndex > 0)
            {
                --startIndex;
            }
            if (startIndex >= contigFeatures.size())
            {
                startIndex = contigFeatures.size() - 1;
            }


            // Below finds the first feature whose lower is greater
            // than query's upper + maxWidth.
            count = contigFeatures.size();
            size_t last = count - 1;
            while (count > 0)
            {
                loc = last;
                step = count/2;
                if (loc >= step)
                {
                    loc -= step;
                }
                if ((queryUp + maxWidth) <
                        get<1>(this->theData)[contigFeatures[loc]].position().lower())
                {
                    if (loc > 0)
                    {
                        --loc;
                    }
                    last = loc;
                    count -= step + 1;
                }
                else
                {
                    count = step;
                }
            }
            endIndex = last + 2;
            if (endIndex > contigFeatures.size())
            {
                endIndex = contigFeatures.size();
            }
        }
        else
        {
            startIndex = 0;
            endIndex = contigFeatures.size();
        }

        for (size_t i = startIndex; i < endIndex; ++i)
        {
            if ((isSubset &&
                    subset(queryInterval,
                           get<1>(this->theData)[contigFeatures[i]].position(),
                           true)
                ) ||
                    (isOverlap &&
                     overlap(queryInterval,
                             get<1>(this->theData)[contigFeatures[i]].position(),
                             true)
                    ) ||
                    (isSuperset &&
                     subset(get<1>(this->theData)[contigFeatures[i]].position(),
                            queryInterval, true)
                    )
               )
            {
                // Below done thusly to simplify the logic and structure.
                bool doPush = false;
                if (ignoreStrand)
                {
                    doPush = true;
                }
                else
                {
                    // Must check strand.
                    if (sense)
                    {
                        if (queryStrand ==
                                get<1>(this->theData)[contigFeatures[i]].strand())
                        {
                            doPush = true;
                        }
                    }
                    else
                    {
                        if (queryStrand !=
                                get<1>(this->theData)[contigFeatures[i]].strand())
                        {
                            doPush = true;
                        }
                    }
                }
                if (doPush)
                {
                    theFeaturesIndices.push_back(contigFeatures[i]);
                }
            }
        }
    }


    IndexTable featuresWithLocation(location_type& theLoc,
                                    bool ignoreStrand = false,
                                    bool sense = true,
                                    bool isSubset = false,
                                    bool isOverlap = true,
                                    bool isSuperset = false)
    {
#ifdef DEBUG
        assert(!disableContigMaps);
#endif
        IndexTable theResult;
        featuresWithLocation(theLoc, theResult, ignoreStrand, sense,
                             isSubset, isOverlap, isSuperset);
        return theResult;
    }


    void setSorted(bool isSorted)
    {
        sorted = isSorted;
    }

    void setDisableContigMaps()
    {
        disableContigMaps = true;
    }
    void setDisableQualifierMaps()
    {
        disableContigMaps = true;
    }
    void setDisableMaps()
    {
        disableContigMaps = true;
        disableQualifierMaps = true;
    }

protected:
    // contigMapFeatures map the contig name to the feature indices which are
    // associated with that contig.
    unordered_map<string, IndexTable> contigMapFeatures;
    // qualifierMapFeatures maps a key and its value to the indices which are
    // associated with that key and value.
    // Going with map instead of unordered_map, since unknown what is the
    // effect of hashing on pair<string, string>
    map<pair<string, string>, IndexTable> qualifierMapFeatures;
    bool sorted;
    bool disableContigMaps, disableQualifierMaps;

private:
    // The below is used in sorted searches since we know the maximum width of
    // the feature intervals, we can use this to give us a lower and upper bound
    // on the first and last features for a particular query interval.
    element_type maxWidth;

};


/** \file Features.hpp
 * \var typedef Features<SequenceLocation> SequenceFeatures
 * Typedef'd since this is used as the standard location container in MolBioLib.
 */
typedef Features<SequenceLocation> SequenceFeatures;



#endif

