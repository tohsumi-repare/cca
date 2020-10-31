#ifndef MOLBIOLIB_RIDGELINEPEAKS_H
#define MOLBIOLIB_RIDGELINEPEAKS_H 1

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


#include "src/Objects/Peaks.hpp"



/** \file RidgeLinePeaks.hpp
 * \fn template<typename T, bool USE_POSITION = true> void checkIfClosestPeak(T& currPos, T& pos, size_t& i, bool& found, size_t& result, T& dist)
 * Function to test if the current peak at index i is really the closest
 * and if so, change found, result, and dist to match that.
 */
template<typename T, bool USE_POSITION = true>
void checkIfClosestPeak(T& currPos, T& pos, size_t& i,
                        bool& found, size_t& result, T& dist)
{
    T tempDist = numeric_limits<T>::max();
    if (!USE_POSITION)
    {
        // Do below so do not have wraparound 0 for size_t.
        if (pos < currPos)
        {
            tempDist = currPos - pos;
        }
        else
        {
            tempDist = pos - currPos;
        }
    }
    else
    {
        tempDist = abs(pos - currPos);
    }
    if (dist > tempDist)
    {
        found = true;
        result = i;
        dist = tempDist;
    }
}



/** \file RidgeLinePeaks.hpp
 *\fn template<typename P, typename T, bool USE_POSITION = true> size_t closestPeak(Peaks<P>& currPeaks, vector<bool>& peakUsed, T pos, T window)
 * Given a set of peaks and peaksUsed and position, find the closest
 * peak, such that it is not used and within +/- 2*window of pos.
 * T below is either size_t or long.
 */
template<typename P, typename T, bool USE_POSITION = true>
size_t closestPeak(Peaks<P>& currPeaks, vector<bool>& peakUsed,
                   T pos, T window)
{

#ifdef PROG_DEBUG
    if (USE_POSITION)
    {
        assert(typeid(T) == typeid(long));
    }
    else
    {
        assert(typeid(T) == typeid(size_t));
    }
#endif

    T dist = numeric_limits<T>::max();
    size_t result = numeric_limits<size_t>::max();
    if (currPeaks.size() == 0)
    {
        return result;
    }
    bool found = false;
    P tempPeak;
    size_t low = numeric_limits<size_t>::max(),
           up = numeric_limits<size_t>::max();
    if (!USE_POSITION)
    {
        tempPeak.posIndex = pos;
        low = lower_boundIndex(currPeaks.peaks, 0, numeric_limits<size_t>::max(),
                               compareTablePeaksPosIndex<P>);
        up = upper_boundIndex(currPeaks.peaks, 0, numeric_limits<size_t>::max(),
                              compareTablePeaksPosIndex<P>);
    }
    else
    {
        tempPeak.pos.position().assign(pos);
        low = lower_boundIndex(currPeaks.peaks, 0, numeric_limits<size_t>::max(),
                               compareTablePeaksPos<P>);
        up = upper_boundIndex(currPeaks.peaks, 0, numeric_limits<size_t>::max(),
                              compareTablePeaksPosIndex<P>);
    }
    // size_t i is used later so define now and set to low for below's loop.
    size_t i = low;
    T currPos = 0;
    // Check for peaks from low to high
    for ( ; i < up; ++i)
    {
        if (peakUsed[i])
        {
            continue;
        }
        currPos = USE_POSITION ?  currPeaks[i].pos.position().lower() :
                  currPeaks[i].posIndex;
        checkIfClosestPeak<T, USE_POSITION>(currPos, pos, i, found, result, dist);
    }
    // Check for peaks on left [negative] side of pos, up to 2*window away
    i = low;
    currPos = USE_POSITION ?  currPeaks[i].pos.position().lower() :
              currPeaks[i].posIndex;
    while (!found && (i > 0) && ((currPos + 2*window) >= pos))
    {
        --i;
        if (peakUsed[i])
        {
            continue;
        }
        currPos = USE_POSITION ?  currPeaks[i].pos.position().lower() :
                  currPeaks[i].posIndex;
        checkIfClosestPeak<T, USE_POSITION>(currPos, pos, i, found, result, dist);
    }
    // Check for peaks on right [positive] side of pos, up to 2*window away
    i = up;
    if (i < currPeaks.size())
        currPos = USE_POSITION ?  currPeaks[i].pos.position().lower() :
                  currPeaks[i].posIndex;
    while (!found && (i < currPeaks.size()) && ((pos + 2*window) <= currPos))
    {
        if (peakUsed[i])
        {
            continue;
        }
        currPos = USE_POSITION ?  currPeaks[i].pos.position().lower() :
                  currPeaks[i].posIndex;
        checkIfClosestPeak<T, USE_POSITION>(currPos, pos, i, found, result, dist);
        ++i;
    }

    return result;
}



/** \file RidgeLinePeaks.hpp
 * \fn template<typename P, bool USE_POSITION = true> bool findRidge(size_t& currPos, size_t& maxLevel, vector< Peaks<P> >& multiscalePeaks, vector< vector<bool> >& usedPeaks, size_t window, size_t minLen, size_t gap, size_t& levelLength, size_t& numGaps, vector<size_t>& ridgeIndices)
 * This routine will not return true if it does find a ridge that satisfies
 * the gap parameter, but does *not* satisfy the minLen parameter.  However,
 * in such a case, the usedPeaks will be marked as used.
 */
template<typename P, bool USE_POSITION = true>
bool findRidge(size_t& currPos, size_t& maxLevel,
               vector< Peaks<P> >& multiscalePeaks,
               vector< vector<bool> >& usedPeaks, size_t window,
               size_t minLen, size_t gap,
               size_t& levelLength, size_t& numGaps,
               vector<size_t>& ridgeIndices)
{

    ridgeIndices.clear();

    if (maxLevel == 0)
    {
        ridgeIndices.push_back(currPos);
        if (usedPeaks[0][currPos])
        {
            return false;
        }
        usedPeaks[0][currPos] = true;
        if (gap > 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    numGaps = 0;
    vector<size_t> goodIndices;
    vector<size_t> theLevels;

    size_t numLevels;

    // Either a start of a new ridge or part of one - either way, used is true.
    usedPeaks[maxLevel][currPos] = true;
    size_t& thePosIndex = multiscalePeaks[maxLevel][currPos].posIndex;
    size_t& thePosLong = multiscalePeaks[maxLevel][currPos].pos.position().lower();
    for (long currLevel = maxLevel - 1; currLevel >= 0; --currLevel)
    {
        size_t currLevelIndex = static_cast<size_t>(currLevel);
        size_t result = numeric_limits<size_t>::max();
        if (USE_POSITION)
        {
            result = closestPeak<P, long>(multiscalePeaks[currLevelIndex], usedPeaks[currLevelIndex], thePosLong, window);
        }
        else
        {
            result = closestPeak<P, size_t, false>(multiscalePeaks[currLevelIndex], usedPeaks[currLevelIndex], thePosIndex, window);
        }

        goodIndices.push_back(result);
        if (result != numeric_limits<size_t>::max())
        {
            theLevels.push_back(currLevelIndex);
        }
        else
        {
            theLevels.push_back(numeric_limits<size_t>::max());
        }
    }
    numLevels = theLevels.size();
    for (size_t i = 0; i < numLevels; ++i)
    {
        if (goodIndices[i] == numeric_limits<size_t>::max())
        {
            ++numGaps;
        }
        else
        {
            if (USE_POSITION)
            {
                if (static_cast<long>(window) < abs(thePosLong - multiscalePeaks[theLevels[i]][goodIndices[i]].pos.position().lower()))
                {
                    ++numGaps;
                }
            }
            else
            {
                size_t testPos = multiscalePeaks[theLevels[i]][goodIndices[i]].posindex;
                if (testPos < thePosIndex)
                {
                    if (window < (thePosIndex - testPos))
                    {
                        ++numGaps;
                    }
                }
                else
                {
                    if (window < (testPos - thePosIndex))
                    {
                        ++numGaps;
                    }
                }
            }
        }
    }

    if (numGaps <= gap)
    {
        for (size_t i = 0; i < numLevels; ++i)
            // Below check needed since a -1 means there is an accepted gap.
            if (goodIndices[i] != numeric_limits<size_t>::max())
            {
                usedPeaks[theLevels[i]][goodIndices[i]] = true;
            }
        if (numLevels >= minLen)
        {
            levelLength = numLevels;
            ridgeIndices.resize(maxLevel+1);
            for (size_t j = 0; j < maxLevel; ++j)
            {
                ridgeIndices[j] = 0;
            }
            ridgeIndices[maxLevel] = currPos;
            for (size_t i = 0; i < numLevels; ++i)
                if (theLevels[i] != numeric_limits<size_t>::max())
                {
                    ridgeIndices[theLevels[i]] = goodIndices[i];
                }
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}






/** \file RidgeLinePeaks.hpp
 * Find maximums defined by peaks at multiple scales.
 * \fn template<typename P, bool USE_POSTION = true> void ridgeLinePeaks(vector< Peaks<P> >& multiscalePeaks, vector<size_t>& window, Peaks<PeakRidgeLines>& outPeaks, size_t minLen = 2, size_t gap = 1)
 * Finds ridge lines peaks from peaks sampled at various scalings.
 * The only point of USE_POSTION is to use pos or posIndex.
 */
template<typename P, bool USE_POSITION = true>
void ridgeLinePeaks(vector< Peaks<P> >& multiscalePeaks,
                    vector<size_t>& window,
                    Peaks<PeakRidgeLines>& outPeaks,
                    size_t minLen = 2, size_t gap = 1)
{
    if (minLen > 0)
    {
        --minLen;
    }

    outPeaks.clear();
    size_t numLevels = multiscalePeaks.size();

#ifdef DEBUG
    if (numLevels == 0)
    {
        cerr << "Error in ridgeLinePeaks!  Number of levels in multiscalePeaks is zero.  Exiting." << endl;
        exit(1);
    }
    if (numLevels < minLen)
    {
        cerr  << "Error in ridgeLinePeaks!  Number of levels < minLen.  Exiting." << endl;
        exit(1);
    }
#endif

    // First sort each level's peaks.
    for (size_t i = 0; i < numLevels; ++i)
    {
#ifdef PROG_DEBUG
        cerr << "Num peaks at level " << i << " = " << multiscalePeaks[i].size()
             << endl;
#endif
        multiscalePeaks[i].sortPos();
    }
    vector< vector<bool> > usedPeak(numLevels);
    for (size_t i = 0; i < numLevels; ++i)
    {
        usedPeak.resize(multiscalePeaks[i].size(), false);
    }

    // Using long so --levelIndex when levelIndex = 0 will give -1 and
    // exit the loop properly.
    long lowestLevelLong = static_cast<long>(minLen);
    for (long levelLong = static_cast<long>(numLevels - 1);
            levelLong >= lowestLevelLong; --levelLong)
    {
        size_t level = static_cast<size_t>(levelLong);
        Peaks<P>& currLevel = multiscalePeaks[level];
        size_t numPeaks = currLevel.size();

        size_t levelLength = 0,
               numGaps = numeric_limits<size_t>::max();
        vector<size_t> ridgeIndices;

        for (size_t topPeak = 0; topPeak < numPeaks; ++topPeak)
        {
            if (usedPeak[level][topPeak])
            {
                continue;
            }
            if (findRidge<P, USE_POSITION>(topPeak, level, multiscalePeaks,
                                           usedPeak, window[level], minLen, gap,
                                           levelLength, numGaps, ridgeIndices))
            {
                PeakRidgeLines newPeak = currLevel[topPeak];
                newPeak.maxLevel = level + 1;
                newPeak.numGaps = numGaps;
                newPeak.levelLength = levelLength + 1;
                newPeak.ridgeIndices = ridgeIndices;
                outPeaks.push_back(newPeak);
            }
        }
    }
}

#endif

