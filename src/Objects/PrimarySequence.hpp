#ifndef MOLBIOLIB_PRIMARYSEQUENCE_H
#define MOLBIOLIB_PRIMARYSEQUENCE_H 1

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



/** \file PrimarySequence.hpp
 * Contains only the #PrimarySequence class.
 */

/** An extended <code>string</code> class to hold a base sequence.
 * This is an extension to the <code>string</code> class
 * to hold sequences and do
 * operations on them.  The usage is:
 * \code
 * PrimarySequence myDNA("ACTGACGTA");
 *    // NOT PrimarySequence myDNA = "ACTGACGTA";
 *    // But PrimarySequence myDNA; myDNA = "ACTGACGTA"; is acceptable.
 * \endcode
 * or
 * \code
 * string mySeq = "ACTGACGTA";
 * PrimarySequence myDNA(mySeq);
 * \endcode
 * or
 * \code
 * PrimarySequence myDNA;
 * myDNA.assign(mySeq);  // or myDNA.assign("ACTGACGTA");
 * \endcode
 * One can do all the things one can do
 * on a <code>string</code> one can do on a #PrimarySequence,
 * such as printing:
 * \code
 * cout << myDNA << endl;
 * \endcode
 * and reading:
 * \code
 * cin >> myDNA;
 * \endcode
 * Get the size/length of the sequence:
 * \code
 * cout << myDNA.size() << endl;
 * \endcode
 * or
 * \code
 * cout << myDNA.length() << endl;
 * \endcode
 * where the size type is given by size_t
 * (which are of the same type of variable as
 * a <code>string</code>'s size.  Typically one can just use
 * <code>size_t</code> to represent a length of a string., e.g.
 * \code
 * size_t m;
 * m = myDNA.size();
 * \endcode
 * For concatenation:
 * \code
 * PrimarySequence mySecondDNA("TTA");
 * myDNA = myDNA + mySecondDNA;
 * \endcode
 * One can find the first and last position within a sequence:
 * \code
 * myDNA.find("TTA");
 * \endcode
 * and
 * \code
 * myDNA.rfind("TTA", 2);
 * \endcode
 * respectively.  One may add the starting position as a second optional
 * parameter, as shown in rfind.  For find, it will then consider between
 * the pos and the end.  For rfind, it will then consider between the start
 * and the pos.  A sequence can have parts of
 * it replaced.  For example, to
 * replace the second through fifth character by "XYZ", one does
 * \code
 * myDNA.replace(1, 4, "XYZ");
 * \endcode
 * where 1 is the second position (recall everything is 0-based indexing)
 * and the 4 is the length of the bit to be replaced (positions 2,3,4,5 =
 * length 4.  One can delete the second through fifth characters by
 * using an empty string:
 * \code
 * myDNA.replace(1, 4, "");
 * \endcode
 * and one can insert without replacing by specifying a position, but
 * making the length 0:
 * \code
 * myDNA.replace(1, 0, "XYZ");
 * \endcode
 * One can use a combination of find and replace to replace an instance
 * of something in a sequence:
 * \code
 * myDNA.replace(myDNA.find("AACGT"), 5, "TTGCATTGCA");
 * \endcode
 * where the 5 is the length of what's in the find.  Similarly, one
 * could do
 * \code
 * PrimarySequence sourceSeq, targetSeq;
 * sourceSeq = "AACGT";
 * targetSeq = "TTGCATTGCA";
 * myDNA.replace(myDNA.find(sourceSeq), sourceSeq.size(), targetSeq);
 * \endcode
 * A sequence can be cleared/erased via:
 * \code
 * myDNA.clear();
 * \endcode
 * A sequence can be resized via:
 * \code
 * myDNA.resize(25);
 * \endcode
 * One can get subsequence
 * \code
 * mySecondDNA = myDNA.subsequence(2, 8);   // Start at the third position and
 *                                          // get 8 characters inclusive of
 *                                          // the third characters
 * \endcode
 * Alternatively,
 * \code
 * myDNA.subsequence(2, 8, mySecondDNA);
 * \endcode
 * We will not show the second variant on the remaining subusequence
 * examples, but assume they exist.
 * \code
 * cout << myDNA.trim5prime(3);   // Get myDNA, trimmed of 3 bases in the front.
 * cout << myDNA.trim3prime(5);   // Get myDNA, trimmed of 5 bases in the end.
 * cout << myDNA.subsequenceAfter("CCT");   // Get all characters after the
 *                                          // first occurance of CCT.
 * cout << myDNA.subsequenceBefore("CCT");
 * cout << myDNA.subsequenceAfterLast("AA");   // Get all characters after the
 *                                             // last occurance of CCT.
 * cout << myDNA.subsequenceBeforeLast("AA");
 * cout << myDNA.subsequenceBetween("CCT", "AA");
 * \endcode
 * Above is same as combining <code>After("CCT")</code>
 * with <code>BeforeLast("AA")</code>.
 * One can find the location of all instances of a regular expression.
 * <b>Warning: this can be slow!</b>  <i>Don't even think of
 * using this for alignments.</i>  :-)
 * \code
 * IndexTable theLocations, theLengths;
 * string theRegex = "UUA";   // or something more complicated like
 *                            // string theRegex = "U.UA[1-9]5x";
 * myDNA.findRegex(theRegex, theLocations);   // Puts all the locations
 *                                            // (0-based indices) of the
 *                                            // pattern found in myDNA in
 *                                            // theLocations.
 * myDNA.findRegex(theRegex, theLocations, theLengths);   // Same as above, but
 *                                                        // also gives the
 *                                                        // lengths.
 * \endcode
 * For example, if <code>myDNA</code> contained the
 * characters xla-425xlaXX5xla, then
 * finding xla would result in locations 0, 7, and 13 with match lengths
 * of 3, 3, and 3.  If instead one tried to find 5xl.*, then this results in
 * position 6 and length 10 since the .* uses up the rest of the sequence.
 *
 * If for some reason, one wishes to use #PrimarySequence in a C-type
 * function, like atof (characters [ASCII] to float), then one uses the
 * <code>c_str()</code> function to pass to the
 * function the C-version of what is in the
 * PrimarySequence object, <i>e.g.</i>
 * \code
 * double x = atof(myDNA.c_str());
 * \endcode
 * Hopefully, one will not ever need to use this functionality.
 * There is a more elegant way to convert sequences to
 * doubles using string streams, <i>e.g.</i>
 * \code
 *  PrimarySequence myNum("3.14");
 *  double x = convertFromString<double>(myNum);
 * \endcode
 * The above is really for your information only.  It is hoped in
 * the course of coding, one will never need to use this functionality.
 */
class PrimarySequence : public string
{
public:
    typedef size_t size_type;

    // Constructors below.
    // By default, we'll assume this is a DNA sequence...
    PrimarySequence() { }
    PrimarySequence(const PrimarySequence& newSeq)
    {
        assign(newSeq.c_str());
    }
    PrimarySequence(string newString)
    {
        assign(newString);
    }

    PrimarySequence& operator=(const string& theString)
    {
        assign(theString);
        return *this;
    }

    // size() and length() [from string] is available.

    PrimarySequence subsequence(size_t pos, size_t length)
    {
#ifdef PROG_DEBUG
        assert(size() >= pos + length);
#endif
        // Since PrimarySequence is derived from string, can do below.
        string tempStr = this->substr(pos, length);
        PrimarySequence tempSeq(tempStr);
        return tempSeq;
    }
    void subsequence(size_t pos, size_t length,
                     PrimarySequence& theSequence)
    {
#ifdef PROG_DEBUG
        assert(size() >= pos + length);
#endif
        theSequence.assign(substr(pos,length));
    }

    PrimarySequence trim5prime(size_t length)
    {
        size_t currSize = size();
#ifdef PROG_DEBUG
        assert(currSize >= length);
#endif
        return subsequence(length, currSize-length);
    }
    void trim5prime(size_t length, PrimarySequence& theSequence)
    {
        size_t currSize = size();
#ifdef PROG_DEBUG
        assert(currSize >= length);
#endif
        theSequence.assign(substr(length, currSize-length));
    }

    PrimarySequence trim3prime(size_t length)
    {
        size_t currSize = size();
#ifdef PROG_DEBUG
        assert(currSize >= length);
#endif
        return subsequence(0, currSize-length);
    }
    void trim3prime(size_t length, PrimarySequence& theSequence)
    {
        size_t currSize = size();
#ifdef PROG_DEBUG
        assert(currSize >= length);
#endif
        theSequence.assign(substr(0, currSize-length));
    }

    PrimarySequence subsequenceAfter(string delim)
    {
        size_t pos = this->find(delim) + delim.size();
        string tempStr = "";
        if (pos != string::npos)
        {
            tempStr = substr(pos, size() - pos);
        }
        PrimarySequence tempSeq(tempStr);
        return tempSeq;
    }
    void subsequenceAfter(string delim, PrimarySequence& theSequence)
    {
        size_t pos = this->find(delim) + delim.size();
        theSequence = "";
        if (pos != string::npos)
        {
            theSequence = substr(pos, size() - pos);
        }
    }
    PrimarySequence subsequenceBefore(string delim)
    {
        size_t pos = find(delim);
        string tempStr = "";
        if (pos != string::npos)
        {
            tempStr = substr(0, pos);
        }
        PrimarySequence tempSeq(tempStr);
        return tempSeq;
    }
    void subsequenceBefore(string delim, PrimarySequence& theSequence)
    {
        size_t pos = find(delim);
        theSequence = "";
        if (pos != string::npos)
        {
            theSequence = substr(0, pos);
        }
    }

    PrimarySequence subsequenceAfterLast(string delim)
    {
        size_t pos = rfind(delim) + delim.size();
        string tempStr = "";
        if (pos != string::npos)
        {
            tempStr = substr(pos, size() - pos);
        }
        PrimarySequence tempSeq(tempStr);
        return tempSeq;
    }
    void subsequenceAfterLast(string delim, PrimarySequence& theSequence)
    {
        size_t pos = rfind(delim) + delim.size();
        theSequence = "";
        if (pos != string::npos)
        {
            theSequence = substr(pos, size() - pos);
        }
    }
    PrimarySequence subsequenceBeforeLast(string delim)
    {
        size_t pos = rfind(delim);
        string tempStr = "";
        if (pos != string::npos)
        {
            tempStr = substr(0, pos);
        }
        PrimarySequence tempSeq(tempStr);
        return tempSeq;
    }
    void subsequenceBeforeLast(string delim, PrimarySequence& theSequence)
    {
        size_t pos = rfind(delim);
        theSequence = "";
        if (pos != string::npos)
        {
            theSequence = substr(0, pos);
        }
    }


    PrimarySequence subsequenceBetween(string firstDelim, string lastDelim)
    {
        PrimarySequence tempSeq, tempSeq2;
        subsequenceAfter(firstDelim, tempSeq);
        tempSeq.subsequenceBeforeLast(lastDelim, tempSeq2);
        return tempSeq2;
    }
    void subsequenceBetween(string firstDelim, string lastDelim,
                            PrimarySequence& theSequence)
    {
        PrimarySequence tempSeq;
        subsequenceAfter(firstDelim, tempSeq);
        tempSeq.subsequenceBeforeLast(lastDelim, theSequence);
    }



    void findRegex(string theRegex, IndexTable& theLocs,
                   IndexTable& theLengths)
    {
        regex_t preg;
#ifdef DEBUG
        if (regcomp(&preg, theRegex.c_str(), 0))
        {
            cerr << "Error in PrimarySequence::findRegex!  "
                 << theRegex << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&preg, theRegex.c_str(), 0);
#endif
        regmatch_t pmatch;
        size_t offset = 0;

        theLocs.clear();
        theLengths.clear();

        // While we keep finding a valid match (i.e. regexec returns 0),
        // push back location and lengths and move forward by pmatch.rm_eo,
        // which is the match length.
        while(regexec(&preg, (*this).c_str() + offset, 1, &pmatch, 0) == 0)
        {
#ifdef DEBUG
            assert(pmatch.rm_so >= 0);
            assert(pmatch.rm_eo > pmatch.rm_so);
#endif
            theLocs.push_back(offset + pmatch.rm_so);
            theLengths.push_back(pmatch.rm_eo - pmatch.rm_so);
            offset += pmatch.rm_eo;
        }
    }
    void findRegex(string theRegex, IndexTable& theLocs)
    {
        IndexTable theLengths;
        findRegex(theRegex, theLocs, theLengths);
    }


};


#endif

