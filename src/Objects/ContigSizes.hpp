#ifndef MOLBIOLIB_CONTIGSIZES_H
#define MOLBIOLIB_CONTIGSIZES_H 1

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


#include "src/Objects/Fasta.hpp"
#include "src/Objects/ReadOnlyDelimitedFile.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"


/** \file ContigSizes.hpp
 * Simple class to hold contig names and sizes.
 * Typically, this is used to hold output of
 * Application/Annotation/Reads/readsHeaderSizes.cpp specifying the
 * HEADER_SIZES_OUT parameter.  Usage:
 * \code
 * ContigSizes theSizes("ref.headerSizes", 1);
 * // Or ContigSizes theSizes;
 * //    theSizes.read("ref.headerSizes", 1);
 * // where the 1 (default is numeric_limits<size_t>::max()) is the column in
 * // which the contig size is given.  It is assumed the contig name is in the
 * // first column (i.e. column 0).
 * //
 * // Or Fasta ref;
 * //    // init ref.
 * //    ContigSizes theSizes(ref);
 * //
 * // Can also set sizes manually:
 * theSizes.setSize("chr1", 247249719);
 * // One may need to ignore the fact that chr1 is a new contig.  If so,
 * // pass true after the contig size in the above example.
 * //
 * // Can then get sizes (type size_t):
 * cout << theSizes["chr1"] << endl;
 * // or
 * cout << theSizes.getSize("chr1") << endl;
 * //
 * // Can get all the headers
 * StringTable headers = theSizes.getContigHeaders();
 * //
 * // Can write out headerSizes
 * theSizes.write("ref.headerSizes");
 * //
 * // Can clear everything:
 * theSizes.clear();
 * \endcode
 */
class ContigSizes
{
public:
    void read(string filename,
              size_t lengthColumn = numeric_limits<size_t>::max())
    {
        vector<string> tokens;
        bool haveCol = false;
        if (lengthColumn != numeric_limits<size_t>::max())
        {
            haveCol = true;
        }

        ReadOnlyTSVFile ifp(filename, true, true);  // Seq mode, no size.
        for (size_t i = 0; !ifp.fail(); ++i)
        {
            tokens = ifp[i];
            if (!haveCol)
            {
                contigSizes[tokens[0]] = convertFromString<size_t>(tokens[tokens.size() - 1]);
            }
            else
            {
                contigSizes[tokens[0]] = convertFromString<size_t>(tokens[lengthColumn]);
            }
        }
        ifp.close();
    }

    void write(string filename)
    {
        Ofstream ofp(filename);
        for (map<string, size_t>::iterator i = contigSizes.begin();
                i != contigSizes.end(); ++i)
        {
            ofp << i->first << "\t" << i->second << endl;
        }
        ofp.close();
    }

    ContigSizes() { }

    ContigSizes(string filename,
                size_t lengthColumn = numeric_limits<size_t>::max())
    {
        read(filename, lengthColumn);
    }

    ContigSizes(Fasta& initFasta)
    {
        size_t numFastaEntries = initFasta.size();
        for (size_t i = 0; i < numFastaEntries; ++i)
            // Ignore missing contigs in below calls.
        {
            setSize(initFasta.getHeader(i), initFasta.seqSize(i), true);
        }
    }

    size_t operator[](string contig)
    {
#ifdef DEBUG
        if (contigSizes.find(contig) == contigSizes.end())
        {
            cerr << "Error in contigSizes[].  " << contig << " not found.  Returning zero." << endl;
            return 0;
        }
#endif
        return contigSizes[contig];
    }

    size_t size()
    {
        return contigSizes.size();
    }

    size_t getSize(string contig)
    {
#ifdef DEBUG
        if (contigSizes.find(contig) == contigSizes.end())
        {
            cerr << "Error in contigSizes.getSize().  " << contig << " not found.  Returning zero." << endl;
            return 0;
        }
#endif
        return contigSizes[contig];
    }

    void setSize(string contig, size_t theSize,
                 bool ignoreMissingContig = false)
    {
#ifdef DEBUG
        if (!ignoreMissingContig &&
                (contigSizes.find(contig) == contigSizes.end()))
        {
            cerr << "Warning in contigSizes.setSize().  " << contig << " not found.  Creating new contig entry and continuing execution." << endl;
        }
#endif
        contigSizes[contig] = theSize;
    }

    StringTable getContigHeaders()
    {
        StringTable theResult;
        for (map<string, size_t>::iterator i = contigSizes.begin();
                i != contigSizes.end(); ++i)
        {
            theResult.push_back(i->first);
        }
        return theResult;
    }

    void clear()
    {
        contigSizes.clear();
    }

    ~ContigSizes()
    {
        clear();
    }

private:
    map<string, size_t> contigSizes;

};

#endif

