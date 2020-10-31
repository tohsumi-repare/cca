#ifndef MOLBIOLIB_FASTATABLEREADERWRITER_H
#define MOLBIOLIB_FASTATABLEREADERWRITER_H 1

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
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"



/** \file FastaTableReaderWriter.hpp
 * This file defines two functions to read and write FASTA files.
 * \fn template<size_t N = 0, typename... T> void readFastaTable(Table<T...>& theTable, string filename, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0, bool readFullHeader = false, bool checkUnique = false)
 * The target columns (3 consecutive ones starting at N [default = 0]) must be
 * of the form (or implicitly castable to) <size_t, string, Sequence>.  The
 * leading column
 * of numbers is the index in case a subset the table is used and compared
 * against some aligners' output who use the index as the identifier instead
 * of the name.  Two functions are
 * defined, #readFastaTable and #writeFastaTable with
 * string filename arguments.  Optionally, one may specify the minimum
 * number of rows to allocate (thus making reads go faster) and then
 * afterwards optionaly specify the length of string to reserve (often
 * making reads of long lines a *lot* faster).  This is similar to the
 * method #readTSVTable.  Rows are filled from row 0.  (Can change this
 * by pushing back empty rows in the #Table.)
 * \section Usage Usage:
 * \code
 * Table<string, FASTA_ROW_TYPE, int, int> myTable;
 * \endcode
 * Below reserves 100000 lines in the table and 10000000 for the
 * input string line.  Since 0 and numeric_limits<size_t>::max() are specified, it reads from the
 * first line to the last (numeric_limits<size_t>::max() = last).  Obviously, 0 and numeric_limits<size_t>::max() could be
 * changed.  0 and numeric_limits<size_t>::max() are the default values.  The
 * <code>false</code> means not to read the full header.  Reading in the full
 * header (instead of the first non-whitespace word) is a rare occurance for
 * when one needs to do regex search on the headers.  The
 * last <code>false</code> below is to check for unique headers.  If true,
 * then will check, otherwise will not.
 * \code
 * readFastaTable<1>(myTable, "myGenome.fasta", 0, numeric_limits<size_t>::max(), 100000, 10000000, false, false);
 * \endcode
 * Note that one could use fileSize to overestimate the size of the line
 * length to use.  One departure from a correct fasta file is allowed by
 * readFastaTable: blank lines may appear anywhere in the body.  They are
 * ignored.
 */
template<size_t N = 0, typename... T>
void readFastaTable(Table<T...>& theTable, string filename,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    size_t numAlloc = 0,
                    size_t minStringReserve = 0,
                    bool readFullHeader = false,
                    bool checkUnique = false)
{

    assert(tuple_size<typename Table<T...>::row_type>::value >= N+2);

    if (filename == "")
    {
#ifdef DEBUG
        cerr << "Warning from readFastaTable.  An empty filename was passed.  This will return an empty Fasta Table." << endl;
#endif
        return;
    }

    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error!  The file " << filename << " cannot be opened." << endl;
        assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
    }
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    set<string> uniqueNamesCheck;

    theTable.labels[N] = "Header";
    theTable.labels[N+1] = "Sequence";


    size_t currentRead = 0;
    if (numAlloc > 0)
    {
        theTable.reserve(numAlloc);
    }
    string line, tempLine;
    vector<string> tokens;
    if (minStringReserve > 0)
    {
        line.reserve(minStringReserve);
        tempLine.reserve(minStringReserve);
    }
    Ifstream fin(filename);
    getline(fin, line);
    bool proceed = true, stop = false;
    while (!fin.fail() && !stop)
    {
        assert(line[0] == '>');  // Check for validity of FASTA file.
        proceed = true;
        if (currentRead < first)
        {
            proceed = false;
        }
        if (last != numeric_limits<size_t>::max() &&
                currentRead > last)
        {
            proceed = false;
            stop = true;
        }
        if (proceed)
        {
            if (!readFullHeader)
            {
                splitString(line, "\t", tokens);
                tokens[0].erase(0, 1);   // Get rid of leading >
                vector<string> tempTokens;
                splitString(tokens[0], " ", tempTokens);
                bool breakNow = false;
                size_t goodToken = 0;
                for (size_t i = 0; !breakNow && (i < tempTokens.size()); ++i)
                {
                    if (tempTokens[i] != "")
                    {
                        goodToken = i;
                        breakNow = true;
                    }
                }
                tokens[0] = tempTokens[goodToken];
            }
            else
            {
                tokens.resize(1, "");
                tokens[0] = line;  // Make the whole thing the first "token".
                tokens[0].erase(0, 1);   // Get rid of leading >
            }
            // Want unique names for FASTA sequences, else can give rise to
            // errors in other routines which depend on uniqueness of names.
            if (checkUnique &&
                    (uniqueNamesCheck.find(tokens[0]) != uniqueNamesCheck.end()))
            {
                cerr << "Error in readFastaTable.  A non-unique name - " << tokens[0] << " - was found.  Exiting." << endl;
            }
            assert(!checkUnique ||
                   (uniqueNamesCheck.find(tokens[0]) == uniqueNamesCheck.end()));
            if (checkUnique)
            {
                uniqueNamesCheck.insert(tokens[0]);
            }
            get<N>(theTable.theData).push_back(tokens[0]);  // Push back name without >.
        }
        ++currentRead;
        size_t currRow = 0;
        if (proceed)
        {
            string nullString = "";
            get<N+1>(theTable.theData).push_back(nullString);
            currRow = get<N+1>(theTable.theData).size() - 1;
            if (minStringReserve > 0)
            {
                get<N+1>(theTable.theData)[currRow].reserve(minStringReserve);
            }
            get<N+1>(theTable.theData)[currRow].clear();
        }
        do
        {
            getline(fin, line);
            if (fin.fail())
            {
                break;
            }
            tempLine = trimSpacesString(line);
            if (tempLine == "")
            {
                continue;    // blank line.
            }
            if (tempLine[0] == ';')
            {
                continue;    // If start with ;, means comment
            }
            if (tempLine[0] != '>' && proceed)
            {
                get<N+1>(theTable.theData)[currRow] += tempLine;
            }
        }
        while (tempLine[0] != '>');
    }
    fin.close();
}



/** \file FastaTableReaderWriter.hpp
 * This file defines two functions to read and write FASTA files.
 * \fn template<size_t N = 0, typename... T> void writeFastaTable(Table<T...>& theTable, string filename, size_t first = 0, size_t last = numeric_limits<size_t>::max())
 * \section Usage Usage:
 * \code
 * writeFastaTable(myTable, "theirGenome.fasta");
 * \endcode
 */
template<size_t N = 0, typename... T>
void writeFastaTable(Table<T...>& theTable, string filename,
                     size_t first = 0, size_t last = numeric_limits<size_t>::max())
{

#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
    assert(tuple_size<typename Table<T...>::row_type>::value >= N+2);
#endif

    size_t numRows = theTable.size();
    if (last != numeric_limits<size_t>::max())
    {
        numRows = last + 1;
    }

    Ofstream fout(filename);
    for (size_t i = first; i < numRows; ++i)
    {
        fout << '>' << get<N>(theTable.theData)[i] << endl;
        fout << get<N+1>(theTable.theData)[i] << endl;
    }
    fout.close();
}


#endif


