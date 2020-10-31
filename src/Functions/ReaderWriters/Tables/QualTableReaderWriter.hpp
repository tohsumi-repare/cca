#ifndef MOLBIOLIB_QUALTABLEREADERWRITER_H
#define MOLBIOLIB_QUALTABLEREADERWRITER_H 1

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



/** \file QualTableReaderWriter.hpp
 * This file defines two functions to read and write Qual files.
 * \fn template<size_t N = 0, typename... T> void readQualTable(Table<T...>& theTable, string filename, int OFFSET = 33, int MAX = 93, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0, bool readFullHeader = false, bool checkUnique = false)
 * The target columns (3 consecutive ones starting at N [default = 0]) must be
 * of the form (or implicitly castable to) <size_t, string, string>.  The
 * leading column
 * of numbers is the index in case a subset the table is used and compared
 * against some aligners' output who use the index as the identifier instead
 * of the name.  Two functions are defined, #readQualTable and
 * #writeQualTable with
 * string filename arguments.  Optionally, one may specify the minimum
 * number of rows to allocate (thus making reads go faster) and then
 * afterwards optionaly specify the length of string to reserve (often
 * making reads of long lines a *lot* faster).  This is similar to the
 * method #readTSVTable.  Rows are filled from row 0.  The table
 * is filled using Phred conventions
 * by default (OFFSET = 33).  The maximum value is 93, but this can be
 * changed to any number.  numeric_limits<size_t>::max() means no maximum value, though the encoding is
 * limited to the maximum character value of 255.
 * \section Usage Usage:
 * \code
 * Table<string, QUAL_ROW_TYPE, int, int> myTable;
 * \endcode
 * Below reserves 100000 lines in the table and 10000000 for the
 * input string line.  Since 0 and numeric_limits<size_t>::max() are specified, it reads from the
 * first line to the last (numeric_limits<size_t>::max() = last).  Obviously, 0 and numeric_limits<size_t>::max() could be
 * changed.  0 and numeric_limits<size_t>::max() are the default values.  The 33 is the ASCII offset
 * and 93 is the maximum value.  The first <code>false</code> means to
 * read only the
 * first word and not the entire line for the headers.  The second
 * <code>false</code> is not to check for uniqueness of the headers.
 * \code
 * readQualTable<1>(myTable, "myGenome.fasta", 33, 93, 0, numeric_limits<size_t>::max(), 100000, 10000000, false, false);
 * \endcode
 * Note that one could use fileSize to overestimate the size of the line
 * length to use.  One departure from a correct qual file is allowed by
 * readQualTable: blank lines may appear anywhere in the body.  They are
 * ignored.
 * \todo NOT TESTED YET
 */
template<size_t N = 0, typename... T>
void readQualTable(Table<T...>& theTable, string filename,
                   int OFFSET = 33, int MAX = 93,
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
        cerr << "Warning from readQualTable.  An empty filename was passed.  This will return an empty Qual Table." << endl;
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
    theTable.labels[N+1] = "Quality";


    size_t currentRead = 0;
    if (numAlloc > 0)
    {
        theTable.reserve(numAlloc);
    }
    string line, tempLine;
    vector<int> theNums;
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
        assert(line[0] == '>');  // Check for validity of Qual file.
        proceed = true;
        if (currentRead < first)
        {
            proceed = false;
        }
        if ((last != numeric_limits<size_t>::max()) &&
                (currentRead > last))
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
                tokens[0] = line;
                tokens[0].erase(0, 1);   // Get rid of leading >
            }
            // Want unique names for Qual values, else can give rise to
            // errors in other routines which depend on uniqueness of names.
            assert(!checkUnique ||
                   (uniqueNamesCheck.find(tokens[0]) == uniqueNamesCheck.end()));
            if (checkUnique)
            {
                uniqueNamesCheck.insert(tokens[0]);
            }
            get<N>(theTable.theData).push_back(tokens[0]);  // Push back name without >.
        }
        ++currentRead;
        theNums.clear();
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
            if (tempLine[0] != '>')
            {
                splitString(tempLine, " ", tokens);
                size_t numTokens = tokens.size();
                for (size_t i = 0; i < numTokens; ++i)
                {
                    theNums.push_back(convertFromString<int>(tokens[i]));
                }
            }
        }
        while (tempLine[0] != '>');
        if (proceed)
        {
            string nullString = "";
            get<N+1>(theTable.theData).push_back(nullString);
            size_t currRow = get<N+1>(theTable.theData).size() - 1;
            size_t numValues = theNums.size();
            get<N+1>(theTable.theData)[currRow].resize(numValues);
            for (size_t i = 0; i < numValues; ++i)
            {
                // Sanger Phred scores: ASCII char - 33 = phred score.
                int theNum = theNums[i];
                // If > 93, then specification is to make it 93.
                if (MAX >= 0)
                {
                    if (theNum > MAX)
                    {
                        theNum = MAX;
                    }
                }
                get<N+1>(theTable.theData)[currRow][i] = static_cast<char>(theNum + OFFSET);
            }
        }
    }
    fin.close();
}



/** \file QualTableReaderWriter.hpp
 * This file defines two functions to read and write FASTA files.
 * \fn template<size_t N = 0, int OFFSET, typename... T> void writeQualTable(Table<T...>& theTable, string filename, int OFFSET = 33, size_t first = 0, size_t last = numeric_limits<size_t>::max())
 * \section Usage Usage:
 * \code
 * writeQualTable(myTable, "theirGenome.qual", 33);
 * \endcode
 * where 33 is the ASCII offset.
 */
template<size_t N = 0, typename... T>
void writeQualTable(Table<T...>& theTable, string filename,
                    int OFFSET = 33,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max())
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
        string& theQuals = get<N+1>(theTable.theData)[i];
        size_t numQuals = theQuals.size();
        for (size_t j = 0; j < numQuals; ++j)
        {
            fout << static_cast<int>(theQuals[j]) - OFFSET;
            if (j != numQuals - 1)
            {
                fout << " ";
            }
        }
        fout << endl;
    }
    fout.close();
}


#endif


