#ifndef MOLBIOLIB_FASTQTABLEREADERWRITER_H
#define MOLBIOLIB_FASTQTABLEREADERWRITER_H 1

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


/** \file FastqTableReaderWriter.hpp
 * This file defines two functions to read and write FASTQ files.
 * \fn template<size_t N = 0, size_t M = 0, typename... T, typename...U> void readFastqTable(Table<T...>& theFastaTable, Table<U...>& theQualTable, string filename, int OFFSET = 0, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool skipFasta = false, bool skipQual = false, size_t numAlloc = 0, size_t minStringReserve = 0, bool readFullHeader = false, bool checkUnique = false, bool preserveEmptyQualHeaders = false)
 * The target columns (3 consecutive ones starting at N [default = 0]) must be
 * of the form (or implicitly castable to) <size_t, string, Sequence>.  The
 * OFFSET default is set for the Sanger standard of 0 - that is, no
 * offset is added to the character value.  This is used when converting other
 * formats, e.g. Solexa/Illumina's.  The leading column
 * of numbers is the index in case a subset the table is used and compared
 * against some aligners' output who use the index as the identifier instead
 * of the name.  Two functions are defined, #readFastqTable and
 * #writeFastqTable with
 * string filename arguments.  Optionally, one may specify the minimum
 * number of rows to allocate (thus making reads go faster) and then
 * afterwards optionaly specify the length of string to reserve (often
 * making reads of long lines a *lot* faster).  This is similar to the
 * method #readTSVTable.  Rows are filled from row 0.  (Can change this
 * by pushing back empty rows in the #Table.)
 * \section Usage Usage:
 * \code
 * Table<string, FASTA_ROW_TYPE, int, int> myFastaTable;
 * Table<string, int, QUAL_ROW_TYPE> myQualTable;
 * \endcode
 * Below reserves 100000 lines in the table and 10000000 for the
 * input string line.  Since 0 and numeric_limits<size_t>::max() are specified, it reads from the
 * first line to the last (numeric_limits<size_t>::max() = last).  Obviously, 0 and numeric_limits<size_t>::max() could be
 * changed.  0 and numeric_limits<size_t>::max() are the default values.  The two <code>false</code>
 * below mean not to ignore reading in the fasta and qual part respectively.
 * The third false means only read in first word of header, else read in whole
 * header.  The fourth false is not to check for the uniqueness of the header.
 * The last false is to preserve the quality headers.  This is if the headers
 * themselves are empty and are to be preserved as such.
 * \code
 * readFastqTable<1>(myFastaTable, myQualTable, "myGenome.fasta", 0, numeric_limits<size_t>::max(), false, false, 100000, 10000000, false, false, false);
 * \endcode
 * Note that one could use fileSize to overestimate the size of the line
 * length to use.  One departure from a correct fastq file is allowed by
 * readFastqTable: blank lines may appear anywhere in the body.  They are
 * ignored.
 * \todo NOT TESTED YET.  Someday, resolve why g++ has this bug to force
 * needing a foundAt Boolean variable.
 */
template<size_t N = 0, size_t M = 0, typename... T, typename... U>
void readFastqTable(Table<T...>& theFastaTable, Table<U...>& theQualTable,
                    string filename,
                    int OFFSET = 0,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    bool skipFasta = false, bool skipQual = false,
                    size_t numAlloc = 0,
                    size_t minStringReserve = 0,
                    bool readFullHeader = false,
                    bool checkUnique = false,
                    bool preserveEmptyQualHeaders = false)
{

    assert(skipFasta || tuple_size<typename Table<T...>::row_type>::value >= N+2);
    assert(skipFasta || tuple_size<typename Table<U...>::row_type>::value >= N+2);

    if (filename == "")
    {
#ifdef DEBUG
        cerr << "Warning from readFastqTable.  An empty filename was passed.  This will return an empty Fastq Table." << endl;
#endif
        return;
    }

    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error!  The file " << filename << " cannot be opened." << endl;
        assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
    }
#ifdef PROG_DEBUG
    assert(skipFasta || theFastaTable.sameRowSizes() == true);
    assert(skipQual || theQualTable.sameRowSizes() == true);
#endif
    set<string> uniqueNamesCheck;

    if (!skipFasta)
    {
        theFastaTable.labels[N] = "Header";
        theFastaTable.labels[N+1] = "Sequence";
    }
    if (!skipQual)
    {
        theQualTable.labels[N] = "Header";
        theQualTable.labels[N+1] = "Quality";
    }


    size_t currentRead = 0;
    if (numAlloc > 0)
    {
        if (!skipFasta)
        {
            theFastaTable.reserve(numAlloc);
        }
        if (!skipQual)
        {
            theQualTable.reserve(numAlloc);
        }
    }
    string line, tempLine, fastaHeader = "";
    vector<string> tokens;
    if (minStringReserve > 0)
    {
        line.reserve(minStringReserve);
        tempLine.reserve(minStringReserve);
    }
    Ifstream fin(filename);
    // Get first sequence header.
    getline(fin, line);
    bool proceed = true, stop = false;
    string currHeader;
    while (!fin.fail() && !stop)
    {
        assert(line[0] == '@');  // Check for validity of FASTQ file.
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
            break;
        }
        if (proceed)
        {
            if (!readFullHeader)
            {
                currHeader = line;
                splitString(line, "\t", tokens);
                tokens[0].erase(0, 1);   // Get rid of leading @
                currHeader.erase(0, 1);
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
                tokens[0].erase(0, 1);   // Get rid of leading @
            }
            // Want unique names for FASTA sequences, else can give rise to
            // errors in other routines which depend on uniqueness of names.
            if (checkUnique &&
                    (uniqueNamesCheck.find(tokens[0]) != uniqueNamesCheck.end()))
            {
                cerr << "Error in readFastqTable.  A non-unique name - " << tokens[0] << " - was found.  Exiting." << endl;
            }
            assert(!checkUnique ||
                   (uniqueNamesCheck.find(tokens[0]) == uniqueNamesCheck.end()));
            if (checkUnique)
            {
                uniqueNamesCheck.insert(tokens[0]);
            }
            if (!skipFasta)
            {
                // Push back name without >.
                fastaHeader = tokens[0];
                get<N>(theFastaTable.theData).push_back(fastaHeader);
            }
        }
        // Allocate space for sequence.
        size_t currRow = 0;
        if (proceed && !skipFasta)
        {
            string nullString = "";
            get<N+1>(theFastaTable.theData).push_back(nullString);
            currRow = get<N+1>(theFastaTable.theData).size() - 1;
            if (minStringReserve > 0)
            {
                get<N+1>(theFastaTable.theData)[currRow].reserve(minStringReserve);
            }
            get<N+1>(theFastaTable.theData)[currRow].clear();
        }
        // Get sequence.
        getline(fin, line);
#ifdef DEBUG
        if (fin.fail())
        {
            cerr << "Error in readFastqTable!  Premature ending of fastq file " << filename << ".  Exiting." << endl;
            assert(true == false);
        }
#endif
        tempLine = trimSpacesString(line);
        if (proceed && !skipFasta)
        {
            get<N+1>(theFastaTable.theData)[currRow] += tempLine;
        }
        // Now, get quality score header.
        // Should be same as sequence header so skip.
        getline(fin, line);
#ifdef DEBUG
        if (fin.fail())
        {
            cerr << "Error in readFastqTable!  Premature ending of fastq file " << filename << ".  Exiting." << endl;
            assert(true == false);
        }
        assert(line[0] == '+');
#endif
        tempLine = trimSpacesString(line);
        if (proceed)
        {
            // Now get the qual part...
            // Need not assert(line[0] == '+') since by above, is it.
            if (!readFullHeader)
            {
                splitString(tempLine, "\t", tokens);
                tokens[0].erase(0, 1);   // Get rid of leading +
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
                tokens[0] = tempLine;
                tokens[0].erase(0, 1);   // Get rid of leading +
            }
            // DEBUG: We are not checking if the quality names are unique since
            // it may be that the quality names are omitted for space conservation.
            if (!skipQual)
            {
                // Below is so that we can get the qual headers, if need be.
                // If there is no fasta, then fastaHeader would be empty anyway.
                if (!preserveEmptyQualHeaders && (tokens[0] == ""))
                {
                    tokens[0] = fastaHeader;
                }
                // Push back name without >.
                get<M>(theQualTable.theData).push_back(tokens[0]);
                string nullString = "";
                get<M+1>(theQualTable.theData).push_back(nullString);
                currRow = get<M+1>(theQualTable.theData).size() - 1;
                if (minStringReserve > 0)
                {
                    get<M+1>(theQualTable.theData)[currRow].reserve(minStringReserve);
                }
                get<M+1>(theQualTable.theData)[currRow].clear();
            }
        }


        // Finally, get quality line.
        getline(fin, line);
#ifdef DEBUG
        if (fin.fail())
        {
            cerr << "Error in readFastqTable!  Premature ending of fastq file " << filename << ".  Exiting." << endl;
            assert(true == false);
        }
#endif
        tempLine = trimSpacesString(line);
        if (proceed && !skipQual)
        {
            if (OFFSET != 0)
            {
                size_t numBases = tempLine.size();
                for (size_t i = 0; i < numBases; ++i)
                {
                    char theBase = tempLine[i];
                    int theValue = static_cast<int>(theBase) + OFFSET;
                    tempLine[i] = static_cast<char>(theValue);
                }
            }
            get<M+1>(theQualTable.theData)[currRow] += tempLine;
        }

        // Get next sequence header.
        getline(fin, line);
        ++currentRead;
    }
    fin.close();
}







/** \file FastqTableReaderWriter.hpp
 * This file defines two functions to read and write FASTQ files.
 * \fn template<size_t N = 0, size_t M = 0, typename... T, typename...U> void readFastqTableWraparound(Table<T...>& theFastaTable, Table<U...>& theQualTable, string filename, int OFFSET = 0, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool skipFasta = false, bool skipQual = false, size_t numAlloc = 0, size_t minStringReserve = 0, bool readFullHeader = false, bool checkUnique = false)
 * The target columns (3 consecutive ones starting at N [default = 0]) must be
 * of the form (or implicitly castable to) <size_t, string, Sequence>.  The
 * OFFSET default is set for the Sanger standard of 0 - that is, no
 * offset is added to the character value.  This is used when converting other
 * formats, e.g. Solexa/Illumina's.  The leading column
 * of numbers is the index in case a subset the table is used and compared
 * against some aligners' output who use the index as the identifier instead
 * of the name.  Two functions are defined, #readFastqTable and
 * #writeFastqTable with
 * string filename arguments.  Optionally, one may specify the minimum
 * number of rows to allocate (thus making reads go faster) and then
 * afterwards optionaly specify the length of string to reserve (often
 * making reads of long lines a *lot* faster).  This is similar to the
 * method #readTSVTable.  Rows are filled from row 0.  (Can change this
 * by pushing back empty rows in the #Table.)
 * \section Usage Usage:
 * \code
 * Table<string, FASTA_ROW_TYPE, int, int> myFastaTable;
 * Table<string, int, QUAL_ROW_TYPE> myQualTable;
 * \endcode
 * Below reserves 100000 lines in the table and 10000000 for the
 * input string line.  Since 0 and numeric_limits<size_t>::max() are specified, it reads from the
 * first line to the last (numeric_limits<size_t>::max() = last).  Obviously, 0 and numeric_limits<size_t>::max() could be
 * changed.  0 and numeric_limits<size_t>::max() are the default values.  The two <code>false</code>
 * below mean not to ignore reading in the fasta and qual part respectively.
 * The third false means only read in first word of header, else read in whole
 * header.  The last false means not to check for the uniqueness of the headers.
 * \code
 * readFastqTableWraparound<1>(myFastaTable, myQualTable, "myGenome.fasta", 0, numeric_limits<size_t>::max(), false, false, 100000, 10000000, false, false);
 * \endcode
 * Note that one could use fileSize to overestimate the size of the line
 * length to use.  One departure from a correct fastq file is allowed by
 * readFastqTable: blank lines may appear anywhere in the body.  They are
 * ignored.
 *
 * This one differs from readFastqTable in that it accepts multiline quality
 * and sequence strings.  However, it will incorrectly parse any quality string
 * that starts with an ampersand.
 *
 * \todo NOT TESTED YET.  Someday, resolve why g++ has this bug to force
 * needing a foundAt Boolean variable.
 */
template<size_t N = 0, size_t M = 0, typename... T, typename... U>
void readFastqTableWraparound(Table<T...>& theFastaTable,
                              Table<U...>& theQualTable,
                              string filename,
                              int OFFSET = 0,
                              size_t first = 0,
                              size_t last = numeric_limits<size_t>::max(),
                              bool skipFasta = false, bool skipQual = false,
                              size_t numAlloc = 0,
                              size_t minStringReserve = 0,
                              bool readFullHeader = false,
                              bool checkUnique = false)
{

    assert(skipFasta || tuple_size<typename Table<T...>::row_type>::value >= N+2);
    assert(skipFasta || tuple_size<typename Table<U...>::row_type>::value >= N+2);

    if (filename == "")
    {
#ifdef DEBUG
        cerr << "Warning from readFastqTableWraparound.  An empty filename was passed.  This will return an empty Fastq Table." << endl;
#endif
        return;
    }

    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error!  The file " << filename << " cannot be opened." << endl;
        assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
    }
#ifdef PROG_DEBUG
    assert(skipFasta || theFastaTable.sameRowSizes() == true);
    assert(skipQual || theQualTable.sameRowSizes() == true);
#endif
    set<string> uniqueNamesCheck;

    if (!skipFasta)
    {
        theFastaTable.labels[N] = "Header";
        theFastaTable.labels[N+1] = "Sequence";
    }
    if (!skipQual)
    {
        theQualTable.labels[N] = "Header";
        theQualTable.labels[N+1] = "Quality";
    }


    size_t currentRead = 0;
    if (numAlloc > 0)
    {
        if (!skipFasta)
        {
            theFastaTable.reserve(numAlloc);
        }
        if (!skipQual)
        {
            theQualTable.reserve(numAlloc);
        }
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
        assert(line[0] == '@');  // Check for validity of FASTQ file.
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
            break;
        }
        if (proceed)
        {
            if (!readFullHeader)
            {
                splitString(line, "\t", tokens);
                tokens[0].erase(0, 1);   // Get rid of leading @
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
                tokens[0].erase(0, 1);   // Get rid of leading @
            }
            // Want unique names for FASTA sequences, else can give rise to
            // errors in other routines which depend on uniqueness of names.
            if (checkUnique &&
                    (uniqueNamesCheck.find(tokens[0]) != uniqueNamesCheck.end()))
            {
                cerr << "Error in readFastqTableWraparound.  A non-unique name - " << tokens[0] << " - was found.  Exiting." << endl;
            }
            assert(!checkUnique ||
                   (uniqueNamesCheck.find(tokens[0]) == uniqueNamesCheck.end()));
            if (checkUnique)
            {
                uniqueNamesCheck.insert(tokens[0]);
            }
            if (!skipFasta)
            {
                get<N>(theFastaTable.theData).push_back(tokens[0]);    // Push back name without >.
            }
        }
        size_t currRow = 0;
        if (proceed && !skipFasta)
        {
            string nullString = "";
            get<N+1>(theFastaTable.theData).push_back(nullString);
            currRow = get<N+1>(theFastaTable.theData).size() - 1;
            if (minStringReserve > 0)
            {
                get<N+1>(theFastaTable.theData)[currRow].reserve(minStringReserve);
            }
            get<N+1>(theFastaTable.theData)[currRow].clear();
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
            if (tempLine[0] != '+' && proceed && !skipFasta)
            {
                get<N+1>(theFastaTable.theData)[currRow] += tempLine;
            }
        }
        while (tempLine[0] != '+');
        if (fin.fail())
        {
            cerr << "Error!  File " << filename << " has a sequence entry, but terminates before a quality entry.  Exiting.\n";
            assert(true == false);
        }
        // Now get the qual part...
        // Need not assert(line[0] == '+') since by above, is it.
        if (!readFullHeader)
        {
            splitString(tempLine, "\t", tokens);
            tokens[0].erase(0, 1);   // Get rid of leading +
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
            tokens[0] = tempLine;
            tokens[0].erase(0, 1);   // Get rid of leading +
        }
        // DEBUG: We are not checking if the quality names are unique since
        // it may be that the quality names are omitted for space conservation.
        if (proceed && !skipQual)
        {
            get<M>(theQualTable.theData).push_back(tokens[0]);  // Push back name without >.
            string nullString = "";
            get<M+1>(theQualTable.theData).push_back(nullString);
            currRow = get<M+1>(theQualTable.theData).size() - 1;
            if (minStringReserve > 0)
            {
                get<M+1>(theQualTable.theData)[currRow].reserve(minStringReserve);
            }
            get<M+1>(theQualTable.theData)[currRow].clear();
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
            // Looking for next sequence entry.
            if (proceed && !skipQual && tempLine[0] != '@')
            {
                if (OFFSET != 0)
                {
                    size_t numBases = tempLine.size();
                    for (size_t i = 0; i < numBases; ++i)
                    {
                        char theBase = tempLine[i];
                        int theValue = static_cast<int>(theBase) + OFFSET;
                        tempLine[i] = static_cast<char>(theValue);
                    }
                }
                get<M+1>(theQualTable.theData)[currRow] += tempLine;
            }
        }
        while (tempLine[0] != '@');
        ++currentRead;
    }
    fin.close();
}




















/** \file FastqTableReaderWriter.hpp
 * This file defines two functions to read and write FASTA files.
 * \fn template<size_t N = 0, size_t M = 0, typename... T, typename... U> void writeFastqTable(Table<T...>& theFastaTable, Table<U...>& theQualTable, string filename, int OFFSET = 0, size_t first = 0, size_t last = numeric_limits<size_t>::max())
 * \section Usage Usage:
 * \code
 * Fasta myFasta;
 * Qual myQual;
 * // ...
 * writeFastqTable(myFasta, myQual, "theirGenome.fastq", 0, 0, numeric_limits<size_t>::max());
 * \endcode
 */
template<size_t N = 0, size_t M = 0, typename... T, typename... U>
void writeFastqTable(Table<T...>& theFastaTable,
                     Table<U...>& theQualTable,
                     string filename,
                     int OFFSET = 0,
                     size_t first = 0,
                     size_t last = numeric_limits<size_t>::max())
{

#ifdef PROG_DEBUG
    assert(theFastaTable.sameRowSizes() == true);
    assert(theQualTable.sameRowSizes() == true);
    assert(tuple_size<typename Table<T...>::row_type>::value >= N+2);
    assert(tuple_size<typename Table<U...>::row_type>::value >= N+2);
#endif

    size_t numRows = theFastaTable.size();
    if (last != numeric_limits<size_t>::max())
    {
        numRows = last + 1;
    }

    Ofstream fout(filename);
    for (size_t i = first; i < numRows; ++i)
    {
        fout << '@' << get<N>(theFastaTable.theData)[i] << endl;
        fout << get<N+1>(theFastaTable.theData)[i] << endl;
        fout << '+' << get<M>(theQualTable.theData)[i] << endl;
        string& theQuals = get<M+1>(theQualTable.theData)[i];
        size_t numQuals = theQuals.size();
        if (OFFSET == 0)
        {
            fout << get<N+1>(theQualTable.theData)[i] << endl;
        }
        else
        {
            for (size_t j = 0; j < numQuals; ++j)
            {
                char theQual = theQuals[j];
                int qualValue = static_cast<int>(theQual) + OFFSET;
                fout << static_cast<char>(qualValue);
            }
            fout << endl;
        }
    }
    fout.close();
}


#endif


