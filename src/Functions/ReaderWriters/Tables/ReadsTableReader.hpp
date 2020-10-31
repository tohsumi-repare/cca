#ifndef MOLBIOLIB_READSTABLEREADER_H
#define MOLBIOLIB_READSTABLEREADER_H 1

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


#include "src/Objects/Sequence.hpp"
#include "src/Functions/ReaderWriters/Tables/FastaTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/FastqTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/QualTableReaderWriter.hpp"


/** \file ReadsTableReader.hpp
 * Universal reader for fasta or fastq files.
 * \fn bool isFastqFile(string filename)
 * Determines whether a file is a fastq file or not.
 */
bool isFastqFile(string filename)
{
    Ifstream ifp(filename);
    string line = "";
    getline(ifp, line);
    // Skip blank lines
    while (!ifp.fail() && (line == ""))
    {
        getline(ifp, line);
    }
    if (!ifp.fail())
    {
        if (line[0] == '@')
        {
            return true;
        }
        else if (line[0] == '>')
        {
            return false;
        }
    }
    else
    {
        return true;    // Empty fastq file.
    }
    cerr << "Not a valid Fasta, Qual, or Fastq file.  Exiting." << endl;
    exit(1);
}


/** \file ReadsTableReader.hpp
 * Universal reader for fasta or fastq files.
 * \fn template<size_t N = 0, size_t M = 0, typename... T, typename... U> void readReadsTable(Table<T...>& theFastaTable, Table<U...>& theQualTable, string filename, int OFFSET = 0, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0, bool readFullHeader = false, bool checkUnique = false)
 * Reads in a Fastq file into two tables. This is used when it is known that
 * a Fastq file is being read.
 */
template<size_t N = 0, size_t M = 0, typename... T, typename... U>
void readReadsTable(Table<T...>& theFastaTable, Table<U...>& theQualTable,
                    string filename,
                    int OFFSET = 0,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    size_t numAlloc = 0,
                    size_t minStringReserve = 0,
                    bool readFullHeader = false,
                    bool checkUnique = false)
{
    readFastqTable<N, M>(theFastaTable, theQualTable, filename, OFFSET, first, last, false, false, numAlloc, minStringReserve, readFullHeader, checkUnique);
}


/** \file ReadsTableReader.hpp
 * Universal reader for fasta or fastq files.
 * \fn template<size_t N = 0, typename... T> void readReadsTable(Table<T...>& theFastaQualTable, string filename, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0, int OFFSET_FASTQ = 0, int OFFSET_QUAL = 33, int MAX = 93, bool fastaTrue = true, bool readFullHeader = false, bool checkUnique = false)
 * Reads in either a fasta, qual, or fastq file.  Depending on
 * if <code>fastaTrue</code> is
 * <code>true</code> (read into a fasta table) or <code>false</code> (read into
 * a qual table), the appropriate reader and parameters are called.  By
 * default, a fasta file is read in.  Two different offsets are needed since
 * offset values are different between the two formats.
 */
template<size_t N = 0, typename... T>
void readReadsTable(Table<T...>& theFastaQualTable,
                    string filename,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    size_t numAlloc = 0,
                    size_t minStringReserve = 0,
                    int OFFSET_FASTQ = 0, int OFFSET_QUAL = 33, int MAX = 93,
                    bool fastaTrue = true, bool readFullHeader = false,
                    bool checkUnique = false)
{
    if (isFastqFile(filename))
    {
        if (fastaTrue)
        {
            Table<string, string> tempQ;
            readFastqTable<N>(theFastaQualTable, tempQ, filename, OFFSET_FASTQ, first, last, false, true, numAlloc, minStringReserve, readFullHeader, checkUnique);
        }
        else
        {
            Table<string, Sequence> tempF;
            readFastqTable<N>(tempF, theFastaQualTable, filename, OFFSET_FASTQ, first, last, true, false, numAlloc, minStringReserve, readFullHeader, checkUnique);
        }
    }
    else
    {
        if (fastaTrue)
        {
            readFastaTable<N>(theFastaQualTable, filename, first, last, numAlloc, minStringReserve, readFullHeader, checkUnique);
        }
        else
        {
            readQualTable<N>(theFastaQualTable, filename, OFFSET_QUAL, MAX, first, last, numAlloc, minStringReserve, readFullHeader, checkUnique);
        }
    }
}


#endif


