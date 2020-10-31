#ifndef MOLBIOLIB_QUAL_H
#define MOLBIOLIB_QUAL_H 1

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
#include "src/Functions/ReaderWriters/Tables/QualTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/ReadsTableReader.hpp"



#ifndef QUAL_ROW_TYPE
/** \file Qual.hpp
 * Contains only the #Qual class and methods and associated DEFINE.
 * \def QUAL_ROW_TYPE
 * We define the row_type used in the #Table so
 * it can be used elsewhere as part of a different table.  <i>E.g.</i>
 *
 * <code>
 * Table<int, string, QUAL_ROW_TYPE, int, int, bool> myDifferentTable;
 * </code>
 */
#define QUAL_ROW_TYPE string, string
#endif



/** Holds Quality information in a specialized #Table.
 * The below is a #Qual class, derived from a #Table.  (This class was not
 * called Quality, since it is too close to the #Qualifiers class' name.)
 * The leading column of numbers (of type size_t) is there
 * in case subsets of this table is used and compared against some
 * aligners' output who use the index as the identifier instead
 * of the name.
 *
 * \section Usage Usage:
 * \code
 * string tempStr;
 * Qual myQual("myGenome.qual");
 * \endcode
 * or
 * \code
 * Qual myQual;
 * myQual.read("myGenome.qual");
 * \endcode
 * Note that the above may also be fastq files.  Alternatively,
 * if one has a Fasta file and a default quality value, one
 * can do
 * \code
 * Fasta myFasta;
 * // ...
 * Qual myQual(myFasta, 20, 66);
 * // Or alternatively,
 * // Qual myQual;
 * // myQual.initializeFromFasta(myFasta, 20, 66);
 * \endcode
 * which will copy the headers from myFasta and fill myQual with quality
 * values of 20 with an ASCII offset of 66.
 *
 * Can then write a Qual file:
 * \code
 * myQual.write("myOtherGenome.qual");
 * \endcode
 * The above reads and writes a qual sequence file.  See the documentation
 * for the FASTQ reader and writer.
 *
 * Can get the fifth header via:
 * \code
 * tempStr = myQual.getHeader(4);
 * // Or, if one is sure of the uniqueness of a header, one can do:
 * tempStr = myQual.getHeader("Sample Header");
 * // and similarly below (in place of the index value of 4).
 * \endcode
 * or
 * \code
 * myQual.getHeader(4, tempStr);
 * \endcode
 * Can set the header via
 * \code
 * myQual.setHeader(4, tempStr);
 * \endcode
 * Can similarly get and set sequences with getQual and setQual, which is
 * of type <code>string</code>.  Qualities are coded in Sanger Phred, which is
 * from 0 to 93 using ASCII code 33 to 126.  Reader and writers of qualities
 * are expected to translate to and from this.  One may get a specific quality
 * character value (of type <code>char</code>) or quality score (of type
 * <code>int</code>) by doing
 * \code
 * cout << myQual.getQualChar(4, 2) << endl;
 * cout << myQual.getQualChar("Sample Header", 2) << endl;
 * cout << myQual.getQualScore(4, 2, 0) << endl;
 * cout << myQual.getQualScore("Sample Header", 2, 0) << endl;
 * \endcode
 * Furthermore, see <code>getQualScoreTable</code> if one wishes to get the
 * scores in a <code>Table</code> of <code>int</code>s form:
 * \code
 * // Below gets the quality score with an offset of 0.
 * Table<int> theQuals = myQual.getQualScoreTable("Sample Header", 0);
 * \endcode
 * A <code>Table</code> is used instead of a <code>vector</code> since then
 * the statistics on a <code>Table</code> column can be used (such as mean).
 *
 * *** <b> Note that the quality scores are by default offset by 0 (the last
 * parameter is optional), but may be changed to subtract the ASCII character
 * value by some integer. </b> ***
 *
 * One may add new entries via
 * \code
 * string myHeader = "New header";
 * string myQualValues = ">>>>>>>";
 * myQual.push_back(myHeader, myQualValues);
 * \endcode
 *
 * Can get all headers by doing
 * \code
 * StringTable myHeaders = myQual.getHeaders();
 * \endcode
 * after which, one may wish to do <code>myQual.clear(true)</code> to delete
 * the Qual information.  (This is when one needs to just extract the headers.)
 * order that the sequence appears in the FASTA file.
 *
 * One can also determine which indices correspond to which names headers.
 * The mapping between indices and headers is automatically done when a FASTA
 * file is loaded, but if the headers are changed or new sequences are added,
 * one must do
 * \code
 * myQual.regenerateHeaderMaps();
 * \endcode
 * One can then get the list of indices:
 * \code
 * IndexTable theIndices;
 * theIndices = myQual.getIndices("Sample Header");
 * // Or, if one is sure there is only one header...
 * size_t theIndex = myQual.getIndex("Sample Header");
 * \endcode
 * One may check to see if an index is present via
 * \code
 * bool isPresent = myQual.hasIndex("Sample Header");
 * \endcode
 *
 * There is also a parameter, maxLoadFactor that comes after all the other
 * parameters, which may also be set via the function
 * <code>setMaxLoadFactor(float maxLoadFactor)</code>.)  See the documentation
 * for this in Fasta.hpp.  The default is 1.0.
 *
 * \todo Someday, change headerMapping's type from map to unordered_map.  Due
 * to bug in gcc 4.3.1, gives weird performance hit when doing regenerateMaps  .
 */
class Qual : public Table<QUAL_ROW_TYPE>
{
public:
    void setMaxLoadFactor(float maxLoadFactor)
    {
        // headerMappings.max_load_factor(maxLoadFactor);
    }

    Qual()
    {
        labels.resize(2);
        labels[0] = "Header";
        labels[1] = "Quality";
        setMaxLoadFactor(1.0);
    }

    void clear()
    {
        Table<QUAL_ROW_TYPE>::clear();
        headerMappings.clear();
    }

    void regenerateHeaderMaps()
    {
        headerMappings.clear();
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            headerMappings[getHeader(i)].push_back(i);
        }
    }

    // Regarding readFullHeader, see comments above Fasta(...) constructor.
    Qual(string filename, size_t first = 0,
         size_t last = numeric_limits<size_t>::max(),
         size_t numAlloc = 0, size_t minStringReserve = 0,
         int offset = 33, int max = 93,
         float maxLoadFactor = 1.0,
         bool readFullHeader = false)
    {
        labels.resize(2);
        labels[0] = "Header";
        labels[1] = "Quality";
        setMaxLoadFactor(maxLoadFactor);
        readReadsTable(*this, filename, first, last, numAlloc, minStringReserve,
                       0, offset, max, false, readFullHeader);
        regenerateHeaderMaps();
    }

    void initializeFromFasta(Fasta& fastaTable, int defaultQuality, int offset = 33)
    {
        string theHeader;
        string theQual;
        size_t theSize;
        size_t numRows = fastaTable.size();
        char theQualValue = static_cast<char>(defaultQuality + offset);
        for (size_t i = 0; i < numRows; ++i)
        {
            fastaTable.getHeader(i, theHeader);
            theSize = fastaTable.seqSize(i);
            theQual.resize(theSize);
            for (size_t j = 0; j < theSize; ++j)
            {
                theQual[j] = theQualValue;
            }
            push_back(theHeader, theQual);
        }
    }

    Qual(Fasta& fastaTable, int defaultQuality, int offset = 33)
    {
        labels.resize(2);
        labels[0] = "Header";
        labels[1] = "Quality";
        initializeFromFasta(fastaTable, defaultQuality, offset);
    }

    string getHeader(size_t i)
    {
        row_type tempRow;
        getRow(i, tempRow);
        return get<0>(tempRow);
    }
    void getHeader(size_t i, string& theHeader)
    {
        row_type tempRow;
        getRow(i, tempRow);
        theHeader = get<0>(tempRow);
    }

    StringTable getHeaders()
    {
        StringTable result;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            result.push_back(getHeader(i));
        }
        return result;
    }


    void setHeader(size_t i, string& theHeader)
    {
        row_type tempRow;
        getRow(i, tempRow);
        get<0>(tempRow) = theHeader;
        setRow(i, tempRow);
    }


    void getQual(size_t i, string& theQuality)
    {
        row_type tempRow;
        getRow(i, tempRow);
        theQuality = get<1>(tempRow);
    }
    string getQual(size_t i)
    {
        row_type tempRow;
        getRow(i, tempRow);
        return get<1>(tempRow);
    }

    void getQual(string header, string& theQuality)
    {
        row_type tempRow;
        getRow(getIndex(header), tempRow);
        theQuality = get<1>(tempRow);
    }
    string getQual(string header)
    {
        row_type tempRow;
        getRow(getIndex(header), tempRow);
        return get<1>(tempRow);
    }

    void setQual(size_t i, string& theQuality)
    {
        row_type tempRow;
        getRow(i, tempRow);
        get<1>(tempRow) = theQuality;
        setRow(i, tempRow);
    }
    void setQual(string header, string& theQuality)
    {
        row_type tempRow;
        size_t i = getIndex(header);
        getRow(i, tempRow);
        get<1>(tempRow) = theQuality;
        setRow(i, tempRow);
    }


    char getQualChar(size_t contig, size_t pos)
    {
        return ((get<1>(this->theData))[contig])[pos];
    }
    char getQualChar(string header, size_t pos)
    {
        return ((get<1>(this->theData))[getIndex(header)])[pos];
    }

    int getQualScore(size_t contig, size_t pos, int offset = 0)
    {
        char theQual = getQualChar(contig, pos);
        return (static_cast<int>(theQual) - offset);
    }
    int getQualScore(string header, size_t pos, int offset = 0)
    {
        char theQual = getQualChar(header, pos);
        return (static_cast<int>(theQual) - offset);
    }

    Table<int> getQualScoreTable(size_t contig, int offset = 0)
    {
        size_t length = ((get<1>(this->theData))[contig]).size();
        Table<int> result(length, 0);
        for (size_t i = 0; i < length; ++i)
        {
            result[i] = getQualScore(contig, i, offset);
        }
        return result;
    }

    Table<int> getQualScoreTable(string header, int offset = 0)
    {
        return getQualScoreTable(getIndex(header), offset);
    }

    void push_back(string header, string theQuality)
    {
        get<0>(this->theData).push_back(header);
        get<1>(this->theData).push_back(theQuality);
    }



    IndexTable& getIndices(string header)
    {
        return headerMappings[header];
    }


    size_t getIndex(string header)
    {
#ifdef DEBUG
        // If using this function, should only be if the header is unique.
        assert(headerMappings[header].size() == 1);
#endif
        return headerMappings[header][0];
    }


    bool hasIndex(string header)
    {
        if (headerMappings.find(header) == headerMappings.end())
        {
            return false;
        }
        else
        {
            return true;
        }
    }



    void read(string filename,
              size_t first = 0, size_t last = numeric_limits<size_t>::max(),
              int offset = 33, int max = 93,
              size_t numAlloc = 0, size_t minStringReserve = 0)
    {
        readReadsTable(*this, filename, first, last, numAlloc, minStringReserve,
                       0, offset, max, false);
        regenerateHeaderMaps();
    }

    void write(string filename, int offset = 0)
    {
        writeQualTable(*this, filename, offset);
    }


protected:
    map<string, IndexTable> headerMappings;
    // unordered_map<string, IndexTable> headerMappings;

};

#endif

