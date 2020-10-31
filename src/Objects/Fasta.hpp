#ifndef MOLBIOLIB_FASTA_H
#define MOLBIOLIB_FASTA_H 1

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
#include "src/Objects/Table.hpp"
#include "src/Objects/Sequence.hpp"
#include "src/Functions/ReaderWriters/Tables/FastaTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/ReadsTableReader.hpp"



/** \file Fasta.hpp
 * Contains only the #Fasta class and methods and associated DEFINE.
 * \def FASTA_ROW_TYPE
 * We define the row_type used in the #Table so
 * it can be used elsewhere as part of a different table.  <i>E.g.</i>
 *
 * <code>
 * Table<int, string, FASTA_ROW_TYPE, int, int, bool> myDifferentTable;
 * </code>
 */
#ifndef FASTA_ROW_TYPE
#define FASTA_ROW_TYPE string, Sequence
#endif



/** Holds FASTA information in a specialized #Table.
 * The below is a #Fasta class, derived from a #Table.  The
 * leading column of numbers (of type size_t) is there
 * in case subsets of this table is used and compared against some
 * aligners' output who use the index as the identifier instead
 * of the name.
 *
 * \section Usage Usage:
 * \code
 * string tempStr;
 * Fasta myFasta("myGenome.fasta");
 * \endcode
 * or
 * \code
 * Fasta myFasta;
 * myFasta.read("myGenome.fasta");
 * \endcode
 * Note that the above may also be fastq files.  Can then write a FASTA file:
 * \code
 * myFasta.write("myOtherGenome.fasta");
 * \endcode
 * Can get the fifth header via:
 * \code
 * tempStr = myFasta.getHeader(4);
 * // Or, if one is sure of the uniqueness of a header, one can do:
 * tempStr = myFasta.getHeader("Sample Header");
 * // and similarly below (in place of the index value of 4).
 * \endcode
 * or
 * \code
 * myFasta.getHeader(4, tempStr);
 * \endcode
 * Can set the header via
 * \code
 * myFasta.setHeader(4, tempStr);
 * \endcode
 * <b> It is vital that if headers are set, that one does the following </b>:
 * \code
 * myFasta.regenerateHeaderMaps();
 * \endcode
 * so that the <code>getIndex(string& header)</code> and other similar
 * operators work.  One can similarly get and set sequences with
 * <code>getSequence</code> and <code>setSequence</code>.
 * One may add new entries via
 * \code
 * string myHeader = "New header";
 * Sequence mySeq = "ACGGCTA";
 * myFasta.push_back(myHeader, mySeq);
 * \endcode
 * One may need to call <code>regenerateHeaderMaps</code> after the above is
 * done (or the sequence of above is done).
 *
 * <b>*** Since getSequence gives a return value of type #Sequence, all of
 * the operations associated with #Sequence (and therefore with
 * #PrimarySequence and string) apply to the result, such as size().  Please
 * refer to the documentation of string and #PrimarySequence for more
 * details. ***</b>
 *
 * One can also print out the sequence length, if desired:
 * \code
 * cout << myFasta.seqSize(tempStr) << endl;
 * \\ or
 * cout << myFasta.seqSize(4) << endl;
 * \endcode
 *
 * Can get all headers by doing
 * \code
 * StringTable myHeaders = myFasta.getHeaders();
 * \endcode
 * after which, one may wish to do <code>myFasta.clear(true)</code> to delete
 * the Fasta information.  (This is when one needs to just extract the headers.)
 * One may also write the headers to a file:
 * \code
 * Ofstream ofp("headers.txt");
 * myFasta.writeHeaders(ofp);
 * ofp.close();
 * \endcode
 *
 * The #Fasta class also has the ability to write out the sizes of each
 * sequence:
 * \code
 * Ofstream ofp("sizes.txt");
 * myFasta.writeSizes(ofp);
 * // Or myFasta.writeHeaderSizes(ofp) to write header [tab] size.
 * // This can be read in and used by CoverageReader's #ContigSizes class.
 * ofp.close();
 * \endcode
 * with the size of each sequence on one line followed by a newline in the
 * order that the sequence appears in the FASTA file.
 *
 * One can also determine which indices correspond to which names headers.
 * The mapping between indices and headers is automatically done when a FASTA
 * file is loaded, but if the headers are changed or new sequences are added,
 * one must do
 * \code
 * myFasta.regenerateHeaderMaps();
 * \endcode
 * One can then get the list of indices:
 * \code
 * IndexTable theIndices;
 * theIndices = myFasta.getIndices("Sample Header");
 * // Or, if one is sure there is only one header...
 * size_t theIndex = myFasta.getIndex("Sample Header");
 * \endcode
 * One may check to see if an index is present via
 * \code
 * bool isPresent = myFasta.hasIndex("Sample Header");
 * \endcode
 *
 * One may get a specific base (of type <code>char</code>) via contig number or
 * header and specific position (0-indexed), i.e.
 * \code
 * cout << myFasta.getBase(5, 2) << endl;
 * cout << myFasta.getBase("Sample Header", 32) << endl;
 * // Below get bases 5 though 10 inclusive from "Sample Header"
 *   // First number is the postion (0-index based) and the
 *   // second number is the length.
 * cout << myFasta.getBases("Sample Header", 4, 6) << endl;
 * // One can also pass the index instead, if one wishes.
 * \endcode
 * Similarly, one may set a specific base:
 * \code
 * myFasta.setBase(5, 2, 'C');
 * myFasta.setBase("Sample Header", 32, 'T');
 * \endcode
 *
 * Given a user-defined function of the form
 * \code
 * void f(Sequence& the Seq) {...}
 * \endcode
 * or
 * \code
 * void f(string& theHead, Sequence& the Seq) {...}
 * \endcode
 * or
 * \code
 * void f(size_t& seqNum, string& theHead, Sequence& the Seq) {...}
 * \endcode
 * Then one can apply these functions via:
 * \code
 * myFasta.applyToHeader(tempStr, f);   // Apply f to all sequences whose
 *                                      // header matches tempStr.
 * myFasta.applyToHeaderRegex(tempStr, f, true);
 *   // Apply f to all sequences whose header causes a hit in the regular
 *   // expression given in tempStr.  If false instead of true, then apply to
 *   // all sequences whose does not cause a hit.
 * // The above trailing parameter applies to all Regex methods below.
 * myFasta.applyToSequenceRegex(tempStr, f);   // Apply f to all sequences
 *                                             // that causes a hit in the
 *                                             // regular expression tempStr.
 * \endcode
 * One may need to call <code>regenerateHeaderMaps</code> after the above is
 * done if the function being applied alters the headers.
 *
 * If one wishes to just get the indices of the FASTA entries for a particular
 * regular expression applied on a header or sequence, then do:
 * \code
 * // string tempStr contains the regular expression
 * IndexTable myIndices;
 * myFasta.getIndiciesByHeaderRegex(tempStr, myIndices);
 *   // or the less-memory-and-cpu-efficient
 *   // myIndices = myFasta.getIndiciesByHeaderRegex(tempStr);
 * // Similarly,
 * myFasta.getIndiciesBySequenceRegex(tempStr, myIndices);
 *   // or the less-memory-and-cpu-efficient
 *   // myIndices = myFasta.getIndiciesBySequenceRegex(tempStr);
 * \endcode
 *
 * If one wishes to generate a mapping of kmers to contigs and positions, one
 * does (let's say for kmers of size 12):
 * \code
 * myFasta.generateKmerMappings(12);
 * // ...
 * cout << myFasta.getKmerValue() << endl;
 * map< string, Table<size_t, size_t> >& kmerMappings = myFasta.getKmerMappings();
 * // where the table has the Fasta index in the first column and the
 * // position of the kmer in the second column.
 * \endcode
 *
 * There is also a parameter, maxLoadFactor that comes after all the other
 * parameters, which is a floating point number.  (One may also set this via
 * the function <code>setMaxLoadFactor(float maxLoadFactor)</code>.)  First
 * note that the alignments are index by query name.  To get the index, an
 * unordered_map data structure is used.  (Think of it as a hash
 * table, though technically, it is not so.)  Then, the maximum_load_factor is
 * the largest number of entries in a bucket, before a new bucket is created.
 * Thus, this function is used to set the average maximum number of entries.
 * Making the number bigger makes the program go slower, but may save memory.
 * The default is 1.0.
 *
 *
 * \section ProgNote Programming note:
 * <ul>
 *    <li> We could have used the Table::applyToAllRows() functionality.
 *         However, what is coded below is more efficient since we know
 *         what columns everything is in.
 *    </li>
 *    <li> Originally, there was a function createDistinctFastaTable, but this
 *         is made obsolete/redundant by the function #distinctTable.
 *    </li>
 * </ul>
 *
 */
class Fasta : public Table<FASTA_ROW_TYPE>
{
public:
    void setMaxLoadFactor(float maxLoadFactor)
    {
        headerMappings.max_load_factor(maxLoadFactor);
        kmerMappings.max_load_factor(maxLoadFactor);
    }

    Fasta()
    {
        labels.resize(2);
        labels[0] = "Header";
        labels[1] = "Sequence";
        kValue = 0;
        inRegenHeaders = false;
        setMaxLoadFactor(1.0);
    }

    void clear()
    {
        Table<FASTA_ROW_TYPE>::clear();
        kmerMappings.clear();
        headerMappings.clear();
    }

    void regenerateHeaderMaps()
    {
        headerMappings.clear();
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
#ifdef PROG_DEBUG
            if (headerMappings.find(getHeader(i)) != headerMappings.end())
            {
                cerr << "Warning from Fasta::regenerateHeaderMaps.  Headers are not unique.  Continuing execution." << endl;
            }
#endif
            headerMappings[getHeader(i)].push_back(i);
        }
    }

    // readFullHeader may be set to true for those rare instances when the full
    // header (not just the first non-whitespace word) must be read in, so one
    // can do regex search on the headers.
    Fasta(string filename, size_t first = 0,
          size_t last = numeric_limits<size_t>::max(),
          size_t numAlloc = 0, size_t minStringReserve = 0,
          float maxLoadFactor = 1.0,
          bool readFullHeader = false,
          bool inRegenHeadersInput = true)
    {
        kValue = 0;
        labels.resize(2);
        labels[0] = "Header";
        labels[1] = "Sequence";
        inRegenHeaders = inRegenHeadersInput;
        setMaxLoadFactor(maxLoadFactor);
        readReadsTable(*this, filename, first, last, numAlloc, minStringReserve,
                       0, 33, 93, true, readFullHeader);
        if (inRegenHeaders)
        {
            regenerateHeaderMaps();
        }
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

    void writeHeaders(ostream& out)
    {
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            out << getHeader(i) << endl;
        }

    }

    void setHeader(size_t i, string& theHeader)
    {
        row_type tempRow;
        getRow(i, tempRow);
        get<0>(tempRow) = theHeader;
        setRow(i, tempRow);
    }


    /* \file Fasta.hpp
     * \fn void getSequence(size_t i, Sequence& theSequence)
     */
    void getSequence(size_t i, Sequence& theSequence)
    {
        row_type tempRow;
        getRow(i, tempRow);
        theSequence = get<1>(tempRow);
    }
    /* \file Fasta.hpp
     * \fn void getSequence(size_t i)
     */
    Sequence getSequence(size_t i)
    {
        row_type tempRow;
        getRow(i, tempRow);
        return get<1>(tempRow);
    }

    /* \file Fasta.hpp
     * \fn void getSequence(string header, Sequence& theSequence)
     */
    void getSequence(string header, Sequence& theSequence)
    {
        row_type tempRow;
        getRow(getIndex(header), tempRow);
        theSequence = get<1>(tempRow);
    }
    /* \file Fasta.hpp
     * \fn void getSequence(string header);
     */
    Sequence getSequence(string header)
    {
        row_type tempRow;
        getRow(getIndex(header), tempRow);
        return get<1>(tempRow);
    }


    void push_back(string header, Sequence theSequence)
    {
        get<0>(this->theData).push_back(header);
        get<1>(this->theData).push_back(theSequence);
    }


    // Need this since above push_back confuses compiler into thinking there
    // is only one push_back.  Really want to call base class' push_back.
    void push_back(Table<FASTA_ROW_TYPE>::row_type& theRow)
    {
        Table<FASTA_ROW_TYPE>::push_back(theRow);
    }


    char getBase(size_t contig, size_t pos)
    {
        return ((get<1>(this->theData))[contig])[pos];
    }
    char getBase(string header, size_t pos)
    {
        return ((get<1>(this->theData))[getIndex(header)])[pos];
    }


    void setBase(size_t contig, size_t pos, char newBase)
    {
        ((get<1>(this->theData))[contig])[pos] = newBase;
    }
    void setBase(string header, size_t pos, char newBase)
    {
        ((get<1>(this->theData))[getIndex(header)])[pos] = newBase;
    }


    Sequence getBases(size_t contig, size_t pos, size_t length)
    {
        return (((get<1>(this->theData))[contig]).substr(pos, length));
    }

    Sequence getBases(string header, size_t pos, size_t length)
    {
        return (((get<1>(this->theData))[getIndex(header)]).substr(pos, length));
    }


    size_t seqSize(size_t i)
    {
        return ((get<1>(this->theData))[i]).size();
    }


    void writeSizes(ostream& out)
    {
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            out << seqSize(i) << endl;
        }
    }

    void writeHeaderSizes(ostream& out)
    {
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            out << getHeader(i) << "\t" << seqSize(i) << endl;
        }
    }


    /* \file Fasta.hpp
     * \fn void setSequence(size_t i, Sequence& theSequence)
     */
    void setSequence(size_t i, Sequence& theSequence)
    {
        row_type tempRow;
        getRow(i, tempRow);
        get<1>(tempRow) = theSequence;
        setRow(i, tempRow);
    }
    /* \file Fasta.hpp
     * \fn void setSequence(string header, Sequence& theSequence)
     */
    void setSequence(string header, Sequence& theSequence)
    {
        row_type tempRow;
        size_t i = getIndex(header);
        getRow(i, tempRow);
        get<1>(tempRow) = theSequence;
        setRow(i, tempRow);
    }


    void applyToHeader(string theHeader, void f(Sequence& theSeq))
    {
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            if (get<0>(theData)[i] == theHeader)
            {
                f(get<1>(theData)[i]);
            }
        }
    }
    void applyToHeader(string theHeader, void f(string& theHead,
                       Sequence& theSeq))
    {
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            if (get<0>(theData)[i] == theHeader)
            {
                f(get<0>(theData)[i], get<1>(theData)[i]);
            }
        }
    }


    void applyToHeaderRegex(string regexStr, void f(Sequence& theSeq),
                            bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
        // Compile the regular expression stored in regexStr to
        // regex_t form stored in theRegex.
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToHeaderRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            // Apply the regular expression (stored in theRegex) on
            // header i. regexec returns 0 for a successful match (thus !reti
            // below to indicate a match).  The last three parameters are ignored.
            reti = regexec(&theRegex, get<0>(theData)[i].c_str(), 0, NULL, 0);
            if (!reti)
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(get<1>(theData)[i]);
            }
        }
    }
    void applyToHeaderRegex(string regexStr,
                            void f(string& theHead, Sequence& theSeq), 
                            bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToHeaderRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<0>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(get<0>(theData)[i], get<1>(theData)[i]);
            }
        }
    }
    void applyToHeaderRegex(string regexStr,
                            void f(size_t& i, string& theHead, Sequence& theSeq),
                            bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToHeaderRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<0>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(i, get<0>(theData)[i], get<1>(theData)[i]);
            }
        }
    }

    void applyToSequenceRegex(string regexStr,
                              void f(Sequence& theSeq), bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToSequenceRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<1>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(get<1>(theData)[i]);
            }
        }
    }
    void applyToSequenceRegex(string regexStr,
                              void f(string& theHead, Sequence& theSeq),
                              bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToSequenceRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<1>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(get<0>(theData)[i], get<1>(theData)[i]);
            }
        }
    }
    void applyToSequenceRegex(string regexStr,
                              void f(size_t& i, string& theHead, 
                                     Sequence& theSeq),
                              bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::applyToSequenceRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<1>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
                // if (regex_match(get<1>(theData)[i].c_str(), what, theRegex))
            {
                f(i, get<0>(theData)[i], get<1>(theData)[i]);
            }
        }
    }


    IndexTable& getIndices(string header)
    {
#ifdef DEBUG
        // Cannot do below if not intended to create the header mappings.
        assert(inRegenHeaders);
#endif
        return headerMappings[header];
    }


    void getIndicesByHeaderRegex(string regexStr,
                                 IndexTable& theResult, bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::getIndicesByHeaderRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<0>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
            {
                theResult.push_back(i);
            }
        }
    }

    IndexTable getIndicesByHeaderRegex(string regexStr, bool doHit = true)
    {
        IndexTable theResult;
        getIndicesByHeaderRegex(regexStr, theResult, doHit);
        return theResult;
    }

    void getIndicesBySequenceRegex(string regexStr,
                                   IndexTable& theResult, bool doHit = true)
    {
        // G++'s regex does not work, so we use the POSIX one.
        // regex theRegex(regexStr);
        // cmatch what;
        regex_t theRegex;
#ifdef PROG_DEBUG
        if (regcomp(&theRegex, regexStr.c_str(), 0))
        {
            cerr << "Error in Fasta::getIndicesBySequenceRegex!  "
                 << regexStr << " failed to compile.  Exiting." << endl;
            exit(1);
        }
#else
        regcomp(&theRegex, regexStr.c_str(), 0);
#endif
        int reti;
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            reti = regexec(&theRegex, get<1>(theData)[i].c_str(), 0, NULL, 0);
            if ( (doHit && !reti) || (!doHit && reti) )
            {
                theResult.push_back(i);
            }
        }
    }

    IndexTable getIndicesBySequenceRegex(string regexStr, bool doHit = true)
    {
        IndexTable theResult;
        getIndicesBySequenceRegex(regexStr, theResult, doHit);
        return theResult;
    }


    bool hasIndex(string header)
    {
#ifdef DEBUG
        // Cannot do below if not intended to create the header mappings.
        assert(inRegenHeaders);
#endif
        if (headerMappings.find(header) == headerMappings.end())
        {
            return false;
        }
        else
        {
            return true;
        }
    }


    size_t getIndex(string header)
    {
#ifdef DEBUG
        // Cannot do below if not intended to create the header mappings.
        assert(inRegenHeaders);
        // If using this function, should only be if the header is unique.
        if (headerMappings[header].size() != 1)
        {
            string hMTypeErr = "unique";
            if (headerMappings[header].size() == 0)
            {
                hMTypeErr = "present";
            }
            cerr << "Error in Fasta::getIndex.  Header " << header << " is not " << hMTypeErr << ".  Exiting." << endl;
            assert(true == false);
        }
#endif
        return headerMappings[header][0];
    }


    size_t seqSize(string header)
    {
        return ((get<1>(this->theData))[getIndex(header)]).size();
    }


    void generateKmerMappings(size_t theKValue)
    {
        // For each sequence in Fasta, go through the sequence using
        // theKValue window and push back the sequence index in the
        // mapping associated with the current kmer.
        kValue = theKValue;
        kmerMappings.clear();
        tuple<size_t, size_t> theRow;
        size_t numRows = size();
        Sequence currSeq;
        string currKmer;
        for (size_t row = 0; row < numRows; ++row)
        {
            get<0>(theRow) = row;
            currSeq.clear();
            getSequence(row, currSeq);
            size_t numPos = row - kValue + 1;
            for (size_t pos = 0; pos < numPos; ++pos)
            {
                currKmer = currSeq.subsequence(pos, kValue);
                get<1>(theRow) = pos;
                kmerMappings[currKmer].push_back(theRow);
            }
        }
    }

    size_t getKmerValue()
    {
        return kValue;
    }

    unordered_map< string, Table<size_t, size_t> >& getKmerMappings()
    {
        return kmerMappings;
    }


    void read(string filename,
              size_t first = 0, size_t last = numeric_limits<size_t>::max(),
              size_t numAlloc = 0, size_t minStringReserve = 0)
    {
        readReadsTable(*this, filename, first, last, numAlloc, minStringReserve);
        regenerateHeaderMaps();
    }

    void write(string filename)
    {
        writeFastaTable(*this, filename);
    }


protected:
    unordered_map<string, IndexTable> headerMappings;
    size_t kValue;
    unordered_map< string, Table<size_t, size_t> > kmerMappings;

private:
    bool inRegenHeaders;

};

#endif

