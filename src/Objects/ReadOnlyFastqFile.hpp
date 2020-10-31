#ifndef MOLBIOLIB_READONLYFASTQFILE_H
#define MOLBIOLIB_READONLYFASTQFILE_H 1

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
#include "src/Objects/Sequence.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"


/** \file ReadOnlySequentialFastqFile.hpp
 * Contains only the #ReadOnlyFastqFile class.
 */

/** Object that read-only fastq files with a vector interface.
 * Usage is similar to that of the #Fasta object.
 * #ReadOnlyTSVFile.  The #Fasta methods that do not write to the sequence or
 * headers are available <b>in random access mode</b>.  Most of these methods
 * are not available in sequential mode.  Typical usage is:
 * \code
 * string FILE = ...;
 * ReadOnlyFastqFile myReads(FILE);
 * size_t numReads = myReads.size();
 * // May access randomly, though getNextRead() is available too.
 * for (size_t i = 0; i < numReads; ++i)
 * {
 *     // Header and Qual are of type string
 *     // while Sequence is of type Sequence.
 *     // All of the below with an argument perform by calling
 *     // myReads.getRead(i) first and then returning currentSeq, currentHeader,
 *     // currentQual, or currentQualHeader, all public variables of myReads.
 *     // headerChar() returns a _string_ of either > or @ depending on
 *     // whether the file is a fasta or fastq.
 *     cout << myReads.headerChar() << myReads.getHeader(i) << endl;
 *     cout << myReads.getSequence(i) << endl;
 *     // There is also getQualHeader(i) and qualHeaderChar()
 *     cout << myReads.getQual(i) << endl;
 * }
 * // Can do below if header maps (default true) are generated.
 * size_t numSeqWithHeader = myReads.size("chr2");
 * // And similarly with other operations above.
 * for (size_t i = 0; i < numSeqWithHeader; ++i)
 * {
 *     // If know unique, can omit i below and will always get first one.
 *     // Like above, below does myReads.getRead("chr2") first and similarly to
 *     // above.
 *     Sequence tempSeq = myReads.getSequence("chr2", i);
 *     // Process tempSeq here...
 * }
 * myReads.close();
 * \endcode
 *
 * There is also a sequential mode that uses almost no memory, unlike the above:
 * \code
 * ReadOnlyFastqFile myReads(FILE, true);
 * // Header maps are never generated in sequential mode.
 * // myReads will initially read the first read in automatically.
 * // myReads.size() and myReads.get...(...) are not available in this mode.
 * bool goodRead = myReads.hasReads();
 * while (goodRead)
 * {
 *     // In case of a continue, one must get the next read, e.g.
 *     if (myReads.getCurrentSequence() == "")
 *     {
 *         // Skip current blank sequence
 *         goodRead = myReads.getNextRead();
 *         continue;
 *     }
 *     // Notice how all of the gets do not have an argument!
 *     cout << myReads.getCurrentHeader() << endl;
 *     cout << myReads.getCurrentSequence() << endl;
 *     cout << myReads.getCurrentQual() << endl;
 *     // There is also getCurrentQualHeader(), but this is probably
 *     // not so useful, as many times, the qual header is blank.
 *     goodRead = myReads.getNextRead();
 * }
 * myReads.close();
 * \endcode
 *
 */
class ReadOnlyFastqFile
{
public:

    string headerChar()
    {
        return "@";
    }

    string qualHeaderChar()
    {
        return "+";
    }

    bool hasReads()
    {
#ifdef PROG_DEBUG
        assert(ifp.isOpen() || ifpstream.is_open());
#endif
        return !isEmpty;
    }

    void getRead(size_t readNum)
    {
#ifdef DEBUG
        assert(ifp.isOpen() && !isEmpty && !sequentialMode);
#endif
#ifdef PROG_DEBUG
        assert(readNum < numReads);
#endif
        if (readNum != currentRead)
        {
            currentRead = readNum;
            string tempLine = "";
            tuple<size_t, size_t, size_t, size_t, size_t, size_t> tempRow;
            fastqIndex.getRow(readNum, tempRow);
            ifp[get<0>(tempRow)];
            currentHeader = ifp.line.substr(1);
            size_t start = get<1>(tempRow),
                   end = get<1>(tempRow) + get<2>(tempRow);
            currentSeq = "";
            for (size_t i = start; i < end; ++i)
            {
                currentSeq += ifp[i];
            }
            ifp[get<3>(tempRow)];
            currentQualHeader = ifp.line.substr(1);
            start = get<4>(tempRow),
            end = get<4>(tempRow) + get<5>(tempRow);
            currentQual = "";
            for (size_t i = start; i < end; ++i)
            {
                currentQual += ifp[i];
            }
        }
    }

    size_t size(string header = "")
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        if (header == "")
        {
            return numReads;
        }
        else
        {
#ifdef PROG_DEBUG
            assert(generateHeaderMap);
#endif
            return headerMapping[header].size();
        }
    }

    void getRead(string header, size_t i = 0)
    {
#ifdef PROG_DEBUG
        assert(generateHeaderMap);
#endif
#ifdef DEBUG
        assert(!sequentialMode &&
               (headerMapping.find(header) != headerMapping.end()) &&
               (headerMapping[header].size() > 0));
#endif
        getRead(headerMapping[header][i]);
    }


    IndexTable& getIndices(string header)
    {
#ifdef PROG_DEBUG
        // Cannot do below if not intended to create the header mappings.
        assert(generateHeaderMap);
#endif
        return headerMapping[header];
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
            cerr << "Error in ReadOnlySequencesFile::getIndicesByHeaderRegex!  "
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
            getRead(i);
            string tempStr(currentHeader);
            reti = regexec(&theRegex, tempStr.c_str(), 0, NULL, 0);
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
            cerr << "Error in ReadOnlySequencesFile::getIndicesBySequenceRegex!  "
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
            getRead(i);
            string tempStr(currentSeq);
            reti = regexec(&theRegex, tempStr.c_str(), 0, NULL, 0);
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
        assert(generateHeaderMap);
#endif
        if (headerMapping.find(header) == headerMapping.end())
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
        assert(generateHeaderMap);
        // If using this function, should only be if the header is unique.
        if (headerMapping[header].size() != 1)
        {
            string hMTypeErr = "unique";
            if (headerMapping[header].size() == 0)
            {
                hMTypeErr = "present";
            }
            cerr << "Error in ReadOnlySequencesFile::getIndex.  Header " << header << " is not " << hMTypeErr << ".  Exiting." << endl;
            assert(true == false);
        }
#endif
        return headerMapping[header][0];
    }


    // Gets the sequence.
    Sequence& operator[](size_t readNum)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(readNum);
        return currentSeq;
    }


    // Gets the sequence.
    Sequence& operator[](string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentSeq;
    }


    string& getHeader(size_t readNum)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(readNum);
        return currentHeader;
    }
    void getHeader(size_t readNum, string& theHeader)
    {
        theHeader = getHeader(readNum);
    }

    string& getHeader(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentHeader;
    }
    void getHeader(string header, string& theHeader)
    {
        theHeader = getHeader(header);
    }

    StringTable getHeaders()
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
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
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        size_t numRows = size();
        for (size_t i = 0; i < numRows; ++i)
        {
            out << getHeader(i) << endl;
        }

    }

    Sequence& getSequence(size_t readNum)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(readNum);
        return currentSeq;
    }
    void getSequence(size_t readNum, Sequence& theSequence)
    {
        theSequence = getSequence(readNum);
    }

    Sequence& getSequence(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentSeq;
    }
    void getSequence(string header, Sequence& theSequence)
    {
        theSequence = getSequence(header);
    }

    char getBase(size_t contig, size_t pos)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(contig);
        return currentSeq[pos];
    }
    char getBase(string header, size_t pos)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentSeq[pos];
    }

    Sequence getBases(size_t contig, size_t pos, size_t length)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(contig);
        return currentSeq.substr(pos, length);
    }

    Sequence getBases(string header, size_t pos, size_t length)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentSeq.substr(pos, length);
    }

    size_t seqSize(size_t contig)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(contig);
        return currentSeq.size();
    }

    size_t seqSize(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentSeq.size();
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



    void generateKmerMappings(size_t theKValue)
    {
        // For each sequence in the file, go through the sequence using
        // theKValue window and push back the sequence index in the
        // mapping associated with the current kmer.
        kValue = theKValue;
        kmerMappings.clear();
        tuple<size_t, size_t> theRow;
        size_t numRows = size();
        string currKmer;
        for (size_t row = 0; row < numRows; ++row)
        {
            get<0>(theRow) = row;
            getRead(row);
            size_t numPos = row - kValue + 1;
            for (size_t pos = 0; pos < numPos; ++pos)
            {
                currKmer = currentSeq.subsequence(pos, kValue);
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



    size_t qualSize(size_t contig)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(contig);
        return currentQual.size();
    }

    size_t qualSize(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        return currentQual.size();
    }

    string& getQualHeader(size_t readNum)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(readNum);
        return currentQualHeader;
    }
    void getQualHeader(size_t readNum, string& theHeader)
    {
        theHeader = getQualHeader(readNum);
    }

    string& getQualHeader(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentQualHeader;
    }
    void getQualHeader(string header, string& theHeader)
    {
        theHeader = getQualHeader(header);
    }


    string& getQual(size_t readNum)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(readNum);
        return currentQual;
    }

    string& getQual(string header)
    {
#ifdef DEBUG
        assert(!sequentialMode);
#endif
        getRead(getIndex(header));
        return currentQual;
    }
    void getQual(string header, string& theQuality)
    {
        theQuality = getQual(header);
    }

    char getQualChar(size_t contig, size_t pos)
    {
        getRead(contig);
        return currentQual[pos];
    }
    char getQualChar(string header, size_t pos)
    {
        getRead(getIndex(header));
        return currentQual[pos];
    }

    int getQualScore(size_t contig, size_t pos, int offset = 0)
    {
        char theQual = getQualChar(contig, pos);
        int theQualScore = static_cast<int>(theQual) - offset;
        if (theQualScore > qualMax)
        {
            theQualScore = qualMax;
        }
        return theQualScore;
    }

    int getQualScore(string header, size_t pos, int offset = 0)
    {
        char theQual = getQualChar(header, pos);
        int theQualScore = static_cast<int>(theQual) - offset;
        if (theQualScore > qualMax)
        {
            theQualScore = qualMax;
        }
        return theQualScore;
    }

    Table<int> getQualScoreTable(size_t contig, int offset = 0)
    {
        size_t length = currentQual.size();
        Table<int> result(length, 0);
        for (size_t i = 0; i < length; ++i)
        {
            result[i] = getQualScore(contig, i) - offset;
        }
        return result;
    }

    Table<int> getQualScoreTable(string header, int offset = 0)
    {
        return getQualScoreTable(getIndex(header), offset);
    }


XXXX GOT TO HERE XXXX
    bool getNextRead()
    {
        if (!sequentialMode)
        {
            if (currentRead >= (numReads-1))
            {
                return false;
            }
            ++currentRead;
            getRead(currentRead);
        }
        else
        {
#ifdef DEBUG
            assert(ifpstream.is_open() && !isEmpty);
#endif
            string line = "";

            // Get header = next non-blank line.
            getline(ifpstream, line);
            while (!ifpstream.fail() && (trimSpacesString(line) == ""))
            {
                if (ifpstream.fail())
                {
                    return false;
                }
                getline(ifpstream, line);
            }
            if (ifpstream.fail())
            {
                return false;
            }
#ifdef DEBUG
            assert(line[0] == headerIndicator);
#endif
            currentHeader = line.substr(1);

            // Now get sequence - continue until either fail or nextIndicator.
            getline(ifpstream, line);
            while (!ifpstream.fail() && (trimSpacesString(line) == ""))
            {
                if (ifpstream.fail())
                {
                    return false;
                }
                getline(ifpstream, line);
            }
            if (ifpstream.fail())
            {
                return false;
            }
            currentSeq = trimSpacesString(line);
#ifdef DEBUG
            assert(currentSeq[0] != nextIndicator);
#endif
            char c = static_cast<char>(ifpstream.peek());
            while ((c != nextIndicator) && !ifpstream.fail())
            {
                getline(ifpstream, line);
                if (!ifpstream.fail())
                {
                    currentSeq += trimSpacesString(line);
                    c = static_cast<char>(ifpstream.peek());
                }
            }
            if (!fastqMode)
            {
                return true;
            }
            else
            {
                // fastq - so need qual header and quals
                if (ifpstream.fail())
                {
                    return false;
                }
                // From above, should be at qual header (nextIndicator)
                getline(ifpstream, line);
                if (ifpstream.fail())
                {
                    return false;
                }
                currentQualHeader = line.substr(1);
                // Now get quals - continue until either fail or nextIndicator.
                getline(ifpstream, line);
                while (!ifpstream.fail() && (trimSpacesString(line) == ""))
                {
                    if (ifpstream.fail())
                    {
                        return false;
                    }
                    getline(ifpstream, line);
                }
                if (ifpstream.fail())
                {
                    return false;
                }

                currentQual = trimSpacesString(line);
                char c = static_cast<char>(ifpstream.peek());
                while ((c != headerIndicator) && !ifpstream.fail())
                {
                    getline(ifpstream, line);
                    if (!ifpstream.fail())
                    {
                        currentQual += trimSpacesString(line);
                        c = static_cast<char>(ifpstream.peek());
                    }
                }
            }
        }
        return true;
    }



    void open(string filename,
              bool sequentialModeInput = false,
              bool generateHeaderMapInput = false,
              string indexExtInput = "")
    {
        sequentialMode = sequentialModeInput;
        generateHeaderMap = generateHeaderMapInput;
        isEmpty = false;
        currentRead = numeric_limits<size_t>::max();
        // Need random access version of ReadOnlyStringFile.

        // Determine if this is a Fasta or Fastq file.
        fastqMode = true;
        // Do this using a stream.
        Ifstream fastqTestIfp(filename);
        bool goodFlag = !fastqTestIfp.fail();
        if (goodFlag)
        {
            string tempLine = "";
            getline(fastqTestIfp, tempLine);
            if (fastqTestIfp.fail())
            {
                goodFlag = false;
            }
            else
            {
                while (goodFlag && (trimSpacesString(tempLine) == ""))
                {
                    getline(fastqTestIfp, tempLine);
                    if (fastqTestIfp.fail())
                    {
                        goodFlag = false;
                    }
                }
                if (goodFlag)
                {
                    if (tempLine[0] == '>')
                    {
                        // This is a Fasta file.
                        fastqMode = false;
                    }
#ifdef DEBUG
                    else
                    {
                        // Make sure this is really a Fastq file.
                        if (tempLine[0] != '@')
                        {
                            cerr << "Error in ReadOnlySequencesFile!  " << filename
                                 << " is not a valid Fasta or Fastq file.  Exiting."
                                 << endl;
                            exit(1);
                        }
                    }
#endif

                }
            }
        }
        fastqTestIfp.close();
        if (!goodFlag)
        {
            isEmpty = true;
            numReads = 0;
#ifdef DEBUG
            cerr << "Warning from ReadOnlySequencesFile!  " << filename
                 << " is empty.  Continuing execution." << endl;
#endif
        }

        if (fastqMode)
        {
            headerIndicator = '@';
            nextIndicator = '+';
        }
        else
        {
            headerIndicator = '>';
            nextIndicator = '>';
        }


        if (!sequentialMode)
        {
            // Open in random access mode
            ifp.openResetParams(filename);
            size_t numLines = ifp.size();

            setReadsIndexFileExtension(indexExtInput);
            if (!isEmpty)
            {
                numReads = 0;
                string readsIndexFilename = filename + "." +
                                            readsIndexExtension;
                if ((fileSize(readsIndexFilename) ==
                        static_cast<ifstream::pos_type>(-1)) ||
                        isFileOlderThan(readsIndexFilename, filename))
                {
                    if (fileSize(readsIndexFilename) !=
                            static_cast<ifstream::pos_type>(-1))
                    {
                        // If got to here, then index file exists and
                        // is older than file.
                        int removeStat = remove(readsIndexFilename);
                        if (removeStat != 0)
                        {
                            cerr << "Warning in ReadOnlySequencesFile::open!  removeStat is non-zero.  Something went wrong with the file deletion.  Continuing execution." << endl;
                        }
                    }
                    // Generate index
                    string currLine = "";
                    size_t i = findNextNonblankLineInIfp(0, goodFlag);
                    while (goodFlag && (i < numLines))
                    {
                        // Currently i should point to the header
#ifdef DEBUG
                        // Checking for expected line.
                        assert(goodFlag);
#endif
                        currLine = ifp[i];
#ifdef DEBUG
                        assert(currLine[0] == headerIndicator);
#endif
                        // The header is 1 line long.
                        string header = currLine.substr(1);
                        size_t headerPos = i;

                        // Next get sequence.
                        ++i;
                        i = findNextNonblankLineInIfp(i, goodFlag);
#ifdef DEBUG
                        // Checking for expected line.
                        assert(goodFlag);
#endif
                        string currSeq = ifp[i];
                        size_t seqPos = i, seqLength = 1;
#ifdef DEBUG
                        // We assume that the next line is the [start of a]
                        // sequence as the Fasta/Fastq doesn't permit
                        // blank entries.  Below checks if it on to
                        // the next thing, which is bad.
                        assert(currSeq[0] != headerIndicator);
                        if (fastqMode)
                        {
                            assert(currSeq[0] != '+');
                        }
#endif
                        ++i;
                        while ((i < numLines) && (ifp[i][0] != nextIndicator))
                        {
                            ++seqLength;
                            ++i;
                        }

                        // Now get quals.
                        size_t qualHeaderPos = i,
                               qualPos = numeric_limits<size_t>::max(),
                               qualLength = 0;
                        if (fastqMode)
                        {
                            ++i;
                            qualPos = findNextNonblankLineInIfp(i, goodFlag);
#ifdef DEBUG
                            // Check we have the quals.
                            // Note that we do NOT check if this is the
                            // next header since it is possible that the first
                            // character of the qual is an ampersand.
                            assert(goodFlag);
#endif
                            qualLength = 1;
                            ++i;
                            while ((i < numLines) && (ifp[i][0] != '@'))
                            {
                                ++qualLength;
                                ++i;
                            }
                            // Record this in table:
                            fastqIndex.push_back(headerPos, seqPos, seqLength, qualHeaderPos, qualPos, qualLength);
                            if (generateHeaderMap)
                            {
                                headerMapping[header].push_back(fastqIndex.size() - 1);
                            }
                        }
                        else
                        {
                            // fasta mode
                            fastaIndex.push_back(headerPos, seqPos, seqLength);
                            if (generateHeaderMap)
                            {
                                headerMapping[header].push_back(fastaIndex.size() - 1);
                            }
                        }
                    }
                    // Write index
                    Ofstream ofp(readsIndexFilename);
                    if (fastqMode)
                    {
                        numReads = fastqIndex.size();
                        tuple<size_t, size_t, size_t, size_t, size_t, size_t> tempRow;
                        for (size_t i = 0; i < numReads; ++i)
                        {
                            fastqIndex.getRow(i, tempRow);
                            ofp << get<0>(tempRow) << "\t"
                                << get<1>(tempRow) << "\t"
                                << get<2>(tempRow) << "\t"
                                << get<3>(tempRow) << "\t"
                                << get<4>(tempRow) << "\t"
                                << get<5>(tempRow) << endl;
                        }
                    }
                    else
                    {
                        numReads = fastaIndex.size();
                        tuple<size_t, size_t, size_t> tempRow;
                        for (size_t i = 0; i < numReads; ++i)
                        {
                            fastaIndex.getRow(i, tempRow);
                            ofp << get<0>(tempRow) << "\t"
                                << get<1>(tempRow) << "\t"
                                << get<2>(tempRow) << endl;
                        }
                    }
                    // Below will differentiate the Table rows
                    // from the rows storing the header mappings.
                    // This has seven columns, not four or six like above.
                    ofp << "\t\t\t\t\t\t" << endl;
                    if (generateHeaderMap)
                    {
                        for (unordered_map<string, IndexTable>::iterator i = headerMapping.begin(); i != headerMapping.end(); ++i)
                        {
                            ofp << i->first << endl;
                            IndexTable& currTable = i->second;
                            size_t numCurrTableRows = currTable.size();
#ifdef DEBUG
                            if (numCurrTableRows > 1)
                            {
                                cerr << "Warning from ReadOnlySequencesFile!  In "
                                     << filename << " read " << i->first
                                     << " is not a unique header.  Operations "
                                     << "with a header argument will not work. "
                                     << "Continuing execution." << endl;
                            }
#endif
                            ofp << numCurrTableRows << endl;
                            for (size_t j = 0; j < numCurrTableRows; ++j)
                            {
                                ofp << currTable[j] << endl;
                            }
                        }
                        ofp << "\t\t\t\t\t\t" << endl;
                    }
                    ofp.close();
                }
                else
                {
                    // Read index
                    ReadOnlyTSVFile indexFp(readsIndexFilename, true, true);
                    size_t i = 0;
                    vector<string> indexTokens = indexFp[0];
                    while (indexTokens.size() != 7)
                    {
                        ++numReads;
                        if (fastqMode)
                        {
                            tuple<size_t, size_t, size_t, size_t, size_t, size_t> tempRow;
                            tuplizeTokens(indexTokens, tempRow);
                            fastqIndex.push_back(tempRow);
                        }
                        else
                        {
                            tuple<size_t, size_t, size_t> tempRow;
                            tuplizeTokens(indexTokens, tempRow);
                            fastaIndex.push_back(tempRow);
                        }
                        ++i;
                        indexTokens = indexFp[i];
                    }
                    if (generateHeaderMap)
                    {
                        ++i;
                        indexTokens = indexFp[i];
#ifdef DEBUG
                        bool foundHeaderLines = false;
#endif
                        while (indexTokens.size() != 7)
                        {
#ifdef DEBUG
                            foundHeaderLines = true;
#endif
                            string inputHeader = indexFp.line;
                            ++i;
                            indexTokens = indexFp[i];
                            size_t numHeaderLines = convertFromString<size_t>(indexTokens[0]);
#ifdef DEBUG
                            if (numHeaderLines > 1)
                            {
                                cerr << "Warning from ReadOnlySequencesFile!  In "
                                     << filename << " read " << inputHeader
                                     << " is not a unique header.  Operations "
                                     << "with a header argument will not work. "
                                     << "Continuing execution." << endl;
                            }
#endif
                            for (size_t j = 0; j < numHeaderLines; ++j)
                            {
                                ++i;
                                indexTokens = indexFp[i];
                                headerMapping[inputHeader].push_back(convertFromString<size_t>(indexTokens[0]));
                            }
                            ++i;
                            indexTokens = indexFp[i];
                        }
#ifdef DEBUG
                        if (!foundHeaderLines)
                        {
                            cerr << "Warning from ReadOnlySequencesFile!  Did not find any header map lines.  Perhaps the index must be deleted and regenerated?  This may be if the index is originally a non-header-map index.  Continuing execution." << endl;
                        }
#endif
                    }
                    indexFp.close();
                }
                getRead(0);
            }
        }
        else
        {
            // Open in sequential mode
            ifpstream.open(filename);
            generateHeaderMap = false;
            if (goodFlag)
            {
                getNextRead();
            }
        }
    }

    ReadOnlySequencesFile()
    {
        qualOffset = 0;
        qualMax = 255 - qualOffset;
        qualDefault = 190;
#ifdef PROG_DEBUG
        assert((qualOffset + qualDefault) <= qualMax);
#endif
        qualDefaultChar = static_cast<char>(qualDefault + qualOffset);
        currentSeq = "";
        currentHeader = "";
        currentQualHeader = "";
        currentQual = "";
    }
    ReadOnlySequencesFile(string filename,
                          bool sequentialModeInput = false,
                          bool generateHeaderMapInput = true,
                          int qualOffsetInput = 0,
                          int qualMaxInput = 255,
                          int qualDefaultInput = 190,
                          string indexExtInput = "")
    {
        currentSeq = "";
        currentHeader = "";
        currentQualHeader = "";
        currentQual = "";
        qualOffset = qualOffsetInput;
        qualMax = qualMaxInput;
        if (qualMax - qualOffset > 255)
        {
#ifdef DEBUG
            cerr << "In ReadOnlySequencesFile, setting qualMax = 255 - qualOffset." << endl;
#endif
            qualMax = 255 - qualOffset;
        }
        qualDefault = qualDefaultInput;
        if (qualDefault > qualMax)
        {
#ifdef DEBUG
            cerr << "In ReadOnlySequencesFile, setting qualDefault = qualMax - qualOffset." << endl;
#endif
        }
        qualDefaultChar = static_cast<char>(qualDefault + qualOffset);
        open(filename,
             sequentialModeInput, generateHeaderMapInput,
             indexExtInput);
    }

    void close()
    {
        if (!sequentialMode)
        {
            ifp.close();
        }
        else
        {
            ifpstream.close();
        }
        numReads = 0;
    }

    // Alias close as clear too.
    void clear()
    {
        close();
    }

    ~ReadOnlySequencesFile()
    {
        close();
    }


    void setReadsIndexFileExtension(string theExt = "")
    {
        if (theExt == "")
        {
            string headMapStr = "";
            if (generateHeaderMap)
            {
                headMapStr = ".generateHeaderMap";
            }
            readsIndexExtension = "MolBioLib.ReadOnlySequencesFile" +
                                  headMapStr + ".index";
        }
        else
        {
            readsIndexExtension = theExt;
        }
    }


    class iterator
    {
    public:
        iterator(ReadOnlySequencesFile* inputRORF, size_t startIndex) :
            theReadOnlySequencesFile(inputRORF), currentIndex(startIndex) { }
        ~iterator() { }

        iterator& operator=(const iterator& other)
        {
            theReadOnlySequencesFile = other.theReadOnlySequencesFile;
            currentIndex = other.currentIndex;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            return ((theReadOnlySequencesFile ==
                     other.theReadOnlySequencesFile) &&
                    (currentIndex == other.currentIndex));
        }

        bool operator!=(const iterator& other)
        {
            return ((theReadOnlySequencesFile !=
                     other.theReadOnlySequencesFile) ||
                    (currentIndex != other.currentIndex));
        }

        iterator& operator++()
        {
            if (currentIndex != theReadOnlySequencesFile->size())
            {
                ++currentIndex;
            }
            return (*this);
        }

        iterator operator++(int)
        {
            iterator tmp(*this);
            ++(*this);
            return tmp;
        }

        iterator& operator--()
        {
            if (currentIndex != 0)
            {
                --currentIndex;
            }
            return (*this);
        }

        iterator operator--(int)
        {
            iterator tmp(*this);
            --(*this);
            return tmp;
        }

        // operator*() not defined.  Does not make sense in this context.

    private:
        ReadOnlySequencesFile* theReadOnlySequencesFile;
        size_t currentIndex;
    };

    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
        return iterator(this, this->size());
    }

    // Below is because STL convention is lower case while in the framework,
    // all classes start with an upper case.
    typedef iterator Iterator;


    size_t& getCurrentRead()
    {
        return currentRead;
    }

    string& getCurrentHeader()
    {
        return currentHeader;
    }

    Sequence& getCurrentSequence()
    {
        return currentSeq;
    }

    string& getCurrentQualHeader()
    {
        return currentQualHeader;
    }

    string& getCurrentQual()
    {
        return currentQual;
    }



protected:
    size_t numReads;
    bool generateHeaderMap, sequentialMode;
    string readsIndexExtension;

    size_t currentRead;
    string currentHeader, currentQualHeader, currentQual;
    Sequence currentSeq;

    unordered_map<string, IndexTable> headerMapping;
    // header pos, seq pos, seq # lines
    Table<size_t, size_t, size_t> fastaIndex;
    // header pos, seq pos, seq # lines, qual header pos, qual pos, qual # lines
    Table<size_t, size_t, size_t, size_t, size_t, size_t> fastqIndex;

    size_t kValue;
    unordered_map< string, Table<size_t, size_t> > kmerMappings;

    int qualOffset, qualMax, qualDefault;
    char qualDefaultChar;


private:
    Ifstream ifpstream;
    ReadOnlyStringFile ifp;

    bool isEmpty;
    char headerIndicator, nextIndicator;

    size_t findNextNonblankLineInIfp(size_t startLine, bool& goodFlag)
    {
        goodFlag = true;
        if (isEmpty || (startLine >= ifp.size()))
        {
            goodFlag = false;
            return numeric_limits<size_t>::max();
        }
        else
        {
            for (size_t i = startLine; i < ifp.size(); ++i)
            {
                if (trimSpacesString(ifp[i]) != "")
                {
                    return i;
                }
            }
            // Got here means got to the end of the file.
            goodFlag = false;
            return numeric_limits<size_t>::max();
        }
    }

};


#endif

