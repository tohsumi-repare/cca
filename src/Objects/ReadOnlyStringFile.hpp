#ifndef MOLBIOLIB_READONLYSTRINGFILE_H
#define MOLBIOLIB_READONLYSTRINGFILE_H 1

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
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/System.hpp"
// Below for remove
#include "src/Functions/SystemUtilities/Cstdio.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"


/** \file ReadOnlyStringFile.hpp
 * Contains only the #ReadOnlyStringFile class.
 */

/** Object that read-only files with a vector interface.  Does not have to be a
 * TSV file.  Any file with strings in it can be read using this class.
 * \section Usage Usage:
 * \code
 * ReadOnlyStringFile myROTable("myFile.txt", false, false, false, "index", 100000);
 * \endcode
 * or alternatively
 * \code
 * ReadOnlyStringFile myROTable;
 * myROTable.reserve(100000);
 * myROTable.setSequentialMode(false, false);
 * myROTable.open("myFile.txt", false, "index");
 * \endcode
 * An index file, myFile.txt.index is created if one is not found.  If the
 * index extension is not given, then MolBioLib.ReadOnlyStringFile.index, is
 * used.  The first <code>false</code>
 * means that the file may be access in random.  If <code>true</code>, then
 * other than by <code>resetToBegin</code>, access of a line is strictly in
 * increasing order of index, typically sequentially.  However, this means no
 * index is generated, but rather a file to indicate the number of lines
 * (including the header) and maximum length of a line.  This saves memory and
 * time, unless access is done in some other order.  Default
 * is <code>false</code>.  The second <code>false</code> means that the
 * size() method is used.  Otherwise, the maximum length of a line nor the
 * number of lines in the file is not computed.  Typically, then the stopping
 * condition is <code>!myROTable.fail()</code> in the for-loop.  One does this
 * when one knows one will go through the input file only once sequentially.
 * The third <code>false</code> means that if there
 * is no index file, the first line of myFile.txt does <b>not</b> have
 * a header line.  Default is <code>false</code>.
 * The <code>100000</code> option indicates
 * that a <code>reserve</code> will be called on the string variable
 * that holds the current line with a size of 100000.  This is important for
 * performance when creating indexes of files with very long lines.  The
 * index file also has the maximum line length and is set appropriately.
 * \code
 * myROTable.close();  // Cannot access data after this, but then can open
 *                     // a new file.  The index array is cleared too.
 * \endcode
 * The size [number of rows in the file] can
 * be obtained as in the below example.
 * \code
 * size_t fileSize = myROTable.size();
 *    // Alternatively, one can use SizeType instead of size() above.
 * cout << myROTable.size() << endl;   // Prints the number of rows in the file.
 * \endcode
 * Suppose we have
 * \code
 * string myRow;
 * \endcode
 * Then we can get a row of the file by doing one of the three below:
 * \code
 * myROTable.getRow(3, myRow);    // Rows indices start from 0
 *                                // so this is the forth line.
 * \endcode
 * or
 * \code
 * myRow = myROTable.getRow(3);
 * \endcode
 * or
 * \code
 * myRow = myROTable[3];
 * \endcode
 * If one wishes to forgo the implicit copy from line to myRow, one could work
 * on the input line directly by first reading in the row then acting on the
 * variable line:
 * \code
 * myROTable.readRow(3);
 * cout << myROTable.line.size() << endl;
 * \endcode
 * Note that it is really reading a line, but since subsequently derived 
 * classes have the notion of a row, it was decided to make this base class
 * method name the same as the derived.  One may iterator through the file via
 * \code
 * for (ReadOnlyStringFile::iterator i = myROTable.begin(); i != myROTable.end(); ++i) {
 *    cout << *i << endl;
 * }
 * \endcode
 * <code>Iterator</code> may be used in place of <code>iterator</code>.  The
 * iterator may only be used when the size is computed.
 *
 * If there are labels, then can get them via
 * \code
 *    cout << myROTable.getLabels() << endl;
 * \endcode
 *
 * For developement, the following function returns the index extension strings.
 * \code
 * string indexExt = "", seqIndexExt = "";
 * // The below places the file index extension into indexExt and the
 * // sequential index extension into seqIndexExt.
 * myROTable.getIndexExtensions(indexExt, seqIndexExt);
 * \endcode
 *
 * When DEBUG is defined, this will check for DOS/Mac type carriage return in
 * in the file and flag a warning.  Typically, a \\r will produce odd results
 * either in reading in a file and/or printing.  These can be removed using
 * tools such as dos2unix.
 */
class ReadOnlyStringFile
{
public:
    typedef size_t size_type;
    string line;

    void setSequentialMode(bool theMode, bool theNoSize = false)
    {
#ifdef PROG_DEBUG
        // Can only have theNoSize be true if in sequential mode.
        assert(theMode || !theNoSize);
#endif
        sequentialMode = theMode;
        noSize = theNoSize;
    }

    // Constructor
    ReadOnlyStringFile()
    {
        resetParams();
    }
    ReadOnlyStringFile(string filename,
                       bool sequentialMode = false,
                       bool noSize = false,
                       bool assumeExistLabels = false,
                       string fileExtension = "",
                       size_t resLength = 0)
    {
        // resLength only matters when initially creating the index file.  Once
        // that is written, the optimal resLength is retrieved from the
        // index file.
        openResetParams(filename, sequentialMode, noSize, assumeExistLabels, fileExtension, resLength);
    }


    void reserve(size_t resLength)
    {
        line.reserve(resLength);
    }


    virtual void resetToBegin()
    {
#ifdef DEBUG
        assert(fileOpen == true);
#endif
        atEnd = false;
        currentLine = 0;
        fin.clear();   // Resets EOF pointer.
        fin.seekg(0, ios_base::end);   // Go to end.
        endPos = fin.tellg();          // Determine end position.
        fin.clear();   // Resets EOF pointer.
        fin.seekg(0, ios_base::beg);   // Go to beginning of file.
        if (hasLabels)
        {
            getline(fin, labels);
        }
        if (sequentialMode)
        {
            if (fin.fail() || (endPos == static_cast<ifstream::pos_type>(0)))
            {
                atEnd = true;
                line = "";
            }
            if (!noSize)
            {
                if (size() > 0)
                {
                    getline(fin, line);
                }
                else
                {
                    line = "";
                }
            }
            else
            {
                // Special case for noSize.  Must allow for first line read in
                // getLine else fail() will flag true.
                notReadFirst = true;
            }
        }
    }


    virtual void open(string filename,
                      bool assumeExistLabels = false,
                      string fileExtension = "",
                      size_t reserveLength = 0)
    {
        hasLabels = assumeExistLabels;

        if (reserveLength > 0)
        {
            reserve(reserveLength);
        }

        ifstream::pos_type theFileSize = fileSize(filename);
        // Check to see if the file exists.
        assert(fileOpen == false);
        fileOpen = true;
        if (theFileSize == static_cast<ifstream::pos_type>(-1))
        {
            cerr << "Error!  The file " << filename << " cannot be opened." << endl;
            assert(theFileSize != static_cast<ifstream::pos_type>(-1));
        }

        if (!sequentialMode)
        {
            size_t maxLength = 0;
            string indexFilename;
            numEntries = 0;
            setIndexFileExtension(fileExtension);
            indexFilename = filename + "." + indexExtension;
            ifstream::pos_type theIndexSize = fileSize(indexFilename);
            findex.open(indexFilename);
            if (!findex || isFileOlderThan(indexFilename, filename))
            {
                // Index not found or is older than original file, so create.
                if (findex)
                {
                    // If got in here, then index file exists and is older than file.
                    findex.close();
                    int removeStat = remove(indexFilename);
                    if (removeStat != 0)
                    {
                        cerr << "Warning in ReadOnlyStringFile::open!  removeStat is non-zero.  Something went wrong with the file deletion.  Continuing execution." << endl;
                    }
                }
                fin.open(filename);
                Ofstream findexOut(indexFilename);
                // The idea is that first, figure out the maximum number of lines.
                // This will determine the number of digits (indexWidth) per line,
                // minimizing the index's file size.
                if (theFileSize == 0)
                {
                    indexWidth = -1;
                }
                else
                {
                    indexWidth = static_cast<long>(ceil(log10(static_cast<double>(theFileSize)))) + 1;
                }
                findexOut << setfill(' ') << setw(indexWidth);
                findexOut << indexWidth << endl;
                string line;
                if (hasLabels)
                {
                    getline(fin, labels);
                }
                streampos currentPos;
#ifdef DEBUG
                bool alreadyWarnedDos = false;
#endif
                while (!fin.fail())
                {
                    // This gives the current position of the line on the disk file.
                    currentPos = fin.tellg();
                    getline(fin, line);
                    if (!fin.fail())
                    {
#ifdef DEBUG
                        if (!alreadyWarnedDos && (line.find("\r") != string::npos))
                        {
                            cerr << "Warning in ReadOnlyStringFile.hpp on file " << filename << "!  DOS-type or Mac-type newline (escape-r) found in line" << endl << line << endl << "This may result in an unintentional read or write.  Suggest fixing with a DOS to Unix converter.  Continuing execution..." << endl;
                            alreadyWarnedDos = true;
                        }
#endif
                        if (line.size() > maxLength)
                        {
                            maxLength = line.size();
                        }
                        findexOut << setfill(' ') << setw(indexWidth);
                        findexOut << currentPos << endl;
                        ++numEntries;
                    }
                }
                findexOut << setfill(' ') << setw(indexWidth);
                // Write out the maximum length of all the strings
                // in the original file.
                findexOut << maxLength << endl;
                fin.close();
                findexOut.close();
                // Reload the index file size below now that it has been created.
                theIndexSize = fileSize(indexFilename);
                findex.open(indexFilename);
            }
            line.clear();
            findex >> indexWidth;
            if (indexWidth >= 0)
            {
                findex.seekg(static_cast<ifstream::pos_type>(theIndexSize) -
                             static_cast<ifstream::pos_type>(indexWidth+1));
                getline(findex, line);
                streamoff tempLineLength = convertFromString<streamoff>(line);
                // Last value is the maximum line length.
                reserve(static_cast<size_t>(tempLineLength));
                // -2 below
                numEntries = theIndexSize/(indexWidth+1) - 2;
            }
            else
            {
                line = "";
                numEntries = 0;
            }
        }
        else
        {
            // Sequential mode
            if (!noSize)
            {
                setSequentialIndexFileExtension(fileExtension);
                string indexFilename = filename + "." + sequentialIndexExtension;
                findex.open(indexFilename);
                if (!findex || isFileOlderThan(indexFilename, filename))
                {
                    // Index not found or is older than original file, so create.
                    if (findex)
                    {
                        // If got in here, then index file exists and is older than file.
                        findex.close();
                        int removeStat = remove(indexFilename);
                        if (removeStat != 0)
                        {
                            cerr << "Warning in ReadOnlyStringFile::open!  removeStat is non-zero.  Something went wrong with the file deletion.  Continuing execution." << endl;
                        }
                    }
                    makeFileSize(filename, indexFilename);
                    findex.open(indexFilename);
                }
                getline(findex, line);
                vector<string> tokens;
                splitString(line, "\t", tokens);
                numEntries = convertFromString<size_t>(tokens[0]);
                if (hasLabels && (numEntries > 0))
                {
                    --numEntries;
                }
                if (tokens.size() > 1)
                {
                    size_t lineLength = convertFromString<size_t>(tokens[1]);
                    reserve(lineLength + 1);
                }
                line.clear();
            }
        }
        fin.open(filename);
        resetToBegin();
    }

    void openResetParams(string filename,
                         bool sequentialMode = false,
                         bool noSize = false,
                         bool assumeExistLabels = false,
                         string fileExtension = "",
                         size_t resLength = 0)
    {
        // setIndexFileExtension happens in the resetParams() function.
        // setSequentialIndexFileExtension happens in the resetParams()
        // function.
        resetParams(sequentialMode, noSize, assumeExistLabels);
        if (resLength > 0)
        {
            reserve(resLength);
        }
        open(filename, assumeExistLabels, fileExtension);
    }


    void close()
    {
        // Below has indicies size == 0 since destructor will call this even
        // after the user already called close.
        fileOpen = false;
        atEnd = true;
        endPos = static_cast<ifstream::pos_type>(-1);
        fin.close();
        findex.close();
        line.clear();
    }


    size_t size()
    {
#ifdef PROG_DEBUG
        // Cannot use this when have no size.
        assert(!noSize);
#endif
        return numEntries;
    }


    bool fail()
    {
#ifdef PROG_DEBUG
        if (!sequentialMode)
        {
            cerr << "Error in ReadOnlyStringFile::fail().  This should not be used outside of sequential mode.  Exiting." << endl;
            assert(true == false);
        }
#endif
        if (atEnd ||
                (fin.tellg() == endPos) ||
                (fin.tellg() == static_cast<ifstream::pos_type>(-1)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    virtual void readRow(size_t index)
    {
#ifdef DEBUG
        assert(fileOpen == true);
#endif
#ifdef PROG_DEBUG
        assert((!noSize && index < size()) || (noSize && !atEnd));
#endif
        if (!sequentialMode)
        {
            // index+1 below since must skip over the first line which is the width.
            findex.clear();
            findex.seekg(static_cast<ifstream::pos_type>((indexWidth + 1)*(index+1)));
            getline(findex, line);
            ifstream::pos_type tempIndex = convertFromString<ifstream::pos_type>(line);

            fin.clear();
            fin.seekg(tempIndex);
            // Due to the indexing, the below is guaranteed never to fail
            // (unless the file was changed after indexing).
            getline(fin, line);
        }
        else
        {
            // Sequential mode
            if (!noSize && (index == currentLine))
            {
                return;    // already read.
            }
            if (noSize && !notReadFirst && (index == currentLine) && (currentLine > 0))
            {
                return;
            }
            if (!(noSize && notReadFirst) && (index > currentLine))
            {
                size_t numLinesDiff = index - currentLine;
                for (size_t i = 0; !atEnd && (i < numLinesDiff); ++i)
                {
                    getline(fin, line);
                    if (fin.fail())
                    {
                        atEnd = true;
                        line = "";
                    }
                }
            }
            else
            {
                resetToBegin();
                if (noSize && (endPos > static_cast<ifstream::pos_type>(0)))
                {
                    getline(fin, line);
                }
                notReadFirst = false;
                for (size_t i = 0; !atEnd && (i < index); ++i)
                {
                    getline(fin, line);
                    if (fin.fail())
                    {
                        atEnd = true;
                        line = "";
                    }
                }
            }
            currentLine = index;
        }
    }
    void getRow(size_t index, string& theRow)
    {
        readRow(index);
        theRow = line;
    }
    string getRow(size_t index)
    {
        readRow(index);
        return line;
    }
    string operator[](size_t index)
    {
        readRow(index);
        return line;
    }

    string getLabels()
    {
        return labels;
    }

    bool isOpen()
    {
        return fileOpen;
    }

    ~ReadOnlyStringFile()
    {
        close();
    }




    class iterator
    {
    public:
        iterator(ReadOnlyStringFile* inputROSF, size_t startIndex) :
            theReadOnlyStringFile(inputROSF), currentIndex(startIndex) { }
        ~iterator() { }

        iterator& operator=(const iterator& other)
        {
            theReadOnlyStringFile = other.theReadOnlyStringFile;
            currentIndex = other.currentIndex;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            return ((theReadOnlyStringFile == other.theReadOnlyStringFile) &&
                    (currentIndex == other.currentIndex));
        }

        bool operator!=(const iterator& other)
        {
            return ((theReadOnlyStringFile != other.theReadOnlyStringFile) ||
                    (currentIndex != other.currentIndex));
        }

        iterator& operator++()
        {
            if (currentIndex != theReadOnlyStringFile->size())
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

        string operator*()
        {
#ifdef PROG_DEBUG
            // Cannot use iterators in noSize mode.
            assert(!(theReadOnlyStringFile->noSize));
#endif
            return (theReadOnlyStringFile->getRow(currentIndex));
        }

    private:
        ReadOnlyStringFile* theReadOnlyStringFile;
        size_t currentIndex;
    };


    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
#ifdef PROG_DEBUG
        // Cannot use this when have no size.
        assert(!noSize);
#endif
        return iterator(this, this->size());
    }

    // Below is because STL convention is lower case while in the framework,
    // all classes start with an upper case.
    typedef iterator Iterator;

    void getIndexExtensions(string& outIndexExtension,
                            string& outSequentialIndexExtension)
    {
        outIndexExtension = indexExtension;
        outSequentialIndexExtension = sequentialIndexExtension;
    }

protected:
    Ifstream fin, findex;
    ifstream::pos_type endPos;
    bool atEnd, notReadFirst;
    size_t numEntries, currentLine;
    long indexWidth;
    bool fileOpen, sequentialMode, noSize, hasLabels;
    string indexExtension, sequentialIndexExtension;
    string labels;

    void resetParams(bool sequentialMode = false,
                     bool noSize = false,
                     bool assumeExistLabels = false)
    {
        fileOpen = false;
        atEnd = true;
        notReadFirst = true;
        endPos = static_cast<ifstream::pos_type>(-1);
        hasLabels = assumeExistLabels;
        labels = "";
        setIndexFileExtension();
        setSequentialIndexFileExtension();
        setSequentialMode(sequentialMode, noSize);
        line = "";
        currentLine = 0;
    }

    void makeFileSize(string inputFile, string indexFile)
    {
#ifdef HAVE_WC_AWK
        string theComm = "gwc -l -L " + inputFile + " | awk '{print $1\"\\t\"$2}'>& " + indexFile;
        system(theComm);
        // Where the first number $1 is the number of lines.  The second $2 is
        // the max length of a line in the file.
#else
        size_t numLines = 0, maxLen = 0;
        string line;
        Ifstream ifp(inputFile);
        while (!ifp.fail())
        {
            getline(ifp, line);
            if (ifp.fail())
            {
                break;
            }
            ++numLines;
            if (line.size() > maxLen)
            {
                maxLen = line.size();
            }
        }
        ifp.close();

        Ofstream ofp(indexFile);
        ofp << numLines << "\t" << maxLen << endl;
        ofp.close();
#endif
    }

    virtual void setIndexFileExtension(string theExt = "")
    {
        if (theExt == "")
        {
            indexExtension = "MolBioLib.ReadOnlyStringFile.index";
        }
        else
        {
            indexExtension = theExt;
        }
    }

    virtual void setSequentialIndexFileExtension(string theExt = "")
    {
        if (theExt == "")
        {
            sequentialIndexExtension = "MolBioLib.ReadOnlyStringFile.Sequential.index";
        }
        else
        {
            sequentialIndexExtension = theExt;
        }
    }

};


#endif

