#ifndef MOLBIOLIB_READONLYDELIMITEDFILE_H
#define MOLBIOLIB_READONLYDELIMITEDFILE_H 1

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


#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Objects/TextFileTypes.hpp"


/** \file ReadOnlyDelimitedFile.hpp
 * Contains only the #ReadOnlyDelimitedFile class.
 */

/** Object that read-only files with a vector interface.  Has to be a
 * type of text file associated with a text type.  Usage is the same as
 * a #ReadOnlyStringFile (as it is derived from it) but that there are
 * new operators for getRow (both ways) and operator[], which return a
 * <code>vector<string></code>, depending on the type of text type passed
 * in.  CSV and TSV types have been aliased to #ReadOnlyCSVFile and
 * #ReadOnlyTSVFile.  Typical usage is:
 * \code
 * string FILE = ...;
 * ReadOnlyTSVFile myFile(FILE);
 * size_t numRows = myFile.size();
 * for (size_t i = 0; i < numRows; ++i)
 *    cout << (myFile[i])[4] << endl;   // Print out 5th column of file.
 *            // The operator[] returns a vector of string here.
 *    // More efficient (no copy) but less elegant is:
 *    { myFile.readRow(i); cout << myFile.tokens[4] << endl; }
 * myFile.close();
 * \endcode
 */
template<typename Q>
class ReadOnlyDelimitedFile : public ReadOnlyStringFile
{
public:
    vector<string> tokens;
    Q splitter;
    ReadOnlyDelimitedFile() : ReadOnlyStringFile() { }

    ReadOnlyDelimitedFile(string filename,
                          bool sequentialMode = false,
                          bool noSize = false,
                          bool assumeExistLabels = false,
                          string fileExtension = "",
                          size_t resLength = 0) :
        ReadOnlyStringFile(filename, sequentialMode, noSize, assumeExistLabels,
                           fileExtension, resLength) { }


    void readRow(size_t index)
    {
        ReadOnlyStringFile::readRow(index);
        splitter.splitLine(line, tokens);
    }

    void getRow(size_t index, vector<string>& theRow)
    {
        readRow(index);
        theRow = tokens;
    }
    vector<string> getRow(size_t index)
    {
        readRow(index);
        return tokens;
    }
    vector<string> operator[](size_t index)
    {
        readRow(index);
        return tokens;
    }

    vector<string> getLabels()
    {
        vector<string> labelTokens;
        splitter.splitLine(labels, labelTokens);
        return labelTokens;
    }

    // Do not need a destructor.  ReadOnlyStringFile's destructor is
    // implicitly called.


    class iterator
    {
    public:
        iterator(ReadOnlyDelimitedFile<Q>* inputROSF, size_t startIndex) :
            theReadOnlyDelimitedFile(inputROSF), currentIndex(startIndex) { }
        ~iterator() { }

        iterator& operator=(const iterator& other)
        {
            theReadOnlyDelimitedFile = other.theReadOnlyDelimitedFile;
            currentIndex = other.currentIndex;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            return ((theReadOnlyDelimitedFile ==
                     other.theReadOnlyDelimitedFile) &&
                    (currentIndex == other.currentIndex));
        }

        bool operator!=(const iterator& other)
        {
            return ((theReadOnlyDelimitedFile !=
                     other.theReadOnlyDelimitedFile) ||
                    (currentIndex != other.currentIndex));
        }

        iterator& operator++()
        {
            if (currentIndex != theReadOnlyDelimitedFile->size())
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

        vector<string> operator*()
        {
            return (theReadOnlyDelimitedFile->getRow(currentIndex));
        }

    private:
        ReadOnlyDelimitedFile<Q>* theReadOnlyDelimitedFile;
        size_t currentIndex;
    };


    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
#ifdef PROG_DEBUG
        assert(!(this->noSize));
#endif
        return iterator(this, this->size());
    }

    // Below is because STL convention is lower case while in the framework,
    // all classes start with an upper case.
    typedef iterator Iterator;
};


/** \file ReadOnlyDelimitedFile.hpp
 * \var typedef ReadOnlyDelimitedFile<CSVType> ReadOnlyCSVFile
 * Define the #ReadOnlyCSVFile.
 */
typedef ReadOnlyDelimitedFile<CSVType> ReadOnlyCSVFile;



/** \file ReadOnlyDelimitedFile.hpp
 * \var typedef ReadOnlyDelimitedFile<TSVType> ReadOnlyTSVFile
 * Define the #ReadOnlyTSVFile.
 */
typedef ReadOnlyDelimitedFile<TSVType> ReadOnlyTSVFile;



/** \file ReadOnlyDelimitedFile.hpp
 * \var typedef ReadOnlyDelimitedFile<SSVType> ReadOnlySSVFile
 * Define the #ReadOnlySSVFile.
 */
typedef ReadOnlyDelimitedFile<SSVType> ReadOnlySSVFile;



/** \file ReadOnlyDelimitedFile.hpp
 * \var typedef ReadOnlyDelimitedFile<SpacesType> ReadOnlySpacesFile
 * Define the #ReadOnlySpacesFile.
 */
typedef ReadOnlyDelimitedFile<SpacesType> ReadOnlySpacesFile;





#endif

