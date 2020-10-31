#ifndef MOLBIOLIB_READONLYTEXTFILETABLE_H
#define MOLBIOLIB_READONLYTEXTFILETABLE_H 1

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
#include "src/Objects/Table.hpp"
#include "src/Objects/TextFileTypes.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/TupleOperations/TokenizeTuple.hpp"
#include "src/Functions/Transformers/TupleOperations/JoinTuples.hpp"
#include "src/Functions/Transformers/TupleOperations/GetTupleValues.hpp"


/** \file ReadOnlyTextFileTable.hpp
 * Contains only the #ReadOnlyTextFileTable, #ReadOnlyCSVFileTable,
 * #ReadOnlyTSVFileTable, #ReadOnlySSVFileTable, and #ReadOnlySpacesFileTable.
 * the former, using classes defined in TextFileTypes.hpp.
 */

/** Object that read-only text files with an interface like #Table
 * \section Usage Usage:
 * \code
 * typedef ReadOnlyTextFileTable<TSVType, string, string, unsigned long> myFastaTable;
 * myFastaTable myROTable("myFile.tsv", "index", true, false, 10000);
 * \endcode
 * or alternatively
 * \code
 * myFastaTable myROTable;
 * myROTable.setSequentialMode(false);
 * myROTable.open("myFile.tsv", true, "index", 10000);
 * \endcode
 * The above #TSVType is so that one may read in tab-separated-valued (TSV)
 * tables.  Alternatively, one of the other types defined in
 * TextFileTypes.hpp may be used, such as #CSVType (used to read in
 * comma-separated-valued (CSV) tables.  Alternatively to using
 * #ReadOnlyTextFileTable is to use #ReadOnlyCSVFileTable and
 * #ReadOnlyTSVFileTable both of which have the same template parameters
 * except missing the #CSVType and #TSVType parameters.  (The constructor
 * parameters are the same.)  Similarly for SSV, Spaces, and String.  The
 * String version does not allow passage of table column type, since there is
 * only one - string.
 * An index file, myFile.tsv.index, is created if one is not found.
 * If the index extension is not given, then MolBioLib.ReadOnlyStringFile.index
 * is used.  The false means that if there is no index file, the first line of
 * myFile.tsv does *not* have a header line.  Default is <code>true</code>.
 * The <code>false</code> argument is if the table is accessed in sequential
 * mode, in which case no index file is created or loaded, but is much slower
 * at accessing rows that are not right after the current one.
 * The 10000 optional argument is
 * if there is no index table, first do a reserve of 10000 on the
 * string variable reading the file.  For performance sakes, it is important
 * to make this number larger than the length of the longest line in the file.
 * \code
 * myROTable.close();  // Cannot access data after this, but then can open
 *                     // a new file.  The index array is cleared too.
 * Table<string, string, unsigned long, string, string, unsigned long> myTable;
 * ...
 * myFastaTable::row_type myRow;   // Hold a row of data.
 *    // Alternatively, can use RowType instead of row_type above.
 * size_t myRowNumber;   // The indexing of a vector uses this type.
 * \endcode
 * To get the column type, one can do:
 * \code
 * COL_TYPE_TABLE(4, myFastaTable) myColumn;
 * \endcode
 * where <code>COL_TYPE_TABLE</code> is
 * a defined macro.  (Templated typedefs are not yet
 * implemented in g++.)  Yes, this is the one place where the column
 * numbering appears within parenthesis instead of between angle brackets.
 * Alternatively, one could do:
 * \code
 * tuple_element<4, myFastaTable::cols_type>::type myColumn;
 *    // Can use ColsType instead of cols_type above.
 * \endcode
 * Alternatively, if one knew the type, such as int for column 4, one
 * could simply define
 * \code
 * vector<int> myColumn;
 * \endcode
 * Note that row_type, cols_type, and size_type are lower cased so that the
 * naming scheme conforms to the STL container convention.  Continuing, we have
 * \code
 * cout << myROTable.size();   // Returns the number of rows in the table,
 * \endcode
 * and
 * \code
 * myROTable.getRow(3, myRow);
 * \endcode
 * or alternatively
 * \code
 * myRow = myROTable.getRow(3);
 * \endcode
 * or
 * \code
 * myRow = myROTable[3];
 * \endcode
 * One may also get the rows by passing in variables to <code>getRow</code>:
 * \code
 * string s1 = "", s2 = "";
 * unsigned long t1 = 0;
 * myROTable.getRow(3, s1, s2, t1);
 * \endcode
 * In the above example, it is important to order the variables in the same
 * order as they appear in the declaration of <code>myROTable</code>.
 *
 * To access a specific element, one may do
 * \code
 * unsigned long aLong = myROTable.getElement<2>(3);
 * \endcode
 * where the 3 designates the 4th row and 2 is the 3rd column (unsigned ints).
 *
 * For columns <b>(only available in random access mode)</b>,
 * \code
 * myROTable.getCol<4>(myColumn);   // 4 is the last column, i.e. column 5.
 * \endcode
 * \code
 * myROTable.getCol<4>(myColumn, 24, 99, 2);   // Get rows 25 through
 *                                             // 100 (inclusive) only and
 *                                             // fill from third row.
 * \endcode
 * Alternatively (but more expensive and uses more memory), one could do the
 * syntatically simpler
 * \code
 * myColumn = myROTable.getCol<4>();
 * \endcode
 * and
 * \code
 * myColumn = myROTable.getCol<4>(24, 99, 2);
 * \endcode
 * <b>myTable below must be resized to the appropriate size or else
 * a segmentation error will occur! </b>
 * \code
 * myROTable.copyToTable(myTable);   // Copies all of myROTable to myTable
 *                                   // starting from the first column in
 *                                   // myTable.
 * \endcode
 * and
 * \code
 * myROTable.copyToTable<3>(myTable, 24, 99, 4, true);
 * \endcode
 * which copies data to <code>myTable</code>, but starts at the 4th column in
 * <code>myTable</code>.
 * \section ProgNotes
 * <ul>
 *    <li>
 *    We would have liked to have implemented
 *    \code
 *    myFastaTable::cols_type<4> myColumn;
 *    \endcode
 *    but unfortunately g++ does not yet implement templated typedefs.
 *    </li>
 * </ul>
 *
 * Finally, we can iterate through the table:
 * \code
 * for (myFastaTable::iterator i = myROTable.begin(); i != myROTable.end(); ++i) {
 *    myFastaTable::row_type theRow = *i;
 *    cout << theRow << endl;   // This uses the tuple streams
 * }
 * \endcode
 * <code>Iterator</code> may be used in place of <code>iterator</code>.
 */
template<typename Q, typename... T>
class ReadOnlyTextFileTable : public ReadOnlyStringFile
{
public:
    typedef tuple< vector<T>... > cols_type;
    cols_type theData;  // makes tuple<vector<T1>,vector<T2>,...>
    // Below is the size_type (and the index variable type) of the vector.
    typedef size_t size_type;
    typedef tuple<T...> row_type;
    vector<string> labels;
    vector<string> tokens;
    row_type row;

    Q splitter;

    // Constructor
    ReadOnlyTextFileTable<Q, T...>()
    {
        ReadOnlyStringFile();
        labels.resize(tuple_size<row_type>::value, "");
    }
    // No choice - must use random access mode due to headers.
    ReadOnlyTextFileTable<Q, T...>(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        labels.resize(tuple_size<row_type>::value, "");
        openResetParams(filename, false, false, assumeExistLabels, fileExtension, reserveLength);
    }

    virtual void open(string filename, bool assumeExistLabels = true,
                      string fileExtension = "", size_t reserveLength = 0)
    {

        // Check if file exists
        if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
        {
            cerr << "Error!  The filename " << filename << " cannot be opened." << endl;
            assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
        }
        if (reserveLength > 0)
        {
            line.reserve(reserveLength);
        }
        // Open in random access mode.  Needed in case of headers.
        if (fileExtension == "")
        {
            ReadOnlyStringFile::open(filename, assumeExistLabels, this->indexExtension, assumeExistLabels);
        }
        else
        {
            ReadOnlyStringFile::open(filename, assumeExistLabels, fileExtension);
        }
        // Since fin.open() in ReadOnlyStringFile opens at the beginning, no need
        // to seek there.
        ifstream::pos_type tempIndex;
        if (!sequentialMode)
        {
            findex.clear();
            findex.seekg(indexWidth + 1);
            getline(findex, line);
            tempIndex = convertFromString<ifstream::pos_type>(line);
        }
        else
        {
            if (assumeExistLabels)
            {
                tempIndex = static_cast<ifstream::pos_type>(1);
            }
            else
            {
                tempIndex = static_cast<ifstream::pos_type>(0);
            }
        }
        if (tempIndex > 0)
        {
            // There are labels
            fin.clear();
            fin.seekg(0, ios_base::beg);
            getline(fin, line);
            splitter.splitLine(line, tokens);
            for (size_t i = 0; i < tuple_size<row_type>::value; ++i)
            {
                labels[i] = tokens[i];
            }
        }
        resetToBegin();
    }


    void getRow(size_t index, row_type& theRow)
    {
        readRow(index);
        splitter.splitLine(line, tokens);
        tuplizeTokens(tokens, theRow);
        tuplizeTokens(tokens, row);
    }
    void getRow(size_t index, T&... args)
    {
        row_type tempRow;
        getRow(index, tempRow);
        getTupleValues(tempRow, args...);
    }
    row_type getRow(size_t index)
    {
        getRow(index, row);
        return row;
    }
    row_type operator[](size_t index)
    {
        return getRow(index);
    }

    template<size_t N>
    typename tuple_element<N,row_type>::type& getElement(size_t index)
    {
        getRow(index);
        return get<N>(row);
    }


    template<size_t N>
    void getCol(typename tuple_element<N,cols_type>::type& theCol,
                size_t min = 0, size_t max = numeric_limits<size_t>::max(),
                size_t offset = 0)
    {
#ifdef DEBUG
        assert(fileOpen == true);
        assert(!sequentialMode);
#endif
        size_t numRows = size() - 1;
        if (max != numeric_limits<size_t>::max())
        {
            numRows = max;
        }
        string line;
        findex.clear();
        if (hasLabels && (min == 0))
        {
            findex.seekg((indexWidth+1)*(min+1));
        }
        else
        {
            findex.seekg((indexWidth+1)*min);
        }
        getline(findex, line);
        ifstream::pos_type tempIndex = convertFromString<ifstream::pos_type>(line);
        fin.clear();
        fin.seekg(tempIndex);
        if (theCol.size() <= (numRows-min+offset))
        {
            theCol.resize(numRows-min+offset+1);
        }
        for (size_t i = min; i <= numRows; ++i)
        {
            getline(fin, line);
            splitter.splitLine(line, tokens);
            istringstream is(tokens[N]);
            is >> theCol[i-min+offset];
        }
    }


    template<size_t N>
    typename tuple_element<N,cols_type>::type getCol(size_t min = 0,
            size_t max = numeric_limits<size_t>::max(),
            size_t offset = 0)
    {
        typename tuple_element<N,cols_type>::type theCol;
        getCol(theCol, min, max, offset);
        return theCol;
    }




    template<size_t N = 0, typename... U>
    void copyToTable(Table<U...>& theTable,
                     size_t first = 0,
                     size_t last = numeric_limits<size_t>::max(),
                     size_t offset = 0, bool copyLabels = true)
    {
        size_t numRows = size() - 1;
        if (last != numeric_limits<size_t>::max())
        {
            numRows = last;
        }

        row_type thisRow;
        typename Table<U...>::row_type tableRow;

        if (copyLabels)
        {
            for (size_t i = 0; i < labels.size(); ++i)
            {
                theTable.labels[i+N] = labels[i];
            }
        }

        if (theTable.size() <= (numRows+offset))
        {
            theTable.resize(numRows+offset+1);
        }
        for (size_t i = first; i <= numRows; ++i)
        {
            getRow(i, thisRow);
            theTable.getRow(i+offset, tableRow);
            joinTuples<N>(thisRow, tableRow);
            theTable.setRow(i+offset, tableRow);
        }
    }


    // Do not need a destructor.  ReadOnlyStringFile's destructor is
    // called implicitly.


    class iterator
    {
    public:
        iterator(ReadOnlyTextFileTable* inputROTFT, size_t startIndex) :
            theReadOnlyTextFileTable(inputROTFT), currentIndex(startIndex) { }
        ~iterator() { }

        iterator& operator=(const iterator& other)
        {
            theReadOnlyTextFileTable = other.theReadOnlyTextFileTable;
            currentIndex = other.currentIndex;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            return ((theReadOnlyTextFileTable ==
                     other.theReadOnlyTextFileTable) &&
                    (currentIndex == other.currentIndex));
        }

        bool operator!=(const iterator& other)
        {
            return ((theReadOnlyTextFileTable !=
                     other.theReadOnlyTextFileTable) ||
                    (currentIndex != other.currentIndex));
        }

        iterator& operator++()
        {
            if (currentIndex != theReadOnlyTextFileTable->size())
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

        row_type operator*()
        {
            return (theReadOnlyTextFileTable->getRow(currentIndex));
        }

    private:
        ReadOnlyTextFileTable* theReadOnlyTextFileTable;
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

};






/** Object that read-only text files with an interface like #Table
 * #ReadOnlyCSVFileTable is a specialization of #ReadOnlyTextFileTable
 * such that it reads comma-separated-valued (CSV) tables.
 */
template<typename... T>
class ReadOnlyCSVFileTable : public ReadOnlyTextFileTable<CSVType, T...>
{
public:
    ReadOnlyCSVFileTable()
    {
        ReadOnlyTextFileTable<CSVType, T...>();
    }
    ReadOnlyCSVFileTable(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        ReadOnlyTextFileTable<CSVType, T...>::labels.resize(tuple_size< tuple<T...> >::value, "");
        if (reserveLength > 0)
        {
            ReadOnlyTextFileTable<CSVType, T...>::line.reserve(reserveLength);
        }
        if (fileExtension == "")
        {
            ReadOnlyTextFileTable<CSVType, T...>::open(filename, assumeExistLabels, this->indexExtension, reserveLength);
        }
        else
        {
            ReadOnlyTextFileTable<CSVType, T...>::open(filename, assumeExistLabels, fileExtension, reserveLength);
        }
    }
    ~ReadOnlyCSVFileTable()
    {
        ReadOnlyTextFileTable<CSVType, T...>::close();
    }
};



/** Object that read-only text files with an interface like #Table
 * #ReadOnlyTSVFileTable is a specialization of #ReadOnlyTextFileTable
 * such that it reads tab-separated-valued (TSV) tables.
 */
template<typename... T>
class ReadOnlyTSVFileTable : public ReadOnlyTextFileTable<TSVType, T...>
{
public:
    ReadOnlyTSVFileTable()
    {
        ReadOnlyTextFileTable<TSVType, T...>();
    }
    ReadOnlyTSVFileTable(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        ReadOnlyTextFileTable<TSVType, T...>::labels.resize(tuple_size< tuple<T...> >::value, "");
        if (reserveLength > 0)
        {
            ReadOnlyTextFileTable<TSVType, T...>::line.reserve(reserveLength);
        }
        if (fileExtension == "")
        {
            ReadOnlyTextFileTable<TSVType, T...>::open(filename, assumeExistLabels, this->indexExtension, reserveLength);
        }
        else
        {
            ReadOnlyTextFileTable<TSVType, T...>::open(filename, assumeExistLabels, fileExtension, reserveLength);
        }
    }
    ~ReadOnlyTSVFileTable()
    {
        ReadOnlyTextFileTable<TSVType, T...>::close();
    }
};



/** Object that read-only text files with an interface like #Table
 * #ReadOnlySSVFileTable is a specialization of #ReadOnlyTextFileTable
 * such that it reads space-separated-valued (SSV) tables.
 */
template<typename... T>
class ReadOnlySSVFileTable : public ReadOnlyTextFileTable<SSVType, T...>
{
public:
    ReadOnlySSVFileTable()
    {
        ReadOnlyTextFileTable<SSVType, T...>();
    }
    ReadOnlySSVFileTable(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        ReadOnlyTextFileTable<SSVType, T...>::labels.resize(tuple_size< tuple<T...> >::value, "");
        if (reserveLength > 0)
        {
            ReadOnlyTextFileTable<SSVType, T...>::line.reserve(reserveLength);
        }
        if (fileExtension == "")
        {
            ReadOnlyTextFileTable<SSVType, T...>::open(filename, assumeExistLabels, this->indexExtension, reserveLength);
        }
        else
        {
            ReadOnlyTextFileTable<SSVType, T...>::open(filename, assumeExistLabels, fileExtension, reserveLength);
        }
    }
    ~ReadOnlySSVFileTable()
    {
        ReadOnlyTextFileTable<SSVType, T...>::close();
    }
};



/** Object that read-only text files with an interface like #Table
 * #ReadOnlySpacesFileTable is a specialization of #ReadOnlyTextFileTable
 * such that it reads space separated (Spaces) tables.
 */
template<typename... T>
class ReadOnlySpacesFileTable : public ReadOnlyTextFileTable<SpacesType, T...>
{
public:
    ReadOnlySpacesFileTable()
    {
        ReadOnlyTextFileTable<SpacesType, T...>();
    }
    ReadOnlySpacesFileTable(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        ReadOnlyTextFileTable<SpacesType, T...>::labels.resize(tuple_size< tuple<T...> >::value, "");
        if (reserveLength > 0)
        {
            ReadOnlyTextFileTable<SpacesType, T...>::line.reserve(reserveLength);
        }
        if (fileExtension == "")
        {
            ReadOnlyTextFileTable<SpacesType, T...>::open(filename, assumeExistLabels, this->indexExtension, reserveLength);
        }
        else
        {
            ReadOnlyTextFileTable<SpacesType, T...>::open(filename, assumeExistLabels, fileExtension, reserveLength);
        }
    }
    ~ReadOnlySpacesFileTable()
    {
        ReadOnlyTextFileTable<SpacesType, T...>::close();
    }
};



/** Object that read-only text files with an interface like #Table
 * #ReadOnlySpacesFileTable is a specialization of #ReadOnlyTextFileTable
 * such that it reads space separated (Spaces) tables.
 */
class ReadOnlyStringFileTable : public ReadOnlyTextFileTable<StringType, string>
{
public:
    ReadOnlyStringFileTable()
    {
        ReadOnlyTextFileTable<StringType, string>();
    }
    ReadOnlyStringFileTable(string filename, bool assumeExistLabels = true, string fileExtension = "", size_t reserveLength = 0)
    {
        ReadOnlyTextFileTable<StringType, string>::labels.resize(1, "");
        if (reserveLength > 0)
        {
            ReadOnlyTextFileTable<StringType, string>::line.reserve(reserveLength);
        }
        if (fileExtension == "")
        {
            ReadOnlyTextFileTable<StringType, string>::open(filename, assumeExistLabels, this->indexExtension, reserveLength);
        }
        else
        {
            ReadOnlyTextFileTable<StringType, string>::open(filename, assumeExistLabels, fileExtension, reserveLength);
        }
    }
    ~ReadOnlyStringFileTable()
    {
        ReadOnlyTextFileTable<StringType, string>::close();
    }
};


#endif

