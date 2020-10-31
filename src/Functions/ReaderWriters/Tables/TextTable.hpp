#ifndef MOLBIOLIB_TEXTTABLE_H
#define MOLBIOLIB_TEXTTABLE_H 1

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
#include "src/Objects/TextFileTypes.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Functions/Transformers/TupleOperations/TokenizeTuple.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"



/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename S, typename... T> void readTextTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read text tables.
 * \section Usage Usage examples:
 * \code
 * Table<string, string, unsigned int, double> myTable;
 * \endcode
 * Below fails there is an unequal number of rows in myTable between columns.
 * \code
 * readTextTable<CSVType>(myTable, "newfile.csv", true);
 * \endcode
 * where <code>CSVType</code> is a predefined class that describes how
 * comma-separated-valued tables are to be read in.  The
 * predefined classes currently available are:
 * - CSVType
 *   - Comma separated valued tables.
 * - TSVType
 *   - Tab separated valued tables.
 * - StringType
 *   - Single string valued table (StringTable).
 * These types are defined in TextFileTypes.hpp.
 * The <code>true</code> in the above means there is no header line.  (Default
 * is false.)  If one wishes to read only a subset of the tables, one can do
 * \code
 * readTextTable<TSVType>(myTable, "newfile.tsv", false, 4, 99);
 *    // To read in the 5th through 100th lines.  If one wishes to read to the
 *    // end of the file, specify numeric_limits<size_t>::max() as the second argument.
 *    // Also, can replace readTextTable<TSVType> by readTSVTable.
 * \endcode
 * For a faster read, one can also do
 * \code
 * readTextTable<TSVType>(myTable, "newfile.tsv", false, 0, numeric_limits<size_t>::max(), 100);
 *   // where the 100 means that this routine will reserve 100 entries
 *   // thus making the reads go faster.  Optimally, one should
 *   // reserve more than enough.
 * \endcode
 * For long lines, one may need to specify a string reserve size:
 * \code
 * readTextTable<CSVType>(myTable, "newfile.csv", false, 0, numeric_limits<size_t>::max(), 100, 200000);
 *    // Can replace readTextTable<CSVType> by readCSVTable.
 * \endcode
 * where the above 200000 is the reserve size.  In this case, one must
 * specify the number of entries to reserve beforehand.  There is also
 * \code
 * readTextTable<SSVType>(myTable, "newfile.csv", false, 0, numeric_limits<size_t>::max(), 100, 200000);
 *    // Can replace readTextTable<SSVType> by readSSVTable.
 * \endcode
 * which reads in space-delimted files.  Note that there is the special case:
 * \code
 * StringTable sTable;
 * readTextTable<StringType>(sTable, "newfile.txt");
 * // or readStringTable(sTable, "newfile.txt");
 * \endcode
 * where the above parameters also apply.
 * \section ProgNotes Programming Notes:
 * In the below, int's are used to traverse columns, since the
 * tuple's indexing is also of type int.
 */
template<typename Q, typename... T>
void readTextTable(Table<T...>& theTable, string filename, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   size_t numAlloc = 0,
                   size_t minStringReserve = 0)
{

    if (filename == "")
    {
#ifdef DEBUG
        cerr << "Warning from readTextTable.  An empty filename was passed.  This will return an empty Table." << endl;
#endif
        return;
    }

    // Check to see the file exists.
    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error!  The file " << filename << " cannot be opened." << endl;
        assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
    }
    Ifstream fin(filename);
    typedef typename Table<T...>::row_type row_type;
    string line;
    Q splitter;
    if (minStringReserve > 0)
    {
        line.reserve(minStringReserve);
    }
    if (theTable.labels.size() < static_cast<size_t>(tuple_size<row_type>::value))
    {
        theTable.labels.resize(static_cast<size_t>(tuple_size<row_type>::value));
    }
    if (!noLabels)
    {
        vector<string> tokens;
        getline(fin, line);
        splitter.splitLine(line, tokens);
        for (size_t i = 0; i < tuple_size<row_type>::value; ++i)
        {
            theTable.labels[i] = tokens[i];
        }
    }
    if (numAlloc > 0)
    {
        theTable.reserve(numAlloc);
    }
    size_t currentRow = 0;
    row_type tempRow;
    // Need proceed variable so can terminate while early.
    bool proceed = true;
    while(!fin.fail() && proceed)
    {
        getline(fin, line);
        if (!fin.fail())
        {
            if (currentRow >= first)
            {
                proceed = true;
                if ((last != numeric_limits<size_t>::max()) &&
                        (currentRow > last))
                {
                    proceed = false;
                }
                if (proceed)
                {
                    vector<string> tokens;
                    splitter.splitLine(line, tokens);
                    tuplizeTokens(tokens, tempRow);
                    theTable.push_back(tempRow);
                }
            }
            ++currentRow;
        }
    }
    fin.close();
}

/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void readCSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read comma-separated-valued (CSV) tables.  This is just #readTextTable with
 * type #CSVType passed.
 */
template<typename... T>
void readCSVTable(Table<T...>& theTable, string filename,
                  bool noLabels = false,
                  size_t first = 0, size_t last = numeric_limits<size_t>::max(),
                  size_t numAlloc = 0,
                  size_t minStringReserve = 0)
{
    readTextTable<CSVType>(theTable, filename, noLabels,
                           first, last, numAlloc, minStringReserve);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void readTSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read tab-separated-valued (TSV) tables.  This is just #readTextTable with
 * type #TSVType passed.
 */
template<typename... T>
void readTSVTable(Table<T...>& theTable, string filename,
                  bool noLabels = false,
                  size_t first = 0,
                  size_t last = numeric_limits<size_t>::max(),
                  size_t numAlloc = 0,
                  size_t minStringReserve = 0)
{
    readTextTable<TSVType>(theTable, filename, noLabels,
                           first, last, numAlloc, minStringReserve);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void readSSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read space-separated-valued (SSV) tables.  This is just #readTextTable with
 * type #SSVType passed.
 */
template<typename... T>
void readSSVTable(Table<T...>& theTable, string filename,
                  bool noLabels = false,
                  size_t first = 0,
                  size_t last = numeric_limits<size_t>::max(),
                  size_t numAlloc = 0,
                  size_t minStringReserve = 0)
{
    readTextTable<SSVType>(theTable, filename, noLabels,
                           first, last, numAlloc, minStringReserve);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void readSpacesTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read spaces-separated-valued tables.  This is just #readTextTable with
 * type #SpacesType passed.
 */
template<typename... T>
void readSpacesTable(Table<T...>& theTable, string filename,
                     bool noLabels = false,
                     size_t first = 0,
                     size_t last = numeric_limits<size_t>::max(),
                     size_t numAlloc = 0,
                     size_t minStringReserve = 0)
{
    readTextTable<SpacesType>(theTable, filename, noLabels,
                              first, last, numAlloc, minStringReserve);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void readStringTable(StringTable& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t numAlloc = 0, size_t minStringReserve = 0)
 * Read a StringTable.  This is just #readTextTable with
 * type #StringType passed.
 */
void readStringTable(StringTable& theTable, string filename,
                     bool noLabels = false,
                     size_t first = 0,
                     size_t last = numeric_limits<size_t>::max(),
                     size_t numAlloc = 0,
                     size_t minStringReserve = 0)
{
    readTextTable<StringType>(theTable, filename, noLabels,
                              first, last, numAlloc, minStringReserve);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename Q, typename... T> void writeTextTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool generateIndex = false)
 * Write text tables.
 * \section Usage Usage examples:
 * \code
 * Table<string, string, unsigned int, double> myTable;
 * \endcode
 * Below fails if there is an unequal number of rows in myTable between columns.
 * \code
 * writeTextTable<CSVType>(myTable, "newfile.csv");
 *                                         // Print by default header line.
 *                                         // Info stored in labels variable.
 *    // Can replace writeTextTable<CSVType> by writeCSVTable.
 * \endcode
 * Note that one may also specify a first and last line for the above.  See
 * <code>readTextTable</code> for details on the type of classes (such as
 * <code>CSVType</code> above and <code>SSVType</code>) one may
 * pass as a parameter.
 * \code
 * writeTextTable<TSVType>(myTable, "newfile.tsv", false, 0, numeric_limits<size_t>::max(), true);
 *    // Can replace writeTextTable<TSVType> by writeTSVTable.
 * \endcode
 * and
 * \code
 * writeTextTable<SSVType>(myTable, "newfile.tsv", false, 0, numeric_limits<size_t>::max(), true);
 *    // Can replace writeTextTable<SSVType> by writeSSVTable.
 * \endcode
 * For the above, an associated index file, newfile.tsv.index, is
 * created because of the <code>true</code> parameter.  The index
 * file contains [std::]streamoff (stream offset)
 * data as a column of numbers.  One may then do random access of the
 * file with this index.  The column of numbers only point to the data.
 * It does not point to the labels.  One can read back in streamoff
 * then convert to streampos.
 * As before, one also has
 * \code
 * StringTable sTable;
 * // sTable filled in...
 * writeTextTable<StringType>(myTable, "newfile.txt");
 *    // or replace with writeStringTable(myTable, "newfile.txt");
 * \endcode
 * with the aforementioned optional parameters available.
 * \section ProgNotes Programming Notes:
 * In the below, int's are used to traverse columns, since the
 * tuple's indexing is also of type int.
 */
template<typename Q, typename... T>
void writeTextTable(Table<T...>& theTable,
                    string filename, bool noLabels = false,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    bool generateIndex = false)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    Ofstream fout(filename);
    Ofstream findex;
    Q joiner;
    string line;
    if (generateIndex)
    {
        string indexFilename = filename + ".index";
        findex.open(indexFilename);
    }
    typedef typename Table<T...>::row_type row_type;
    if (!noLabels)
    {
        joiner.joinTokens(theTable.labels, line);
        fout << line << endl;
    }
    size_t numRows = theTable.size();
    if (last != numeric_limits<size_t>::max())
    {
        numRows = last + 1;
    }
    row_type tempRow;
    for (size_t i = first; i < numRows; ++i)
    {
        theTable.getRow(i, tempRow);
        vector<string> tokens;
        tokenizeTuple(tempRow, tokens);
        if (generateIndex)
        {
            findex << fout.tellp() << endl;
        }
        joiner.joinTokens(tokens, line);
        fout << line << endl;
    }
    fout.close();
    if (generateIndex)
    {
        findex.close();
    }
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void writeCSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool generateIndex = false)
 * Write comma-separated-valued (CSV) tables.  This is just #writeTextTable with
 * type #CSVType passed.
 */
template<typename... T>
void writeCSVTable(Table<T...>& theTable,
                   string filename, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   bool generateIndex = false)
{
    writeTextTable<CSVType>(theTable, filename, noLabels,
                            first, last, generateIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void writeSSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool generateIndex = false)
 * Write space-separated-valued (SSV) tables.  This is just #writeTextTable with
 * type #SSVType passed.
 */
template<typename... T>
void writeSSVTable(Table<T...>& theTable,
                   string filename, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   bool generateIndex = false)
{
    writeTextTable<SSVType>(theTable, filename, noLabels,
                            first, last, generateIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void writeTSVTable(Table<T...>& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool generateIndex = false)
 * Write tab-separated-valued (TSV) tables.  This is just #writeTextTable with
 * type #TSVType passed.
 */
template<typename... T>
void writeTSVTable(Table<T...>& theTable,
                   string filename, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   bool generateIndex = false)
{
    writeTextTable<TSVType>(theTable, filename, noLabels,
                            first, last, generateIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void writeStringTable(StringTable& theTable, string filename, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), bool generateIndex = false)
 * Write string tables.  This is just #writeTextTable with
 * type #StringType passed.
 */
void writeStringTable(StringTable& theTable,
                      string filename, bool noLabels = false,
                      size_t first = 0,
                      size_t last = numeric_limits<size_t>::max(),
                      bool generateIndex = false)
{
    writeTextTable<StringType>(theTable, filename, noLabels,
                               first, last, generateIndex);
}



/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename Q, typename... T> void printTextTable(Table<T...>& theTable, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), ostream& ofp = cout, bool printIndex = false)
 * Print text tables.
 * \section Usage Usage examples:
 * \code
 * Table<string, string, unsigned int, double> myTable;
 * \endcode
 * Below fails there is an unequal number of rows in myTable between columns.
 * \code
 * printTextTable<TSVType>(myTable, 0, numeric_limits<size_t>::max(), cout, false);
 *                          // Like writeTSV, except does not take in a
 *                          // filename and instead prints to the screen.
 * // Start from 0 to last (max() = last), outputting to cout, and not
 * // adding row number as first column (false).
 * \endcode
 * Similar to above, one may replace TSV by CSV or SSV.
 * Note that one may also specify a first and last line for the above.  See
 * <code>readTextTable</code> for text table types that may be passed.
 * \section ProgNotes Programming Notes:
 * In the below, int's are used to traverse columns, since the
 * tuple's indexing is also of type int.
 */
template<typename Q, typename... T>
void printTextTable(Table<T...>& theTable, bool noLabels = false,
                    size_t first = 0,
                    size_t last = numeric_limits<size_t>::max(),
                    ostream& ofp = cout,
                    bool printIndex = false)
{
#ifdef PROG_DEBUG
    assert(theTable.sameRowSizes() == true);
#endif
    typedef typename Table<T...>::row_type row_type;
    Q joiner;
    string line;
    if (!noLabels)
    {
        joiner.joinTokens(theTable.labels, line);
        if (!printIndex)
        {
            ofp << line << endl;
        }
        else
        {
            ofp << "Row number" << joiner.delimiter() << line << endl;
        }
    }
    size_t numRows = theTable.size();
    if (last != numeric_limits<size_t>::max())
    {
        numRows = last;
    }
    row_type tempRow;
    for (size_t i = first; i < numRows; ++i)
    {
        theTable.getRow(i, tempRow);
        vector<string> tokens;
        tokenizeTuple(tempRow, tokens);
        joiner.joinTokens(tokens, line);
        if (!printIndex)
        {
            ofp << line << endl;
        }
        else
        {
            ofp << convertToString<size_t>(i) << joiner.delimiter() << line << endl;
        }
    }
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void printCSVTable(Table<T...>& theTable, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), ostream& ofp = cout, bool printIndex = false)
 * Print comma-separated-valued (CSV) tables.  This is just #printTextTable with
 * type #CSVType passed.
 */
template<typename... T>
void printCSVTable(Table<T...>& theTable, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   ostream& ofp = cout, bool printIndex = false)
{
    printTextTable<CSVType>(theTable, noLabels, first, last, ofp, printIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void printSSVTable(Table<T...>& theTable, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), ostream& ofp = cout, bool printIndex = false)
 * Print space-separated-valued (SSV) tables.  This is just #printTextTable with
 * type #SSVType passed.
 */
template<typename... T>
void printSSVTable(Table<T...>& theTable, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   ostream& ofp = cout, bool printIndex = false)
{
    printTextTable<SSVType>(theTable, noLabels, first, last, ofp, printIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void printTSVTable(Table<T...>& theTable, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), ostream& ofp = cout, bool printIndex = false)
 * Print tab-separated-valued (TSV) tables.  This is just #printTextTable with
 * type #TSVType passed.
 */
template<typename... T>
void printTSVTable(Table<T...>& theTable, bool noLabels = false,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   ostream& ofp = cout, bool printIndex = false)
{
    printTextTable<TSVType>(theTable, noLabels, first, last, ofp, printIndex);
}


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void printStringTable(StringTable& theTable, bool noLabels = false, size_t first = 0, size_t last = numeric_limits<size_t>::max(), ostream& ofp = cout, bool printIndex = false)
 * Print string tables.  This is just #printTextTable with
 * type #StringType passed.
 */
void printStringTable(StringTable& theTable, bool noLabels = false,
                      size_t first = 0,
                      size_t last = numeric_limits<size_t>::max(),
                      ostream& ofp = cout, bool printIndex = false)
{
    printTextTable<StringType>(theTable, noLabels, first, last, ofp, printIndex);
}



#ifdef DEBUG

template<size_t N, typename C, typename S>
class TableWrite
{
public:
    static void print(C& theData, string filename)
    {
        TableWrite<N-1,C,S>::print(theData, filename);
        Ofstream fout;
        fout.open(filename, ios::app);
        S numRows = get<N>(theData).size();
        for (S i = 0; i < numRows; ++i)
        {
            fout << N << "\t" << i << "\t" << get<N>(theData)[i] << endl;
        }
        fout.close();
    }
};
template<typename C, typename S>
class TableWrite<0,C,S>
{
public:
    static void print(C& theData, string filename)
    {
        Ofstream fout;
        fout.open(filename, ios::app);
        S numRows = get<0>(theData).size();
        for (S i = 0; i < numRows; ++i)
        {
            fout << "0\t" << i << "\t" << get<0>(theData)[i] << endl;
        }
        fout.close();
    }
};


/** \file TextTable.hpp
 * Read and write text tables.
 * \fn template<typename... T> void writeTable(Table<T...>& theTable, string filename)
 * Write #Table in tab-separated form that do not have the same number of rows
 * in each column.  These are not strictly TSV tables.  This is included for
 * purely debugging purposes to visually inspect tables with unequal rows.
 * \section Usage Usage examples:
 * \code
 * Table<string, string, unsigned int, double> myTable;
 * \endcode
 * If DEBUG is defined, then the below is available.
 * \code
 * writeTable(myTable, "newfile.txt");  // Prints a textual output of the table
 *                                      // with first line being labels and
 *                                      // subsequent lines being of the form
 *                                      // vector # [tab] row # [tab] entry.
 *                                      // Can have unequal rows.
 * \endcode
 * Note that one may *NOT* specify a first and last line for the above.
 * The above is because #writeTable can write unevenrow'd #Tables.
 * \section ProgNotes Programming Notes:
 * In the below, int's are used to traverse columns, since the
 * tuple's indexing is also of type int.
 */
template<typename... T>
void writeTable(Table<T...>& theTable, string filename)
{
    typedef typename Table<T...>::row_type row_type;
    typedef typename Table<T...>::cols_type cols_type;
    Ofstream fout;
    fout.open(filename);
    for (size_t i = 0; i < tuple_size<row_type>::value - 1; ++i)
    {
        fout << theTable.labels[i] << "\t";
    }
    fout << theTable.labels[tuple_size<row_type>::value - 1] << endl;
    fout.close();
    TableWrite<tuple_size<row_type>::value - 1,
               cols_type, size_t>::print(theTable.theData, filename);
}
// Below is endif of DEBUG
#endif

#endif


