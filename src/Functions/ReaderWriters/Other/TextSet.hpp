#ifndef MOLBIOLIB_TEXTSET_H
#define MOLBIOLIB_TEXTSET_H 1

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
#include "src/Objects/TextFileTypes.hpp"



/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename S, typename T, size_t N = 0> void readTextSet(set<T>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = -1, size_t minStringReserve = 0)
 * Read text set.
 * \section Usage Usage examples:
 * \code
 * set<double> mySet;
 * readTextSet<CSVType, double, 1>(mySet, "newfile.csv", true);
 * \endcode
 * where <code>CSVType</code> is a predefined class that describes how
 * comma-separated-valued sets are to be read in.  The
 * predefined classes currently available are:
 * - CSVType
 *   - Comma separated valued sets.
 * - TSVType
 *   - Tab separated valued files.
 * - StringType
 *   - Single string valued table (StringType).
 * These types are defined in TextFileTypes.hpp.  The 1 says to use the
 * second column of the file.  If no number is specified, then the first column
 * (0) is used.  The <code>true</code> in the above
 * means there is no header line.  (Default
 * is true.)  If one wishes to read only a subset of the tables, one can do
 * \code
 * readTextSet<TSVType, double, 1>(mySet, "newfile.tsv", false, 4, 99);
 *    // To read in the 5th through 100th lines.  If one wishes to read to the
 *    // end of the file, specify -1 as the second argument.
 *    // Also, can replace readTextSet<TSVType, double, 1> by
 *    // readTSVSet<double, 1>.
 * \endcode
 * If the trailing column (above 1) is omitted, it is assumed to be read
 * from the first column.  For long lines, one may need to specify a
 * string reserve size:
 * \code
 * readTextSet<CSVType, double>(mySet, "newfile.csv", true, 0, -1, 200000);
 *    // Can replace readTextSet<CSVType, double> by readCSVSet<double>.
 * \endcode
 * where the above 200000 is the reserve size.  In this case, one must
 * specify the number of entries to reserve beforehand.  There is also the
 * type
 * \code
 * readTextSet<SSVType, double>(mySet, "newfile.csv", true, 0, -1, 200000);
 *    // Can replace readTextSet<SSVType, double> by readSSVSet<double>.
 * \endcode
 * which reads in space-delimeted files.  Note that there is the special case:
 * \code
 * set<string> sSet;
 * readTextSet<StringType>(sSet, "newfile.txt");
 * // or readStringSet(sSet, "newfile.txt");
 * \endcode
 * where the above parameters also apply.  Here, no column number can be
 * specified.
 */
template<typename Q, typename T, size_t N = 0>
void readTextSet(set<T>& theSet, string filename, bool noLabels = true,
                 size_t first = 0,
                 size_t last = numeric_limits<size_t>::max(),
                 size_t minStringReserve = 0)
{

    // Check to see the file exists.
    if (fileSize(filename) == static_cast<ifstream::pos_type>(-1))
    {
        cerr << "Error!  The file " << filename << " cannot be opened." << endl;
        assert(fileSize(filename) != static_cast<ifstream::pos_type>(-1));
    }
    Ifstream fin(filename);
    string line;
    Q splitter;
    if (minStringReserve > 0)
    {
        line.reserve(minStringReserve);
    }
    if (!noLabels)
    {
        getline(fin, line);
    }
    size_t currentRow = 0;
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
                if (last != numeric_limits<size_t>::max() &&
                        currentRow > last)
                {
                    proceed = false;
                }
                if (proceed)
                {
                    vector<string> tokens;
                    splitter.splitLine(line, tokens);
                    theSet.insert(convertFromString<T>(tokens[N]));
                }
            }
            ++currentRow;
        }
    }
    fin.close();
}

/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename T, size_t N = 0> void readCSVSet(set<T>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t minStringReserve = 0)
 * Read comma-separated-valued (CSV) sets.  This is just #readTextSet with
 * type #CSVType passed.
 */
template<typename T, size_t N = 0>
void readCSVSet(set<T>& theSet, string filename,
                bool noLabels = true,
                size_t first = 0,
                size_t last = numeric_limits<size_t>::max(),
                size_t minStringReserve = 0)
{
    readTextSet<CSVType, T, N>(theSet, filename, noLabels,
                               first, last, minStringReserve);
}


/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename... T> void readTSVSet(set<T>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t minStringReserve = 0)
 * Read tab-separated-valued (TSV) file set.  This is just #readTextSet with
 * type #TSVType passed.
 */
template<typename T, size_t N = 0>
void readTSVSet(set<T>& theSet, string filename,
                bool noLabels = true,
                size_t first = 0,
                size_t last = numeric_limits<size_t>::max(),
                size_t minStringReserve = 0)
{
    readTextSet<TSVType, T, N>(theSet, filename, noLabels,
                               first, last, minStringReserve);
}



/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename... T> void readSSVSet(set<T>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t minStringReserve = 0)
 * Read space-separated-valued (SSV) file set.  This is just #readTextSet with
 * type #SSVType passed.
 */
template<typename T, size_t N = 0>
void readSSVSet(set<T>& theSet, string filename,
                bool noLabels = true,
                size_t first = 0,
                size_t last = numeric_limits<size_t>::max(),
                size_t minStringReserve = 0)
{
    readTextSet<SSVType, T, N>(theSet, filename, noLabels,
                               first, last, minStringReserve);
}




/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename... T> void readSpacesSet(set<T>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t minStringReserve = 0)
 * Read spaces-separated-valued file set.  This is just #readTextSet with
 * type #SpacesType passed.
 */
template<typename T, size_t N = 0>
void readSpacesSet(set<T>& theSet, string filename,
                   bool noLabels = true,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   size_t minStringReserve = 0)
{
    readTextSet<SpacesType, T, N>(theSet, filename, noLabels,
                                  first, last, minStringReserve);
}




/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename... T> void readStringSet(set<string>& theSet, string filename, bool noLabels = true, size_t first = 0, size_t last = numeric_limits<size_t>::max(), size_t minStringReserve = 0)
 * Read a string set.  This is just #readTextSet with
 * type #StringType passed.
 */
void readStringSet(set<string>& theSet, string filename,
                   bool noLabels = true,
                   size_t first = 0,
                   size_t last = numeric_limits<size_t>::max(),
                   size_t minStringReserve = 0)
{
    readTextSet<StringType, string, 0>(theSet, filename, noLabels,
                                       first, last, minStringReserve);
}


/** \file TextSet.hpp
 * Read and write text sets.
 * \fn template<typename T> void writeTextSet(set<T>& theSet, string filename)
 * Write set into a file.
 * \section Usage Usage examples:
 * \code
 * set<double> mySet;
 * \endcode
 * \code
 * writeTextSet(mySet, "newfile.txt");
 * \endcode
 */
template<typename T>
void writeTextSet(set<T>& theSet, string filename)
{
    Ofstream fout(filename);
    for (typename set<T>::iterator i = theSet.begin(); i != theSet.end(); ++i)
    {
        fout << convertToString<T>(*i) << endl;
    }
    fout.close();
}



/** \file TextSet.hpp
 * Read and write sets into text files.
 * \fn template<typename T> void printTextSet(set<T>& theSet, ostream& ofp = cout)
 * Print sets.
 * \section Usage Usage examples:
 * \code
 * set<size_t>
 * printTextSet(mySet);
 * // or ofstream ofp(...);
 * //    printTextSet(mySet, ofp);
 * \endcode
 */
template<typename T>
void printTextSet(set<T>& theSet, ostream& ofp = cout)
{
    for (typename set<T>::iterator i = theSet.begin(); i != theSet.end(); ++i)
    {
        ofp << convertToString<T>(*i) << endl;
    }
}



#endif
