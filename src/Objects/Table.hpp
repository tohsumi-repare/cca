#ifndef MOLBIOLIB_TABLE_H
#define MOLBIOLIB_TABLE_H 1

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


#include "src/include/PrimitiveTypes.hpp"
// Below is for getting a Table row by a variadic set of arguments.
#include "src/Functions/Transformers/TupleOperations/GetTupleValues.hpp"



#ifndef COL_TYPE_TABLE_TYPENAME
/** \file Table.hpp
 * Contains the #Table class and associated typedefs and DEFINEs.
 * \def COL_TYPE_TABLE_TYPENAME(N, TABLE)
 * Define macro to get the table column type with typenames for use in code.
 * Typically, one uses COL_TYPE_TABLE for application code and
 * COL_TYPE_TABLE_TYPENAME for within methods.
 * \todo Change to templated typedefs once it is implemented in g++.
 */
#define COL_TYPE_TABLE_TYPENAME(N, TABLE) typename tuple_element<N, typename TABLE::cols_type>::type
#endif

#ifndef COL_TYPE_TABLE
/** \file Table.hpp
 * Contains the #Table class and associated typedefs and DEFINEs.
 * \def COL_TYPE_TABLE(N, TABLE)
 * Define macro to get the table column type.
 * \todo Change to templated typedefs once it is implemented in g++.
 */
#define COL_TYPE_TABLE(N, TABLE) tuple_element<N, TABLE::cols_type>::type
#endif



/** The #Table is the main class on which everything in MolBioLib is based.
 * \section Usage Usage examples:
 * \code
 * typedef Table<string, string, unsigned long> myFastaTable;
 * myFastaTable myTable, myTable2;
 * Table<unsigned long, string> myNewTable;
 * ...
 * myFastaTable myTable3 = myTable2;   // copy constructor
 * ...
 * myTable3 = myTable;   // assignment copy constructor
 * \endcode
 * Alternatively,
 * \code
 * myTable3.assign(myTable);
 * myFastaTable::row_type myRow;
 * TableIndexType myRowNumber;   // the indexing of a vector
 *                               // uses this type.
 * \endcode
 * Alternatively for the above, one may also use
 * \code
 * myFastaTable::size_type myRowNumber;
 * // Can use instead of myFastaTable::size_type, size_t (of std namespace).
 * \endcode
 * One can initialize the #Table to a scalar value:
 * \code
 * // Create a table of length 25 and initialize each row to 1.25, blank, and 2.
 * Table<double, string, int> myTable4(25, 1.25, "blank", 2);
 * \endcode
 * If the #Table is of a single column, then one can assign a vector to the
 * #Table:
 * \code
 * vector<int> intVec;
 * ...
 * Table<int> intTable = intVec;
 * \endcode
 * <b> Because of a C++ standards inconsistency, Table<bool>, while can be
 * instantiated cannot use the operator[], because of a reference to the
 * bit vector (vector<bool> specialization).  Use getElement, if needed.</b>
 *
 * To set the labels of the table, create a string vector of size at least the
 * number of columns and pass it as such:
 * \code
 * vector<string> myNewLabels(3);
 * myNewLabels[0] = "first col";
 * myNewLabels[1] = "second col";
 * myNewLabels[2] = "third col";
 * myFastaTable.setLabels(myNewLabels);
 * \endcode
 *
 * To get the column type, one can do:
 * \code
 * COL_TYPE_TABLE(4, myFastaTable) myColumn;
 * \endcode
 * where #COL_TYPE_TABLE is a defined macro.  (Templated typedefs are not yet
 * implemented in g++.)  Yes, this is the one place where the column
 * numbering appears within parenthesis instead of between angle brackets.
 * Alternatively, one could do:
 * \code
 * tuple_element<4, myFastaTable::cols_type>::type myColumn;
 * \endcode
 * Alternatively, if one knew the type, such
 * as <code>int</code> for column 4, one could simply define
 * \code
 * vector<int> myColumn;
 * \endcode
 * <ul>
 *    <li>
 *       We would have liked to have
 *       implemented <code>myFastaTable::cols_type<4> myColumn;</code> but
 *       but unfortunately g++ does not yet implement templated typedefs.
 *    </li>
 * </ul>
 * Note that row_type, cols_type, and size_type are lower cased so that the
 * naming scheme conforms to the STL container convention.  Continuing,
 * \code
 * myTable.reserve(100);   // Reserve all vectors to this length.
 * myTable.resize(100);   // Resize all vectors to this length.
 * // If myTable has only one column, say of type string, then one may do
 * myTable.resize(100, "test");
 *    // where each element is set to "test".
 * myTable.size();   // Returns the number of rows in the table,
 *                   // assuming all rows of the same column.
 *    // One also use myTable.numRows();
 * \endcode
 * There are four ways to get a row:
 * \code
 * myTable.getRow(3, myRow);   // May not always be possible if have
 *                             // unequal number of rows.
 * \endcode
 * or the slower (due to a copy)
 * \code
 * myRow = myTable.getRow(3);
 * \endcode
 * One may pass the variable arguments directly:
 * \code
 * string outStr1 = "", outStr2 = "";
 * unsigned long i1 = 0;
 * // There is a setRow analog for the below where instead of getting values,
 * // it sets values.  The table must have allocated the appropriate number
 * // of rows.
 * myTable.getRow(3, outStr1, outStr2, i1);
 * \endcode
 * As noted in the comments above, there is an analogous <code>setRow</code>.
 * If one wishes to get just the first column's value (usually used when Table
 * is one column), then one may use the square brackets operator.  For example,
 * , if the Table's first column is a column of integers then
 * \code
 * int myNum = myTable[3];
 * \endcode
 * To get the number of columns (rarely needed, since one knows this from
 * the declaration), one can use the return value (type int) of
 * \code
 * myTable.numCols();
 * \endcode
 * To get a reference to a specific entry, one can do
 * \code
 *                  // 6th element of 1st column below
 * string aString = myFastaTable.getElement<0>(5);
 *                      // 10th element of 3rd column
 * unsigned long aNum = myFastaTable.getElement<2>(9);
 * //
 * aString = "test string";
 * myFastaTable.setElement<0>(3, aString);   // set 4th element of 1st column.
 * aNum = 50;
 * myFastaTable.setElement<2>(8, aNum);      // set 9th element of 3rd column.
 * \endcode
 *
 * The #Table class does have an iterator class within it.  It has the
 * assignment (=), equality (==), inequality (!=), and prefix and postfix
 * increments (++) and decrements (--).  The dereferencing operator (*)
 * returns a tuple (return by copy) of the current row.  Thus, one can
 * iterator through the table by:
 * \code
 * for (myFastaTable::iterator i = myTable.begin(); i != myTable.end(); ++i) {
 *    myFastaTable::row_type currentRow = *i;
 *    // One cannot do cout << *i << endl; below since g++ cannot
 *    // see that *i is of tuple type that has a stream associated with it.  :-(
 *    cout << currentRow << endl;
 * }
 * \endcode
 *
 * One can get the reference to a column via something like if one knew
 * the second column is of type <code>vector<int></code>, one could do
 * \code
 * vector<int>& myColumnReference = myTable.getColRef<1>();
 * \endcode
 * Alteratively, one could use
 * \code
 * COL_TYPE_TABLE(1, myFastaTable)& myColumnReference = myTable.getColRef<1>();
 * \endcode
 * <ul>
 *    <li>
 *       Alternatively to getColRef above, one can also use
 *       \code
 *       get<1>(myTable.theData)
 *       \endcode
 *       <b> Often, one will have to use this since the compiler may not be
 *       able to parse the code for the method getColRef.</b>
 *    </li>
 * </ul>
 * For getting columns,
 * \code
 * myTable.getCol<4>(myColumn);   // 4 is the last column, i.e. column 5.
 * \endcode
 * Similarly with setRow and setCol.
 * <ul>
 *    <li>
 *       However, one cannot do <code>myTable[3] = myRow</code> or
 *       <code>myTable3.setRow(3) = myRow</code> due to
 *       how the table is set up.  One must
 *       do <code>myTable3.setRow(3, myRow)</code>.
 *    </li>
 * </ul>
 * Optionally, with getCol, one may do:
 * \code
 * myTable.getCol<4>(myColumn, 24, 99, 2);   // Get rows 25 through
 *                                           // 100 (inclusive) only and
 *                                           // fill myColumn from the third
 *                                           // row.
 * \endcode
 * Alternatively (albeit slower and uses more memory due to a copy),
 * \code
 * myColumn = myTable.getCol<4>();
 * \endcode
 * and
 * \code
 * myColumn = myTable.getCol<4>(24, 99, 2);
 * \endcode
 * respectively.
 * \code
 * myTable.fillCol<4>(42, 9, 99);   // Sets the 5th column starting from row 10
 *                                  // through 100 inclusive to the value 42.
 *   // It may be necessary to do instead of above
 *   //    myTable.template fillCol<4>(42, 9, 99)
 *   // since the compiler may not be able to figure out the templating
 *   // if one is using a templated class derived from Table, giving an
 *   // unresolved overloaded function type error.  See
 *   // http://www.informit.com/guides/content.aspx?g=cplusplus&seqNum=388
 * myTable.setCol<4>(myColumn, 24);   // Sets 5th column, but start at row 25
 *                                    // in the #Table (myColumn is read from
 *                                    // the first row to the end).
 *                                    // There is no stop parameter for setCol.
 *                                    // All of myColumn is copied to myTable.
 * myTable.push_back(myRow);   // push new row at the end.
 * \endcode
 * For pushing back a set of arguments, one can do for example:
 * \code
 * myTable.push_back("string 1", "string 2", 3);
 * \endcode
 * If one wishes to push_back a value to the first column only, one may use
 * the push_back with a value.  For example, if the first column is an integer,
 * \code
 * myTable.push_back(3);
 * int n = 42;
 * myTable.push_back(n);
 * \endcode
 * If one wishes to push_back another table of the same type onto the current
 * table, one may do
 * \code
 * myTable.push_back(myTable2, 0, numeric_limits<size_t>::max());
 * \endcode
 * where the last two are optional with the defaults shown.  0 means start from
 * the first row of myTable2 and numeric_limits<size_t>::max() means go to the end of myTable2.
 * (Otherwise, the second parameter must be equal or larger to the second.)
 *
 * One may insert either a Table, a row, or an element similarly to push_back
 * by using the insert operator, as shown below.  If no index is provided, the
 * argument is inserted at the beginning of the Table.
 * \code
 * // Below inserts 7 below the third row.
 * myTable.insert(7, 2);
 * myTable.insert(myRow);   // inserts myRow at the beginning.
 * myTable.insert(myTable2, 3);  // inserts myTable2 into myTable before row 4.
 * \endcode
 * For a set of arguments, the insert_args passes the row first.  It is
 * <b>not</b> an optional argument in the below case.
 * \code
 * // Below inserts the row with the last three arguments before row 4.
 * myTable.insert_args(3, "string 3", "string 4", 42);
 * \endcode
 *
 * Continuing with the operators:
 * \code
 * myTable.erase(3);   // Row index 3 for all columns is deleted and all
 *                     // lower entries are moved up.
 * myTable.erase(3, 5);   // Row index 3 to 5 inclusive for all columns
 *                        // is deleted and all lower entries are moved up.
 * myTable.clearCol<4>();   // Column 4's vector is cleared (0-sized now).
 * myTable.clear();   // Clear all columns, but not the labels.
 *                    // Instead of clear(),
 *                    // use clear(true) to clear the labels.
 * \endcode
 * If one has a function
 * \code
 * template<typename... T>
 * void myRowFunct(size_t& theNum, tuple<T...>& theRow) { ... }
 * \endcode
 * or alternatively
 * \code
 * void myRowFunct(size_t& theNum, myFastaTable::row_type& theRow) { ... }
 * \endcode
 * then one could apply this function to all rows by doing
 * \code
 * myTable.applyToAllRows(myRowFunct);
 * \endcode
 * The above only works when all columns in the table have
 * the same number of rows.
 *
 *
 * If one wishes to apply some operation to all vectors, one must do it
 * manually due to how C++'s templates work, <i>e.g.</i>
 * \code
 * COL_TYPE_TABLE(0, myFastaTable)& myColumn0 = myTable.getColRef<0>();
 * \endcode
 * Alternatively, if one knows for example that the first column is a column
 * of strings, then one could do
 * \code
 * vector<string>& myColumn0 = myTable.getColRef<0>();
 * COL_TYPE_TABLE(1, myFastaTable)& myColumn1 = myTable.getColRef<1>();
 *   ...  // Do whatever operators to vectors myColumn0 and myColumn1.
 *        // Since they are references, no need to do a setCol-type thing.
 * \endcode
 * Alternatively, one could do something like
 * \code
 * COL_TYPE_TABLE(0, myFastaTable) myColumn0;
 * COL_TYPE_TABLE(1, myFastaTable) myColumn1;
 * myTable.getCol<0>(myColumn0);
 * myTable.getCol<1>(myColumn1);
 *   ...  // Do whatever operators to vectors myColumn0 and myColumn1.
 * myTable.setCol<0>(myColumn0);
 * myTable.setCol<1>(myColumn1);
 * \endcode
 * but this is much more computationally expensive and uses more memory.
 * We would have liked to written a function to do this, but at this time,
 * one cannot have a templated function as a parameter to a function.
 *
 * A type of table that is commonly used is
 * \code
 * Table<size_t> indexTable;
 * \endcode
 * Because of this, it is typedef'd to
 * \code
 * IndexTable myIndexTable;
 * \endcode
 * Another type of table that occurs often is
 * \code
 * Table<string> myStringTable;
 * \endcode
 * which is typedef'd to
 * \code
 * StringTable myStringTable;
 * \endcode
 *
 * \todo
 * Write functions to handle labels of tables better, like setting
 * individual ones and getting one or all and setting all.
 * Write a function to get and set columns once templated functions can be
 * passed as a parameter to a function.
 */
template<typename... T>
class Table
{
public:
    typedef tuple<vector<T>...> cols_type;
    cols_type theData;  // makes tuple<vector<T1>,vector<T2>,...>
    // Below is the size_type (and the index variable type) of the vector.
    typedef size_t size_type;
    typedef tuple<T...> row_type;
    vector<string> labels;

    // Not const below since assertion will not allow it.
    size_t size()
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
#endif
        return get<0>(theData).size();
    }

    size_t numRows() const
    {
        return size();
    }

    size_t numCols() const
    {
        return tuple_size<row_type>::value;
    }


    void setLabels(vector<string>& newLabels)
    {
#ifdef PROG_DEBUG
        assert(newLabels.size() >= labels.size());
#endif
        for (size_t i = 0; i < labels.size(); ++i)
        {
            labels[i] = newLabels[i];
        }
    }


    void getRow(size_t i, row_type& theRow)
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
        assert(i < get<0>(theData).size());
#endif
        TableGetRow<tuple_size<row_type>::value - 1,
                    row_type, cols_type>::getRow(theRow, theData, i);
    }
    row_type getRow(size_t i)
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
        assert(i < get<0>(theData).size());
#endif
        row_type theRow;
        getRow(i, theRow);
        return theRow;
    }
    void getRow(size_t i, T&... args)
    {
        row_type theRow;
        getRow(i, theRow);
        getTupleValues(theRow, args...);
    }

    // Below is lhs (of equals sign) operator[], when one has to set something.
    typename tuple_element<0,row_type>::type& operator[](size_t i)
    {
#ifdef PROG_DEBUG
        if (tuple_size<row_type>::value != 1)
        {
            cerr << "Warning.  Table::operator[] used on a table with multiple columns.\n";
        }
        assert(i < (get<0>(theData)).size());
#endif
        return ((get<0>(theData))[i]);
    }


    template<size_t N>
    typename tuple_element<N,row_type>::type getElement(size_t i)
    {
#ifdef PROG_DEBUG
        assert(i < (get<N>(theData)).size());
#endif
        return ((get<N>(theData))[i]);
    }


    template<size_t N>
    void setElement(size_t i, typename tuple_element<N, row_type>::type x)
    {
#ifdef PROG_DEBUG
        assert(i < (get<N>(theData)).size());
#endif
        (get<N>(theData))[i] = x;
    }


    void setRow(size_t i, row_type& theRow)
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
        assert(i < get<0>(theData).size());
#endif
        TableSetRow<tuple_size<row_type>::value - 1,
                    row_type, cols_type>::setRow(theRow, theData, i);
    }

    // Below is when we have all variables, so can do a direct reference.
    void setRow(size_t i, T&... args)
    {
        row_type tempRow = make_tuple(forward<T>(args)...);
        setRow(i, tempRow);
    }

    // Below is when we some or no variables, so do pass by copy.
    template<typename... U>
    void setRow(size_t i, U... args)
    {
        row_type tempRow = make_tuple(forward<U>(args)...);
        setRow(i, tempRow);
    }

    void push_back(row_type& theRow)
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
#endif
        TablePushRow<tuple_size<row_type>::value - 1,
                     row_type, cols_type>::pushRow(theRow, theData);
    }

    // Below is when we have all variables, so can do a direct reference.
    void push_back(T&... args)
    {
        row_type tempRow = make_tuple(forward<T>(args)...);
        push_back(tempRow);
    }

    // Below is when we some or no variables, so do pass by copy.
    template<typename... U>
    void push_back(U... args)
    {
#ifdef PROG_DEBUG
        if (tuple_size< tuple<U...> >::value != tuple_size<row_type>::value)
        {
            cerr << "Warning.  Table::push_back(U... args) used with arguments of different length than the number of columns in the table.  Continuing." << endl;
        }
#endif
        row_type tempRow = make_tuple(forward<U>(args)...);
        push_back(tempRow);
    }


    void push_back(Table<T...>& theTable,
                   size_t start = 0,
                   size_t end = numeric_limits<size_t>::max())
    {
        if (theTable.size() == 0)
        {
            return;
        }
        if (end == numeric_limits<size_t>::max())
        {
            end = theTable.size() - 1;
        }
#ifdef PROG_DEBUG
        assert(end >= start);
#endif
        row_type tempRow;
        for (size_t i = start; i <= end; ++i)
        {
            theTable.getRow(i, tempRow);
            push_back(tempRow);
        }
    }

    void insert(row_type& theRow, size_t insertBeforeRow = 0)
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
#endif
        TableInsertRow<tuple_size<row_type>::value - 1,
                       row_type, cols_type>::insertRow(theRow, theData, insertBeforeRow);
    }

    // Below is when we have all variables, so can do a direct reference.
    void insert_args(size_t insertBeforeRow, T&... args)
    {
        row_type tempRow = make_tuple(forward<T>(args)...);
        insert(tempRow, insertBeforeRow);
    }

    // Below is when we some or no variables, so do pass by copy.
    template<typename... U>
    void insert_args(size_t insertBeforeRow, U... args)
    {
        row_type tempRow = make_tuple(forward<U>(args)...);
        insert(tempRow, insertBeforeRow);
    }

    void insert(typename tuple_element<0,row_type>::type theElement,
                size_t insertBeforeRow = 0)
    {
#ifdef PROG_DEBUG
        if (tuple_size<row_type>::value != 1)
        {
            cerr << "Warning.  Table::insert(typename tuple_element<0,row_type>::type theElement, size_t insertBeforeRow = 0) used on a table with multiple columns.\n";
        }
#endif
        // Below is an iterator to the first element of the vector.
        typename tuple_element<0,cols_type>::type::iterator iter = get<0>(theData).begin();
        get<0>(theData).insert(iter + insertBeforeRow, theElement);
    }
    void insert(Table<T...>& theTable,
                size_t insertBeforeRow = 0,
                size_t start = 0,
                size_t end = numeric_limits<size_t>::max())
    {
        if (theTable.size() == 0)
        {
            return;
        }
        if (end == numeric_limits<size_t>::max())
        {
            end = theTable.size() - 1;
        }
        size_t origThisSize = size();
        size_t theTableSize = end - start + 1;
        resize(origThisSize + theTableSize);
        row_type tempRow;
        // First, copy old row down.
        if (size() > 0)
        {
            bool zeroExit = false;   // Needed since --i for i=0 will wraparound.
            for (size_t i = (origThisSize - 1); !zeroExit && (i >= insertBeforeRow); --i)
            {
                getRow(i, tempRow);
                setRow(i + theTableSize, tempRow);
                if (i == 0)
                {
                    zeroExit = true;
                }
            }
        }
        // Now, copy newTable into created space.
        size_t currTableIndex = 0;
        for (size_t i = start; i <= end; ++i, ++currTableIndex)
        {
            theTable.getRow(i, tempRow);
            setRow(insertBeforeRow + currTableIndex, tempRow);
        }
    }

    void erase(size_t i, size_t j = numeric_limits<size_t>::max())
    {
        if (j == numeric_limits<size_t>::max())
        {
            j = i;
        }
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
        assert(i < get<0>(theData).size());
        assert(j <= get<0>(theData).size());
#endif
        TableDelRows<tuple_size<row_type>::value - 1,
                     cols_type>::delRows(theData, i, j + 1);
    }

    template<size_t N>
    typename tuple_element<N,cols_type>::type& getColRef()
    {
        return get<N>(theData);
    }

    template<size_t N>
    void getCol(typename tuple_element<N,cols_type>::type& theCol,
                size_t min = 0,
                size_t max = numeric_limits<size_t>::max(),
                size_t offset = 0)
    {
        size_t colSize = 0;
        if (max == numeric_limits<size_t>::max())
        {
            colSize = get<N>(theData).size() - 1;
        }
        else
        {
            colSize = max;
        }
        if (theCol.size() < (offset + colSize + 1))
        {
            theCol.reserve(offset + colSize + 1);
        }
        for (size_t i = min; i <= colSize; ++i)
        {
            theCol[i-min+offset] = get<N>(theData)[i];
        }
    }
    template<size_t N>
    typename tuple_element<N,cols_type>::type getCol(size_t min = 0,
            size_t max = numeric_limits<size_t>::max(),
            size_t offset = 0)
    {
        typename tuple_element<N,cols_type>::type theCol;
        getCol<N>(theCol, min, max, offset);
        return theCol;
    }



    template<size_t N>
    void setCol(typename tuple_element<N,cols_type>::type& theCol,
                size_t start = 0)
    {
        size_t colSize = theCol.size();
        if (get<N>(theData).size() < colSize + start)
        {
            get<N>(theData).reserve(colSize + start);
        }
        for (size_t i = 0; i < colSize; ++i)
        {
            get<N>(theData)[i+start] = theCol[i];
        }
    }
    template<size_t N>
    void fillCol(typename tuple_element<N,row_type>::type fillValue,
                 size_t start = 0, size_t end = numeric_limits<size_t>::max())
    {

        size_t colSize = get<N>(theData).size();
        if (end != numeric_limits<size_t>::max())
        {
            colSize = end + 1;
        }
        for (size_t i = start; i < colSize; ++i)
        {
            get<N>(theData)[i] = fillValue;
        }
    }
    template<size_t N>
    void clearCol()
    {
        get<N>(theData).clear();
    }



    void reserve(size_t theSize)
    {
        TableReserve<tuple_size<row_type>::value - 1,
                     cols_type>::reserve(theData, theSize);
    }



    void resize(size_t theSize)
    {
        TableResize<tuple_size<row_type>::value - 1,
                    cols_type>::resize(theData, theSize);
    }


    void resize(size_t theSize, const typename tuple_element<0, Table<T...>::row_type>::type& elem)
    {
#ifdef PROG_DEBUG
        // Only use this if Table is acting like a vector.
        assert(tuple_size<row_type>::value == 1);
#endif
        get<0>(theData).resize(theSize);
        for (size_t i = 0; i < theSize; ++i)
        {
            get<0>(theData)[i] = elem;
        }
    }




    void applyToAllRows(void f(size_t& currentRow, row_type& theRow))
    {
#ifdef PROG_DEBUG
        assert(sameRowSizes() == true);
#endif
        size_t numRows = get<0>(theData).size();
        row_type tempRow;
        for (size_t i = 0; i < numRows; ++i)
        {
            getRow(i, tempRow);
            f(i, tempRow);
            setRow(i, tempRow);
        }
    }


    // Equality and inequality comparator - all other types not well defined.
    // Below assumes vector has operator==
    bool operator==(const Table<T...>& rhs)
    {
        if (labels != rhs.labels)
        {
            return false;
        }
        bool result = true;
        TableCompareData<tuple_size<row_type>::value - 1,
                         cols_type>::compareData(theData, rhs.theData,
                                 result);
        return result;
    }

    bool operator!=(const Table<T...>& rhs)
    {
        return (!(operator==(rhs)));
    }




    // Constructor
    Table<T...>()
    {
        labels.resize(tuple_size<row_type>::value, "");
    }


    // Constructor with size
    Table<T...>(size_t newSize)
    {
        labels.resize(tuple_size<row_type>::value, "");
        resize(newSize);
    }


    // Constructor with a scalar type and size
    Table<T...>(size_t newSize, T... args)
    {
        labels.resize(tuple_size<row_type>::value, "");
        resize(newSize);
        row_type tempRow = make_tuple(forward<T>(args)...);
        for (size_t i = 0; i < newSize; ++i)
        {
            setRow(i, tempRow);
        }
    }


    // Assignment
    void assign(const Table<T...>& oldTable)
    {
        labels.resize(tuple_size<row_type>::value);
        TableCopyData<tuple_size<row_type>::value - 1,
                      cols_type>::copyData(theData, oldTable.theData,
                                           labels, oldTable.labels);
        // for (size_t i = 0; i < tuple_size<row_type>::value; ++i)
        //    labels[i] = oldTable.labels[i];
        // theData = oldTable.theData;
    }


    // Copy constructor
    Table<T...>(const Table<T...>& copy)
    {
        assign(copy);
    }


    // Copy constructor in case of vector
    Table<T...>(const vector<typename tuple_element<0,Table<T...>::row_type>::type>& rhs)
    {
#ifdef PROG_DEBUG
        if (tuple_size<row_type>::value != 1)
        {
            cerr << "Warning.  Table::Table(const vector<typename tuple_element<0,Table<T...>::row_type>::type>& rhs) used on a table with multiple columns.\n";
        }
#endif
        size_t rhsSize = rhs.size();
        get<0>(theData).resize(rhsSize);
        for (size_t i = 0; i < rhsSize; ++i)
        {
            get<0>(theData)[i] = rhs[i];
        }
    }


    // Copy assignment
    // Need both const and non-const versions to satisfy the compiler.
    // copyData really does just const for the rhs, so no problems.
    Table<T...>& operator=(Table<T...>& rhs)
    {
        TableCopyData<tuple_size<row_type>::value - 1,
                      cols_type>::copyData(theData, rhs.theData,
                                           labels, rhs.labels);
        return (*this);
    }
    Table<T...>& operator=(const Table<T...>& rhs)
    {
        TableCopyData<tuple_size<row_type>::value - 1,
                      cols_type>::copyData(theData, rhs.theData,
                                           labels, rhs.labels);
        return (*this);
    }



    // Copy assignment in case of vector
    Table<T...>& operator=(vector<typename tuple_element<0,Table<T...>::row_type>::type>& rhs)
    {
#ifdef PROG_DEBUG
        if (tuple_size<row_type>::value != 1)
        {
            cerr << "Warning.  Table::operator=(vector<typename tuple_element<0,Table<T...>::row_type>::type>& rhs) used on a table with multiple columns.\n";
        }
#endif
        size_t rhsSize = rhs.size();
        theData.resize(rhsSize);
        // Use vector's operator=
        get<0>(theData) = rhs;
        if (labels.size() == 0)
        {
            labels.resize(tuple_size<row_type>::value, "");
        }
        return (*this);
    }



    // clear all
    void clear(bool clearLabels = false)
    {
        if (clearLabels)
        {
            size_t numLabels = labels.size();
            for (size_t i = 0; i < numLabels; ++i)
            {
                labels[i] = "";
            }
        }
        TableClearer<tuple_size<row_type>::value - 1,
                     cols_type>::clearTable(theData);
    }


    // Destructor
    ~Table()
    {
        clear();
        labels.clear();
    }


    class iterator
    {
    public:
        iterator(Table<T...>* startTable, size_t startIndex) : theTable(startTable), currentIndex(startIndex) { }
        ~iterator() { }

        iterator& operator=(const iterator& other)
        {
            theTable = other.theTable;
            currentIndex = other.currentIndex;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            return ((theTable == other.theTable) &&
                    (currentIndex == other.currentIndex));
        }

        bool operator!=(const iterator& other)
        {
            return ((theTable != other.theTable) ||
                    (currentIndex != other.currentIndex));
        }

        iterator& operator++()
        {
            if (currentIndex != theTable->size())
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
            return (theTable->getRow(currentIndex));
        }


    private:
        Table<T...>* theTable;
        size_t currentIndex;
    };


    iterator begin()
    {
        return iterator(this, 0);
    }


    iterator end()
    {
        return iterator(this, size());
    }


#ifdef PROG_DEBUG
    bool sameRowSizes()
    {
        if (tuple_size<row_type>::value == 0)
        {
            return false;
        }
        if (tuple_size<row_type>::value == 1)
        {
            return true;
        }
        vector<size_t> rowSizes(tuple_size<row_type>::value, 0);
        TableRowSizes<tuple_size<row_type>::value - 1,
                      cols_type>::rowSizes(theData, rowSizes);
        // Need to create below variable because otherwise clang++ complains
        // about a warning: comparison of unsigned expression < 0.
        // This is obviously a bug in clang++'s parsing.
        size_t currTupleSize = tuple_size<row_type>::value;
        for (size_t i = 0; i < (currTupleSize - 1); ++i)
            // for (size_t i = 0; i < (tuple_size<row_type>::value - 1); ++i)
        {
            if (rowSizes[i] != rowSizes[i+1])
            {
                return false;
            }
        }
        return true;
    }
#endif


private:
    // Tuple templates have parameters of type size_t


    /** \file Table.hpp
     *\todo We would have liked to use partial template specialization for
     * functions, but they do not exist yet.
     */

    // Unfortunately, these definitions have to go here since the
    // export [template ...] keyword needed to separate the below
    // is implemented by only a few compilers.
    template<size_t N, typename R, typename C>
    class TableGetRow
    {
    public:
        static void getRow(R& theRow, C& theData, size_t& i)
        {
            TableGetRow<N-1,R,C>::getRow(theRow, theData, i);
            get<N>(theRow) = get<N>(theData)[i];
        }
    };
    template<typename R, typename C>
    class TableGetRow<0,R,C>
    {
    public:
        static void getRow(R& theRow, C& theData, size_t& i)
        {
            get<0>(theRow) = get<0>(theData)[i];
        }
    };

    template<size_t N, typename R, typename C>
    class TableSetRow
    {
    public:
        static void setRow(R& theRow, C& theData, size_t& i)
        {
            TableSetRow<N-1,R,C>::setRow(theRow, theData, i);
            get<N>(theData)[i] = get<N>(theRow);
        }
    };
    template<typename R, typename C>
    class TableSetRow<0,R,C>
    {
    public:
        static void setRow(R& theRow, C& theData, size_t& i)
        {
            get<0>(theData)[i] = get<0>(theRow);
        }
    };

    template<size_t N, typename R, typename C>
    class TablePushRow
    {
    public:
        static void pushRow(R& theRow, C& theData)
        {
            TablePushRow<N-1,R,C>::pushRow(theRow, theData);
            get<N>(theData).push_back(get<N>(theRow));
        }
    };
    template<typename R, typename C>
    class TablePushRow<0,R,C>
    {
    public:
        static void pushRow(R& theRow, C& theData)
        {
            get<0>(theData).push_back(get<0>(theRow));
        }
    };

    template<size_t N, typename R, typename C>
    class TableInsertRow
    {
    public:
        static void insertRow(R& theRow, C& theData, size_t insertBeforeRow)
        {
            TableInsertRow<N-1,R,C>::insertRow(theRow, theData, insertBeforeRow);
            typename tuple_element<N,cols_type>::type::iterator iter = get<N>(theData).begin();
            get<N>(theData).insert(iter + insertBeforeRow, get<N>(theRow));
        }
    };
    template<typename R, typename C>
    class TableInsertRow<0,R,C>
    {
    public:
        static void insertRow(R& theRow, C& theData, size_t insertBeforeRow)
        {
            typename tuple_element<0,cols_type>::type::iterator iter = get<0>(theData).begin();
            get<0>(theData).insert(iter + insertBeforeRow, get<0>(theRow));
        }
    };

    template<size_t N, typename C>
    class TableDelRows
    {
    public:
        static void delRows(C& theData, size_t i, size_t j)
        {
            TableDelRows<N-1,C>::delRows(theData, i, j);
            if (j < get<N>(theData).size())
                get<N>(theData).erase(get<N>(theData).begin() + i,
                                      get<N>(theData).begin() + j);
            else
                get<N>(theData).erase(get<N>(theData).begin() + i,
                                      get<N>(theData).end());
        }
    };
    template<typename C>
    class TableDelRows<0,C>
    {
    public:
        static void delRows(C& theData, size_t i, size_t j)
        {
            if (j < get<0>(theData).size())
                get<0>(theData).erase(get<0>(theData).begin() + i,
                                      get<0>(theData).begin() + j);
            else
                get<0>(theData).erase(get<0>(theData).begin() + i,
                                      get<0>(theData).end());
        }
    };




    template<size_t N, typename C>
    class TableReserve
    {
    public:
        static void reserve(C& theData, size_t& theSize)
        {
            TableReserve<N-1,C>::reserve(theData, theSize);
            get<N>(theData).reserve(theSize);
        }
    };
    template<typename C>
    class TableReserve<0,C>
    {
    public:
        static void reserve(C& theData, size_t& theSize)
        {
            get<0>(theData).reserve(theSize);
        }
    };


    template<size_t N, typename C>
    class TableResize
    {
    public:
        static void resize(C& theData, size_t& theSize)
        {
            TableResize<N-1,C>::resize(theData, theSize);
            get<N>(theData).resize(theSize);
        }
    };
    template<typename C>
    class TableResize<0,C>
    {
    public:
        static void resize(C& theData, size_t& theSize)
        {
            get<0>(theData).resize(theSize);
        }
    };



    template<size_t N, typename C>
    class TableCopyData
    {
    public:
        static void copyData(C& theData, const C& rhsData, vector<string>& labels, const vector<string>& rhslabels)
        {
            get<N>(theData) = get<N>(rhsData);
            labels[N] = rhslabels[N];
            TableCopyData<N-1,C>::copyData(theData, rhsData, labels, rhslabels);
        }
    };
    template<typename C>
    class TableCopyData<0,C>
    {
    public:
        static void copyData(C& theData, const C& rhsData, vector<string>& labels, const vector<string>& rhslabels)
        {
            get<0>(theData) = get<0>(rhsData);
            labels[0] = rhslabels[0];
        }
    };



    template<size_t N, typename C>
    class TableCompareData
    {
    public:
        // It is assumed that result is initially true on the first call to
        // compareData.  We only set the result to false if something
        // is not equal.
        static void compareData(C& theData, const C& rhsData, bool& result)
        {
            if (get<N>(theData) == get<N>(rhsData))
            {
                TableCompareData<N-1,C>::compareData(theData, rhsData, result);
            }
            else
            {
                result = false;
            }
        }
    };
    template<typename C>
    class TableCompareData<0,C>
    {
    public:
        static void compareData(C& theData, const C& rhsData, bool& result)
        {
            if (!(get<0>(theData) == get<0>(rhsData)))
            {
                result = false;
            }
        }
    };



    template<size_t N, typename C>
    class TableClearer
    {
    public:
        static void clearTable(C& theData)
        {
            get<N>(theData).clear();
            TableClearer<N-1,C>::clearTable(theData);
        }
    };
    template<typename C>
    class TableClearer<0,C>
    {
    public:
        static void clearTable(C& theData)
        {
            get<0>(theData).clear();
        }
    };




#ifdef PROG_DEBUG
    template<size_t N, typename C>
    class TableRowSizes
    {
    public:
        static void rowSizes(C& theData, vector<size_t>& theSizes)
        {
            theSizes[N] = get<N>(theData).size();
            TableRowSizes<N-1,C>::rowSizes(theData, theSizes);
        }
    };
    template<typename C>
    class TableRowSizes<0,C>
    {
    public:
        static void rowSizes(C& theData, vector<size_t>& theSizes)
        {
            theSizes[0] = get<0>(theData).size();
        }
    };
#endif

};



/** \file Table.hpp
 * \var typedef Table<size_t> IndexTable;
 * Define the #IndexTable, since it is used so often.
 */
typedef Table<size_t> IndexTable;

/** \file Table.hpp
 * \var typedef Table<string> StringTable;
 * Define the #StringTable, since it is used so often.
 */
typedef Table<string> StringTable;

#endif

