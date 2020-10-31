#ifndef MOLBIOLIB_TABLEINDEXMAP_H
#define MOLBIOLIB_TABLEINDEXMAP_H 1

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


// Below, we do not use an unordered_map, since we do not know
// what type of key is being hashed.


/** \file TableIndexMap.hpp
 * Contains the #TableIndexMap and #TableColumnIndexMap classes
 * that maps keys to row indices (returned as a one-column #Table).
 */

/** \file TableIndexMap.hpp
 * \fn template<typename... T> bool tableIndexMapAlwaysGoodRow(typename Table<T...>::row_type& theRow)
 * The function #tableIndexMapAlwaysGoodRow always returns true.
 */
template<typename... T>
bool tableIndexMapAlwaysGoodRow(typename Table<T...>::row_type& theRow)
{
    return true;
}


/** The class #TableIndexMap is used to map keys to rows (via indices) in a
 * Table object.  The template parameters must be passed in.  The compiler
 * does not try to figure out the types implicitly for classes.
 *
 * \section Usage Usage examples:
 * Suppose one has the table
 * \code
 * #define TABLEROW string, string, unsigned long, string, string, unsigned long
 * Table<TABLEROW> exTable;
 * \endcode
 * and the user-defined functions
 * \code
 * string myKeyGenerator(Table<TABLEROW>::row_type& theRow) {
 * ...
 * }
 * bool goodRow(Table<TABLEROW>::row_type& theRow) {
 * ...
 * }
 * \endcode
 * where <code>myKeyGenerator</code> returns a key based on a row's contents.
 * Note that one can extract a particular entry (say the fourth) from a row by
 * doing <code>get<5>(theRow)</code>.  <code>goodRow</code> returns
 * <code>true</code> or <code>false</code> depending on whether the row is to
 * be added (<code>true</code>) to the mappings
 * or not.  (For example, one may wish to
 * exclude rows where the string has rRNA in it.)  <code>goodRow</code> is
 * optional.  One need not define it and pass it below, in which case all rows
 * are considered good.
 *
 * Now, one may construct the mapping via
 * \code
 * TableIndexMap<string, TABLEROW> exTableMap(exTable, myKeyGenerator, goodRow);
 * \endcode
 * The map is automatically generated.  One may similar to above by replacing
 * <code>string</code> with any other type or class, as long as it has an
 * less-than operator associated with it.  (For this reason, the
 * <code>all</code> data type cannot be supported.  It does not even have an
 * equality comparison operator.)
 *
 * If one wishes to use a specific column as the key, one may use the
 * class as shown below (using the <em>third</em> column):
 * \code
 * TableColumnIndexMap<2, TABLEROW> exTableColMap(exTable);
 * \endcode
 *
 * To retrieve the #Table of rows associated with a particular key, one may
 * do something like the following (again, one can replace <code>string</code>
 * with some other data type):
 * \code
 * string specificKey = "...";    // put specific string in place of ellipses
 * IndexTable myIndices = exTableMap[specificKey];
 * \endcode
 * or alternatively
 * \code
 * IndexTable myIndices = exTableMap.getMapping(specificKey);
 * \endcode
 * If for performance-sake, direct access to the mappings is needed (really
 * only for advanced programmers), the STL C++ map that resides within
 * the class is
 * \code
 * map<K, IndexTable> theMapping;
 * \endcode
 * where <code>K</code> is the template type passed in the declaration of the
 * class.
 *
 * One also define a function which specifies if a particular key should be
 * retrieved (<code>true</code>) or not.  For example (again, we are using
 * a <code>string</code> key, but this may be replaced by any other
 * comparable type):
 * \code
 * bool theseKeys(string theKey) {
 *    ...
 * }
 * ...
 * IndexTable mySpecifickeys = exTableMap.getMapping(theseKeys);
 * \endcode
 * Note that the above is exactly the same as the alternate form of getMapping
 * above but we pass in a function rather than a key.  <em>Warning: this
 * function is slow since it must go through each map and check!</em>
 *
 * One may clear the mapping by doing
 * \code
 * exTableMap.clear();
 * \endcode
 *
 * If the table has changed, the mapping must be regenerated as such.
 * \code
 * exTableMap.regenerateMap();
 * \endcode
 */
template<typename K, typename... T>
class TableIndexMap
{
public:

    // Constructor
    TableIndexMap<K, T...>(Table<T...>& initTable,
                           K (*theKeyGenerator)(typename Table<T...>::row_type& theRow),
                           bool (*theGoodRow)(typename Table<T...>::row_type& theRow),
                           bool generateInitialMap = true) :
        theTable(initTable),
        keyGenerator(theKeyGenerator),
        goodRow(theGoodRow)
    {
        if (generateInitialMap)
        {
            regenerateMap();
        }
    }


    // Constructor with no goodRow function passed.  Must do it this way
    // since the compiler complains if a default parameter [function] is used.
    TableIndexMap<K, T...>(Table<T...>& initTable,
                           K (*theKeyGenerator)(typename Table<T...>::row_type& theRow),
                           bool generateInitialMap = true) :
        theTable(initTable),
        keyGenerator(theKeyGenerator)
    {
        goodRow = tableIndexMapAlwaysGoodRow<T...>;
        if (generateInitialMap)
        {
            regenerateMap();
        }
    }





    IndexTable& getMapping(K theKey)
    {
        return theMapping[theKey];
    }



    IndexTable getMapping(bool (*goodKey)(K theKey))
    {
        IndexTable theIndices;
        for (typename map<K, IndexTable>::iterator i = theMapping.begin();
                i != theMapping.end(); ++i)
        {
            if (goodKey(i->first))
            {
                IndexTable& currentMap = theMapping[i->first];
                size_t numMaps = currentMap.size();
                for (size_t j = 0; j < numMaps; ++j)
                {
                    theIndices.push_back(currentMap[j]);
                }
            }
        }
        return theIndices;

    }


    IndexTable& operator[](K theKey)
    {
        return getMapping(theKey);
    }


    void regenerateMap()
    {
#ifdef PROG_DEBUG
        // In a PROG_DEBUG since sameRowSizes is only defined
        // if PROG_DEBUG is defined.
        assert(theTable.sameRowSizes() == true);
#endif
        theMapping.clear();
        size_t numRows = get<0>(theTable.theData).size();
        typename Table<T...>::row_type tempRow;
        for (size_t i = 0; i < numRows; ++i)
        {
            theTable.getRow(i, tempRow);
            if (goodRow(tempRow))
            {
                theMapping[keyGenerator(tempRow)].push_back(i);
            }
        }
    }



    Table<K> getKeys()
    {
        Table<K> theKeys;
        for(typename map<K, IndexTable>::iterator i = theMapping.begin(); i != theMapping.end(); ++i)
        {
            theKeys.push_back(i->first);
        }
        return theKeys;
    }



    void clear()
    {
        theMapping.clear();
    }



    map<K, IndexTable> theMapping;

protected:
    Table<T...>& theTable;
    K (*keyGenerator)(typename Table<T...>::row_type& theRow);
    bool (*goodRow)(typename Table<T...>::row_type& theRow);

};




template<size_t N, typename... T>
typename tuple_element<N, typename Table<T...>::row_type>::type tableColumnIndexMapKeyGenerator(typename Table<T...>::row_type& theRow)
{
    return get<N>(theRow);
}
template<typename... T>
bool tableColumnIndexMapGoodRow(typename Table<T...>::row_type& theRow)
{
    return true;
}
/** This is a special case of #TableIndexMap, where one specifies a particular
 * column to use as the key.  For usage, see #TableIndexMap.
 */
template<size_t N, typename... T>
class TableColumnIndexMap :
    public TableIndexMap<typename tuple_element<N, typename Table<T...>::row_type>::type, T...>
{
public:
    // Constructor - just passes tableColumnIndexMapKeyGenerator and
    // tableColumnIndexMapGoodRow to the constructor of TableIndexMap.
    TableColumnIndexMap<N,T...>(Table<T...>& initTable) : TableIndexMap<typename tuple_element<N, typename Table<T...>::row_type>::type, T...>(initTable, tableColumnIndexMapKeyGenerator<N,T...>, tableColumnIndexMapGoodRow<T...>)
    {
    }
};







#endif

