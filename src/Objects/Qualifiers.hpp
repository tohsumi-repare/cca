#ifndef MOLBIOLIB_QUALIFERS_H
#define MOLBIOLIB_QUALIFERS_H 1

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


// Going with a map below since we do not know the type that is being hashed.
// An unordered_map's iterator also does not support the --operator.

/** \file Qualifiers.hpp
 * Contains only the #Qualifiers class.
 */

/** Holds a vector of pair of strings - the key and value.
 * This class simply contains a vector of keys and values.
 * The template parameters are type key type then the value type.
 * The values are stored in a vector.
 *
 * \section Usage Usage:
 * \code
 * StringQualifiers theQualifiers;
 *    // StringQualifiers is just Qualifiers<string, string> typedef'd
 * theQualifiers.push_back("myKey", "Some value");
 * theQualifiers.push_back("myKey", "Some other value");
 * // Can check for a key
 * bool isThere = theQualifiers.hasKey("myKey");
 * StringTable keyValues = theQualifiers["myKey"];
 * string firstKeyValue = theQualifiers.getQualifier("myKey");
 *    // The above gets just the first qualifer value.
 * StringTable theKeys = theQualifiers.getKeys();
 * theQualifiers.clear();
 * \endcode
 *
 * One may iterator through the qualifers as such (again <code>string</code> may
 * be replaced by types of one's choice assuming the #Qualifiers is used):
 * \code
 * if (!theQualifiers.empty()) {
 *    for (StringQualifiers::iterator qualiferIterator = theQualifiers.begin();
 *         qualiferIterator != theQualifiers.end(); ++qualiferIterator) {
 *       pair<string, string> currentQualifier = *qualiferIterator;
 *       cout << "key = " << currentQualifier.first
 *            << "\tvalue = " << currentQualifier.second << endl;
 *    }
 * }
 * \endcode
 * Note that we have to check if the #Qualifiers is empty since the #Qualifiers
 * iterator cannot handle an empty qualifer (since it iterates on interleaved
 * data).
 */
template<typename K, typename F>
class Qualifiers
{
public:
    void push_back(K theKey, F theQualifier)
    {
        theQualifiers[theKey].push_back(theQualifier);
    }

    Table<F>& operator[](K theKey)
    {
#ifdef DEBUG
        assert(theQualifiers.find(theKey) != theQualifiers.end());
#endif
        return theQualifiers[theKey];
    }

    F getQualifier(K theKey)
    {
#ifdef DEBUG
        assert(theQualifiers.find(theKey) != theQualifiers.end());
#endif
        return theQualifiers[theKey][0];
    }

    Table<K> getKeys()
    {
        Table<K> theResult;
        for (typename map< K, Table<F> >::iterator i = theQualifiers.begin();
                i != theQualifiers.end(); ++i)
        {
            theResult.push_back(i->first);
        }
        return theResult;
    }


    bool hasKey(K theKey)
    {
        if (theQualifiers.find(theKey) != theQualifiers.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }


    bool empty()
    {
        return (theQualifiers.size() == 0);
    }


    void clear()
    {
        // The bottom should call the destructor on each Table, which
        // is the [Table.]clear(true).
        theQualifiers.clear();
    }


    class iterator
    {
    public:
        iterator(map< K, Table<F> >* theMap,
                 typename map< K, Table<F> >::iterator mapI,
                 typename Table<F>::iterator tableI) :
            theMapping(theMap), mapIterator(mapI), tableIterator(tableI)
        {
        }

        iterator& operator=(const iterator& other)
        {
            theMapping = other.theMapping;
            mapIterator = other.mapIterator;
            tableIterator = other.tableIterator;
            return (*this);
        }

        bool operator==(const iterator& other)
        {
            if (mapIterator == theMapping->end())
            {
                return (mapIterator == other.mapIterator);
            }
            else
            {
                return ((mapIterator == other.mapIterator) &&
                        (tableIterator == other.tableIterator));
            }
        }

        bool operator!=(const iterator& other)
        {
            if (mapIterator == theMapping->end())
            {
                return (mapIterator != other.mapIterator);
            }
            else
            {
                return ((mapIterator != other.mapIterator) ||
                        (tableIterator != other.tableIterator));
            }
        }

        iterator& operator++()
        {
            if (mapIterator != theMapping->end())
            {
                ++tableIterator;
                if (tableIterator == (*theMapping)[mapIterator->first].end())
                {
                    ++mapIterator;
                    if (mapIterator != theMapping->end())
                    {
                        tableIterator = (*theMapping)[mapIterator->first].begin();
                    }
                }
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
            if (tableIterator == (*theMapping)[mapIterator->first].begin())
            {
                --mapIterator;
                if (mapIterator != theMapping->end())     // empty map
                {
                    tableIterator = (*theMapping)[mapIterator->first].end();
                    --tableIterator;
                }
            }
            else
            {
                --tableIterator;
            }
            return (*this);
        }

        iterator operator--(int)
        {
            iterator tmp(*this);
            --(*this);
            return tmp;
        }

        pair<K, F> operator*()
        {
            pair<K, F> theQualifier;
            theQualifier.first = mapIterator->first;
            theQualifier.second = get<0>(*tableIterator);
            return theQualifier;
        }


    private:
        map< K, Table<F> >* theMapping;
        typename map< K, Table<F> >::iterator mapIterator;
        typename Table<F>::iterator tableIterator;
    };


    iterator begin()
    {
#ifdef DEBUG
        // We cannot handle an empty set of qualifers.
        assert(theQualifiers.size() != 0);
#endif
        return iterator(&theQualifiers, theQualifiers.begin(),
                        (theQualifiers.begin())->second.begin());
    }
    iterator end()
    {
#ifdef DEBUG
        // We cannot handle an empty set of qualifers.
        assert(theQualifiers.size() != 0);
#endif
        // We have to give some iterator for the table,
        // so go for the only known valid one.
        return iterator(&theQualifiers, theQualifiers.end(),
                        (theQualifiers.begin())->second.begin());
    }


protected:
    map< K, Table<F> > theQualifiers;

};



/** \file Qualifiers.hpp
 * <code>Qualifiers<string, string></code> is used so often that
 * <code>StringQualifiers</code> is designated.  Functions to get
 * a tab-separated value of keys and values (multiple values comma-separated)
 * are also provided.
 * \code
 * StringQualifiers theQualifiers;
 * // ...
 * string theKeys, theValues;
 * theQualifiers.getStringKeys(theKeys);
 * // or alternatively, theKeys = theQualifiers.getStringKeys();
 * theQualifiers.getStringValues(theValues);
 * // or alternatively, theValues = theQualifiers.getStringValues();
 * \endcode
 * One can also get the reverse of the qualifiers via
 * \code
 * theQualifiers.generateReverseQualifiers();
 * StringQualifiers::reverseQualifier_type& reversedQualifiers = theQualifiers.getReverseQualifiers();
 * \endcode
 * The values are concatenated in the same order as the keys.
 */
class StringQualifiers : public Qualifiers<string, string>
{
public:
    typedef map< string, Table<string> > reverseQualifier_type;

    void getStringKeys(string& theKeys)
    {
        StringTable theResult = this->getKeys();
        size_t numResults = theResult.size();
        theKeys = "";
        for (size_t i = 0; i < numResults; ++i)
        {
            theKeys += theResult[i];
            if (i < numResults - 1)
            {
                theKeys += "\t";
            }
        }
    }
    string getStringKeys()
    {
        string theKeys;
        getStringKeys(theKeys);
        return theKeys;
    }

    void getStringValues(string& theValues)
    {
        StringTable theKeys = this->getKeys();
        size_t numKeys = theKeys.size();
        theValues = "";
        for (size_t i = 0; i < numKeys; ++i)
        {
            StringTable& currValues = theQualifiers[theKeys[i]];
            size_t numValues = currValues.size();
            string tempValues = "";
            for (size_t j = 0; j < numValues; ++j)
            {
                tempValues += currValues[j];
                if (j < numValues - 1)
                {
                    tempValues += ",";
                }
            }
            theValues += tempValues;
            if (i < numKeys - 1)
            {
                theValues += "\t";
            }
        }
    }
    string getStringValues()
    {
        string theValues;
        getStringValues(theValues);
        return theValues;
    }

    void generateReverseQualifiers(reverseQualifier_type& reverseQualifiers)
    {
        if (!empty())
        {
            for (StringQualifiers::iterator qualiferIterator = begin();
                    qualiferIterator != end(); ++qualiferIterator)
            {
                pair<string, string> currentQualifier = *qualiferIterator;
                reverseQualifiers[currentQualifier.second].push_back(currentQualifier.first);
            }
        }
    }

    void generateReverseQualifiers()
    {
        if (!empty())
        {
            generateReverseQualifiers(theReverseQualifiers);
            generatedReverseQualifiers = true;
        }
    }

    reverseQualifier_type& getReverseQualifiers()
    {
        if (!generatedReverseQualifiers)
        {
#ifdef DEBUG
            cerr << "Warning!  In Qualifiers.hpp, forced to generate the reversed qualifiers.\n";
#endif
            generateReverseQualifiers();
        }
        return theReverseQualifiers;
    }

    void clear()
    {
        theQualifiers.clear();
        theReverseQualifiers.clear();
    }

protected:
    reverseQualifier_type theReverseQualifiers;

private:
    bool generatedReverseQualifiers;


};



/** \file Qualifiers.hpp
 * The below is an instantiation of a null string qualifers.  This is so it
 * can be passed as a default #StringQualifiers parameter.  Cannot be a const
 * since passed into pass-by-reference parametes in Broad alignment readers.
 */
StringQualifiers nullStringQualifiers;

#endif

