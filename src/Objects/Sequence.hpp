#ifndef MOLBIOLIB_SEQUENCE_H
#define MOLBIOLIB_SEQUENCE_H 1

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


#include "src/Objects/PrimarySequence.hpp"


/** \file Sequence.hpp
 * Contains the #Sequence class and associated enum.
 * \enum SeqType
 * #SeqType type currently just <code>DNA</code>,
 * <code>RNA</code>, and <code>Protein</code> to indicate the
 * alphabet to be used.
 * \todo Replace by strongly typed enums when they work.  Instead of
 * <code>enum SeqType</code>, do <code>enum class SeqType</code>.
 */
enum SeqType
{
    DNA, RNA, Protein
};



/** Extends #PrimarySequence with type and type-specific methods.
 * This is an extension of a #PrimarySequence where it include additional
 * information, such as sequence
 * type (<code>DNA</code>, <code>RNA</code>, <code>Protein</code>).
 * \section Usage Usage:
 * \code
 * Sequence myDNA("ACTGACGTA"), mySecondDNA;
 * \endcode
 *
 * <b>*** Important functionality: there are many methods, such as find,
 * replace, etc., which are all part of #PrimarySequence, which is derived
 * from the standard string class.  #Sequence is derived from
 * #PrimarySequence.  Thus, one should also see the documentation for
 * #PrimarySequence as well as the standard C++ string class. ***</b>
 *
 * <i>By default, a #Sequence type is designated as <code>DNA</code>.</i>
 * To get the type, one does
 * \code
 * SeqType whichOne = myDNA.getType();
 * \endcode
 * Similarly, one could set a type, <i>e.g.</i>
 * \code
 * myDNA.setType(RNA);
 * \endcode
 * If the sequence is of type <code>DNA</code> or <code>RNA</code>, one
 * can get the reverse complement:
 * \code
 * mySecondDNA = myDNA.reverseComplement();
 * \endcode
 * or
 * \code
 * myDNA.reverseComplement(mySecondDNA);
 * \endcode
 * If one needs to reverse complement itself, one may do
 * \code
 * // Before below, some sequence
 * myDNA.selfReverseComplement();
 * // Now, myDNA is the reverse complement of the original sequence.
 * \endcode
 * Simiarly, one could to the above analogs of <code>changeTtoU</code> and
 * <code>changeUtoT</code>.
 * If the sequence is of type <code>DNA</code> or
 * <code>RNA</code>, one can check for ambiguous bases:
 * \code
 * if (myDNA.isAmbiguous()) cout << "Ambiguous\n";
 * \endcode
 * If the sequence is of type <code>DNA</code> or <code>RNA</code>, one
 * can convert to a protein
 * sequence (the DNA's T is converted to a U for the purpose of translation)
 * and the * character is used to denote a STOP.
 * \code
 * Sequence myProtein = myDNA.convertToProtein();
 * \endcode
 * or
 * \code
 * Sequence myProtein;
 * myDNA.convertToProtein(myProtein);
 * \endcode
 * If one wishes to know the GC content of a sequence, one may do
 * \code
 * double gcPct = myDNA.gc();
 * \endcode
 * where gcPct will be between 0 and 1 (1 = 100%).
 */
class Sequence : public PrimarySequence
{
public:
    // String sizes are of type size_t.
    typedef size_t size_type;

    SeqType seq_type;

    Sequence()
    {
        seq_type = DNA;
    }
    Sequence(const PrimarySequence& newSeq)
    {
        assign(newSeq.c_str());
        seq_type = DNA;
    }
    Sequence(string newString)
    {
        assign(newString);
        seq_type = DNA;
    }
    Sequence(const Sequence& newSeq)
    {
        assign(newSeq.c_str());
        seq_type = newSeq.seq_type;
    }

    Sequence& operator=(const string& theString)
    {
        assign(theString);
        return *this;
    }
    Sequence& operator=(const PrimarySequence& theSeq)
    {
        assign(theSeq.c_str());
        return *this;
    }


    SeqType getType()
    {
        return seq_type;
    }
    void setType(SeqType newType)
    {
        seq_type = newType;
    }


    Sequence substr( size_t pos = 0, size_t n = string::npos )
    {
        Sequence result(string::substr(pos, n));
        result.seq_type = seq_type;
        return result;
    }


    double gc()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        size_t numChar = size();
        if (numChar == 0)
        {
            return 0.0;
        }
        size_t numGC = 0;
        for (size_t i = 0; i < numChar; ++i)
        {
            if (((*this)[i] == 'g') ||
                    ((*this)[i] == 'G') ||
                    ((*this)[i] == 'c') ||
                    ((*this)[i] == 'C'))
            {
                ++numGC;
            }
        }
        return static_cast<double>(numGC)/static_cast<double>(numChar);
    }

    void changeTtoU(Sequence& theSequence)
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        size_t numChar = size();
        theSequence.resize(numChar);
        char thisBase, outBase;
        for (size_t i = 0; i < numChar; ++i)
        {
            thisBase = (*this)[i];
            switch (thisBase)
            {
                case 't' :
                    outBase = 'u';
                    break;
                case 'T' :
                    outBase = 'U';
                    break;
                default:
                    outBase = thisBase;
                    break;
            }
            theSequence[i] = outBase;
        }
    }
    Sequence changeTtoU()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        changeTtoU(tempSeq);
        tempSeq.seq_type = this->seq_type;
        return tempSeq;
    }
    void selfChangeTtoU()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        changeTtoU(tempSeq);
        assign(tempSeq.c_str());
        tempSeq.seq_type = this->seq_type;
    }

    void changeUtoT(Sequence& theSequence)
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        size_t numChar = size();
        theSequence.resize(numChar);
        char thisBase, outBase;
        for (size_t i = 0; i < numChar; ++i)
        {
            thisBase = (*this)[i];
            switch (thisBase)
            {
                case 'u' :
                    outBase = 't';
                    break;
                case 'U' :
                    outBase = 'T';
                    break;
                default:
                    outBase = thisBase;
                    break;
            }
            theSequence[i] = outBase;
        }
    }
    Sequence changeUtoT()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        changeUtoT(tempSeq);
        tempSeq.seq_type = this->seq_type;
        return tempSeq;
    }
    void selfChangeUtoT()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        changeUtoT(tempSeq);
        tempSeq.seq_type = this->seq_type;
        assign(tempSeq.c_str());
    }

    void reverseComplement(Sequence& theSequence)
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        size_t numChar = size();
        theSequence.resize(numChar);
        char thisBase, outBase;
        for (size_t i = 0; i < numChar; ++i)
        {
            thisBase = (*this)[i];
            switch (thisBase)
            {
                case 'a' :
                    if (seq_type == DNA)
                    {
                        outBase = 't';
                    }
                    else   // RNA
                    {
                        outBase = 'u';
                    }
                    break;
                case 'A' :
                    if (seq_type == DNA)
                    {
                        outBase = 'T';
                    }
                    else   // RNA
                    {
                        outBase = 'U';
                    }
                    break;
                case 'c' :
                    outBase = 'g';
                    break;
                case 'C' :
                    outBase = 'G';
                    break;
                case 'g' :
                    outBase = 'c';
                    break;
                case 'G' :
                    outBase = 'C';
                    break;
                case 't' :
                    outBase = 'a';
                    break;
                case 'T' :
                    outBase = 'A';
                    break;
                case 'u' :
                    outBase = 'a';
                    break;
                case 'U' :
                    outBase = 'A';
                    break;
                case 'r' :
                    outBase = 'y';
                    break;   // R = A or G
                case 'R' :
                    outBase = 'Y';
                    break;
                case 'y' :
                    outBase = 'r';
                    break;   // Y = C or T/U
                case 'Y' :
                    outBase = 'R';
                    break;
                case 's' :
                    outBase = 's';
                    break;   // S = C or G
                case 'S' :
                    outBase = 'S';
                    break;
                case 'w' :
                    outBase = 'w';
                    break;   // W = A or T/U
                case 'W' :
                    outBase = 'W';
                    break;
                case 'k' :
                    outBase = 'm';
                    break;   // K = G or T/U
                case 'K' :
                    outBase = 'M';
                    break;
                case 'm' :
                    outBase = 'k';
                    break;   // M = A or C
                case 'M' :
                    outBase = 'K';
                    break;
                case 'b' :
                    outBase = 'v';
                    break;   // B = C, G, or T/U
                case 'B' :
                    outBase = 'V';
                    break;
                case 'd' :
                    outBase = 'h';
                    break;   // D = A, G, or T/U
                case 'D' :
                    outBase = 'H';
                    break;
                case 'h' :
                    outBase = 'd';
                    break;   // H = A, C, or T/U
                case 'H' :
                    outBase = 'D';
                    break;
                case 'v' :
                    outBase = 'b';
                    break;   // V = A, C, or G
                case 'V' :
                    outBase = 'B';
                    break;
                case 'n' :
                    outBase = 'n';
                    break;   // N = A, C, G, or T/U
                case 'N' :
                    outBase = 'N';
                    break;
                case '.' :
                    outBase = '.';
                    break;   // . and - denote a gap
                case '-' :
                    outBase = '-';
                    break;
                case '*' :
                    outBase = '*';
                    break;   // * is a padding symbol from
                    // the SAM/CIGAR format.
                default: // Should never get here.
                    // Unrecognized base.
                    cerr << "Error in reverseComplement.  Unknown base character = " << thisBase << ".  Exiting." << endl;
                    exit(1);   // Stop right away.
                    break;
            }
            theSequence[numChar - 1 - i] = outBase;
        }
    }

    Sequence reverseComplement()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        reverseComplement(tempSeq);
        tempSeq.seq_type = this->seq_type;
        return tempSeq;
    }


    void selfReverseComplement()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        Sequence tempSeq;
        reverseComplement(tempSeq);
        tempSeq.seq_type = this->seq_type;
        assign(tempSeq.c_str());
    }


    bool isAmbiguous()
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
#endif
        char thisBase;
        size_t numChar = size();
        for (size_t i = 0; i < numChar; ++i)
        {
            // For some reason, toupper returns an int that should be then be
            // cast back to a char.
            thisBase = static_cast<char>(toupper((*this)[i]));
            switch (thisBase)
            {
                case 'R' :
                case 'Y' :
                case 'S' :
                case 'W' :
                case 'K' :
                case 'M' :
                case 'B' :
                case 'D' :
                case 'H' :
                case 'V' :
                case 'N' :
                case '.' :
                case '-' :
                    return true;
                    break;
                default :
                    break;
            }
        }
        return false;
    }


    // forceStartM means if the string starts AUG, GUG, UUG, AUU, or CUG, then
    // put M (start codon) at the beginning.
    void convertToProtein(Sequence& newSeq, size_t start = 0, bool forceStartM = false)
    {
#ifdef PROG_DEBUG
        assert(seq_type == DNA || seq_type == RNA);
        // Since using codon table, must have a multiple of 3.
        assert((size() - start) % 3 == 0);
#endif
#ifdef DEBUG
        // Cannot convert an ambiguous sequence.
        assert(!isAmbiguous());
#endif
        size_t numChar = size();
        newSeq.resize(numChar/3);  // numChar is exactly divisible by 3.
        // In below, cB is currentBase.
        char cB[3], newBase;
        string cBstring = "   ";
        for (size_t i = start; i < numChar; i += 3)
        {
            cB[0] = static_cast<char>(toupper((*this)[i]));
            cB[1] = static_cast<char>(toupper((*this)[i+1]));
            cB[2] = static_cast<char>(toupper((*this)[i+2]));
            if (cB[0] == 'T')
            {
                cB[0] = 'U';
            }
            if (cB[1] == 'T')
            {
                cB[1] = 'U';
            }
            if (cB[2] == 'T')
            {
                cB[2] = 'U';
            }
            newBase = '0';  // We'll designate * as Stop,
            // and 1 as error.
            if (cB[0] == 'U')
            {
                if (cB[1] == 'U')
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'F';
                    }
                    else
                    {
                        newBase = 'L';
                    }
                }
                else if (cB[1] == 'C')
                {
                    newBase = 'S';
                }
                else if (cB[1] == 'A')
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'Y';
                    }
                    else
                    {
                        newBase = '*';
                    }
                }
                else       // cB[1] == 'G'
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'C';
                    }
                    else if (cB[2] == 'A')
                    {
                        newBase = '*';
                    }
                    else
                    {
                        newBase = 'W';
                    }
                }
            }
            else if (cB[0] == 'C')
            {
                if (cB[1] == 'U')
                {
                    newBase = 'L';
                }
                else if (cB[1] == 'C')
                {
                    newBase = 'P';
                }
                else if (cB[1] == 'A')
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'H';
                    }
                    else if (cB[2] == 'A')
                    {
                        newBase = 'Q';
                    }
                    else
                    {
                        newBase = 'Q';
                    }
                }
                else       // cB[1] == 'G'
                {
                    newBase = 'R';
                }
            }
            else if (cB[0] == 'A')
            {
                if (cB[1] == 'U')
                {
                    if (cB[2] == 'G')
                    {
                        newBase = 'M';
                    }
                    else
                    {
                        newBase = 'I';
                    }
                }
                else if (cB[1] == 'C')
                {
                    newBase = 'T';
                }
                else if (cB[1] == 'A')
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'N';
                    }
                    else
                    {
                        newBase = 'K';
                    }
                }
                else       // cB[1] == 'G'
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'S';
                    }
                    else
                    {
                        newBase = 'R';
                    }
                }
            }
            else       // cb[0] == 'G'
            {
                if (cB[1] == 'U')
                {
                    newBase = 'V';
                }
                else if (cB[1] == 'C')
                {
                    newBase = 'A';
                }
                else if (cB[1] == 'A')
                {
                    if (cB[2] == 'U' || cB[2] == 'C')
                    {
                        newBase = 'D';
                    }
                    else
                    {
                        newBase = 'E';
                    }
                }
                else       // cB[1] == 'G'
                {
                    newBase = 'G';
                }
            }
#ifdef DEBUG
            assert(newBase != '0');
#endif
            if (forceStartM && (i == start))
            {
                cBstring[0] = cB[0];
                cBstring[1] = cB[1];
                cBstring[2] = cB[2];
                // Below are all the known start codons.
                if ((cBstring == "AUG") || (cBstring == "GUG") ||
                        (cBstring == "UUG") || (cBstring == "AUU") ||
                        (cBstring == "CUG"))
                {
                    newBase = 'M';
                }
            }
            newSeq[(i-start)/3] = newBase;
        }
        newSeq.seq_type = Protein;
    }

    Sequence convertToProtein(size_t start = 0, bool forceStartM = false)
    {
        Sequence tempSeq;
        convertToProtein(tempSeq, start, forceStartM);
        tempSeq.seq_type = Protein;
        return tempSeq;
    }

};

#endif
