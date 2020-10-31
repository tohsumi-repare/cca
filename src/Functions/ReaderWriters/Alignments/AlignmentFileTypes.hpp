#ifndef MOLBIOLIB_ALIGNMENTFILETYPES_H
#define MOLBIOLIB_ALIGNMENTFILETYPES_H 1

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
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Objects/ReadOnlyDelimitedFile.hpp"

/** \file AlignmentFileTypes.hpp
 * Contains just the AlignmentFileType enum.
 * <ul>
 *    <li> Null = Alignment format type not specified. </li>
 *    <li> Qltout = Broad Institute Arachne's
 *         qltout alignment format. </li>
 *    <li> SAM = Sequence Alignment/MAP (SAM) format. </li>
 *    <li> Sms = Helics sms format. </li>
 * </ul>
 * \enum AlignmentFileType
 * \todo Replace <code>enum AlignmentFileType</code> with
 * <code>enum class AlignmentFileType</code> when strongly typed
 * enums work.
 */
enum AlignmentFileType
{
    Null, Qltout, SAM, Sms
};


/** \file AlignmentFileTypes.hpp
 * Determine the type of alignment file.  Use ignoreExist = true when creating
 * a new file, e.g. for alignment file writers.  maxLengthToCheck is the
 * maximum number of lines to check a potential alignment file for a line that
 * has a characteristic specific to an alignment file format.
 */
AlignmentFileType determineAlignmentFileType(string filename,
        bool ignoreExist = false,
        size_t maxLengthToCheck = 10000)
{
    // First, check if the filename exists
    if (!ignoreExist &&
            (fileSize(filename) == static_cast<ifstream::pos_type>(-1)))
    {
        return Null;
    }

    // Check the file extension.  Cases are ignored.
    if (filename.size() > 7)
    {
        string qltoutExt = filename.substr(filename.size()-7, 7);
        string qltoutExtUpper;
        qltoutExtUpper.resize(7);
        for (size_t i = 0; i < 7; ++i)
        {
            qltoutExtUpper[i] = static_cast<char>(toupper(qltoutExt[i]));
        }
        if (qltoutExtUpper == ".QLTOUT")
        {
            return Qltout;
        }
    }
    if (filename.size() > 4)
    {
        string threeLetterExt = filename.substr(filename.size()-4, 4);
        string threeLetterExtUpper;
        threeLetterExtUpper.resize(4);
        for (size_t i = 0; i < 4; ++i)
        {
            threeLetterExtUpper[i] = static_cast<char>(toupper(threeLetterExt[i]));
        }
        if (threeLetterExtUpper == ".SAM")
        {
            return SAM;
        }
        else if ((threeLetterExtUpper == ".SMS") ||
                 (threeLetterExtUpper == ".TXT"))
        {
            return Sms;
        }
    }


    // Try to determine from content the alignment file type.  This is,
    // of course, not always reliable.
    ReadOnlyTSVFile ifp(filename, true);   // Read in sequential mode.
    size_t numRows = ifp.size();
    vector<string> tokens;
    // Check for qltout format by checking for a line that starts with QUERY.
    size_t maxRows = min(numRows, maxLengthToCheck);
    for (size_t i = 0; i < maxRows; ++i)
    {
        tokens = ifp[i];
        if ((tokens.size() > 0) && (tokens[0] == "QUERY"))
        {
            ifp.close();
            return Qltout;
        }
    }
    // Check for SAM format by checking for a valid CIGAR token.
    // ReadOnlyTSVFile's ifp will automatically call resetToBegin() on i=0.
    for (size_t i = 0; i < maxRows; ++i)
    {
        tokens = ifp[i];
        if ((tokens.size() > 0) && (tokens[0][0] == '@'))
        {
            continue;
        }
        if (tokens.size() > 5)
        {
            string CIGAR = tokens[5];
            size_t cigarSize = CIGAR.size();
            size_t currCigarPos = 0;
            string currSize = "";

            bool goodCigar = true;

            while (goodCigar && (currCigarPos < cigarSize))
            {
                char& currChar = CIGAR[currCigarPos];
                if ((currChar != 'M') && (currChar != 'I') && (currChar != 'D') &&
                        (currChar != 'N') && (currChar != 'S') && (currChar != 'H') &&
                        (currChar != 'P') && (currChar != 'X') && (currChar != '='))
                {
                    if (!isdigit(currChar))
                    {
                        goodCigar = false;
                    }
                }
                else
                {
                    // At operation.  Check next is a number if not at end of cigar.
                    if ((currCigarPos < (cigarSize-1)) &&
                            !isdigit(CIGAR[currCigarPos+1]))
                    {
                        goodCigar = false;
                    }
                }
                ++currCigarPos;
            }

            if (goodCigar)
            {
                ifp.close();
                return SAM;
            }
        }
    }
    // Check for Helicos Sms format by checking for sequences in the last two.
    // ReadOnlyTSVFile's ifp will automatically call resetToBegin() on i=0.
    for (size_t i = 0; i < maxRows; ++i)
    {
        tokens = ifp[i];
        if ((tokens.size() > 0) &&
                ((tokens[0][0] == '#') || (tokens[0] == "Reference_ID")))
        {
            continue;
        }
        if (tokens.size() == 14)
        {
            bool goodSms = true;
            string& tagSeq = tokens[12];
            for (size_t i = 0; goodSms && (i < tagSeq.size()); ++i)
            {
                if ((tagSeq[i] != '-') &&
                        (tagSeq[i] != 'A') && (tagSeq[i] != 'a') &&
                        (tagSeq[i] != 'C') && (tagSeq[i] != 'c') &&
                        (tagSeq[i] != 'G') && (tagSeq[i] != 'g') &&
                        (tagSeq[i] != 'T') && (tagSeq[i] != 't') &&
                        (tagSeq[i] != 'N'))
                {
                    goodSms = false;
                }
            }
            string& refSeq = tokens[12];
            for (size_t i = 0; goodSms && (i < refSeq.size()); ++i)
            {
                if ((refSeq[i] != '-') &&
                        (refSeq[i] != 'A') && (refSeq[i] != 'a') &&
                        (refSeq[i] != 'C') && (refSeq[i] != 'c') &&
                        (refSeq[i] != 'G') && (refSeq[i] != 'g') &&
                        (refSeq[i] != 'T') && (refSeq[i] != 't') &&
                        (refSeq[i] != 'N'))
                {
                    goodSms = false;
                }
            }
            if (goodSms)
            {
                ifp.close();
                return Sms;
            }
        }
    }
    ifp.close();

    // Could not figure out alignment file format.
    return Null;

}



#endif

