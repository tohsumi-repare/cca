#ifndef MOLBIOLIB_COMMANDPARSER_H
#define MOLBIOLIB_COMMANDPARSER_H 1

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
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"

using namespace std;

/** \file CommandParser.hpp
 * Definitions to ease parsing of command-line arguments.
 * All command arguments must be of the form X=Y where Y is either a
 * series of characters without spaces or Y is enclosed in double quotes.
 * Example usage is
 * \code
 * BeginCommandArguments
 * // Use below instead of above to always suppress printing of the header.
 * // BeginCommandArgumentsNoHeader
 *    CommandArgumentString("FASTA", FASTA);
 *    CommandArgumentDoubleDoc("MIN_I", minI, "The minimum intensity.");
 *    CommandArgumentLongDefault("MAXALIGN", maxAlign, 10000);
 *    CommandArgumentBoolDefaultDoc("DO_MORE", doMore, false, "More flag.");
 * EndCommandArguments
 * \endcode
 * The above, String, Double, Long, and Bool, are the only four types of
 * parameters handled.  However, string can be converted to a wide variety
 * of data type via #convertFromString.  Comma-separated valued strings can
 * also be first processed using #splitString.
 */


/** \file CommandParser.hpp
 * Definitions to ease parsing of command-line arguments.
 * \def MAX_TCLS2
 * Longest command line possible.  Usually limited by the shell.
 * POSIX standard is a minimum of 4096 bytes.
 */
#define MAX_TCLS2 4096


/** \file CommandParser.hpp
 * Definitions to ease parsing of command-line arguments.
 * \todo Replace <code>enum CommandParserFileType</code> with
 * <code>enum class CommandParserFileType</code> when strongly typed
 * enums work.
 */
enum CommandParserFileType
{
    FileIn, FileOut
};



/** \file CommandParser.hpp
 * \fn template<typename T> int commandParser(int argc, char* argv[], string name, T& theVal, bool optional, T defaultVal, string theDoc)
 * This function tries to find parameters of the form <code>name=theVal</code>
 * in the command arguments.  If found, the value is stored in the variable
 * <code>theVal</code>, otherwise the default value is used.  If the value is
 * not optional and if the name is not found, then a value of 1 is returned.
 * Otherwise, 0 is returned.  This is so we can catch the error at the end.
 * See the code for <code>BeginCommandArguments</code>.
 */
template<typename T>
int commandParser(int argc, char* argv[], string name, T& theVal, bool optional, T defaultVal, string theDoc)
{

    bool found = false;
    int i;
    theVal = defaultVal;
    for (i = 1; i < argc && !found; ++i)
    {
        string currentArg(argv[i]);
        if (currentArg.find("=") == string::npos)
        {
            cerr << "Error!  Command line parameter " << i << ", " << currentArg << ", does not have an equal sign in it.  Exiting." << endl;
            exit(1);
        }
        string argName = currentArg.substr(0, currentArg.find("="));
        string argValue = currentArg.substr(currentArg.find("=") + 1);
        if (name ==  argName)
        {
            found = true;
            theVal = convertFromString<T>(argValue);
        }
    }
    if (!found && !optional)
    {
        cerr << "Error!  Required command line argument " << name << " missing.  Exiting." << endl;
        return 1;
    }
    return 0;
}
/** \file CommandParser.hpp
 * \fn int commandParserFile(int argc, char* argv[], string name, string& theVal, bool optional, string defaultVal, string theDoc, CommandParserFileType fileType)
 * This function tries to find parameters of the form <code>name=theVal</code>
 * in the command arguments.  If found, the value is stored in the variable
 * <code>theVal</code>, otherwise the default value is used.  If the value is
 * not optional and if the name is not found, then a value of 1 is returned.
 * Otherwise, 0 is returned.  This is so we can catch the error at the end.
 * See the code for <code>BeginCommandArguments</code>.
 */
int commandParserFile(int argc, char* argv[], string name, string& theVal, bool optional, string defaultVal, string theDoc, CommandParserFileType fileType)
{

    bool found = false;
    int i;
    theVal = defaultVal;
    for (i = 1; i < argc && !found; ++i)
    {
        string currentArg(argv[i]);
        if (currentArg.find("=") == string::npos)
        {
            cerr << "Error!  Command line parameter " << i << ", " << currentArg << ", does not have an equal sign in it.  Exiting." << endl;
            exit(1);
        }
        string argName = currentArg.substr(0, currentArg.find("="));
        string argValue = currentArg.substr(currentArg.find("=") + 1);
        if (name ==  argName)
        {
            found = true;
            theVal = argValue;
        }
    }
    if (!found && !optional)
    {
        cerr << "Error!  Required command line argument " << name << " missing.  Exiting." << endl;
        return 1;
    }
    if ((fileType == FileIn) && (theVal != "") && (fileSize(theVal) == static_cast<ifstream::pos_type>(-1)))
    {
        cerr << "Error!  File " << theVal << " is missing.  Exiting." << endl;
        return 1;
    }
#ifndef IGNORE_FILEOUT
    if ((fileType == FileOut) && (theVal != "") && (fileSize(theVal) != static_cast<ifstream::pos_type>(-1)))
    {
        cerr << "Error!  File " << theVal << " exists.  Exiting." << endl;
        return 1;
    }
#endif
    return 0;
}
/** \file CommandParser.hpp
 * \fn int commandParserBool(int argc, char* argv[], string name, bool& theVal, bool optional, bool defaultVal, string theDoc)
 * This is the specialization of #commandParser above for
 * type <code>bool</code>.  Needed since
 * <code>convertFromString<bool>("true")</code> gives <code>false</code>.
 */
int commandParserBool(int argc, char* argv[], string name, bool& theVal, bool optional, bool defaultVal, string theDoc)
{

    bool found = false;
    int i;
    theVal = defaultVal;
    for (i = 1; i < argc && !found; ++i)
    {
        string currentArg(argv[i]);
        if (currentArg.find("=") == string::npos)
        {
            cerr << "Error!  Command line parameter " << i << ", " << currentArg << ", does not have an equal sign in it.  Exiting." << endl;
            exit(1);
        }
        string argName = currentArg.substr(0, currentArg.find("="));
        string argValue = currentArg.substr(currentArg.find("=") + 1);
        if (name ==  argName)
        {
            found = true;
            if ((argValue.find("true") != string::npos) ||
                    (argValue.find("True") != string::npos) ||
                    (argValue.find("TRUE") != string::npos) ||
                    (argValue.find("1") != string::npos))
            {
                theVal = true;
            }
            else if ((argValue.find("false") != string::npos) ||
                     (argValue.find("False") != string::npos) ||
                     (argValue.find("FALSE") != string::npos) ||
                     (argValue.find("0") != string::npos))
            {
                theVal = false;
            }
            else
            {
                cerr << "Error!  Boolean values must be true or false.  Exiting." << endl;
                exit(1);
            }
        }
    }
    if (!found && !optional)
    {
        cerr << "Error!  Required command line argument " << name << " missing.  Exiting." << endl;
        return 1;
    }
    return 0;
}




/** \file CommandParser.hpp
 * \fn int commandParserKeysCheck(int argc, char* argv[], set<string> names) {
 * This function tries to find parameters of the form <code>name=theVal</code>
 * in the command arguments and sees if there is a name which is not contained
 * in names.  If a name is not found, then 1 is returned.  Otherwise, return 0.
 */
int commandParserKeysCheck(int argc, char* argv[], set<string> names)
{
    int result = 0;
    for (int i = 1; i < argc; ++i)
    {
        string currentArg(argv[i]);
        if (currentArg.find("=") == string::npos)
        {
            cerr << "Error!  Command line parameter " << i << ", " << currentArg << ", does not have an equal sign in it.  Exiting.\n";
            exit(1);
        }
        string argName = currentArg.substr(0, currentArg.find("="));
        if (names.find(argName) == names.end())
        {
            cerr << "Error!  " << argName << " is not a known parameter.\n";
            result = 1;
        }
    }
    return result;
}




/** \file CommandParser.hpp
 * \fn void commandParserUsage(vector<string>& commandUsage)
 * This just prints the command-line usage.
 */
void commandParserUsage(vector<string>& commandUsage)
{
    cerr << "All command-line arguments must be of the form PARAMETER=VALUE.  Usage:\n";
    for (vector<string>::size_type i = 0; i < commandUsage.size(); ++i)
    {
        cerr << commandUsage[i] << endl;
    }
    exit(1);
}



/** \file CommandParser.hpp
 * \def BeginCommandArguments
 * Macro signifying the beginning of section where the CommandArgument macros
 * below may be used. If HELP=true, then print out the parameter lists.  If
 * HEADER=true (default - value stored in checkComArgHeader after first call to
 * commandParserBool), print header.  The tCLS2[0] = ' '; is there just so that
 * the compiler warning about tCLS2 never being used for an empty set of
 * command arguments (BeginCommandArguments, then on next line
 * EndCommandArguments) will not occur, even with the -Wall option.
 */
#define BeginCommandArguments bool checkComArgHeader = true; int commandLineErrorNum = commandParserBool(argc, argv, "HEADER", checkComArgHeader, true, true, "Specify as false if the initial command-line printing is to be suppressed."); if (checkComArgHeader) { cout << "Program compiled on " << __DATE__ << " at " << __TIME__ << endl << "Command-line = "; for (int argc_counter = 0; argc_counter < argc; ++argc_counter) cout << argv[argc_counter] << " "; cout << endl; } vector<string> commandLineUsage; bool checkComArgHelp = false; commandLineErrorNum += commandParserBool(argc, argv, "HELP", checkComArgHelp, true, false, "Specify as true if want to see parameters usage list."); if (checkComArgHelp) ++commandLineErrorNum; string tempCommandLineString = ""; char tCLS2[MAX_TCLS2]; tCLS2[0] = ' '; set<string> commandLineParameterNames; commandLineParameterNames.insert("HEADER"); commandLineParameterNames.insert("HELP");

/** \file CommandParser.hpp
 * \def BeginCommandArgumentsStderr
 * Macro signifying the beginning of section where the CommandArgument macros
 * below may be used. If HELP=true, then print out the parameter lists.  If
 * HEADER=true (default - value stored in checkComArgHeader after first call to
 * commandParserBool), print header.  The tCLS2[0] = ' '; is there just so that
 * the compiler warning about tCLS2 never being used for an empty set of
 * command arguments (BeginCommandArguments, then on next line
 * EndCommandArguments) will not occur, even with the -Wall option.  Differs
 * from BeginCommandArguments in that outputs command line to cerr instead of
 * cout.
 */
#define BeginCommandArgumentsStderr bool checkComArgHeader = true; int commandLineErrorNum = commandParserBool(argc, argv, "HEADER", checkComArgHeader, true, true, "Specify as false if the initial command-line printing is to be suppressed."); if (checkComArgHeader) { cerr << "Program compiled on " << __DATE__ << " at " << __TIME__ << endl << "Command-line = "; for (int argc_counter = 0; argc_counter < argc; ++argc_counter) cerr << argv[argc_counter] << " "; cerr << endl; } vector<string> commandLineUsage; bool checkComArgHelp = false; commandLineErrorNum += commandParserBool(argc, argv, "HELP", checkComArgHelp, true, false, "Specify as true if want to see parameters usage list."); if (checkComArgHelp) ++commandLineErrorNum; string tempCommandLineString = ""; char tCLS2[MAX_TCLS2]; tCLS2[0] = ' '; set<string> commandLineParameterNames; commandLineParameterNames.insert("HEADER"); commandLineParameterNames.insert("HELP");

/** \file CommandParser.hpp
 * \def BeginCommandArgumentsNoHeader
 * Macro signifying the beginning of section where the CommandArgument macros
 * below may be used.  This is used when the printing of the header is always
 * suppressed.  The HEADER parameter's value will be ignored.
 */
#define BeginCommandArgumentsNoHeader vector<string> commandLineUsage; bool checkComArgHelp = false; int commandLineErrorNum = commandParserBool(argc, argv, "HELP", checkComArgHelp, true, false, "Specify as true if want to see parameters usage list."); if (checkComArgHelp) ++commandLineErrorNum; string tempCommandLineString = ""; char tCLS2[MAX_TCLS2]; tCLS2[0] = ' '; set<string> commandLineParameterNames; commandLineParameterNames.insert("HEADER"); commandLineParameterNames.insert("HELP");


#define CommandArgumentFirstPart(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, CA_type, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) commandLineErrorNum += commandParser<CA_type>(argc, argv, MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC); if (MOLBIOLIB_CA_OPT) tempCommandLineString = "Optional"; else tempCommandLineString = "Required";  tempCommandLineString += " parameter: "; sprintf(tCLS2, "%s", MOLBIOLIB_CA_NAME); tempCommandLineString += tCLS2 ; if (MOLBIOLIB_CA_DEFAULT) tempCommandLineString += " (Default value = " + convertToString<CA_type>(MOLBIOLIB_CA_VALUE) + ")"; commandLineParameterNames.insert(MOLBIOLIB_CA_NAME);
// Special case needed for file in and out
#define CommandArgumentFirstPartFile(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, CA_type, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT, MOLBIOLIB_CA_FILETYPE) commandLineErrorNum += commandParserFile(argc, argv, MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_FILETYPE); if (MOLBIOLIB_CA_OPT) tempCommandLineString = "Optional"; else tempCommandLineString = "Required";  tempCommandLineString += " parameter: "; sprintf(tCLS2, "%s", MOLBIOLIB_CA_NAME); tempCommandLineString += tCLS2 ; if (MOLBIOLIB_CA_DEFAULT) tempCommandLineString += " (Default value = " + convertToString<CA_type>(MOLBIOLIB_CA_VALUE) + ")"; commandLineParameterNames.insert(MOLBIOLIB_CA_NAME);
// Special case needed for bool - must call commandParseBool
#define CommandArgumentFirstPartBool(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, CA_type, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) commandLineErrorNum += commandParserBool(argc, argv, MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC); if (MOLBIOLIB_CA_OPT) tempCommandLineString = "Optional"; else tempCommandLineString = "Required";  tempCommandLineString += " parameter: "; sprintf(tCLS2, "%s", MOLBIOLIB_CA_NAME); tempCommandLineString += tCLS2 ; if (MOLBIOLIB_CA_DEFAULT) tempCommandLineString += " (Default value = " + convertToString<bool>(MOLBIOLIB_CA_VALUE) + ")"; commandLineParameterNames.insert(MOLBIOLIB_CA_NAME);

#define CommandArgumentLastPart(MOLBIOLIB_CA_DOC) sprintf(tCLS2, "%s", MOLBIOLIB_CA_DOC); tempCommandLineString += tCLS2; commandLineUsage.push_back(tempCommandLineString);




/** \file CommandParser.hpp
 * \def CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT)
 * This is the macro that requires all of the arguments to be present.  The
 * macros #CommandArgumentLongDefaultDoc, #CommandArgumentLongDefault,
 * #CommandArgumentLongDoc, and #CommandArgumentLong all are special cases of
 * this macro.
 * \section Usage Usage example:
 * \code
 * BeginCommandArguments
 * CommandArgumentLongFull("TestLong", tester, true, 314, "Sample optional argument of type long.", true);
 * EndCommandArguments
 * \endcode
 * The above defines the variable <code>tester</code>, then checks to see if in
 * the command-line arguments there is something like
 * \code
 * TestLong=271
 * \endcode
 * and if so then initializes <code>tester</code> to 271.  If no such
 * argument exists, then <code>tester</code> is set to 314.  If there is a
 * usage printed, then string at the end is printed.  The <code>true</code>
 * means that this is an optional argument.  If instead <code>false</code> was
 * passed, then this is a required argument.  In such cases, any value may be
 * passed as the default argument, as it is never used.
 */
#define CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) long MOLBIOLIB_CA_VAR; CommandArgumentFirstPart(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, long, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT); tempCommandLineString += " of type long.  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentLongDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * \section Usage Usage:
 * \code
 * CommandArgumentLongDefaultDoc("TestLong", tester, 314, "Sample optional argument.");
 * \endcode
 */
#define CommandArgumentLongDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentLongDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * \section Usage Usage:
 * \code
 * CommandArgumentLongDefault("TestLong", tester, 314);
 * \endcode
 */
#define CommandArgumentLongDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentLongDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * \section Usage Usage:
 * \code
 * CommandArgumentLongDoc("TestLong", tester, "Sample required argument.");
 * \endcode
 */
#define CommandArgumentLongDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0, MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentLong(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * \section Usage Usage:
 * \code
 * CommandArgumentLong("TestLong", tester);
 * \endcode
 */
#define CommandArgumentLong(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentLongFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0, "", false);


/** \file CommandParser.hpp
 * \def CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>size_t</code>.
 */
#define CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) size_t MOLBIOLIB_CA_VAR; CommandArgumentFirstPart(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, size_t, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT); tempCommandLineString += " of type size_t.  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentSize_tDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>size_t</code>.
 */
#define CommandArgumentSize_tDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentSize_tDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>size_t</code>.
 */
#define CommandArgumentSize_tDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentSize_tDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>size_t</code>.
 */
#define CommandArgumentSize_tDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0, MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentSize_t(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>size_t</code>.
 */
#define CommandArgumentSize_t(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentSize_tFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0, "", false);


/** \file CommandParser.hpp
 * \def CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>double</code>.
 */
#define CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) double MOLBIOLIB_CA_VAR; CommandArgumentFirstPart(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, double, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT); tempCommandLineString += " of type double.  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentDoubleDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>double</code>.
 */
#define CommandArgumentDoubleDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentDoubleDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>double</code>.
 */
#define CommandArgumentDoubleDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentDoubleDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>double</code>.
 */
#define CommandArgumentDoubleDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0.0, MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentDouble(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>double</code>.
 */
#define CommandArgumentDouble(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentDoubleFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, 0.0, "", false);


/** \file CommandParser.hpp
 * \def CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) string MOLBIOLIB_CA_VAR; CommandArgumentFirstPart(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, string, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT); tempCommandLineString += " of type string.  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentStringDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentStringDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentStringDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentStringDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentStringDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentStringDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentString(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentString(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentStringFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", "", false);




/** \file CommandParser.hpp
 * \def CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) string MOLBIOLIB_CA_VAR; CommandArgumentFirstPartFile(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, string, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT, FileIn); tempCommandLineString += " of type input file path and name (a string, file should be present).  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentFileInDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileInDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentFileInDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileInDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentFileInDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileInDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentFileIn(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileIn(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentFileInFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", "", false);




/** \file CommandParser.hpp
 * \def CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) string MOLBIOLIB_CA_VAR; CommandArgumentFirstPartFile(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, string, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT, FileOut); tempCommandLineString += " of type output file path and name (a string, file should not be present).  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentFileOutDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileOutDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentFileOutDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileOutDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentFileOutDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileOutDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentFileOut(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>string</code>.
 */
#define CommandArgumentFileOut(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentFileOutFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", "", false);




/** \file CommandParser.hpp
 * \def CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongFull.  Same except defines
 * type <code>bool</code>.  <em>Important:</em> one must pass either the value
 * <code>True</code> or <code>False</code> (case-sensitive) to this parameter
 * for this to work.
 */
#define CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT) bool MOLBIOLIB_CA_VAR; CommandArgumentFirstPartBool(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, bool, MOLBIOLIB_CA_OPT, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, MOLBIOLIB_CA_DEFAULT); tempCommandLineString += " of type bool.  "; CommandArgumentLastPart(MOLBIOLIB_CA_DOC);
/** \file CommandParser.hpp
 * \def CommandArgumentBoolDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDefaultDoc.  Same except defines
 * type <code>bool</code>.
 */
#define CommandArgumentBoolDefaultDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC) CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, MOLBIOLIB_CA_DOC, true);
/** \file CommandParser.hpp
 * \def CommandArgumentBoolDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE)
 * See documentation for #CommandArgumentLongDefault.  Same except defines
 * type <code>bool</code>.
 */
#define CommandArgumentBoolDefault(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_VALUE) CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, true, MOLBIOLIB_CA_VALUE, "", true);
/** \file CommandParser.hpp
 * \def CommandArgumentBoolDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC)
 * See documentation for #CommandArgumentLongDoc.  Same except defines
 * type <code>bool</code>.
 */
#define CommandArgumentBoolDoc(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, MOLBIOLIB_CA_DOC) CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", MOLBIOLIB_CA_DOC, false);
/** \file CommandParser.hpp
 * \def CommandArgumentBool(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR)
 * See documentation for #CommandArgumentLong.  Same except defines
 * type <code>bool</code>.
 */
#define CommandArgumentBool(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR) CommandArgumentBoolFull(MOLBIOLIB_CA_NAME, MOLBIOLIB_CA_VAR, false, "", "", false);




/** \file CommandParser.hpp
 * \def EndCommandArguments
 * Macro signifying the end of the section where the CommandArgument macros
 * below may be used.  Must only be used after #BeginCommandArguments.
 * Both <code>cout</code> and <code>cerr</code> are set to
 * <code>MOLBIOLIB_IO_PRECISION</code> (defined in PrimitiveTypes.hpp) here.
 */
#define EndCommandArguments commandLineErrorNum += commandParserKeysCheck(argc, argv, commandLineParameterNames); if (commandLineErrorNum > 0) commandParserUsage(commandLineUsage); cout.precision(MOLBIOLIB_IO_PRECISION); cerr.precision(MOLBIOLIB_IO_PRECISION);



#endif

