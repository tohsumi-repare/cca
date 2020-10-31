#ifndef MOLBIOLIB_H
#define MOLBIOLIB_H 1

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


/** \file MolBioLib.hpp
 * Main include file for MolBioLib applications.
 * Include this, so need not include anything else.
 */


#include "src/include/PrimitiveTypes.hpp"



// Below are utilities.  VectorStream.hpp must be placed before
// ConverString.hpp, which must be placed before CommandParser.hpp.
#include "src/Functions/SystemUtilities/System.hpp"
#include "src/Functions/SystemUtilities/Cstdio.hpp"
#include "src/Functions/SystemUtilities/FileSize.hpp"
#include "src/Objects/TextFileTypes.hpp"
#include "src/Objects/FileStream.hpp"
#include "src/Functions/ReaderWriters/Streams/VectorStream.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertString.hpp"
#include "src/Functions/SystemUtilities/CommandParser.hpp"
#include "src/Functions/SystemUtilities/CurrentTime.hpp"
#include "src/Functions/SystemUtilities/FileModificationTime.hpp"
#include "src/Functions/Transformers/StringOperations/UpperLowerCaseString.hpp"
#include "src/Functions/Transformers/StringOperations/ReplaceString.hpp"
#include "src/Functions/Transformers/StringOperations/TrimSpacesString.hpp"
#include "src/Functions/Transformers/StringOperations/CSVString.hpp"
#include "src/Functions/Transformers/StringOperations/SSVString.hpp"
#include "src/Functions/Transformers/StringOperations/SplitString.hpp"
#include "src/Functions/Transformers/StringOperations/JoinStrings.hpp"
#include "src/Functions/Transformers/VectorOperations/ConvertVectorSet.hpp"
#include "src/Functions/Transformers/Other/ConvertCSVSet.hpp"


// Below is the Table include.
#include "src/Objects/Table.hpp"


// Below are not Table objects, but some objects may need Table.
#include "src/Objects/SparseVector.hpp"
#include "src/Objects/ContigSizes.hpp"
#include "src/Objects/Qualifiers.hpp"
#include "src/Functions/ReaderWriters/Streams/QualifiersStream.hpp"
#include "src/Objects/ReadOnlyStringFile.hpp"
#include "src/Objects/ReadOnlyDelimitedFile.hpp"
#include "src/Objects/ReadOnlyTextFileTable.hpp"
#include "src/Objects/PrimarySequence.hpp"
#include "src/Objects/Sequence.hpp"
#include "src/Objects/Interval.hpp"
#include "src/Functions/Algorithms/Intervals/CreateIntervalsFromPoints.hpp"
#include "src/Functions/Algorithms/Intervals/IntervalBounds.hpp"
#include "src/Functions/ReaderWriters/Streams/IntervalStream.hpp"
#include "src/Functions/Algorithms/Vectors/VectorBounds.hpp"
#include "src/Objects/Location.hpp"
#include "src/Functions/ReaderWriters/Streams/TupleStream.hpp"
#include "src/Functions/ReaderWriters/Tuples/WriteTuple.hpp"
#include "src/Functions/Transformers/TupleOperations/GetTupleValues.hpp"
#include "src/Functions/Transformers/TupleOperations/JoinTuples.hpp"
#include "src/Functions/Transformers/TupleOperations/TokenizeTuple.hpp"
#include "src/Functions/Transformers/StringOperations/ConvertTokens.hpp"
#include "src/Functions/ReaderWriters/Other/TextSet.hpp"
#include "src/Objects/RandomLib.hpp"
#include "src/Objects/SyncSort.hpp"


// Below are Table-related includes.
#include "src/Objects/TableIndexMap.hpp"
#include "src/Functions/Transformers/TableOperations/ConcatenateTables.hpp"
#include "src/Functions/Algorithms/Tables/SortTable.hpp"
#include "src/Functions/Algorithms/Tables/TableBounds.hpp"
#include "src/Functions/Transformers/TableOperations/DistinctTable.hpp"
#include "src/Functions/Transformers/TableOperations/VectorDeleteTable.hpp"
#include "src/Functions/Transformers/TableOperations/FilterTable.hpp"
#include "src/Functions/Transformers/TableOperations/DifferenceTables.hpp"
#include "src/Functions/Transformers/TableOperations/IntersectionTables.hpp"
#include "src/Functions/Transformers/TableOperations/JoinTables.hpp"
#include "src/Functions/Transformers/TableOperations/ConvertTableColumnToSet.hpp"
#include "src/Functions/ReaderWriters/Tables/TextTable.hpp"
#include "src/Functions/Algorithms/Math/Scalar/AppliedStatisticsAlgorithms.hpp"
#include "src/Functions/Algorithms/Math/Scalar/DCDFLIB.hpp"
#include "src/Functions/Algorithms/Math/Scalar/DoubleValueFunctions.hpp"
#include "src/Functions/Algorithms/Math/Scalar/Hypergeometric.hpp"
#include "src/Functions/Algorithms/Math/Scalar/PoissonFindMinNum.hpp"
#include "src/Functions/Algorithms/Math/Tables/TableColumnStats.hpp"
#include "src/Functions/Algorithms/Math/Tables/TableColumnsStats.hpp"
#include "src/Functions/Algorithms/Math/Tables/TableNormalizeColumn.hpp"
#include "src/Functions/Algorithms/Math/Tables/TableHistogram.hpp"
#include "src/Functions/Algorithms/Math/Tables/TableSmoother.hpp"
#include "src/Functions/Algorithms/Math/Tables/CheckConstantTable.hpp"


// Below objects are derived from Table, but may also contain non-Table objects.
#include "src/Objects/Features.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesBEDReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesFileTableReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableRefGeneReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableRmskReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableEnsemblQualifiersReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableEnsemblStructuresReader.hpp"
#include "src/Functions/ReaderWriters/Features/FeaturesTableHumanProteinsReader.hpp"
#include "src/Objects/AlignmentElement.hpp"
#include "src/Objects/AlignmentFragment.hpp"
#include "src/Objects/Alignment.hpp"
#include "src/Objects/SequenceAlignmentFragmentToFeaturesTracker.hpp"
#include "src/Functions/Algorithms/Alignments/SmithWaterman.hpp"
#include "src/Functions/ReaderWriters/Trackers/SequenceAlignmentFragmentToFeatureTrackerWriter.hpp"
#include "src/Objects/Peaks.hpp"
#include "src/Functions/Algorithms/Peaks/LocalMaximumWindow.hpp"
#include "src/Functions/Algorithms/Peaks/RidgeLinePeaks.hpp"
#include "src/Functions/Algorithms/Peaks/CorrectPeakPosition.hpp"
#include "src/Functions/Algorithms/Peaks/PeakConstantSizedWidth.hpp"
#include "src/Functions/Algorithms/Peaks/PeakPercentAmpWidth.hpp"
#include "src/Functions/Algorithms/Peaks/CorrectPeaksWidth.hpp"


// Below are files based on some specification, standard, or implmentation.
#include "src/Objects/Fasta.hpp"
#include "src/Objects/ReadOnlySequencesFile.hpp"
#include "src/Functions/ReaderWriters/Tables/FastaTableReaderWriter.hpp"
#include "src/Objects/Qual.hpp"
#include "src/Functions/ReaderWriters/Tables/QualTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/FastqTableReaderWriter.hpp"
#include "src/Functions/ReaderWriters/Tables/ReadsTableReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/SequenceAlignmentFragmentWriter.hpp"
#include "src/Functions/ReaderWriters/Alignments/BroadQltoutAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/BroadQltoutAlignmentFragmentWriter.hpp"
#include "src/Functions/ReaderWriters/Alignments/HelicosSmsAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/NCBIBlastnTableAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/SAMAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/SAMAlignmentFragmentWriter.hpp"
#include "src/Functions/ReaderWriters/Alignments/FileSequenceAlignmentFragmentReader.hpp"
#include "src/Functions/ReaderWriters/Alignments/FileSequenceAlignmentFragmentWriter.hpp"


// Repare
#include "src/Projects/Repare/repare.hpp"


#endif

