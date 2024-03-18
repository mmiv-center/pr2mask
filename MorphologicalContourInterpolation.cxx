#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMedianImageFilter.h>
#include <itkMorphologicalContourInterpolator.h>

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "gdcmUIDGenerator.h"

#include "itkNumericSeriesFileNames.h"



#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>



void CopyDictionary(itk::MetaDataDictionary& fromDict, itk::MetaDataDictionary& toDict) {
  using DictionaryType = itk::MetaDataDictionary;

  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  using MetaDataStringType = itk::MetaDataObject<std::string>;

  while (itr != end) {
    itk::MetaDataObjectBase::Pointer entry = itr->second;

    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType*>(entry.GetPointer());
    if (entryvalue) {
      std::string tagkey = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
    }
    ++itr;
  }
}


int main(int argc, char* argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");
  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  std::string versionString = std::string("0.0.1.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  command.SetVersion(versionString.c_str());
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("MorphologicalContourInterpolation: Creates an interpolated volume label from individual slice segmentations.");
  command.SetCategory("mask editing");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM image series.", MetaCommand::STRING, true);

  command.SetOption(
    "UIDFixed", "u", false,
    "If enabled identifiers are stable - will not change for a given input. This allows image series to overwrite each other - assuming that the PACS "
    "supports this overwrite mode. By default the SeriesInstanceUID and SOPInstanceUID values are generated again every time the processing is done.");
  command.SetOptionLongTag("UIDFixed", "uid-fixed");

  command.SetOption("MaxNumberOfThreads", "t", 4, "Use at most X (4) threads for computation.");
  command.SetOptionLongTag("MaxNumberOfThreads", "maxnumberofthreads");
  command.AddOptionField("MaxNumberOfThreads", "maxnumberofthreads", MetaCommand::INT, 4);


  // convert a specific series
  std::string convertSpecificSeries = "";

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  // be nice
  int maxThreads = 4;
  if (command.GetOptionWasSet("MaxNumberOfThreads")) {
    maxThreads = command.GetValueAsInt("MaxNumberOfThreads");
    if (maxThreads < 1)
      maxThreads = 1;
  } 
  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(maxThreads);

  bool uidFixedFlag = false;
  if (command.GetOptionWasSet("UIDFixed"))
    uidFixedFlag = true;

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  std::string input_path = command.GetValueAsString("indir");
  std::string output_path = command.GetValueAsString("outdir");
  if (input_path == "" || output_path == "") {
    fprintf(stderr, "Error: no input or output directory specified.\n");
    return 1;
  }


  using MaskImageType = itk::Image<unsigned short, 3>;
  using OutputImageType = itk::Image<unsigned short, 2>;
  typedef itk::ImageSeriesReader<MaskImageType> MaskReaderType;
  MaskReaderType::Pointer reader = MaskReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  dicomIO->LoadPrivateTagsOn();

  ImageIOType::Pointer dicomIOImage = ImageIOType::New();
  dicomIOImage->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO);

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(false); // we want to use the keys as SeriesInstanceUIDs
  nameGenerator->AddSeriesRestriction("0008|0060");
  nameGenerator->SetRecursive(true);
  nameGenerator->SetDirectory(input_path);

  // for each found series instance uid do the following
  try {
    using SeriesIdContainer = std::vector<std::string>;
    const SeriesIdContainer& seriesUID = nameGenerator->GetSeriesUIDs();
    auto                      seriesItr = seriesUID.begin();
    auto                      seriesEnd = seriesUID.end();

    if (seriesItr == seriesEnd) {
      std::cout << "No DICOMs in: " << input_path << std::endl;
      return EXIT_SUCCESS;
    }

    seriesItr = seriesUID.begin();
    while (seriesItr != seriesUID.end()) {
      std::string seriesIdentifier;
      if (convertSpecificSeries != "") {
        seriesIdentifier = convertSpecificSeries;
        seriesItr = seriesUID.end();
      } else {
        seriesIdentifier = seriesItr->c_str();
        seriesItr++;
      }
      if (verbose) {
        std::cout << "Reading: ";
        std::cout << seriesIdentifier << std::endl;
      }
      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames; // for the label series

      // labelSeries could be the SeriesInstanceUID
      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      reader->SetFileNames(fileNames);

      try {
        reader->Update();
      } catch (itk::ExceptionObject& ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }

      using mciType = itk::MorphologicalContourInterpolator<MaskImageType>;
      mciType::Pointer mci = mciType::New();
      mci->SetInput(reader->GetOutput());
      bool UseDistanceTransform = true;
      bool ball = true;
      int axis = -1; // all axis
      int label = 1; // use label 1 (0 would be all labels)
      mci->SetUseDistanceTransform(UseDistanceTransform);
      mci->SetUseBallStructuringElement(ball);
      mci->SetAxis(axis);
      mci->SetLabel(label);

      mci->Update();

      int smoothingRadius = 2;

      using MedianType = itk::MedianImageFilter<MaskImageType, MaskImageType>;
      MedianType::Pointer medF = MedianType::New();
      medF->SetInput(mci->GetOutput());
      medF->SetRadius(smoothingRadius);
      medF->Update();

      // Instead of writing a single file, we want to write out a new DICOM series
      // but keep all the input DICOM tags in place. Or at least make them compatible
      // with '-u'.

      if (0) {
        // if we would want to save a single file output
        using WriterType = itk::ImageFileWriter<MaskImageType>;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(output_path);
        writer->SetInput(medF->GetOutput());
        writer->SetUseCompression(true);
        writer->Update();
      }
      //
      // save as individual DICOM files again
      //
      MaskReaderType::DictionaryRawPointer inputDict = (*(reader->GetMetaDataDictionaryArray()))[0];
      MaskReaderType::DictionaryArrayType  outputArray;

      std::string newSeriesInstanceUID("");
      if (uidFixedFlag) {
        std::string derivedSeriesInstanceUID(seriesIdentifier);
        std::string endString = ".1";
        if (derivedSeriesInstanceUID.substr(derivedSeriesInstanceUID.size() - 2, 2) == ".1")
          endString = ".2";

        // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
        derivedSeriesInstanceUID = derivedSeriesInstanceUID.substr(0, 64 - 3) + endString;
        newSeriesInstanceUID = derivedSeriesInstanceUID;
      } else {
        if (newSeriesInstanceUID == "") {
          gdcm::UIDGenerator uid;
          uid.SetRoot("1.3.6.1.4.1.45037");
          newSeriesInstanceUID = std::string(uid.Generate());
        } // keep reusing else
      }

      gdcm::UIDGenerator suid;
      suid.SetRoot("1.3.6.1.4.1.45037");
      // std::string        seriesUID = suid.Generate();
      //gdcm::UIDGenerator fuid;
      //fuid.SetRoot("1.3.6.1.4.1.45037");
      //std::string        frameOfReferenceUID = fuid.Generate();

      const MaskImageType::RegionType& inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
      const MaskImageType::IndexType    start = inputRegion.GetIndex();
      const MaskImageType::SizeType& inputSize = inputRegion.GetSize();

      std::string studyUID;
      //std::string sopClassUID;
      itk::ExposeMetaData<std::string>(*inputDict, "0020|000d", studyUID);
      //itk::ExposeMetaData<std::string>(*inputDict, "0008|0016", sopClassUID);
      dicomIO->KeepOriginalUIDOn();

      for (unsigned int f = 0; f < inputSize[2]; ++f) { // save one DICOM for each slice
        inputDict = (*(reader->GetMetaDataDictionaryArray()))[f];

        // Create a new dictionary for this slice
        auto dict = new MaskReaderType::DictionaryType;

        // Copy the dictionary from the first slice
        CopyDictionary(*inputDict, *dict);

        // Set the UID's for the study, series, SOP  and frame of reference
        itk::EncapsulateMetaData<std::string>(*dict, "0020|000d", studyUID);
        itk::EncapsulateMetaData<std::string>(*dict, "0020|000e", newSeriesInstanceUID);
        std::string oldSOPInstanceUID("");
        itk::ExposeMetaData<std::string>(*inputDict, "0008|0018", oldSOPInstanceUID);

        std::string newSOPInstanceUID("");
        if (uidFixedFlag) {
          std::string derivedSOPInstanceUID(oldSOPInstanceUID);
          std::string endString = ".4";
          if (derivedSOPInstanceUID.substr(derivedSOPInstanceUID.size() - 2, 2) == ".4")
            endString = ".5";

          // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
          derivedSOPInstanceUID = derivedSOPInstanceUID.substr(0, 64 - 3) + endString;
          newSOPInstanceUID = derivedSOPInstanceUID;
        } else {
          gdcm::UIDGenerator uid;
          uid.SetRoot("1.3.6.1.4.1.45037");
          newSOPInstanceUID = std::string(uid.Generate());
        }


        // std::string sopInstanceUID = suid.Generate();
        itk::EncapsulateMetaData<std::string>(*dict, "0008|0018", newSOPInstanceUID);
        itk::EncapsulateMetaData<std::string>(*dict, "0002|0003", newSOPInstanceUID);

        std::string oldSeriesDesc;
        itk::ExposeMetaData<std::string>(*inputDict, "0008|103e", oldSeriesDesc);

        std::ostringstream value;
        value.str("");
        value << oldSeriesDesc << " (interpolated)";
        unsigned lengthDesc = value.str().length();
        std::string seriesDesc(value.str(), 0, lengthDesc > 64 ? 64 : lengthDesc);
        itk::EncapsulateMetaData<std::string>(*dict, "0008|103e", seriesDesc);

        std::string oldSeriesNumber;
        itk::ExposeMetaData<std::string>(*inputDict, "0020|0011", oldSeriesNumber);
        value.str("");
        value << oldSeriesNumber << "1";
        itk::EncapsulateMetaData<std::string>(*dict, "0020|0011", value.str());

        // add how the image was derived
        value.str("");
        value << "Mask volume generated by MorphologicalContourInterpolation " << versionString;
        lengthDesc = value.str().length();
        std::string derivationDesc(value.str(), 0, lengthDesc > 1024 ? 1024 : lengthDesc);
        itk::EncapsulateMetaData<std::string>(*dict, "0008|2111", derivationDesc);

        // Save the dictionary
        outputArray.push_back(dict);
      }

      // Make the output directory and generate the file names.
      itksys::SystemTools::MakeDirectory(output_path);

      // Generate the file names
      using SeriesWriterType = itk::ImageSeriesWriter<MaskImageType, OutputImageType>;
      using OutputNamesGeneratorType = itk::NumericSeriesFileNames;
      auto        outputNames = OutputNamesGeneratorType::New();
      std::string seriesFormat(output_path);
      seriesFormat = seriesFormat + "/" + "IM%04d.dcm";
      outputNames->SetSeriesFormat(seriesFormat.c_str());
      const unsigned int firstSlice = start[2];
      const unsigned int lastSlice = start[2] + inputSize[2] - 1;
      outputNames->SetStartIndex(firstSlice);
      outputNames->SetEndIndex(lastSlice);
      outputNames->SetIncrementIndex(1);

      auto seriesWriter = SeriesWriterType::New();
      seriesWriter->SetInput(medF->GetOutput());

      seriesWriter->SetImageIO(dicomIO);
      seriesWriter->SetFileNames(outputNames->GetFileNames());
      seriesWriter->SetMetaDataDictionaryArray(&outputArray);
      try {
        seriesWriter->Update();
      } catch (const itk::ExceptionObject& excp) {
        std::cerr << "Exception thrown while writing the series " << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
      }

    }
  } catch (const itk::ExceptionObject& ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
