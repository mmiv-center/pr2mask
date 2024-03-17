#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMedianImageFilter.h>
#include <itkMorphologicalContourInterpolator.h>

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>

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

  // convert a specific series
  std::string convertSpecificSeries = "";

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  bool uidFixedFlag = false;
  if (command.GetOptionWasSet("UIDFixed"))
    uidFixedFlag = true;

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  std::string input_path = command.GetValueAsString("indir");
  std::string output_path = command.GetValueAsString("outdir");

  using MaskImageType = itk::Image<unsigned short, 3>;
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
      } else // otherwise convert everything
      {
        seriesIdentifier = seriesItr->c_str();
        seriesItr++;
      }
      std::cout << "Reading: ";
      std::cout << seriesIdentifier << std::endl;

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

      using WriterType = itk::ImageFileWriter<MaskImageType>;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(output_path);
      writer->SetInput(medF->GetOutput());
      writer->SetUseCompression(true);
      writer->Update();
    }
  } catch (const itk::ExceptionObject& ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
