#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMedianImageFilter.h>
#include <itkMorphologicalContourInterpolator.h>
#include <itkGradientMagnitudeImageFilter.h>
#include "itkImageRegionIterator.h"


#include "itkExtractImageFilter.h"
#include "itkPasteImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "gdcmUIDGenerator.h"

#include "itkNumericSeriesFileNames.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkLaplacianSegmentationLevelSetImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"


#include "mytypes.h"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>

#define isValid(x,y) ((x[0] > 0) && (x[0] < y[0]) && (x[1] > 0) && (x[1] < y[1]) && (x[2] > 0) && (x[2] < y[2]))
using MaskImageType = itk::Image<unsigned short, 3>;

// If fine_tune_mask is for example == 1 we will shrink and grow the mask once each and 
// add/remove any voxel with large image gradients.
MaskImageType::Pointer fineTune(MaskImageType::Pointer mask, std::string image_path, int fine_tune_mask, int N, int verbose) {

  // so the rules of the game are that we use the smallest stencil of 7 sampling points
  // We will want the center point to be inside the mask and at least one other point on the background.
  // For each voxel we will have a gradient value (smoothed derived from the image_path).
  // If a voxel from the mask has the highest gradient value we will remove the center voxel (set to background).
  // If a background voxel has the highest gradient we will add that one point to the mask.
  // If the center voxel has the highest gradient we will do nothing. 

  // make a copy of the mask image, lookup values in that and change the mask you return
  using DuplicatorType = itk::ImageDuplicator<MaskImageType>;
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage(mask);
  duplicator->Update();

  // Read in the image volume
  using ImageType3D = itk::Image<unsigned short, 3>;

  typedef itk::ImageSeriesReader<ImageType3D> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

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
  nameGenerator->SetDirectory(image_path);
  
  try {
    using SeriesIdContainer = std::vector<std::string>;
    const SeriesIdContainer& seriesUID = nameGenerator->GetSeriesUIDs();
    auto                      seriesItr = seriesUID.begin();
    auto                      seriesEnd = seriesUID.end();

    if (seriesItr == seriesEnd) {
      std::cout << "No DICOMs in: " << image_path << std::endl;
      return mask;
    }

    seriesItr = seriesUID.begin();
    while (seriesItr != seriesUID.end()) {
      std::string seriesIdentifier = seriesItr->c_str();
      seriesItr++;
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
        return mask;
      }
    }
  } catch (const itk::ExceptionObject& ex) {
    std::cout << ex << std::endl;
    return mask;
  }

  // now reader->GetOutput() has our image, should have the same size as our mask (assert?)
  ImageType3D::Pointer image = reader->GetOutput();
  const ImageType3D::RegionType& imageRegion = image->GetLargestPossibleRegion();
  const MaskImageType::RegionType& maskRegion = mask->GetLargestPossibleRegion();
  const ImageType3D::SizeType& imageSize = imageRegion.GetSize();
  const MaskImageType::SizeType& maskSize = maskRegion.GetSize();

  if (imageSize[0] != maskSize[0] || imageSize[1] != maskSize[1] || imageSize[2] != maskSize[2]) {
    fprintf(stderr, "Error: could not finetune the mask. Mask volume and provided image volume do not match in dimensions.\n");
    return mask;
  }

  using InternalPixelType = float;
  using InternalImageType = itk::Image< InternalPixelType, 3 >;

  using DiffusionFilterType =
  itk::GradientAnisotropicDiffusionImageFilter< ImageType3D,
                                                InternalImageType >;
  DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
  diffusion->SetNumberOfIterations( N );
  diffusion->SetTimeStep(0.05);
  diffusion->SetConductanceParameter( 2.0 );

  using LaplacianSegmentationLevelSetImageFilterType =
    itk::LaplacianSegmentationLevelSetImageFilter< MaskImageType,
            InternalImageType >;
  LaplacianSegmentationLevelSetImageFilterType::Pointer laplacianSegmentation
            = LaplacianSegmentationLevelSetImageFilterType::New();

  laplacianSegmentation->SetCurvatureScaling( 1.0 );
  laplacianSegmentation->SetPropagationScaling( 1.0 );

  laplacianSegmentation->SetMaximumRMSError( 0.002 );
  laplacianSegmentation->SetNumberOfIterations( N );

  laplacianSegmentation->SetIsoSurfaceValue( 0.5 );

  using ThresholdingFilterType =
      itk::BinaryThresholdImageFilter< InternalImageType, MaskImageType >;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();

  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );
  
  thresholder->SetOutsideValue(  1  );
  thresholder->SetInsideValue(  0 );


  diffusion->SetInput( image );
  laplacianSegmentation->SetInput( mask );
  laplacianSegmentation->SetFeatureImage( diffusion->GetOutput() );
  thresholder->SetInput( laplacianSegmentation->GetOutput() );
  thresholder->Update();

  return thresholder->GetOutput();
}


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
  std::string versionString = std::string("0.0.2.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  if (versionString.find("..") != std::string::npos)
    versionString.replace(versionString.find(".."), 2, ".");
  command.SetVersion(versionString.c_str());
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("MorphologicalContourInterpolation: Creates an interpolated volume label from individual slice segmentations.");
  command.SetCategory("mask editing");
  command.AddField("indir", "Directory with input DICOM image (mask) series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM image series.", MetaCommand::STRING, true);

  command.SetOption(
    "UIDFixed", "u", false,
    "If enabled identifiers are stable - will not change for a given input. This allows image series to overwrite each other - assuming that the PACS "
    "supports this overwrite mode. By default the SeriesInstanceUID and SOPInstanceUID values are generated again every time the processing is done.");
  command.SetOptionLongTag("UIDFixed", "uid-fixed");

  command.SetOption("MaxNumberOfThreads", "t", false, "Use at most X (4) threads for computation.");
  command.SetOptionLongTag("MaxNumberOfThreads", "maxnumberofthreads");
  command.AddOptionField("MaxNumberOfThreads", "maxnumberofthreads", MetaCommand::INT, false);

  command.SetOption("FineTuneMask", "f", false, "Adjust the mask (N=10). This option requires that the image series is also provided (option -i).");
  command.SetOptionLongTag("FineTuneMask", "fine-tune-mask");
  command.AddOptionField("FineTuneMask", "value", MetaCommand::INT, false);

  // in case we want to fine tune the mask we need the image series as well
  command.SetOption("ImageSeries", "i", false, "In case the input masks needs to be fine-tuned (-f) we need to image series as well.");
  command.SetOptionLongTag("ImageSeries", "image-series");
  command.AddOptionField("ImageSeries", "value", MetaCommand::STRING, false);

  // convert a specific series
  std::string convertSpecificSeries = "";

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  command.SetOption("Version", "V", false, "Print version information.");
  command.SetOptionLongTag("Version", "version");

  if (argc == 2 && std::string(argv[1]) == std::string("--version")) {
    fprintf(stdout, "Version: %s\n", versionString.c_str());
    return 0;
  }

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

  std::string image_path = "";
  if (command.GetOptionWasSet("ImageSeries")) {
    image_path = command.GetValueAsString("ImageSeries", "value");
    // check if the path exists
    boost::filesystem::path image_path_p = image_path;
    if (!itksys::SystemTools::FileIsDirectory(image_path_p.c_str())) {
      // create the output directory
      fprintf(stderr, "Error: the provided image path (-i) \"%s\" could not be found.\n", image_path_p.c_str());
      exit(-1);
    }
  }

  int fine_tune_mask = 0;
  if (command.GetOptionWasSet("FineTuneMask")) {
    if (image_path.size() == 0) {
      fprintf(stderr, "set an ImageSeries if you want to finetune.\n");
      exit(-1);
    }
    fine_tune_mask = command.GetValueAsInt("FineTuneMask", "value");
    if (fine_tune_mask < 0)
      fine_tune_mask = - fine_tune_mask;
  }
  if (verbose && fine_tune_mask != 0) {
    fprintf(stdout, "use fine tuning of N=%d for mask post-processing\n", fine_tune_mask);
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

      // create an output folder based on the seriesIdentifier
      boost::filesystem::path p_out = output_path + boost::filesystem::path::preferred_separator + seriesIdentifier;
      if (!itksys::SystemTools::FileIsDirectory(p_out.c_str())) {
        // create the output directory
        create_directories(p_out);
      } 

      // In order to speed up this processing we should extract the sub-region of the label, process on the sub-volume
      // and put the result back into the input image before saving it.
      //   see: https://itk.org/Doxygen/html/Examples_2IO_2ImageReadExtractFilterInsertWrite_8cxx-example.html#_a4

      //
      // start by computing the minimal bounding box around our regions of interest
      //
      std::vector<int> boundingBox{0,0,0,0,0,0};
      bool boundingBoxValid = false;
      using IteratorTypeImage = itk::ImageRegionIteratorWithIndex< MaskImageType >;
      const MaskImageType::RegionType& inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
      itk::ImageRegionIteratorWithIndex<MaskImageType> iter(reader->GetOutput(), inputRegion);
      iter.GoToBegin();
      while( !iter.IsAtEnd() ) {
        MaskImageType::PixelType value = iter.Value();
        if (value > 0) {
          MaskImageType::IndexType idx = iter.GetIndex();
          if (!boundingBoxValid) {
            // first time init
            boundingBox[0] = idx[0];
            boundingBox[1] = idx[1];
            boundingBox[2] = idx[2];
            boundingBox[3] = boundingBox[0];
            boundingBox[4] = boundingBox[1];
            boundingBox[5] = boundingBox[2];
            boundingBoxValid = true;
          }
          if (idx[0] < boundingBox[0])
            boundingBox[0] = idx[0];
          if (idx[1] < boundingBox[1])
            boundingBox[1] = idx[1];
          if (idx[2] < boundingBox[2])
            boundingBox[2] = idx[2];
          if (idx[0] > boundingBox[3])
            boundingBox[3] = idx[0];
          if (idx[1] > boundingBox[4])
            boundingBox[4] = idx[1];
          if (idx[2] > boundingBox[5])
            boundingBox[5] = idx[2];
        }
        ++iter;
      }
      if (verbose) {
        fprintf(stdout, "found minimum enclosing bounding box: %d,%d,%d..%d,%d,%d\n", boundingBox[0], boundingBox[1], boundingBox[2], boundingBox[3], boundingBox[4], boundingBox[5]);
      }
      // extend the bounding box (based on fine_tune_mask)
      boundingBox[0] -= 1 + fine_tune_mask;
      if (boundingBox[0] < 0)
        boundingBox[0] = 0;
      boundingBox[1] -= 1 + fine_tune_mask;
      if (boundingBox[1] < 0)
        boundingBox[1] = 0;
      boundingBox[2] -= 1 + fine_tune_mask;
      if (boundingBox[2] < 0)
        boundingBox[2] = 0;
      boundingBox[3] += 1 + fine_tune_mask;
      if (boundingBox[3] >= inputRegion.GetSize()[0])
        boundingBox[3] = inputRegion.GetSize()[0];
      boundingBox[4] += 1 + fine_tune_mask;
      if (boundingBox[4] >= inputRegion.GetSize()[1])
        boundingBox[4] = inputRegion.GetSize()[1];
      boundingBox[5] += 1 + fine_tune_mask;
      if (boundingBox[5] >= inputRegion.GetSize()[2])
        boundingBox[5] = inputRegion.GetSize()[2];

      MaskImageType::SizeType roi_size = inputRegion.GetSize();
      roi_size[0] = boundingBox[3]-boundingBox[0]; // TODO: do we have to add 1 or 2 here?
      roi_size[1] = boundingBox[4]-boundingBox[1];
      roi_size[2] = boundingBox[5]-boundingBox[2];
      MaskImageType::IndexType roi_start = inputRegion.GetIndex();
      roi_start[0] = boundingBox[0];
      roi_start[1] = boundingBox[1];
      roi_start[2] = boundingBox[2];

      // next step is to copy this region into another (smaller) volume
      // Hope is that in a smaller volume all following computations are faster.
      using ExtractFilterType = itk::ExtractImageFilter<MaskImageType, MaskImageType>;
      auto extractFilter = ExtractFilterType::New();
      //extractFilter->SetDirectionCollapseToSubmatrix();
      const MaskImageType *     inputImage = reader->GetOutput();

      MaskImageType::RegionType desiredRegion;
      desiredRegion.SetSize(roi_size);
      desiredRegion.SetIndex(roi_start);
      extractFilter->SetInput(inputImage);
      extractFilter->SetExtractionRegion(desiredRegion);


      using mciType = itk::MorphologicalContourInterpolator<MaskImageType>;
      mciType::Pointer mci = mciType::New();
      mci->SetInput(extractFilter->GetOutput() /*reader->GetOutput() */);
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


      using PasteFilterType = itk::PasteImageFilter<MaskImageType, MaskImageType>;
      auto pasteFilter = PasteFilterType::New();
      pasteFilter->SetSourceImage(medF->GetOutput());
      pasteFilter->SetDestinationImage(inputImage);
      pasteFilter->SetDestinationIndex(roi_start);
      const MaskImageType * medianImage = medF->GetOutput();
      pasteFilter->SetSourceRegion(medianImage->GetBufferedRegion());
      pasteFilter->Update();

      // In case we want to fine-tune the mask we can do this here. We would need to
      // read the image series and see if it matches with the mask volume.
      MaskImageType::Pointer fine_tuned_mask;
      if (fine_tune_mask > 0) {
        fine_tuned_mask = fineTune(pasteFilter->GetOutput(), image_path, fine_tune_mask, fine_tune_mask, verbose);
      }

      // Instead of writing a single file, we want to write out a new DICOM series
      // but keep all the input DICOM tags in place. Or at least make them compatible
      // with '-u'.

      if (0) {
        // if we would want to save a single file output
        using WriterType = itk::ImageFileWriter<MaskImageType>;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(p_out.c_str());
        writer->SetInput(pasteFilter->GetOutput() /*medF->GetOutput()*/);
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

      const MaskImageType::RegionType& inputRegion2 = reader->GetOutput()->GetLargestPossibleRegion();
      const MaskImageType::IndexType    start = inputRegion2.GetIndex();
      const MaskImageType::SizeType& inputSize = inputRegion2.GetSize();

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
      std::string seriesFormat(p_out.c_str());
      seriesFormat = seriesFormat + "/" + "IM%04d.dcm";
      outputNames->SetSeriesFormat(seriesFormat.c_str());
      const unsigned int firstSlice = start[2];
      const unsigned int lastSlice = start[2] + inputSize[2] - 1;
      outputNames->SetStartIndex(firstSlice);
      outputNames->SetEndIndex(lastSlice);
      outputNames->SetIncrementIndex(1);

      auto seriesWriter = SeriesWriterType::New();
      if (fine_tune_mask > 0)
        seriesWriter->SetInput(fine_tuned_mask);
      else
        seriesWriter->SetInput(pasteFilter->GetOutput() /*medF->GetOutput()*/);

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
