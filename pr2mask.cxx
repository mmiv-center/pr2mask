#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
//#include "itkExtractImageFilter.h"
//#include "itkPasteImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRGBPixel.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
//#include <itkPixelAccessor.h>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDataSetHelper.h"
#include "gdcmFileDerivation.h"
#include "gdcmFileExplicitFilter.h"
#include "gdcmGlobal.h"
#include "gdcmImageApplyLookupTable.h"
#include "gdcmImageChangePlanarConfiguration.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmMediaStorage.h"
#include "gdcmRescaler.h"
#include "gdcmStringFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <map>

#include "mytypes.h"

using json = nlohmann::json;
using namespace boost::filesystem;

// forward declaration
void CopyDictionary(itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

template <typename TFilter> class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() {}

public:
  virtual void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE { Execute((const itk::Object *)caller, event); }

  virtual void Execute(const itk::Object *object, const itk::EventObject &event) ITK_OVERRIDE {
    const TFilter *filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) != typeid(itk::IterationEvent)) {
      return;
    }
    if (filter->GetElapsedIterations() == 1) {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1 << std::endl;
    }
    std::cout << "  Iteration " << filter->GetElapsedIterations() << " (of " << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()] << ").  ";
    std::cout << " Current convergence value = " << filter->GetCurrentConvergenceMeasurement() << " (threshold = " << filter->GetConvergenceThreshold() << ")"
              << std::endl;
  }
};

template <typename TValue> TValue Convert(std::string optionString) {
  TValue value;
  std::istringstream iss(optionString);

  iss >> value;
  return value;
}

template <typename TValue> std::vector<TValue> ConvertVector(std::string optionString) {
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find('x', 0);

  if (crosspos == std::string::npos) {
    values.push_back(Convert<TValue>(optionString));
  } else {
    std::string element = optionString.substr(0, crosspos);
    TValue value;
    std::istringstream iss(element);
    iss >> value;
    values.push_back(value);
    while (crosspos != std::string::npos) {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find('x', crossposfrom + 1);
      if (crosspos == std::string::npos) {
        element = optionString.substr(crossposfrom + 1, optionString.length());
      } else {
        element = optionString.substr(crossposfrom + 1, crosspos);
      }
      std::istringstream iss2(element);
      iss2 >> value;
      values.push_back(value);
    }
  }
  return values;
}

json resultJSON;

int main(int argc, char *argv[]) {
  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("PR2NII: Convert presentation state files with polygons to nii.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output image series.", MetaCommand::STRING, true);

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);

  command.SetOption("Force", "f", false, "Ignore existing directories and force reprocessing. Default is to stop processing if directory already exists.");

  command.SetOption("Verbose", "V", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  bool seriesIdentifierFlag = false;
  std::string input = command.GetValueAsString("indir");
  std::string output = command.GetValueAsString("outdir");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;
  bool force = false;
  if (command.GetOptionWasSet("Force"))
    force = true;

  if (command.GetOptionWasSet("SeriesName"))
    seriesIdentifierFlag = true;
  std::string labelfieldfilename = command.GetValueAsString("SaveLabelfield", "labelfieldfilename");
  // todo: the argument could not be there, in this case the labelfieldfilename might be empty

  std::string seriesName = command.GetValueAsString("SeriesName", "seriesname");

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  //  typedef signed short PixelType;
  typedef float FloatPixelType;
  const unsigned int Dimension = 3;

  //  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::Image<FloatPixelType, Dimension> FloatImageType;

  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::Image<unsigned char, Dimension> MaskImageType;
  typedef itk::Image<unsigned char, 2> MaskSliceImageType;

  using StructuringElementType = itk::BinaryBallStructuringElement<PixelType, Dimension>;
  using ErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;
  using DilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;

  typedef itk::ConnectedComponentImageFilter<MaskImageType, ImageType> ConnectedComponentImageFilterType;

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  dicomIO->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO);

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetRecursive(true);
  nameGenerator->SetDirectory(input);

  try {
    typedef std::vector<std::string> SeriesIdContainer;

    const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while (seriesItr != seriesEnd) {
      std::cout << "Found DICOM Series: ";
      std::cout << std::endl;
      std::cout << "  " << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }

    std::string seriesIdentifier;

    SeriesIdContainer runThese;
    if (seriesIdentifierFlag) { // If no optional series identifier
      // seriesIdentifier = seriesName;
      runThese.push_back(seriesName);
    } else {
      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
      // seriesIdentifier = seriesUID.begin()->c_str();
      while (seriesItr != seriesEnd) {
        runThese.push_back(seriesItr->c_str());
        ++seriesItr;
      }
    }

    // we need to identify for each series if they are a presentation state object and if we can extract some
    // contours from them

    seriesItr = runThese.begin();
    seriesEnd = runThese.end();
    while (seriesItr != seriesEnd) {
      seriesIdentifier = seriesItr->c_str();
      ++seriesItr;

      std::cout << "Processing series: " << std::endl;
      std::cout << "  " << seriesIdentifier << std::endl;

      std::string outputSeries = output + "/" + seriesIdentifier;
      if (!force && itksys::SystemTools::FileIsDirectory(outputSeries.c_str())) {
        fprintf(stdout, "Skip this series %s, output directory exists already...\n", outputSeries.c_str());
        exit(0); // this is no skip, that is giving up...
      }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      if (fileNames.size() < 2) {
        std::cout << "skip processing, not enough images in this series..." << std::endl;
        continue;
      }
      fprintf(stdout, "sufficient number of images (%lu) in this series\n", fileNames.size());
      resultJSON["series_identifier"] = seriesIdentifier;

      reader->SetFileNames(fileNames);
      reader->ForceOrthogonalDirectionOff(); // do we need this?

      try {
        reader->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }

      // read the data dictionary
      ImageType::Pointer inputImage = reader->GetOutput();
      typedef itk::MetaDataDictionary DictionaryType;
      DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();
      fprintf(stdout, "pixel spacing of input is: %f %f %f\n", inputImage->GetSpacing()[0], inputImage->GetSpacing()[1], inputImage->GetSpacing()[2]);

      std::string studyDescription;
      std::string seriesDescription;
      std::string patientID;
      std::string patientName;
      std::string sopClassUID;
      std::string seriesDate;
      std::string seriesTime;
      std::string studyDate;
      std::string studyTime;
      std::string patientSex;
      std::string convolutionKernelGroup;
      std::string modality;
      std::string manufacturer;
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|1030", studyDescription))
        resultJSON["SeriesDescription"] = boost::algorithm::trim_copy(seriesDescription);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|103e", seriesDescription))
        resultJSON["StudyDescription"] = boost::algorithm::trim_copy(studyDescription);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID))
        resultJSON["SOPClassUID"] = boost::algorithm::trim_copy(sopClassUID);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0021", seriesDate))
        resultJSON["StudyDate"] = boost::algorithm::trim_copy(studyDate);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0031", seriesTime))
        resultJSON["SeriesTime"] = boost::algorithm::trim_copy(seriesTime);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0020", patientID))
        resultJSON["PatientID"] = boost::algorithm::trim_copy(patientID);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0010", patientName))
        resultJSON["PatientName"] = boost::algorithm::trim_copy(patientName);
      if (itk::ExposeMetaData<std::string>(dictionary, "0010|0040", patientSex))
        resultJSON["PatientSex"] = boost::algorithm::trim_copy(patientSex);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0030", studyTime))
        resultJSON["StudyTime"] = boost::algorithm::trim_copy(studyTime);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0020", studyDate))
        resultJSON["SeriesDate"] = boost::algorithm::trim_copy(seriesDate);
      if (itk::ExposeMetaData<std::string>(dictionary, "0018|9316", convolutionKernelGroup))
        resultJSON["CTConvolutionKernelGroup"] = boost::algorithm::trim_copy(convolutionKernelGroup); // LUNG, BRAIN, BONE, SOFT_TISSUE
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0060", modality))
        resultJSON["Modality"] = boost::algorithm::trim_copy(modality);
      if (itk::ExposeMetaData<std::string>(dictionary, "0008|0080", manufacturer))
        resultJSON["Manufacturer"] = boost::algorithm::trim_copy(manufacturer);

    } // loop over series
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  std::string si(resultJSON["series_identifier"]);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  o << output << "/" << si << ".json";
  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}

void CopyDictionary(itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict) {
  typedef itk::MetaDataDictionary DictionaryType;

  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  typedef itk::MetaDataObject<std::string> MetaDataStringType;

  while (itr != end) {
    itk::MetaDataObjectBase::Pointer entry = itr->second;

    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    if (entryvalue) {
      std::string tagkey = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
    }
    ++itr;
  }
}
