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
#include "gdcmDirectoryHelper.h"
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
#include "gdcmReader.h"
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

// We need to identify for each series if they are a presentation state object and if we can extract some
// contours from them.
struct Polygon {
  std::vector<float> coords;            // the vector of pixel coordinates extracted from presentation state
  std::string ReferencedSOPInstanceUID; // the referenced SOP instance UID (identifies the image)
  std::string ReferencedSeriesInstanceUID;
  std::string StudyInstanceUID;
  std::string SeriesInstanceUID;
  std::string Filename; // the name of the DICOM file
};

bool parseForPolygons(std::string input, std::vector<Polygon> *storage, std::map<std::string, std::string> *SOPInstanceUID2SeriesInstanceUID) {
  // read from input and create polygon structs in storage
  // a local cache of the Series that might be referenced

  for (boost::filesystem::recursive_directory_iterator end, dir(input); dir != end; ++dir) {
    // std::cout << *dir << "\n";  // full path
    // std::cout << dir->path().filename() << "\n"; // just last bit
    std::string filename = dir->path().string();

    // Instanciate the reader:
    gdcm::Reader reader;
    reader.SetFileName(filename.c_str());
    if (!reader.Read()) {
      std::cerr << "Could not read  \"" << filename << "\" as DICOM, ignore." << std::endl;
      continue;
    }

    // The output of gdcm::Reader is a gdcm::File
    gdcm::File &file = reader.GetFile();

    // the dataset is the the set of element we are interested in:
    gdcm::DataSet &ds = file.GetDataSet();

    // is this Modality PR?
    std::string SeriesInstanceUID;
    std::string StudyInstanceUID;
    std::string SOPInstanceUID;
    std::string Modality;

    const gdcm::Tag graphicType(0x0070, 0x0023);              // GraphicType
    const gdcm::Tag numberOfGraphicPoints(0x0070, 0x0021);    // NumberOfGraphicPoints
    const gdcm::Tag graphicData(0x0070, 0x0022);              // GraphicData - n-tupel of Single
    const gdcm::Tag graphicObjectSequence(0x0070, 0x0009);    // GraphicObjectSequence
    const gdcm::Tag referencedImageSequence(0x0008, 0x1140);  // ReferencedImageSequence its inside 0x0070,0x0001
    const gdcm::Tag referencedSOPInstanceUID(0x0008, 0x1155); // ReferencedSOPInstanceUID

    gdcm::Attribute<0x0008, 0x0060> modalityAttr;
    modalityAttr.Set(ds);
    if (modalityAttr.GetValue() != std::string("PR")) { // PR's might reference these
      // we should create a cache here for the SeriesInstanceUID the ReferencedSOPInstanceUID might point to
      gdcm::Attribute<0x0020, 0x000E> seriesinstanceuidAttr;
      seriesinstanceuidAttr.Set(ds);
      SeriesInstanceUID = seriesinstanceuidAttr.GetValue();
      if (SeriesInstanceUID.back() == '\0')
        SeriesInstanceUID.replace(SeriesInstanceUID.end() - 1, SeriesInstanceUID.end(), "");

      gdcm::Attribute<0x0008, 0x0018> sopInstanceUIDAttr;
      sopInstanceUIDAttr.Set(ds);
      SOPInstanceUID = sopInstanceUIDAttr.GetValue();
      if (SOPInstanceUID.back() == '\0')
        SOPInstanceUID.replace(SOPInstanceUID.end() - 1, SOPInstanceUID.end(), "");

      // fprintf(stdout, " add to cache %s -> %s\n", SOPInstanceUID.c_str(), SeriesInstanceUID.c_str());
      SOPInstanceUID2SeriesInstanceUID->insert(std::pair<std::string, std::string>(SOPInstanceUID, SeriesInstanceUID));
      continue;
    }

    // get the StudyInstanceUID
    gdcm::Attribute<0x0020, 0x000D> studyinstanceuidAttr;
    studyinstanceuidAttr.Set(ds);
    StudyInstanceUID = studyinstanceuidAttr.GetValue();
    // fprintf(stdout, " found StudyInstanceUID: %s\n", StudyInstanceUID.c_str());

    // get the StudyInstanceUID
    gdcm::Attribute<0x0020, 0x000E> seriesinstanceuidAttr;
    seriesinstanceuidAttr.Set(ds);
    SeriesInstanceUID = seriesinstanceuidAttr.GetValue();
    // fprintf(stdout, " found SeriesInstanceUID: %s\n", SeriesInstanceUID.c_str());

    const gdcm::DataElement &de = ds.GetDataElement(gdcm::Tag(0x0070, 0x0001));
    // SequenceOfItems * sqi = (SequenceOfItems*)de.GetSequenceOfItems();
    gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = de.GetValueAsSQ();
    if (sqi) {

      /* (0070,0009) SQ (Sequence with explicit length #=1)      # 180, 1 GraphicObjectSequence
            (fffe,e000) na (Item with explicit length #=6)          # 172, 1 Item
              (0070,0005) CS [PIXEL]                                  #   6, 1 GraphicAnnotationUnits
              (0070,0020) US 2                                        #   2, 1 GraphicDimensions
              (0070,0021) US 13                                       #   2, 1 NumberOfGraphicPoints
              (0070,0022) FL 92\113\99\129\110\142\132\146\153\135\164\118\172\98\162\75\132\67... # 104,26 GraphicData
              (0070,0023) CS [POLYLINE]                               #   8, 1 GraphicType
              (0070,0024) CS [N]                                      #   2, 1 GraphicFilled
            (fffe,e00d) na (ItemDelimitationItem for re-encoding)   #   0, 0 ItemDelimitationItem
          (fffe,e0dd) na (SequenceDelimitationItem for re-encod.) #   0, 0 SequenceDelimitationItem
          */

      gdcm::SequenceOfItems::SizeType nitems = sqi->GetNumberOfItems();
      // fprintf(stdout, "found %lu items\n", nitems);
      for (int itemNr = 1; itemNr <= nitems; itemNr++) { // we should create our polys at this level....

        Polygon poly;
        poly.StudyInstanceUID = boost::algorithm::trim_copy(StudyInstanceUID);
        // the string above might contain a utf-8 version of a null character
        if (poly.StudyInstanceUID.back() == '\0')
          poly.StudyInstanceUID.replace(poly.StudyInstanceUID.end() - 1, poly.StudyInstanceUID.end(), "");
        poly.SeriesInstanceUID = boost::algorithm::trim_copy(SeriesInstanceUID);
        // the string above might contain a utf-8 version of a null character
        if (poly.SeriesInstanceUID.back() == '\0')
          poly.SeriesInstanceUID.replace(poly.SeriesInstanceUID.end() - 1, poly.SeriesInstanceUID.end(), "");
        poly.Filename = filename;

        gdcm::Item &item = sqi->GetItem(itemNr);
        gdcm::DataSet &subds = item.GetNestedDataSet();

        // lookup the ReferencedImageSequence
        /*    (0008,1140) SQ (Sequence with explicit length #=1)      # 114, 1 ReferencedImageSequence
                (fffe,e000) na (Item with explicit length #=2)          # 106, 1 Item
                  (0008,1150) UI =MRImageStorage                          #  26, 1 ReferencedSOPClassUID
                  (0008,1155) UI [1.3.6.1.4.1.45037.091a6b29babb817c5ffcc7199d0c0de4cf9d52c155360] #  64, 1 ReferencedSOPInstanceUID
                (fffe,e00d) na (ItemDelimitationItem for re-encoding)   #   0, 0 ItemDelimitationItem
              (fffe,e0dd) na (SequenceDelimitationItem for re-encod.) #   0, 0 SequenceDelimitationItem
        */
        if (!subds.FindDataElement(referencedImageSequence)) {
          // fprintf(stdout, " %d does not have a referenced image sequence\n", itemNr);
          continue;
        } else {
          // fprintf(stdout, " %d found a ReferencedImageSequence\n", itemNr);
        }
        const gdcm::DataElement &de3 = subds.GetDataElement(referencedImageSequence);
        gdcm::SmartPointer<gdcm::SequenceOfItems> sqiReferencedImageSequence = de3.GetValueAsSQ();
        gdcm::SequenceOfItems::SizeType nitems3 = sqiReferencedImageSequence->GetNumberOfItems();
        // fprintf(stdout, " referenced image sequence with %lu item(-s)\n", nitems3);
        for (int itemNr2 = 1; itemNr2 <= nitems3; itemNr2++) {
          gdcm::Item &item3 = sqiReferencedImageSequence->GetItem(itemNr2);
          gdcm::DataSet &subds3 = item3.GetNestedDataSet();
          if (!subds3.FindDataElement(referencedSOPInstanceUID)) {
            // fprintf(stdout, " %d no referencedSOPInstanceUID\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deReferencedSOPInstanceUID = subds3.GetDataElement(referencedSOPInstanceUID);
          const gdcm::ByteValue *bv = deReferencedSOPInstanceUID.GetByteValue();
          std::string refUID(bv->GetPointer(), bv->GetLength());
          // fprintf(stdout, " found: %s as ReferencedSOPInstanceUID\n", refUID.c_str());
          poly.ReferencedSOPInstanceUID = boost::algorithm::trim_copy(refUID);
          // the string above might contain a utf-8 version of a null character
          if (poly.ReferencedSOPInstanceUID.back() == '\0')
            poly.ReferencedSOPInstanceUID.replace(poly.ReferencedSOPInstanceUID.end() - 1, poly.ReferencedSOPInstanceUID.end(), "");
        }

        // this is for items in 0070,0001, now we need to look for 0070,0009
        if (!subds.FindDataElement(graphicObjectSequence)) {
          // fprintf(stdout, " %d does not have GraphicObjectSequence\n", itemNr);
          continue;
        } else {
          // fprintf(stdout, " %d found a GraphicObjectSequence\n", itemNr);
        }
        const gdcm::DataElement &de2 = subds.GetDataElement(graphicObjectSequence);
        gdcm::SmartPointer<gdcm::SequenceOfItems> sqiGraphicObjectSequence = de2.GetValueAsSQ();
        gdcm::SequenceOfItems::SizeType nitems2 = sqiGraphicObjectSequence->GetNumberOfItems();
        // fprintf(stdout, " graphic object sequence with %lu item(-s)\n", nitems2);
        for (int itemNr2 = 1; itemNr2 <= nitems2; itemNr2++) {
          gdcm::Item &item2 = sqiGraphicObjectSequence->GetItem(itemNr2);
          gdcm::DataSet &subds2 = item2.GetNestedDataSet();
          if (!subds2.FindDataElement(graphicType)) {
            // fprintf(stdout, " %d no graphic type\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deGraphicType = subds2.GetDataElement(graphicType);
          const gdcm::ByteValue *bv = deGraphicType.GetByteValue();
          std::string gT(bv->GetPointer(), bv->GetLength());

          if (!subds2.FindDataElement(numberOfGraphicPoints)) {
            // fprintf(stdout, " %d no number of graphic points\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deNumberOfGraphicPoints = subds2.GetDataElement(numberOfGraphicPoints);
          // std::string nGP = gdcm::DirectoryHelper::GetStringValueFromTag(numberOfGraphicPoints, subds2);
          const gdcm::ByteValue *bv2 = deNumberOfGraphicPoints.GetByteValue();
          std::string nGP(bv2->GetPointer(), bv2->GetLength());
          int numberOfPoints = *(nGP.c_str());

          if (!subds2.FindDataElement(graphicData)) {
            // fprintf(stdout, " %d no graphic data\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deGraphicData = subds2.GetDataElement(graphicData);
          const gdcm::ByteValue *bv3 = deGraphicData.GetByteValue();
          std::string gD(bv3->GetPointer(), bv3->GetLength());

          // we can read the values now based on the value representation (VR::FL is single float)
          gdcm::Element<gdcm::VR::FL, gdcm::VM::VM1_n> elwc;
          gdcm::VR vr = gdcm::VR::FL;
          unsigned int vrsize = vr.GetSizeof();
          unsigned int count = gD.size() / vrsize;
          elwc.SetLength(count * vrsize);
          std::stringstream ss1;
          ss1.str(gD);
          elwc.Read(ss1);
          poly.coords.resize(elwc.GetLength());
          for (unsigned int i = 0; i < elwc.GetLength(); ++i) {
            poly.coords[i] = elwc.GetValue(i);
          }
          // now store the poly
          storage->push_back(poly);
        }
      }

    } else {
      fprintf(stdout, "not sequence?\n");
    }
  }

  return true;
}

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

  // to read the PR images we need another ReaderType, they are not images
  // lets call a function that will give us back the polygon information we need to process
  // the image data
  std::vector<Polygon> storage;
  std::map<std::string, std::string> SOPInstanceUID2SeriesInstanceUID;
  parseForPolygons(input, &storage, &SOPInstanceUID2SeriesInstanceUID);

  auto first_element = SOPInstanceUID2SeriesInstanceUID.begin();
  auto pos = SOPInstanceUID2SeriesInstanceUID.find(first_element->first);
  // we would like to add to storage the ReferencedSeriesInstanceUID's, (we only get the ReferencedSOPInstanceUID above)
  int goodStorage = 0;
  for (int i = 0; i < storage.size(); i++) {
    bool found = false;
    for (auto pos = SOPInstanceUID2SeriesInstanceUID.begin(); pos != SOPInstanceUID2SeriesInstanceUID.end(); pos++) {
      if (pos->first == storage[i].ReferencedSOPInstanceUID) {
        storage[i].ReferencedSeriesInstanceUID = pos->second;
        found = true;
        goodStorage++;
        break;
      }
    }
    if (!found) {
      fprintf(stderr, "Error: We have never seen an MR image with the SOPInstanceUID: \"%s\". It is referenced in \"%s\"\n",
              storage[i].ReferencedSOPInstanceUID.c_str(), storage[i].Filename.c_str());
    }
  }
  fprintf(stdout, "We could identify the referenced series in %d referenced polylines.\n", goodStorage);

  // loop over storage and append to resultJSON
  resultJSON["POLYLINES"] = json::array();
  for (int i = 0; i < storage.size(); i++) {
    auto entry = json::object();
    entry["ReferencedSOPInstanceUID"] = storage[i].ReferencedSOPInstanceUID;
    entry["Coordinates"] = storage[i].coords;
    entry["Filename"] = storage[i].Filename;
    entry["StudyInstanceUID"] = storage[i].StudyInstanceUID;
    entry["SeriesInstanceUID"] = storage[i].SeriesInstanceUID;
    entry["ReferencedSeriesInstanceUID"] = storage[i].ReferencedSeriesInstanceUID;
    resultJSON["POLYLINES"].push_back(entry);
  }

  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::Image<unsigned char, 3> MaskImageType;
  typedef itk::Image<unsigned char, 2> MaskSliceImageType;

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  dicomIO->LoadPrivateTagsOn();

  reader->SetImageIO(dicomIO);

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(false); // we want to use the keys as SeriesInstanceUIDs
  nameGenerator->AddSeriesRestriction("0008|0060");
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

    seriesItr = runThese.begin();
    seriesEnd = runThese.end();
    while (seriesItr != seriesEnd) {
      seriesIdentifier = seriesItr->c_str();
      ++seriesItr;
      // do we have this seriesIdentifier in storage? What polygons are part of that?
      bool doSomething = false;
      for (int i = 0; i < storage.size(); i++) {
        if (storage[i].ReferencedSeriesInstanceUID == seriesIdentifier)
          doSomething = true;
      }
      if (!doSomething)
        continue;

      std::cout << "Processing series: " << std::endl;
      std::cout << "  " << seriesIdentifier << std::endl;

      // std::string outputSeries = output + "/" + seriesIdentifier;
      // if (!force && itksys::SystemTools::FileIsDirectory(outputSeries.c_str())) {
      //   fprintf(stdout, "Skip this series %s, output directory exists already...\n", outputSeries.c_str());
      //   exit(0); // this is no skip, that is giving up...
      // }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      // we should check now if any of these files SOPInstanceUID appears in our array of polylines (storage)
      // read the images one by one
      if (1) {
        // loop over all files in this series
        for (int sliceNr = 0; sliceNr < fileNames.size(); sliceNr++) {
          std::vector<std::string> oneFile;
          oneFile.push_back(fileNames[sliceNr]);
          reader->SetFileNames(oneFile);
          try {
            reader->Update();
          } catch (itk::ExceptionObject &ex) {
            std::cout << ex << std::endl;
            return EXIT_FAILURE;
          }
          ImageType::Pointer inputImage = reader->GetOutput();
          typedef itk::MetaDataDictionary DictionaryType;
          DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();

          // example: ./Modules/Filtering/ImageIntensity/test/itkPolylineMask2DImageFilterTest.cxx
        }
      }

      // BELOW: we read the series as a volume - would make it possible to work with this as nii for example
      // Here we should restrict ourselves to series that are mentioned in the ReferenedSeriesInstanceUIDs in storage
      if (1) {
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
      }
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
