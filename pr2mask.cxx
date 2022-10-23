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

#include "itkPolyLineParametricPath.h"
#include "itkPolylineMask2DImageFilter.h"

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
#include <codecvt>
#include <locale> // wstring_convert
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
  std::string UnformattedTextValue;
  std::string Filename; // the name of the DICOM file
};

bool parseForPolygons(std::string input, std::vector<Polygon> *storage, std::map<std::string, std::string> *SOPInstanceUID2SeriesInstanceUID, bool verbose) {
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
      if (verbose) {
        std::cerr << "Could not read  \"" << filename << "\" as DICOM, ignore." << std::endl;
      }
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
    const gdcm::Tag textObjectSequence(0x0070, 0x0008);       // TextObjectSequence
    const gdcm::Tag unformattedTextValue(0x0070, 0x0006);     // unformattedTextValue
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
    if (verbose) {
      fprintf(stdout, " found SeriesInstanceUID: %s\n", SeriesInstanceUID.c_str());
    }
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
          if (verbose)
            fprintf(stdout, " %d does not have a referenced image sequence\n", itemNr);
          continue;
        } else {
          if (verbose)
            fprintf(stdout, " %d found a ReferencedImageSequence\n", itemNr);
        }
        const gdcm::DataElement &de3 = subds.GetDataElement(referencedImageSequence);
        gdcm::SmartPointer<gdcm::SequenceOfItems> sqiReferencedImageSequence = de3.GetValueAsSQ();
        gdcm::SequenceOfItems::SizeType nitems3 = sqiReferencedImageSequence->GetNumberOfItems();
        // fprintf(stdout, " referenced image sequence with %lu item(-s)\n", nitems3);
        for (int itemNr2 = 1; itemNr2 <= nitems3; itemNr2++) {
          gdcm::Item &item3 = sqiReferencedImageSequence->GetItem(itemNr2);
          gdcm::DataSet &subds3 = item3.GetNestedDataSet();
          if (!subds3.FindDataElement(referencedSOPInstanceUID)) {
            if (verbose)
              fprintf(stdout, " %d no referencedSOPInstanceUID\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deReferencedSOPInstanceUID = subds3.GetDataElement(referencedSOPInstanceUID);
          const gdcm::ByteValue *bv = deReferencedSOPInstanceUID.GetByteValue();
          std::string refUID(bv->GetPointer(), bv->GetLength());
          if (verbose)
            fprintf(stdout, " found: %s as ReferencedSOPInstanceUID\n", refUID.c_str());
          poly.ReferencedSOPInstanceUID = boost::algorithm::trim_copy(refUID);
          // the string above might contain a utf-8 version of a null character
          if (poly.ReferencedSOPInstanceUID.back() == '\0')
            poly.ReferencedSOPInstanceUID.replace(poly.ReferencedSOPInstanceUID.end() - 1, poly.ReferencedSOPInstanceUID.end(), "");
        }

        //
        // this is for items in 0070,0008 textObjectSequence
        //
        if (subds.FindDataElement(textObjectSequence)) { // optional
          /*
              (0070,0008) SQ (Sequence with explicit length #=1)      # 150, 1 TextObjectSequence
                (fffe,e000) na (Item with explicit length #=4)          # 142, 1 Item
                  (0070,0004) CS [PIXEL]                                  #   6, 1 AnchorPointAnnotationUnits
                  (0070,0006) ST [Min/Max: 7 / 349
          Mean: 132, Deviation: 103
          Total: 657?860
          Pixel... #  94, 1 UnformattedTextValue
                  (0070,0014) FL 190.5\108.5                              #   8, 2 AnchorPoint
                  (0070,0015) CS [Y]                                      #   2, 1 AnchorPointVisibility
                (fffe,e00d) na (ItemDelimitationItem for re-encoding)   #   0, 0 ItemDelimitationItem
              (fffe,e0dd) na (SequenceDelimitationItem for re-encod.) #   0, 0 SequenceDelimitationItem
          */
          const gdcm::DataElement &de2 = subds.GetDataElement(textObjectSequence);
          gdcm::SmartPointer<gdcm::SequenceOfItems> sqiTextObjectSequence = de2.GetValueAsSQ();
          gdcm::SequenceOfItems::SizeType nitems2 = sqiTextObjectSequence->GetNumberOfItems();
          if (verbose)
            fprintf(stdout, " graphic object sequence with %lu item(-s)\n", nitems2);
          for (int itemNr2 = 1; itemNr2 <= nitems2; itemNr2++) {
            gdcm::Item &item2 = sqiTextObjectSequence->GetItem(itemNr2);
            gdcm::DataSet &subds2 = item2.GetNestedDataSet();
            if (!subds2.FindDataElement(unformattedTextValue)) {
              if (verbose)
                fprintf(stdout, " %d no unformatted text value\n", itemNr2);
              continue;
            }
            const gdcm::DataElement &deUnformattedTextValue = subds2.GetDataElement(unformattedTextValue);
            const gdcm::ByteValue *bv = deUnformattedTextValue.GetByteValue();
            std::string uTV(bv->GetPointer(), bv->GetLength());
            poly.UnformattedTextValue = boost::algorithm::trim_copy(uTV);
          }
        }

        //
        // this is for items in 0070,0001, now we need to look for 0070,0009
        //
        if (!subds.FindDataElement(graphicObjectSequence)) {
          // fprintf(stdout, " %d does not have GraphicObjectSequence\n", itemNr);
          continue;
        } else {
          // fprintf(stdout, " %d found a GraphicObjectSequence\n", itemNr);
        }
        const gdcm::DataElement &de2 = subds.GetDataElement(graphicObjectSequence);
        gdcm::SmartPointer<gdcm::SequenceOfItems> sqiGraphicObjectSequence = de2.GetValueAsSQ();
        gdcm::SequenceOfItems::SizeType nitems2 = sqiGraphicObjectSequence->GetNumberOfItems();
        if (verbose)
          fprintf(stdout, " graphic object sequence with %lu item(-s)\n", nitems2);
        for (int itemNr2 = 1; itemNr2 <= nitems2; itemNr2++) {
          gdcm::Item &item2 = sqiGraphicObjectSequence->GetItem(itemNr2);
          gdcm::DataSet &subds2 = item2.GetNestedDataSet();
          if (!subds2.FindDataElement(graphicType)) {
            if (verbose)
              fprintf(stdout, " %d no graphic type\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deGraphicType = subds2.GetDataElement(graphicType);
          const gdcm::ByteValue *bv = deGraphicType.GetByteValue();
          std::string gT(bv->GetPointer(), bv->GetLength());

          if (!subds2.FindDataElement(numberOfGraphicPoints)) {
            if (verbose)
              fprintf(stdout, " %d no number of graphic points\n", itemNr);
            continue;
          }
          const gdcm::DataElement &deNumberOfGraphicPoints = subds2.GetDataElement(numberOfGraphicPoints);
          // std::string nGP = gdcm::DirectoryHelper::GetStringValueFromTag(numberOfGraphicPoints, subds2);
          const gdcm::ByteValue *bv2 = deNumberOfGraphicPoints.GetByteValue();
          std::string nGP(bv2->GetPointer(), bv2->GetLength());
          int numberOfPoints = *(nGP.c_str());

          if (!subds2.FindDataElement(graphicData)) {
            if (verbose)
              fprintf(stdout, " %d no graphic data\n", itemNr);
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
      if (verbose)
        fprintf(stdout, "not sequence?\n");
    }
  }

  return true;
}

using ImageType2D = itk::Image<PixelType, 2>;

ImageType2D::Pointer createMaskFromStorage(ImageType2D::Pointer im2change, std::vector<int> polyIds, std::vector<Polygon> storage) {
  // Use a copy of im2change and return the image of all the masks applied to that slice.
  ImageType2D::Pointer mask = ImageType2D::New();
  ImageType2D::RegionType maskRegionInput = im2change->GetLargestPossibleRegion();
  mask->SetRegions(maskRegionInput);
  mask->Allocate();
  mask->FillBuffer(itk::NumericTraits<PixelType>::Zero);
  mask->SetOrigin(im2change->GetOrigin());
  mask->SetSpacing(im2change->GetSpacing());
  mask->SetDirection(im2change->GetDirection());

  ImageType2D::Pointer lmask = ImageType2D::New();
  lmask->SetRegions(maskRegionInput);
  lmask->Allocate();
  lmask->SetOrigin(im2change->GetOrigin());
  lmask->SetSpacing(im2change->GetSpacing());
  lmask->SetDirection(im2change->GetDirection());

  for (int p = 0; p < polyIds.size(); p++) {
    lmask->FillBuffer(itk::NumericTraits<PixelType>::One);

    using InputPolylineType = itk::PolyLineParametricPath<2>;
    InputPolylineType::Pointer inputPolyline = InputPolylineType::New();
    using InputFilterType = itk::PolylineMask2DImageFilter<ImageType2D, InputPolylineType, ImageType2D>;
    InputFilterType::Pointer filter = InputFilterType::New();

    int storageIdx = polyIds[p]; // we just use the first one

    // COPY THE NEW IMAGE from polygon data
    using VertexType = InputPolylineType::VertexType;

    // Add vertices to the polyline
    double spacingx = lmask->GetSpacing()[0];
    double spacingy = lmask->GetSpacing()[1];
    double originx = lmask->GetOrigin()[0];
    double originy = lmask->GetOrigin()[1];

    for (int j = 0; j < storage[storageIdx].coords.size(); j += 2) {
      VertexType v0;
      // coordinates are in pixel, we need coordinates based on the bounding box
      v0[0] = originx + spacingx * storage[storageIdx].coords[j];
      v0[1] = originy + spacingy * storage[storageIdx].coords[j + 1];
      // fprintf(stdout, "POLYLINE with %lu elements %f\n", storage[storageIdx].coords.size(), originx + spacingx * storage[storageIdx].coords[j]);
      inputPolyline->AddVertex(v0);
    }

    // Connect the input image
    filter->SetInput1(lmask);

    // Connect the Polyline
    filter->SetInput2(inputPolyline);
    try {
      filter->Update();
    } catch (itk::ExceptionObject &err) {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
    }
    ImageType2D::Pointer lres = filter->GetOutput();

    // now copy the filter output to mask
    ImageType2D::RegionType maskRegion = mask->GetLargestPossibleRegion();
    ImageType2D::RegionType lmaskRegion = lres->GetLargestPossibleRegion();
    itk::ImageRegionIterator<ImageType2D> maskIterator(mask, maskRegion);
    itk::ImageRegionIterator<ImageType2D> lmaskIterator(lres, lmaskRegion);
    while (!maskIterator.IsAtEnd() && !lmaskIterator.IsAtEnd()) {
      if (lmaskIterator.Get() > 0) {
        maskIterator.Set(1);
      }
      ++maskIterator;
      ++lmaskIterator;
    }
  }
  return mask;
}

bool invalidChar(char c) { return !isprint(static_cast<unsigned char>(c)); }
void stripUnicode(std::string &str) { str.erase(remove_if(str.begin(), str.end(), invalidChar), str.end()); }

int main(int argc, char *argv[]) {
  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetVersion("0.0.1");
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("PR2MASK: Convert presentation state files with polygons to label fields in DICOM format.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for images/ and labels/ folder as DICOM.", MetaCommand::STRING, true);

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.SetOptionLongTag("SeriesName", "seriesname");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  bool seriesIdentifierFlag = false;
  std::string input = command.GetValueAsString("indir");
  std::string output = command.GetValueAsString("outdir");

  if (input.size() == 0 || output.size() == 0) {
    return 1;
  }

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  if (command.GetOptionWasSet("SeriesName"))
    seriesIdentifierFlag = true;

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
  parseForPolygons(input, &storage, &SOPInstanceUID2SeriesInstanceUID, verbose);

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
    if (!found && verbose) {
      fprintf(stderr, "Warning: Unknown MR (SOPInstanceUID): \"%s\" referenced in \"%s\"\n", storage[i].ReferencedSOPInstanceUID.c_str(),
              storage[i].Filename.c_str());
    }
  }
  if (verbose)
    fprintf(stdout, "We could identify the referenced series for %d/%lu polylines.\n", goodStorage, storage.size());

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

    std::string a = storage[i].UnformattedTextValue;
    stripUnicode(a);
    entry["UnformattedTextValue"] = a; // these might contain non UTF-8 characters, so we remove anything unicode
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
    /*while (seriesItr != seriesEnd) {
      std::cout << "Found DICOM Series: ";
      std::cout << std::endl;
      std::cout << "  " << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }*/

    std::string seriesIdentifier;

    SeriesIdContainer runThese;
    if (seriesIdentifierFlag) { // If no optional series identifier
      runThese.push_back(seriesName);
    } else {
      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
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
      for (int i = 0; i < storage.size(); i++) { // not just the first key, we need all keys
        if (storage[i].ReferencedSeriesInstanceUID == seriesIdentifier) {
          doSomething = true;
        }
      }
      if (!doSomething)
        continue;

      if (verbose) {
        std::cout << "Processing series: " << std::endl;
        std::cout << "  " << seriesIdentifier << std::endl;
      }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      // we should check now if any of these files SOPInstanceUID appears in our array of polylines (storage)
      // read the images one by one
      if (1) {
        gdcm::UIDGenerator uid;
        uid.SetRoot("1.3.6.1.4.1.45037");
        const char *newSeriesInstanceUID = uid.Generate();
        // fprintf(stdout, "NEW SERIESINSTANCEUID: \"%s\"\n", newSeriesInstanceUID);
        //  loop over all files in this series
        for (int sliceNr = 0; sliceNr < fileNames.size(); sliceNr++) {
          // using ImageType2D = itk::Image<PixelType, 2>;
          typedef itk::ImageFileReader<ImageType2D> Reader2DType;
          typedef itk::ImageFileWriter<ImageType2D> Writer2DType;
          Reader2DType::Pointer r = Reader2DType::New();
          Writer2DType::Pointer w = Writer2DType::New();
          typedef itk::GDCMImageIO ImageIOType;
          ImageIOType::Pointer dicomIO = ImageIOType::New();
          dicomIO->LoadPrivateTagsOn();
          dicomIO->KeepOriginalUIDOn();
          // we need to find out what for this image the ReferencedSOPInstanceUID is
          // only draw the contour on that image

          r->SetImageIO(dicomIO);
          r->SetFileName(fileNames[sliceNr]);
          try {
            r->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }
          // now changed the slice we are importing
          ImageType2D::Pointer im2change = r->GetOutput();

          // make a copy of this image series in the output/images/ folder
          if (1) {
            w->SetInput(im2change);
            // we should have a folder for each image series
            boost::filesystem::path p(fileNames[sliceNr]);
            boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "images" + boost::filesystem::path::preferred_separator +
                                            seriesIdentifier + boost::filesystem::path::preferred_separator + p.filename().c_str() + ".dcm";
            if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
              // create the output directory
              create_directories(p_out.parent_path());
            }
            w->SetFileName(p_out.c_str());
            w->SetImageIO(dicomIO);

            // fprintf(stdout, "write: %s with %s %s\n", p_out.c_str(), newSOPInstanceUID.c_str(), newSeriesInstanceUID);

            try {
              w->Update();
            } catch (itk::ExceptionObject &err) {
              std::cerr << "ExceptionObject caught !" << std::endl;
              std::cerr << err << std::endl;
              return EXIT_FAILURE;
            }
          }

          // get some meta-data from the opened file (find out what polygons are relevant)
          typedef itk::MetaDataDictionary DictionaryType;
          DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();
          std::string SeriesInstanceUID;
          std::string SOPInstanceUID;
          std::string seriesNumber;
          std::string seriesDescription;
          itk::ExposeMetaData<std::string>(dictionary, "0020|000E", SeriesInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0008|0018", SOPInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0011", seriesNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0008|103E", seriesDescription);
          // lookup the correct polygon, we need a loop over multiple polygons we also need
          // outer and inner rings to represent holes. Not something we can get from PACS :-(
          std::vector<int> polyIds;
          for (int i = 0; i < storage.size(); i++) {
            if (storage[i].ReferencedSOPInstanceUID == SOPInstanceUID) {
              // for this slice we should treat this polygon
              polyIds.push_back(i);
            }
          }
          // using InputPolylineType = itk::PolyLineParametricPath<2>;
          // InputPolylineType::Pointer inputPolyline = InputPolylineType::New();
          // using InputFilterType = itk::PolylineMask2DImageFilter<ImageType2D, InputPolylineType, ImageType2D>;
          // InputFilterType::Pointer filter = InputFilterType::New();
          ImageType2D::Pointer maskFromPolys;
          if (polyIds.size() > 0) {
            maskFromPolys = createMaskFromStorage(im2change, polyIds, storage);

            /*
            int storageIdx = polyIds[0]; // we just use the first one

            // COPY THE NEW IMAGE from polygon data
            using VertexType = InputPolylineType::VertexType;

            // Add vertices to the polyline
                        double spacingx = im2change->GetSpacing()[0];
                        double spacingy = im2change->GetSpacing()[1];
                        double originx = im2change->GetOrigin()[0];
                        double originy = im2change->GetOrigin()[1];

                        for (int j = 0; j < storage[storageIdx].coords.size(); j += 2) {
                          VertexType v0;
                          // coordinates are in pixel, we need coordinates based on the bounding box
                          v0[0] = originx + spacingx * storage[storageIdx].coords[j];
                          v0[1] = originy + spacingy * storage[storageIdx].coords[j + 1];
                          inputPolyline->AddVertex(v0);
                        }

                        // fill the input image with 1 (label)
                        im2change->FillBuffer(itk::NumericTraits<PixelType>::One);

                        // Connect the input image
                        filter->SetInput1(im2change);

                        // Connect the Polyline
                        filter->SetInput2(inputPolyline);
                        try {
                          filter->Update();
                        } catch (itk::ExceptionObject &err) {
                          std::cerr << "ExceptionObject caught !" << std::endl;
                          std::cerr << err << std::endl;
                          return EXIT_FAILURE;
                        } */
          } else { // without polygon just return an empty image
            // fill the input image with 1 (label)
            im2change->FillBuffer(itk::NumericTraits<PixelType>::Zero);
          }

          // InputImageType::Pointer inputImage = reader->GetOutput();
          ImageType2D::RegionType region;
          region = im2change->GetBufferedRegion();
          ImageType2D::SizeType size = region.GetSize();
          // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

          ImageType2D::PixelContainer *container;
          container = im2change->GetPixelContainer();
          container->SetContainerManageMemory(false);
          unsigned int bla = sizeof(ImageType2D::PixelType);
          ImageType2D::PixelType *buffer2 = container->GetBufferPointer();

          ImageType2D::Pointer nImage;
          if (polyIds.size() > 0) {
            // nImage = filter->GetOutput();

            // ImageType2D::Pointer nImage = filter->GetOutput();
            ImageType2D::PixelContainer *container2;
            container2 = maskFromPolys->GetPixelContainer();
            ImageType2D::PixelType *buffer3 = container2->GetBufferPointer();

            // Here we copy all values over, that is 0, 1, 2, 3 but also additional labels
            // that have been selected before (air in intestines for example).
            memcpy(buffer2, &(buffer3[0]), size[0] * size[1] * bla);
            // We can clean the data (remove all other label).
            /*for (int k = 0; k < size[0] * size[1]; k++) {
              if (buffer2[k] > 3) {
                buffer2[k] = 0; // set to background
              }
            }*/
          }

          // dilate and erode the im2change (mask) (would be better if we do this in 3D)
          using StructuringElementType = itk::BinaryBallStructuringElement<PixelType, 2>;
          using ErodeFilterType = itk::BinaryErodeImageFilter<ImageType2D, ImageType2D, StructuringElementType>;
          using DilateFilterType = itk::BinaryDilateImageFilter<ImageType2D, ImageType2D, StructuringElementType>;
          DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
          ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

          StructuringElementType structuringElement;
          structuringElement.SetRadius(1); // 3x3 structuring element
          structuringElement.CreateStructuringElement();
          binaryDilate->SetKernel(structuringElement);
          binaryErode->SetKernel(structuringElement);
          binaryErode->SetForegroundValue(1);
          binaryErode->SetBackgroundValue(0);
          binaryDilate->SetForegroundValue(1);
          binaryDilate->SetBackgroundValue(0);
          binaryDilate->SetInput(im2change);
          binaryErode->SetInput(binaryDilate->GetOutput());
          binaryDilate->SetDilateValue(1);
          binaryErode->SetErodeValue(1);

          try {
            binaryErode->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }

          // copy the values back to the im2change buffer
          ImageType2D::Pointer cleanMask = binaryErode->GetOutput();
          ImageType2D::PixelContainer *container3;
          container3 = cleanMask->GetPixelContainer();
          ImageType2D::PixelType *buffer4 = container3->GetBufferPointer();
          memcpy(buffer2, &(buffer4[0]), size[0] * size[1] * bla);

          // now change something to make a new copy of that file
          int newSeriesNumber = 1000 + atoi(seriesNumber.c_str()) + 1;
          gdcm::UIDGenerator uid;
          uid.SetRoot("1.3.6.1.4.1.45037");
          std::string newSOPInstanceUID(uid.Generate());

          dicomIO->KeepOriginalUIDOn();
          itk::MetaDataDictionary &dictionarySlice = r->GetOutput()->GetMetaDataDictionary();
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0011", std::to_string(newSeriesNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0018", newSOPInstanceUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000E", std::string(newSeriesInstanceUID));

          // set the series description (max 64 characters)
          std::string newSeriesDescription = seriesDescription + " (mask)";
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000E", newSeriesDescription.substr(0, 64));

          w->SetInput(im2change);
          // create the output filename
          // we should have a folder for each image series
          boost::filesystem::path p(fileNames[sliceNr]);
          boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "labels" + boost::filesystem::path::preferred_separator +
                                          newSeriesInstanceUID + boost::filesystem::path::preferred_separator + p.filename().c_str() + ".dcm";
          if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
            create_directories(p_out.parent_path());
          }
          w->SetFileName(p_out.c_str());
          w->SetImageIO(dicomIO);

          if (verbose) {
            fprintf(stdout, "write: %s with %s %s\n", p_out.c_str(), newSOPInstanceUID.c_str(), newSeriesInstanceUID);
          }

          try {
            w->Update();
          } catch (itk::ExceptionObject &err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
          }
          // example: ./Modules/Filtering/ImageIntensity/test/itkPolylineMask2DImageFilterTest.cxx
        }
        // we should remember this mapping in a csv file
        boost::filesystem::path csv_out = output + boost::filesystem::path::preferred_separator + "data.csv";
        if (!itksys::SystemTools::FileIsDirectory(csv_out.parent_path().c_str())) {
          create_directories(csv_out.parent_path());
        }
        if (!boost::filesystem::exists(csv_out)) {
          FILE *fp = fopen(csv_out.c_str(), "a");
          fprintf(fp, "ImageSeriesInstanceUID,LabelSeriesInstanceUID\n");
          fclose(fp);
        }
        FILE *fp = fopen(csv_out.c_str(), "a");
        fprintf(fp, "\"images/%s\",\"labels/%s\"\n", seriesIdentifier.c_str(), newSeriesInstanceUID);
        fclose(fp);

        std::string res = resultJSON.dump(4) + "\n";
        // save the json information to a file as well, use folder names
        boost::filesystem::path json_out = output + boost::filesystem::path::preferred_separator + seriesIdentifier + "_" + newSeriesInstanceUID + ".json";
        std::ofstream out(json_out.c_str());
        out << res;
        out.close();
      }

    } // loop over series
  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  if (verbose) {
    std::string res = resultJSON.dump(4) + "\n";
    fprintf(stdout, "%s", res.c_str());
  }

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
