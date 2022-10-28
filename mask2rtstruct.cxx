#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"

#include "itkImageAdaptor.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

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
#include "itkContourExtractor2DImageFilter.h"
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

const float FLOAT_EPSILON = 0.0001;
const unsigned int precision = 16;
const double contourValue = 100.0;

using json = nlohmann::json;
using namespace boost::filesystem;

json resultJSON;

// We need to identify for each series if they are a presentation state object and if we can extract some
// contours from them.
struct Polygon {
  std::vector<float> coords; // the vector of pixel coordinates extracted from presentation state
  std::string StudyInstanceUID;
  std::string SeriesInstanceUID;
  std::string SOPInstanceUID;
  std::string GeometricType;
  int Intensity;        // the original image's intensity for which this is a contour
  std::string Filename; // the name of the DICOM file
};

using ImageType2D = itk::Image<PixelType, 2>;
typedef itk::ContourExtractor2DImageFilter<ImageType2D> ContourExtractorType;
typedef ContourExtractorType::VertexType VertexType;

using ImageTypeRT = itk::Image<PixelType, 2>;

std::vector<Polygon> storeContours(ContourExtractorType::Pointer contourExtractFilter) {

  std::vector<Polygon> ret;
  // can we get the offset and spacing for the data from the filter?
  const ImageType2D *img = contourExtractFilter->GetInput();
  // ImageType2D::RegionType inputRegion = img->GetLargestPossibleRegion();
  // ImageType2D::SizeType size = inputRegion.GetSize();
  const ImageType2D::SpacingType &spacing = img->GetSpacing();
  double originx = img->GetOrigin()[0];
  double originy = img->GetOrigin()[1];

  unsigned int numOutputs = contourExtractFilter->GetNumberOfOutputs();

  unsigned int numVertices;
  VertexType firstVertex;
  VertexType lastVertex;

  for (unsigned int i = 0; i < numOutputs; i++) {
    ContourExtractorType::VertexListConstPointer vertices = contourExtractFilter->GetOutput(i)->GetVertexList();
    // does the orientation of the polygons matter here?
    // for RTStruct its either keyhole technique (bah) or XOR. For XOR we need the correct order or only large area followed by all the holes
    // What about islands?

    numVertices = vertices->Size();

    firstVertex = vertices->ElementAt(0);
    lastVertex = vertices->ElementAt(numVertices - 1);

    Polygon poly;

    if ((fabs(firstVertex[0] - lastVertex[0]) < FLOAT_EPSILON) && (fabs(firstVertex[1] - lastVertex[1]) < FLOAT_EPSILON)) {
      // It's a closed contour.
      // So, the last vertex won't be written as it is same as the 1st vertex.
      // WriteCommonData(currentSlice, numVertices - 1, CLOSED_PLANAR, file1);
      poly.GeometricType = "CLOSED_PLANAR"; // should this be CLOSEDPLANAR_XOR ?

      for (unsigned int j = 0; j < (numVertices - 2); j++) {
        //   WriteVertexCoordinates(vertices->ElementAt(j), offset_index, spacing, zValue, file1);
        VertexType v = vertices->ElementAt(j);
        poly.coords.push_back((v[0] - originx) * spacing[0]);
        poly.coords.push_back((v[1] - originy) * spacing[1]);
      }
    } else {
      // The contour is a open planar one.
      // WriteCommonData(currentSlice, numVertices, OPEN_PLANAR, file1);
      poly.GeometricType = "OPEN_PLANAR";

      for (unsigned int j = 0; j < (numVertices - 1); j++) {
        //   WriteVertexCoordinates(vertices->ElementAt(j), offset_index, spacing, zValue, file1);
        VertexType v = vertices->ElementAt(j);
        poly.coords.push_back((v[0] - originx) * spacing[0]);
        poly.coords.push_back((v[1] - originy) * spacing[1]);
      }
    }
    ret.push_back(poly);
  }

  return ret;
}

gdcm::DataElement CreateFakeElement(gdcm::Tag const &tag, bool toremove) {
  static const gdcm::Global &g = gdcm::Global::GetInstance();
  static const gdcm::Dicts &dicts = g.GetDicts();
  static const gdcm::Dict &pubdict = dicts.GetPublicDict();
  static size_t countglobal = 0;
  static std::vector<gdcm::Tag> balcptags = gdcm::Anonymizer::GetBasicApplicationLevelConfidentialityProfileAttributes();
  size_t count = countglobal % balcptags.size();

  const gdcm::DictEntry &dictentry = pubdict.GetDictEntry(tag);

  gdcm::DataElement de;
  de.SetTag(tag);
  using gdcm::VR;
  const VR &vr = dictentry.GetVR();
  // if( vr != VR::INVALID )
  if (vr.IsDual()) {
    if (vr == VR::US_SS) {
      de.SetVR(VR::US);
    } else if (vr == VR::US_SS_OW) {
      de.SetVR(VR::OW);
    } else if (vr == VR::OB_OW) {
      de.SetVR(VR::OB);
    }
  } else {
    de.SetVR(vr);
  }
  const char str[] = "BasicApplicationLevelConfidentialityProfileAttributes";
  const char safe[] = "This is safe to keep";
  if (de.GetVR() != VR::SQ) {
    if (toremove)
      de.SetByteValue(str, (uint32_t)strlen(str));
    else
      de.SetByteValue(safe, (uint32_t)strlen(safe));
  } else {
    // Create an item
    gdcm::Item it;
    it.SetVLToUndefined();
    gdcm::DataSet &nds = it.GetNestedDataSet();
    // Insert sequence into data set
    assert(de.GetVR() == gdcm::VR::SQ);
    gdcm::SmartPointer<gdcm::SequenceOfItems> sq = new gdcm::SequenceOfItems();
    sq->SetLengthToUndefined();
    de.SetValue(*sq);
    de.SetVLToUndefined();
    // ds.Insert(de);

    if (!toremove) {
      nds.Insert(CreateFakeElement(balcptags[count], true));
      countglobal++;
    } else {
      gdcm::Attribute<0x0008, 0x0000> at1 = {0}; // This element has no reason to be 'anonymized'...
      nds.Insert(at1.GetAsDataElement());
      gdcm::Attribute<0x000a, 0x0000> at2 = {0};
      nds.Insert(at2.GetAsDataElement());
    }
    sq->AddItem(it);
  }
  return de;
}

void getRT(std::vector<Polygon> storage, ImageType2D::Pointer im2change, gdcm::File *f) {
  // ImageTypeRT::Pointer rt_object = ImageTypeRT::New();

  gdcm::DataSet &ds = f->GetDataSet();

  using gdcm::Tag;
  using gdcm::VR;

  std::vector<gdcm::Tag> balcptags = gdcm::Anonymizer::GetBasicApplicationLevelConfidentialityProfileAttributes();

  // gdcm::Writer w;
  //  gdcm::File &f = w.GetFile();
  // gdcm::DataSet &ds = f->GetDataSet();

  // Add attribute that need to be anonymized:
  std::vector<gdcm::Tag>::const_iterator it = balcptags.begin();
  for (; it != balcptags.end(); ++it) {
    ds.Insert(CreateFakeElement(*it, true));
  }

  // Add attribute that do NOT need to be anonymized:
  static const gdcm::Global &g = gdcm::Global::GetInstance();
  static const gdcm::Dicts &dicts = g.GetDicts();
  static const gdcm::Dict &pubdict = dicts.GetPublicDict();

  using gdcm::Dict;
  Dict::ConstIterator dictit = pubdict.Begin();
  for (; dictit != pubdict.End(); ++dictit) {
    const gdcm::Tag &dicttag = dictit->first;
    if (dicttag == Tag(0x6e65, 0x6146))
      break;
    // const gdcm::DictEntry &dictentry = dictit->second;
    ds.Insert(CreateFakeElement(dicttag, false));
  }
  ds.Remove(gdcm::Tag(0x400, 0x500));
  ds.Remove(gdcm::Tag(0x12, 0x62));
  ds.Remove(gdcm::Tag(0x12, 0x63));

  // Make sure to override any UID stuff
  gdcm::UIDGenerator uid;
  gdcm::DataElement de(Tag(0x8, 0x18)); // SOP Instance UID
  de.SetVR(VR::UI);
  const char *u = uid.Generate();
  de.SetByteValue(u, (uint32_t)strlen(u));
  // ds.Insert( de );
  ds.Replace(de);

  de.SetTag(Tag(0x8, 0x16)); // SOP Class UID
  de.SetVR(VR::UI);
  gdcm::MediaStorage ms(gdcm::MediaStorage::RTStructureSetStorage);
  de.SetByteValue(ms.GetString(), (uint32_t)strlen(ms.GetString()));
  ds.Replace(de); // replace !
                  /*
                    // gdcm::FileMetaInformation &fmi = f->GetHeader();
                    //  fmi.SetDataSetTransferSyntax( gdcm::TransferSyntax::ImplicitVRLittleEndian );
                    // fmi.SetDataSetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
                
                    // gdcm::UIDGenerator uid;
                    gdcm::DataElement de2(gdcm::Tag(0x8, 0x18)); // SOP Instance UID
                    de2.SetVR(gdcm::VR::UI);
                    // const char *u = uid.Generate();
                    de2.SetByteValue(u, (uint32_t)strlen(u));
                    // ds.Insert( de );
                    ds.Replace(de2);
                
                    // add an element
                    gdcm::Attribute<0x0008, 0x0000> at1 = {0};
                    ds.Insert(at1.GetAsDataElement());
                    de = at1.GetAsDataElement();
                
                    // add a sequence
                    gdcm::Item item;
                    item.SetVLToUndefined();
                    gdcm::DataSet &nds = item.GetNestedDataSet();
                    // Insert sequence into data set
                    // assert(de.GetVR() == gdcm::VR::SQ);
                    gdcm::SmartPointer<gdcm::SequenceOfItems> sq = new gdcm::SequenceOfItems();
                    sq->SetLengthToUndefined();
                    de.SetValue(*sq);
                    de.SetVLToUndefined();
                    sq->AddItem(item);
                
                    ds.Insert(de);
                
                    de.SetTag(gdcm::Tag(0x8, 0x16)); // SOP Class UID
                    de.SetVR(gdcm::VR::UI);
                    gdcm::MediaStorage ms2(gdcm::MediaStorage::RTStructureSetStorage);
                    de.SetByteValue(ms2.GetString(), (uint32_t)strlen(ms2.GetString()));
                    ds.Replace(de); // replace !
                  */
  // set the  meta info header
  gdcm::FileMetaInformation &fmi = f->GetHeader();
  fmi.SetDataSetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);

  // return rt_object;
}

std::vector<int> getIntensities(ImageType2D::Pointer im2change) {
  std::vector<int> intensities;
  std::map<int, bool> memory;

  ImageType2D::RegionType maskRegion = im2change->GetLargestPossibleRegion();
  itk::ImageRegionIterator<ImageType2D> maskIterator(im2change, maskRegion);
  while (!maskIterator.IsAtEnd()) {
    memory.insert(std::pair<int, bool>(maskIterator.Get(), true));
    ++maskIterator;
  }
  for (std::map<int, bool>::iterator iter = memory.begin(); iter != memory.end(); ++iter) {
    intensities.push_back(iter->first);
  }

  return intensities;
}

static bool endsWith(std::string_view str, std::string_view suffix) {
  return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}
bool invalidChar(char c) { return !isprint(static_cast<unsigned char>(c)); }
void stripUnicode(std::string &str) { str.erase(remove_if(str.begin(), str.end(), invalidChar), str.end()); }

int main(int argc, char *argv[]) {
  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetVersion("0.0.1");
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("MASK2RTSTRUCT: Convert a mask series to polygons and finally an RT-Struct DICOM file.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM.", MetaCommand::STRING, true);

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

  // loop over storage and append to resultJSON
  resultJSON["POLYLINES"] = json::array();

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

      if (verbose) {
        std::cout << "Processing series: " << std::endl;
        std::cout << "  " << seriesIdentifier << std::endl;
      }

      typedef std::vector<std::string> FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      if (1) {
        gdcm::UIDGenerator uid;
        uid.SetRoot("1.3.6.1.4.1.45037");
        const char *newSeriesInstanceUID = uid.Generate();
        //  loop over all files in this series
        for (int sliceNr = 0; sliceNr < fileNames.size(); sliceNr++) {
          // using ImageType2D = itk::Image<PixelType, 2>;
          typedef itk::ImageFileReader<ImageType2D> Reader2DType;
          Reader2DType::Pointer r = Reader2DType::New();
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

          //
          // extract the contours for this slice
          //
          // what are the intensities we need to extract contours for?
          std::vector<int> intensities = getIntensities(im2change);
          std::vector<Polygon> polys;

          // for each intensity we need to extract the pixel into a mask and run the contour extraction
          // but what about holes? How are they represented in the RTStruct objects?
          for (int iidx = 0; iidx < intensities.size(); iidx++) {
            // create a new slice object and copy the intensities over
            int intensity = intensities[iidx];
            if (intensity == 0)
              continue; // ignore background material

            ImageType2D::Pointer labelField = ImageType2D::New();
            ImageType2D::RegionType im2changeRegion = im2change->GetLargestPossibleRegion();
            labelField->SetRegions(im2changeRegion);
            labelField->Allocate();
            labelField->SetOrigin(im2change->GetOrigin());
            labelField->SetSpacing(im2change->GetSpacing());
            labelField->SetDirection(im2change->GetDirection());
            ImageType2D::RegionType labelFieldRegion = labelField->GetLargestPossibleRegion();
            itk::ImageRegionIterator<ImageType2D> im2changeFieldIterator(im2change, im2changeRegion);
            itk::ImageRegionIterator<ImageType2D> labelFieldIterator(labelField, labelFieldRegion);
            while (!labelFieldIterator.IsAtEnd() && !im2changeFieldIterator.IsAtEnd()) {
              if (im2changeFieldIterator.Get() == intensity) {
                labelFieldIterator.Set(1);
              }
              ++im2changeFieldIterator;
              ++labelFieldIterator;
            }

            ContourExtractorType::Pointer contourExtractFilter = ContourExtractorType::New();
            contourExtractFilter->SetContourValue(0.5);
            contourExtractFilter->ReverseContourOrientationOn();
            contourExtractFilter->SetInput(labelField);
            try {
              contourExtractFilter->Update();
            } catch (itk::ExceptionObject &err) {
              std::cerr << "ExceptionObject caught!" << std::endl;
              std::cerr << err << std::endl;
              return EXIT_FAILURE;
            }
            std::vector<Polygon> tmpPolys = storeContours(contourExtractFilter);

            // add them to polys in the correct order for XOR
            for (int i = 0; i < tmpPolys.size(); i++) {
              // we would need to add the first clock-wise polygon and afterwards add the counter-clockwise in order
              // this is very strange... islands need to be correct...
              tmpPolys[i].Intensity = intensity;
              polys.push_back(tmpPolys[i]);
            }
          }

          // get some meta-data from the opened file (find out what polygons are relevant)
          typedef itk::MetaDataDictionary DictionaryType;
          DictionaryType &dictionary = dicomIO->GetMetaDataDictionary();
          std::string SeriesInstanceUID;
          std::string SOPInstanceUID;
          std::string seriesNumber;
          std::string seriesDescription;
          std::string StudyInstanceUID;
          itk::ExposeMetaData<std::string>(dictionary, "0020|000E", SeriesInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0020|000D", StudyInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0008|0018", SOPInstanceUID);
          itk::ExposeMetaData<std::string>(dictionary, "0020|0011", seriesNumber);
          itk::ExposeMetaData<std::string>(dictionary, "0008|103E", seriesDescription);

          for (int i = 0; i < polys.size(); i++) { // add the missing information and put into storage
            polys[i].Filename = fileNames[sliceNr];
            polys[i].StudyInstanceUID = StudyInstanceUID;
            polys[i].SeriesInstanceUID = SeriesInstanceUID;
            polys[i].SOPInstanceUID = SOPInstanceUID;
            storage.push_back(polys[i]);
          }

          // remember all polygons in the output json
          for (int i = 0; i < storage.size(); i++) {
            auto entry = json::object();
            entry["Coordinates"] = storage[i].coords;
            entry["Filename"] = storage[i].Filename;
            entry["StudyInstanceUID"] = storage[i].StudyInstanceUID;
            entry["SeriesInstanceUID"] = storage[i].SeriesInstanceUID;
            entry["SOPInstanceUID"] = storage[i].SOPInstanceUID;
            entry["GeometricType"] = storage[i].GeometricType;
            entry["Intensity"] = storage[i].Intensity;
            resultJSON["POLYLINES"].push_back(entry);
          }

          // now change something to make a new copy of that file
          int newSeriesNumber = atoi(seriesNumber.c_str()) + 1;
          gdcm::UIDGenerator uid;
          uid.SetRoot("1.3.6.1.4.1.45037");
          std::string newSOPInstanceUID(uid.Generate());

          dicomIO->KeepOriginalUIDOn();
          itk::MetaDataDictionary &dictionarySlice = r->GetOutput()->GetMetaDataDictionary();
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|0011", std::to_string(newSeriesNumber));
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|0018", newSOPInstanceUID);
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0020|000E", std::string(newSeriesInstanceUID));

          // set the series description (max 64 characters)
          std::string newSeriesDescription = seriesDescription + " (RT-STRUCT)";
          itk::EncapsulateMetaData<std::string>(dictionarySlice, "0008|103E", newSeriesDescription.substr(0, 64));

          if (verbose) {
            fprintf(stdout, "process mask: %s (input SeriesInstanceUID: %s)\n", newSOPInstanceUID.c_str(), newSeriesInstanceUID);
          }
          // create an RTStruct object for this slice and save it
          boost::filesystem::path p(fileNames[sliceNr]);
          // assume we have an extension already
          std::string fname = p.filename().c_str();
          if (endsWith(fname.c_str(), ".dcm")) {
            fname.replace(fname.end() - 4, fname.end(), "");
          }
          boost::filesystem::path p_out = output + boost::filesystem::path::preferred_separator + "RTStructs" + boost::filesystem::path::preferred_separator +
                                          newSeriesInstanceUID + boost::filesystem::path::preferred_separator + fname.c_str() + ".dcm";
          if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
            create_directories(p_out.parent_path());
          }

          gdcm::Writer w;
          gdcm::File &f = w.GetFile();
          if (1) {
          }

          // getRT(storage, im2change, &file);
          if (verbose)
            fprintf(stdout, " Writing file: %s\n", p_out.c_str());
          w.SetCheckFileMetaInformation(true);
          w.SetFileName(p_out.c_str());
          if (!w.Write()) {
            std::cerr << "Error writing file!" << std::endl;
            return EXIT_FAILURE;
          }
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

        // compute computational time
        boost::posix_time::ptime timeLocalEnd = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_period tp(timeLocal, timeLocalEnd);
        resultJSON["wall_time"] = boost::posix_time::to_simple_string(timeLocalEnd - timeLocal);
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
