#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkMetaDataObject.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkLinearInterpolateImageFunction.h"

#include "gdcmPhotometricInterpretation.h"
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


#include "itkRGBPixel.h"
#include "itkImageAdaptor.h"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>

#include <codecvt>

#include <ft2build.h>
#include FT_FREETYPE_H

std::string versionString;

typedef signed short PixelType;

// input image and mask
using ImageType3D = itk::Image<PixelType, 3>;
using MaskImageType3D = itk::Image<PixelType, 3>;

// output fused image
using CPixelType = itk::RGBPixel<unsigned char>;
using CImageType = itk::Image<CPixelType, 3>;

using ImageType2D = itk::Image<PixelType, 2>;
using CImageType2D = itk::Image<CPixelType, 2>;

// accumulate the overlay text into a separate kbuffer at this size
unsigned int image_buffer_gen_size[2] = {0,0};
unsigned char **image_buffer_gen = NULL;

typedef struct {
  std::string title;
  std::string subtitle;
  std::string info;
  float brightnessContrastLL;
  float brightnessContrastUL;
  bool  votemapMode;
  float votemapMax; // user provided possible max value for votemap
  float votemapAgree; // blue to yellow border
  float obtained_max_votemap_value; // calculate during fusing
} OverlayInfos_t;

bool verbose = false;
bool stableUIDs = true;

std::vector<std::vector<float>> labelColors2 = {{0, 0, 0}, {166,206,227}, {31,120,180}, {178,223,138}, {51,160,44}, {251,154,153}, {227,26,28}, {253,191,111}, {255,127,0}, {202,178,214}, {106,61,154}, {255,255,153}, {177,89,40}};
// only three colors for background, low agreement and high agreement
std::vector<std::vector<float>> labelColorsVotemap = {{0, 0, 0}, {255,255,0}, {0,0,255} };

// read mask image series
ImageType3D::Pointer readImageSerie(std::string dirName) {
  using NamesGeneratorType = itk::GDCMSeriesFileNames;
  auto nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetGlobalWarningDisplay(false);
  nameGenerator->SetDirectory(dirName);

  try {
    using SeriesIdContainer = std::vector<std::string>;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    auto                      seriesItr = seriesUID.begin();
    auto                      seriesEnd = seriesUID.end();

    if (seriesItr != seriesEnd) {
      std::cout << "The directory: ";
      std::cout << dirName << std::endl;
      std::cout << "Contains the following DICOM Series: ";
      std::cout << std::endl;
    } else {
      std::cout << "No DICOMs in: " << dirName << std::endl;
      return nullptr;
    }

    while (seriesItr != seriesEnd) {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }

    seriesItr = seriesUID.begin();
    while (seriesItr != seriesUID.end()) {
      std::string seriesIdentifier;
      //if (argc > 3) {
      //  seriesIdentifier = argv[3];
      //  seriesItr = seriesUID.end();
     // }
      //else // otherwise convert everything
      //{
        seriesIdentifier = seriesItr->c_str();
        seriesItr++;
      //}
      std::cout << "\nReading: ";
      std::cout << seriesIdentifier << std::endl;
      using FileNamesContainer = std::vector<std::string>;
      FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      using ReaderType = itk::ImageSeriesReader<ImageType3D>;
      auto reader = ReaderType::New();
      using ImageIOType = itk::GDCMImageIO;
      auto dicomIO = ImageIOType::New();
      reader->MetaDataDictionaryArrayUpdateOn();
      //reader->LoadPrivateTagsOn();
      reader->SetImageIO(dicomIO);
      reader->SetFileNames(fileNames);
      reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt
      try {
        reader->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return nullptr;
      }
      return reader->GetOutput();
    }
  } catch (const itk::ExceptionObject & ex) {
    std::cout << ex << std::endl;
    return nullptr;
  }
  return nullptr; 
}

// read mask image series
MaskImageType3D::Pointer readMaskImageSerie(std::string dirName) {
    using NamesGeneratorType = itk::GDCMSeriesFileNames;
  auto nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetGlobalWarningDisplay(false);
  nameGenerator->SetDirectory(dirName);

  try {
    using SeriesIdContainer = std::vector<std::string>;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    auto                      seriesItr = seriesUID.begin();
    auto                      seriesEnd = seriesUID.end();

    if (seriesItr != seriesEnd) {
      std::cout << "The directory: ";
      std::cout << dirName << std::endl;
      std::cout << "Contains the following DICOM Series: ";
      std::cout << std::endl;
    } else {
      std::cout << "No DICOMs in: " << dirName << std::endl;
      return nullptr;
    }

    while (seriesItr != seriesEnd) {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
    }

    seriesItr = seriesUID.begin();
    while (seriesItr != seriesUID.end()) {
      std::string seriesIdentifier;
      //if (argc > 3) {
      //  seriesIdentifier = argv[3];
      //  seriesItr = seriesUID.end();
     // }
      //else // otherwise convert everything
      //{
        seriesIdentifier = seriesItr->c_str();
        seriesItr++;
      //}
      std::cout << "\nReading: ";
      std::cout << seriesIdentifier << std::endl;
      using FileNamesContainer = std::vector<std::string>;
      FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);

      using ReaderType = itk::ImageSeriesReader<MaskImageType3D>;
      auto reader = ReaderType::New();
      using ImageIOType = itk::GDCMImageIO;
      auto dicomIO = ImageIOType::New();
      reader->SetImageIO(dicomIO);
      reader->SetFileNames(fileNames);
      reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt
      try {
        reader->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return nullptr;
      }
      return reader->GetOutput();
    }
  } catch (const itk::ExceptionObject & ex) {
    std::cout << ex << std::endl;
    return nullptr;
  }
  return nullptr; 
}

// add a colored bar over the top of every image
void addBar(CImageType::Pointer img, int height) {
  CImageType::RegionType kregion = img->GetLargestPossibleRegion();
  using ImageSizeType = typename CImageType::SizeType;
  ImageSizeType regionSize;
  regionSize[0] = kregion.GetSize()[0];
  regionSize[1] = height;
  regionSize[2] = kregion.GetSize()[2];

  kregion.SetSize(regionSize);
  itk::ImageRegionIteratorWithIndex<CImageType> kIterator(img, kregion);
  kIterator.GoToBegin();
  while(!kIterator.IsAtEnd()) {
    CPixelType value = kIterator.Value(); 
    value.SetRed((int)(255));
    value.SetGreen((int)(77));
    value.SetBlue((int)(5));
    kIterator.Set(value);
    ++kIterator;
  }

  // a single black line below the red bar
  CImageType::RegionType::IndexType regionStart;
  regionStart = kregion.GetIndex();
  regionStart[1] = height;
  regionSize[0] = kregion.GetSize()[0];
  regionSize[1] = 1;
  regionSize[2] = kregion.GetSize()[2];

  kregion.SetSize(regionSize);
  kregion.SetIndex(regionStart);
  itk::ImageRegionIteratorWithIndex<CImageType> kIterator2(img, kregion);
  kIterator2.GoToBegin();
  while(!kIterator2.IsAtEnd()) {
    CPixelType value = kIterator.Value();

    value.SetRed(0); //(int)(255));
    value.SetGreen(0); //(int)(255));
    value.SetBlue(0); //(int)(255));
    //value.SetGreen((int)(77));
    //value.SetBlue((int)(5));
    kIterator2.Set(value);
    ++kIterator2;
  }
}

void draw_bitmap_gen(FT_Bitmap *bitmap, int width, int height, FT_Int x, FT_Int y) {
  FT_Int i, j, p, q;
  FT_Int x_max = x + bitmap->width;
  FT_Int y_max = y + bitmap->rows;

  for (i = x, p = 0; i < x_max; i++, p++) {
    for (j = y, q = 0; j < y_max; j++, q++) {
      if (i < 0 || j < 0 || i >= width || j >= height)
        continue;

      image_buffer_gen[j][i] |= bitmap->buffer[q * bitmap->width + p];
    }
  }
}

// return one character at position num in wc2
int get_mb(wchar_t *wc2, const char* ptr, int num) {
    std::mbtowc(nullptr, 0, 0); // reset the conversion state
    const char* end = ptr + std::strlen(ptr);
    int ret{};
    int n = 0;
    for (wchar_t wc; (ret = std::mbtowc(&wc, ptr, end - ptr)) > 0; ptr += ret) {
        if (n++ == num) {
           *wc2 = wc;
           break; 
        }
    }
    return ret;
}

void addToReportGen(char *buffer, std::string font_file, int font_size, std::string sstext, int posx, int posy, float radiants) {
  FT_Library library;


  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>, wchar_t> convert;
  std::wstring stext = convert.from_bytes(sstext);

  //bool verbose = 1;
  if (verbose) {
    fprintf(stdout, "  addToReportGen: \"%s\"\n", sstext.c_str());
  }

  double angle;
  int target_height;
  int n, num_chars;

  // int font_length = 20;

  std::string font_path = font_file;
  // int font_size = 42;
  int face_index = 0;

  FT_Face face;
  FT_GlyphSlot slot;
  FT_Matrix matrix; /* transformation matrix */
  FT_Vector pen;    /* untransformed origin  */
  FT_Error error;

  error = FT_Init_FreeType(&library); /* initialize library */

  if (error != 0) {
    fprintf(stderr, "\033[0;31mError\033[0m: The freetype library could not be initialized with this font.\n");
    return;
  }

  int start_px = 40;
  int start_py = 20;
  int text_lines = 1;

  float repeat_spacing = 1.0f;
  int xmax = image_buffer_gen_size[0];
  int ymax = image_buffer_gen_size[1];
  num_chars = stext.size();

  int px = posx;
  // WIDTH - ((num_chars + 2) * font_size - start_px);
  // int py = start_py + (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));
  int py = posy;
  // start_py + 2.0 * (font_size);

  // our image is not of that size but much larger (repeated mosaic tiles in y)
  // reset image_buffer_gen now
  for (int i = 0; i < image_buffer_gen_size[1]; i++) {
    if (image_buffer_gen[i])
      memset(image_buffer_gen[i], 0, sizeof(unsigned char) * image_buffer_gen_size[0]);
  }

  angle = radiants;
  target_height = image_buffer_gen_size[1]; // 512; // shouldn't this be the height of the buffer, e.g. image_buffer_gen_size[1]?

  error = FT_New_Face(library, font_file.c_str(), face_index, &face); /* create face object */

  if (face == NULL) {
    fprintf(stderr, "\033[0;31mError\033[0m: no face found, provide the filename of a ttf file...\n");
    return;
  }

  int font_size_in_pixel = font_size;
  error = FT_Set_Char_Size(face, font_size_in_pixel * 64, 0, 150, 150); // font_size_in_pixel * 64, 0, 96, 0); /* set character size */
  /* error handling omitted */
  if (error != 0) {
    fprintf(stderr, "\033[0;31mError\033[0;31m: FT_Set_Char_Size returned error, could not set size %d.\n", font_size_in_pixel);
    return;
  }

  slot = face->glyph;

  /* set up matrix */
  matrix.xx = (FT_Fixed)(cos(angle) * 0x10000L);
  matrix.xy = (FT_Fixed)(-sin(angle) * 0x10000L);
  matrix.yx = (FT_Fixed)(sin(angle) * 0x10000L);
  matrix.yy = (FT_Fixed)(cos(angle) * 0x10000L);

  /* the pen position in 26.6 cartesian space coordinates; */
  /* start at (300,200) relative to the upper left corner  */
  //pen.x = (num_chars * 1 * 64);
  pen.x = (1 * 64);
  pen.y = (target_height - 40) * 64; // the 60 here is related to the font size!

  int nn = 0; 
  for (std::wstring::iterator it = stext.begin(); it != stext.end(); it++) {
    //wchar_t c = *it;
    wchar_t c;
    int ret = get_mb(&c, sstext.c_str(), nn);
    nn++;
    if (c == '\n') {
      continue; // ignore newlines
    }

    FT_Set_Transform(face, &matrix, &pen);

    FT_UInt glyph_index = FT_Get_Char_Index( face, *it );
    error = FT_Load_Glyph(face, glyph_index, FT_LOAD_RENDER);
    if (error) {
      fprintf(stdout, "\033[0;31mError\033[0m:: [addToReportGen] could not load character: '%ls'\n", &c);
      continue;
    }

    draw_bitmap_gen(&slot->bitmap, image_buffer_gen_size[0], image_buffer_gen_size[1], slot->bitmap_left, target_height - slot->bitmap_top);

    pen.x += slot->advance.x;
    pen.y += slot->advance.y;
  } 

  FT_Done_Face(face);

  float current_image_min_value = 0.0f;
  float current_image_max_value = 255.0f;
  unsigned char *bvals = (unsigned char *)buffer;
  for (int yi = 0; yi < image_buffer_gen_size[1]; yi++) {
    for (int xi = 0; xi < image_buffer_gen_size[0]; xi++) {
      if (image_buffer_gen[yi][xi] == 0)
        continue;
      // I would like to copy the value from image over to
      // the buffer. At some good location...

      int newx = px + xi;
      int newy = py + yi;
      int idx = newy * xmax + newx;
      if (newx < 0 || newx >= xmax || newy < 0 || newy >= ymax)
        continue;
      //if (image_buffer828[yi][xi] == 0)
      //  continue;

      // instead of blending we need to use a fixed overlay color
      // we have image information between current_image_min_value and current_image_max_value
      // we need to scale the image_buffer by those values.
      float f = 1.0;
      float v = 1.0f * image_buffer_gen[yi][xi] / 255.0; // 0 to 1 for color, could be inverted if we have a white background
      float w = f * 1.0f * bvals[idx] / current_image_max_value;
      float alpha_blend = (v + w * (1.0f - v));

      // fprintf(stdout, "%d %d: %d\n", xi, yi, bvals[idx]);
      bvals[idx] = (unsigned char)std::max(
          0.0f, std::min(current_image_max_value, current_image_min_value + (alpha_blend) * (current_image_max_value - current_image_min_value)));
    }
  }
}


// use the reference image folder as a source of DICOM tag values for the fusedImage and save a new RGB image series to outputDir
int saveFusedImageSeries(CImageType::Pointer fusedImage, std::string outputDir, std::string referenceDir, OverlayInfos_t *overlayInfo) {
  // maybe better use https://docs.itk.org/projects/doxygen/en/stable/Examples_2IO_2DicomSeriesReadSeriesWrite_8cxx-example.html

  // we would need to mark this buffer_fusedImage as AI with red bar and text
  std::string font_file = "Menlo.ttf";
  if (const char *env_p = std::getenv("REPORT_FONT_PATH")) {
    font_file = std::string(env_p);
  }
  if (!boost::filesystem::exists(font_file)) {
    fprintf(stderr, "\033[0;31mError\033[0m: no font provided. Set the environment variable REPORT_FONT_PATH to the location of a ttf file.\n");
    fflush(stderr);
    return EXIT_FAILURE;
  }

  CImageType::RegionType kregion = fusedImage->GetLargestPossibleRegion();
  int KWIDTH = kregion.GetSize()[0];
  int KHEIGHT = kregion.GetSize()[1];
  // create a buffer for the text before fusing it ontop of fusedImage
  char *kbuffer = new char[KWIDTH * KHEIGHT];
  memset(&kbuffer[0], 0, sizeof(char)*KWIDTH*KHEIGHT);
  image_buffer_gen_size[0] = KWIDTH;
  image_buffer_gen_size[1] = KHEIGHT;
  image_buffer_gen = (unsigned char**) malloc(KHEIGHT * sizeof(unsigned char*));
  if (image_buffer_gen) { 
    for (int i = 0; i < KHEIGHT; i++) {
      image_buffer_gen[i] = (unsigned char*)malloc(KWIDTH * sizeof(unsigned char));
      if (!image_buffer_gen[i]) {
        fprintf(stderr, "Error allocating memory for fused overlay image\n");
        fflush(stderr);
        exit(-1); // give up
      }
    }
  }

  int kw = kregion.GetSize()[0];
  int barHeight = std::max<int>(26, 0.07 * kw);
  int fontSize = std::max<int>(4, 0.012 * kw);
  addBar(fusedImage, barHeight);
  if (overlayInfo->votemapMode) {
    // in vote map mode we have a possible placeholder for the max votemap value relative to the possible value
    float maxObtainedVotemapValue = 100.0f * overlayInfo->obtained_max_votemap_value / overlayInfo->votemapMax;
    char buff_str[256];
    sprintf(buff_str, "%.0f%%", maxObtainedVotemapValue);
    std::string maxObtainedVotemapValueStr = std::string(buff_str);
    // replace the placeholder in all info strings
    std::string placeholder = "{peak_agreement}";
    size_t pos = overlayInfo->info.find(placeholder);
    while (pos != std::string::npos) {
      overlayInfo->info.replace(pos, placeholder.length(), maxObtainedVotemapValueStr);
      pos = overlayInfo->info.find(placeholder, pos + 1);
    }
    pos = overlayInfo->subtitle.find(placeholder);
    while (pos != std::string::npos) {
      overlayInfo->subtitle.replace(pos, placeholder.length(), maxObtainedVotemapValueStr);
      pos = overlayInfo->subtitle.find(placeholder, pos + 1);
    }
    pos = overlayInfo->title.find(placeholder);
    while (pos != std::string::npos) {
      overlayInfo->title.replace(pos, placeholder.length(), maxObtainedVotemapValueStr);
      pos = overlayInfo->title.find(placeholder, pos + 1);
    }
  }

  addToReportGen(kbuffer, font_file, (fontSize+2)<7?7:fontSize+2, overlayInfo->title, 8, -22, 0);
  addToReportGen(kbuffer, font_file, (fontSize)<7?7:fontSize, overlayInfo->subtitle, 8, 20, 0);
  addToReportGen(kbuffer, font_file, (fontSize)<7?7:fontSize, overlayInfo->info, 8, KHEIGHT-50, 0);

  using ImageSizeType = typename CImageType::SizeType;
  ImageSizeType regionSize;
  int numSlices = kregion.GetSize()[2];
  // now write the kbuffer into the fusedImage (into every single slice) kbuffer is 2D
  for (int sliceNr = 0; sliceNr < numSlices; sliceNr++) {
    kregion = fusedImage->GetLargestPossibleRegion();
    regionSize[0] = kregion.GetSize()[0];
    regionSize[1] = kregion.GetSize()[1];
    regionSize[2] = 1;
    kregion.SetSize(regionSize);
    CImageType::RegionType::IndexType regionStart;
    regionStart = kregion.GetIndex();
    regionStart[2] = sliceNr;
    kregion.SetIndex(regionStart);

    itk::ImageRegionIteratorWithIndex<CImageType> kIterator(fusedImage, kregion);
    kIterator.GoToBegin();
    while(!kIterator.IsAtEnd()) {
      CPixelType value = kIterator.Value(); 
      CImageType::IndexType kidx = kIterator.GetIndex(); 

      // the grayscale input image (fused image)
      int c1 = (int)(value.GetRed());
      int c2 = (int)(value.GetGreen());
      int c3 = (int)(value.GetBlue());

      // the grayscale text image where 0 is transparent and 255 is fully colored text to be fused with fusedImage
      int ttt = (unsigned char)kbuffer[kidx[1]*KWIDTH + kidx[0]]; // 0 to 255

      float alpha_a = (float)(ttt/255.0);
      float alpha_b = 0.5;
      float alpha = alpha_a + alpha_b*(1.0 - alpha_a); // alpha of text is ttt/255, alpha of b is 0
      float C_a_r = 1.0;
      float C_a_g = 237.0/255.0;
      float C_a_b = 160.0/255.0;
      float C_b_r = c1/255.0;
      float C_b_g = c2/255.0;
      float C_b_b = c3/255.0;
      float Cr = (C_a_r * alpha_a + (C_b_r * alpha_b *(1.0 - alpha_a)))/alpha;
      float Cg = (C_a_g * alpha_a + (C_b_g * alpha_b *(1.0 - alpha_a)))/alpha;
      float Cb = (C_a_b * alpha_a + (C_b_b * alpha_b *(1.0 - alpha_a)))/alpha;
      Cr = std::min<float>(1, std::max<float>(0,Cr));
      Cg = std::min<float>(1, std::max<float>(0,Cg));
      Cb = std::min<float>(1, std::max<float>(0,Cb));

      value.SetRed((int)(255.0 * Cr));
      value.SetGreen((int)(255.0 * Cg));
      value.SetBlue((int)(255.0 * Cb));
      kIterator.Set(value);

      ++kIterator;
    }
  }
  free(kbuffer);
  // end putting text into fusedImage

  using ReaderType = itk::ImageSeriesReader<ImageType3D>;

  using ImageIOType = itk::GDCMImageIO;
  using NamesGeneratorType = itk::GDCMSeriesFileNames;
 
  auto gdcmIO = ImageIOType::New();
  auto namesGenerator = NamesGeneratorType::New();

  namesGenerator->SetInputDirectory(referenceDir.c_str());
  
  const ReaderType::FileNamesContainer & filenames = namesGenerator->GetInputFileNames();
    
  size_t numberOfFileNames = filenames.size();

  // this is 3D, we copy into 2D slice
  CImageType::PixelContainer *container_fusedImage = fusedImage->GetPixelContainer();
  CImageType::PixelType *buffer_fusedImage = container_fusedImage->GetBufferPointer();

  std::string newFusedSeriesInstanceUID(""); // done only once
  if (!stableUIDs) {
    gdcm::UIDGenerator uid;
    uid.SetRoot("1.3.6.1.4.1.45037");
    newFusedSeriesInstanceUID = std::string(uid.Generate());
  }

  // the number of files and the last dimension in fusedImage have to match
  if (numberOfFileNames != fusedImage->GetLargestPossibleRegion().GetSize()[2]) {
    std::cerr << "Error: number of files in reference directory (" << numberOfFileNames << ") does not match number of slices in fused image ("
              << fusedImage->GetLargestPossibleRegion().GetSize()[2] << ")" << std::endl;
    return EXIT_FAILURE;
  }

  // for each slice in the list of filenames we need to read the 2D original image, create the output and copy the tag values over that we need
  for (int sliceNr = 0; sliceNr < filenames.size(); sliceNr++) {
    typedef itk::ImageFileReader<ImageType2D> Reader2DType;
    typedef itk::ImageFileWriter<CImageType2D> Writer2DType;
    Reader2DType::Pointer r = Reader2DType::New();

    typedef itk::GDCMImageIO ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    dicomIO->LoadPrivateTagsOn();
    dicomIO->KeepOriginalUIDOn();

    r->SetImageIO(dicomIO);
    r->SetFileName(filenames[sliceNr]);
    try {
      r->Update();
    } catch (itk::ExceptionObject &err) {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
    // now changed the slice we are importing
    ImageType2D::Pointer im2change = r->GetOutput();

    // get some meta-data from the opened file (find out what polygons are relevant)
    typedef itk::MetaDataDictionary DictionaryType;
    itk::MetaDataDictionary &dictionarySlice = im2change->GetMetaDataDictionary();

    std::string frameOfReferenceUID("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0052", frameOfReferenceUID);

    std::string StudyInstanceUID("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|000d", StudyInstanceUID);

    std::string newFusedSOPInstanceUID("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0018", newFusedSOPInstanceUID);
    if (stableUIDs) {
      std::string derivedFusedSOPInstanceUID = newFusedSOPInstanceUID;
      std::string endString = ".4";
      if (derivedFusedSOPInstanceUID.substr(derivedFusedSOPInstanceUID.size() - 2, 2) == ".4")
        endString = ".5";

      // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
      derivedFusedSOPInstanceUID = derivedFusedSOPInstanceUID.substr(0, 64 - 3) + endString;
      newFusedSOPInstanceUID = derivedFusedSOPInstanceUID;
    } else {
      gdcm::UIDGenerator uid;
      uid.SetRoot("1.3.6.1.4.1.45037");
      newFusedSOPInstanceUID = std::string(uid.Generate());
    }

    std::string seriesNumberStr("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0011", seriesNumberStr);

    std::string patientID("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0010|0020", patientID);

    std::string patientName("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0010|0010", patientName);

    std::string studyID("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|000d", studyID);

    std::string accession_number("");
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0050", accession_number);

    std::string imagePositionPatient; // we might be able to get them this way, but can we set them?
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0032", imagePositionPatient);
    // perhaps we have to use the parsed values to write them again further down?
    double origin3D[3];
    sscanf(imagePositionPatient.c_str(), "%lf\\%lf\\%lf", &(origin3D[0]), &(origin3D[1]), &(origin3D[2]));

    std::string imageOrientation;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0037", imageOrientation); // image orientation patient
    double imageOrientationField[6];
    sscanf(imageOrientation.c_str(), "%lf\\%lf\\%lf\\%lf\\%lf\\%lf", &(imageOrientationField[0]), &(imageOrientationField[1]), &(imageOrientationField[2]),
          &(imageOrientationField[3]), &(imageOrientationField[4]), &(imageOrientationField[5]));

    std::string sliceThicknessString;
    double sliceThickness = 0.0;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0018|0050", sliceThicknessString);
    sscanf(sliceThicknessString.c_str(), "%lf", &sliceThickness);

    std::string sliceLocationString;
    float sliceLocation = 0.0f;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|1041", sliceLocationString);
    sscanf(sliceLocationString.c_str(), "%f", &sliceLocation);

    std::string imageInstanceString;
    int imageInstance = 0;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0013", imageInstanceString);
    sscanf(imageInstanceString.c_str(), "%d", &imageInstance);

    std::string imageAcquisitionString;
    int acquisitionNumber = 0;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0012", imageAcquisitionString);
    sscanf(imageAcquisitionString.c_str(), "%d", &acquisitionNumber);

    std::string seriesDescription;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|103e", seriesDescription);

    std::string studyDescription;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|1030", studyDescription);

    std::string referringPhysicianName;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0090", referringPhysicianName);

    std::string institutionName;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0080", institutionName);

    std::string modality;
    itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0060", modality);
    if (modality == "") {
      modality = "OT"; // other
    }

    if (stableUIDs) {
      itk::ExposeMetaData<std::string>(dictionarySlice, "0020|000e", newFusedSeriesInstanceUID);
      std::string derivedFusedSeriesInstanceUID = newFusedSeriesInstanceUID;
      std::string endString = ".3";
      if (derivedFusedSeriesInstanceUID.substr(derivedFusedSeriesInstanceUID.size() - 2, 2) == ".3")
        endString = ".4";
        // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
        derivedFusedSeriesInstanceUID = derivedFusedSeriesInstanceUID.substr(0, 64 - 3) + endString;
        newFusedSeriesInstanceUID = derivedFusedSeriesInstanceUID;
    }

    // how big is the image?
    ImageType2D::RegionType region;
    region = im2change->GetBufferedRegion();
    ImageType2D::SizeType size = region.GetSize();

    // create output image slice
    gdcm::PixelFormat pf = gdcm::PixelFormat::UINT8;
    pf.SetSamplesPerPixel(3);
    gdcm::SmartPointer<gdcm::Image> simage = new gdcm::Image;
    gdcm::Image &image = *simage;
    image.SetNumberOfDimensions(2);
    // typedef itk::Image<PixelType, 2> ImageType2D;
    ImageType2D::RegionType inRegion = im2change->GetLargestPossibleRegion();
    image.SetDimension(0, static_cast<unsigned int>(inRegion.GetSize()[0]));
    image.SetDimension(1, static_cast<unsigned int>(inRegion.GetSize()[1]));
    // image.SetDimension(2, m_Dimensions[2] );
    image.SetSpacing(0, im2change->GetSpacing()[0]);
    image.SetSpacing(1, im2change->GetSpacing()[1]);

    image.SetPixelFormat(pf);
    image.SetPhotometricInterpretation(gdcm::PhotometricInterpretation::RGB);
    image.SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
    // copy the DICOM tags over from inputImage to image
    gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));

    unsigned int size_of_pixel = sizeof(CImageType::PixelType);
    uint32_t len = size[0] * size[1] * size_of_pixel;
    pixeldata.SetByteValue(
      (char *)(((&((char *)buffer_fusedImage)[sliceNr * size[0] * size[1] * size_of_pixel]))), 
      len);
    image.SetDataElement(pixeldata);

    gdcm::SmartPointer<gdcm::File> file = new gdcm::File; // empty file
    gdcm::FileDerivation fd;
    const char ReferencedSOPClassUID[] = "1.2.840.10008.5.1.4.1.1.7"; // Secondary Capture
    // create a new frameOfRefenceUID
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    // do we need a new one here? - only if we don't find one from before
    if (frameOfReferenceUID == "")
      frameOfReferenceUID = fuid.Generate();

    fd.AddReference(ReferencedSOPClassUID, frameOfReferenceUID.c_str());
    fd.SetPurposeOfReferenceCodeSequenceCodeValue(
        121324); // segmentation  (see
               // https://github.com/malaterre/GDCM/blob/master/Source/MediaStorageAndFileFormat/gdcmFileDerivation.cxx)
    // CID 7203 Image Derivation
    // { "DCM",113072,"Multiplanar reformatting" },
    fd.SetDerivationCodeSequenceCodeValue(113076);
    fd.SetFile(*file);
    // If all Code Value are ok the filter will execute properly
    if (!fd.Derive()) {
      std::cerr << "Sorry could not derive using input info" << std::endl;
      return EXIT_FAILURE;
    }

    gdcm::DataSet &ds = fd.GetFile().GetDataSet();
    gdcm::Anonymizer ano;
    ano.SetFile(fd.GetFile());

    // now set all DICOM tags as needed
    using DictionaryType = itk::MetaDataDictionary;
    const DictionaryType &dictionaryIn = gdcmIO->GetMetaDataDictionary();
    using MetaDataStringType = itk::MetaDataObject<std::string>;
    auto itr = dictionaryIn.Begin();
    auto end = dictionaryIn.End();

    int seriesNumber = 0;
    // use the serieNumberStr to set the new series number
    if (seriesNumberStr != "") {
      seriesNumber = atoi(seriesNumberStr.c_str());
    }

    while (itr != end) {
      itk::MetaDataObjectBase::Pointer entry = itr->second;
      MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
      if (entryvalue) {
        std::string tagkey = itr->first;
        std::string labelId;
        bool found = itk::GDCMImageIO::GetLabelFromTag(tagkey, labelId);
        std::string tagvalue = entryvalue->GetMetaDataObjectValue();
        if (strcmp(tagkey.c_str(), "0020|0011") == 0) {
          seriesNumber = atoi(tagvalue.c_str());
        }
        // change window level from -400..600 to 150..180 (why don't we use the computed values? or 0 and 255?)
        if (strcmp(tagkey.c_str(), "0028|1050") == 0) {
          tagvalue = std::string("150");
        }
        if (strcmp(tagkey.c_str(), "0028|1051") == 0) {
          tagvalue = std::string("260");
        }

        unsigned int f1;
        unsigned int f2;
        sscanf(tagkey.c_str(), "%x|%x", &f1, &f2);
        ano.Replace(gdcm::Tag(f1, f2), tagvalue.c_str());
      }
      ++itr;
    }

    // StudyInstanceUID
    gdcm::Attribute<0x0020, 0x000d> at1_2; // Derivative Description
    at1_2.SetValue(StudyInstanceUID);
    ds.Replace(at1_2.GetAsDataElement());

    // StudyDescription
    gdcm::Attribute<0x0008, 0x1030> at1_3; // Derivative Description
    at1_3.SetValue(studyDescription);
    ds.Replace(at1_3.GetAsDataElement());

    //std::string patientID("");
    //itk::ExposeMetaData<std::string>(dictionarySlice, "0010|0020", patientID);
    gdcm::Attribute<0x0010, 0x0020> at1_4; // Derivative Description
    at1_4.SetValue(patientID);
    ds.Replace(at1_4.GetAsDataElement());

    //std::string patientName("");
    //itk::ExposeMetaData<std::string>(dictionarySlice, "0010|0010", patientName);
    gdcm::Attribute<0x0010, 0x0010> at1_5; // Derivative Description
    at1_5.SetValue(patientName);
    ds.Replace(at1_5.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x2111> at1; // Derivative Description
    if (overlayInfo->votemapMode) {
      at1.SetValue("Fused Vote Map, blue agree and yellow disagree");
    } else {
      at1.SetValue("Fused Segmentation");
    }
    ds.Replace(at1.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x0060> at2; // Derivative Description
    at2.SetValue(modality);
    ds.Replace(at2.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000E> at3;
    at3.SetValue(newFusedSeriesInstanceUID);
    ds.Replace(at3.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x103E> at4;
    std::string extension = " (fused segmentation)";
    if (overlayInfo->votemapMode) {
      extension = " (fused vote map)";
    }
    std::ostringstream value;
    value.str("");
    value << seriesDescription;
    // This is a long string and there is a 64 character limit in the
    // standard
    unsigned lengthDesc = value.str().length();
    std::string seriesDesc(value.str(), 0, lengthDesc + extension.length() > 64 ? 64 - extension.length() : lengthDesc + extension.length());
    // itk::EncapsulateMetaData<std::string>(dictionary, "0008|103e", seriesDesc + extension);
    at4.SetValue(seriesDesc + extension);
    ds.Replace(at4.GetAsDataElement());

    // seriesInstance
    gdcm::Attribute<0x0020, 0x0011> at5;
    at5.SetValue(1000 + seriesNumber + 3);
    ds.Replace(at5.GetAsDataElement());

    // use a unique SOPInstanceUID
    gdcm::Attribute<0x0008, 0x0018> at6;
    at6.SetValue(newFusedSOPInstanceUID);
    ds.Replace(at6.GetAsDataElement());

    //ano.Replace(gdcm::Tag(0x0008, 0x0018), newFusedSOPInstanceUID.c_str());

    // set image type to derived
    gdcm::Attribute<0x0008, 0x0008> at_image_type;
    static const gdcm::CSComp values[] = {"DERIVED","SECONDARY","OTHER"};
    at_image_type.SetValues( values, 3, true ); // true => copy data !
    if ( ds.FindDataElement( at_image_type.GetTag() ) ) {
      const gdcm::DataElement &de = ds.GetDataElement( at_image_type.GetTag() );
      //at_image_type.SetFromDataElement( de );
      // Make sure that value #1 is at least 'DERIVED', so override in all cases:
      at_image_type.SetValue( 0, values[0] );
      at_image_type.SetValue( 1, values[1] );
      at_image_type.SetValue( 2, values[2] );
    }
    ds.Replace( at_image_type.GetAsDataElement() );

    // software version
    std::string AIVersion("AI version: ");
    if (const char *env_p = std::getenv("VERSION")) { // part of the Dockerfile for building this application
      AIVersion += std::string(env_p);
    }
    std::string ContainerVersion("container version: ");
    if (const char *env_p = std::getenv("CONTAINER_VERSION")) { // part of the Dockerfile for building this application
      ContainerVersion += std::string(env_p);
    }
    gdcm::Attribute<0x0018, 0x1020> at_software_version;
    static const gdcm::LOComp version_values[] = {(std::string("imageAndMask2Fused ") + versionString).c_str(), ContainerVersion.c_str(), AIVersion.c_str()};
    at_software_version.SetValues( version_values, 3, true ); // true => copy data !
    if ( ds.FindDataElement( at_image_type.GetTag() ) ) {
      const gdcm::DataElement &de = ds.GetDataElement( at_image_type.GetTag() );
      //at_image_type.SetFromDataElement( de );
      // Make sure that value #1 is at least 'DERIVED', so override in all cases:
      at_software_version.SetValue( 0, version_values[0] );
      at_software_version.SetValue( 1, version_values[1] );
      at_software_version.SetValue( 2, version_values[2] );
    }
    ds.Replace( at_software_version.GetAsDataElement() );

    // image position patient from input
    // These values are actually not getting written to the files (RGB has no origin, values are 0\0\0, but see set origin further down)
    gdcm::Attribute<0x0020, 0x0032> at7;
    at7.SetValue(origin3D[0], 0);
    at7.SetValue(origin3D[1], 1);
    at7.SetValue(origin3D[2], 2);
    ds.Replace(at7.GetAsDataElement());
    std::ostringstream value2;
    value2.str("");
    at7.Print(value2);
    // fprintf(stdout, "origin is now supposed to be: %lf\\%lf\\%lf %s\n", origin3D[0], origin3D[1], origin3D[2], value2.str().c_str());
    // For RGB we can set this to make sure they show up at the right location in Horos/OsiriX
    image.SetOrigin(0, origin3D[0]);
    image.SetOrigin(1, origin3D[1]);
    image.SetOrigin(2, origin3D[2]);

    gdcm::Attribute<0x0018, 0x0050> at8;
    at8.SetValue(sliceThickness);
    ds.Replace(at8.GetAsDataElement());

    // TODO: this does not work. We end up with the wrong ImageOrientationPatient information in the output.
    // Either we do not write these attributes for a secondary capture image or we use the correct once 
    // (for linking in Horos for example).
    image.SetDirectionCosines(imageOrientationField);

    gdcm::Attribute<0x0020, 0x1041> at11;
    at11.SetValue(sliceLocation);
    ds.Replace(at11.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x0012> at12;
    at12.SetValue(1000 + acquisitionNumber + 3);
    ds.Replace(at12.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x0013> at13;
    at13.SetValue(imageInstance); // count starts at 1 and increments for all slices
    ds.Replace(at13.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x0052> at14;
    at14.SetValue(frameOfReferenceUID.c_str());
    ds.Replace(at14.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x0010> at16;
    at16.SetValue(studyID.c_str());
    ds.Replace(at16.GetAsDataElement());

    // accession_number
    gdcm::Attribute<0x0008, 0x0050> at16_2;
    at16_2.SetValue(accession_number.c_str());
    ds.Replace(at16_2.GetAsDataElement());

    std::time_t t = std::time(nullptr);
    char mbstr[100];
    gdcm::Attribute<0x0008, 0x0021> at17;
    if (std::strftime(mbstr, sizeof(mbstr), "%Y%m%d", std::localtime(&t))) {
      at17.SetValue(mbstr);
      ds.Replace(at17.GetAsDataElement());
    }

    gdcm::Attribute<0x0008, 0x0031> at18;
    if (std::strftime(mbstr, sizeof(mbstr), "%H%M%S", std::localtime(&t))) {
      at18.SetValue(mbstr);
      ds.Replace(at18.GetAsDataElement());
    }

    gdcm::Attribute<0x0008, 0x0090> at19;
    at19.SetValue(referringPhysicianName.c_str());
    ds.Replace(at19.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x0080> at20;
    at20.SetValue(institutionName.c_str());
    ds.Replace(at20.GetAsDataElement());

    boost::filesystem::path p(filenames[sliceNr]);
    std::string filename_without_extension = (p.filename().string()).substr(0, (p.filename().string()).find_last_of("."));
    std::string output_dir_name = "fused";
    if (overlayInfo->votemapMode) {
      output_dir_name = "fused_vote_map";
    }

    boost::filesystem::path p_out = outputDir + boost::filesystem::path::preferred_separator + output_dir_name + boost::filesystem::path::preferred_separator +
      newFusedSeriesInstanceUID + boost::filesystem::path::preferred_separator + filename_without_extension.c_str() + ".dcm";
    if (!itksys::SystemTools::FileIsDirectory(p_out.parent_path().c_str())) {
      // create the output directory
      create_directories(p_out.parent_path());
    }

    // some tags are not written here. Like the SOPInstanceUID is always the same and the ImageNumber as well.
    gdcm::ImageWriter writer;
    writer.SetImage(image);
    writer.SetFile(fd.GetFile());
    // std::ostringstream o;
    // o << outputSeries << "/dicom" << i << ".dcm";
    writer.SetFileName(p_out.c_str());
    if (!writer.Write()) {
      return EXIT_FAILURE;
    } else {
      if (verbose) {
        std::cout << "Wrote fused image slice: " << p_out.c_str() << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;

}


CImageType::Pointer computeFusedImage(ImageType3D::Pointer inputImage, MaskImageType3D::Pointer maskImage, OverlayInfos_t *overlayInfo) {

  // create the output fused image
  using CPixelType = itk::RGBPixel<unsigned char>;
  using CImageType3D = itk::Image<CPixelType, 3>;
  using CWriterType3D = itk::ImageFileWriter<CImageType3D>;
  CImageType3D::Pointer fused = CImageType3D::New();
  CImageType3D::RegionType fusedRegion = inputImage->GetLargestPossibleRegion();
  fused->SetRegions(fusedRegion);
  fused->Allocate();
  fused->FillBuffer(itk::NumericTraits<CPixelType>::ZeroValue());
  fused->SetOrigin(inputImage->GetOrigin());
  fused->SetSpacing(inputImage->GetSpacing());
  fused->SetDirection(inputImage->GetDirection());

/*  itk::MetaDataDictionary &dictionarySlice = inputImage->GetMetaDataDictionary();

  std::string SeriesNumber("");
  std::string studyID(""); // this is the StudyInstanceUID
  std::string AcquisitionNumber("");
  std::string InstanceNumber("");
  std::string frameOfReferenceUID("");
  std::string InstitutionName("");
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|000d", studyID);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0008|0080", InstitutionName);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0011", SeriesNumber);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0012", AcquisitionNumber);
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0013", InstanceNumber); // keep that number
  itk::ExposeMetaData<std::string>(dictionarySlice, "0020|0052", frameOfReferenceUID);

  int newSeriesNumber = atoi(SeriesNumber.c_str()) + 1003;
*/
  using LabelType = unsigned short;
  using ShapeLabelObjectType = itk::ShapeLabelObject<LabelType, 3>;
  using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;

  // create a labelMap from the mask image
  using I2LType = itk::LabelImageToShapeLabelMapFilter<MaskImageType3D, LabelMapType>;
  auto i2l = I2LType::New();
  i2l->SetInput(maskImage);
  // not needed if we use the votemapmode
  i2l->Update();
  LabelMapType *labelMap = i2l->GetOutput();
  if (verbose) {
    fprintf(stdout, "imageAndMask2Fused: label to shape filter done. Created %lu label%s.\n", labelMap->GetNumberOfLabelObjects(), labelMap->GetNumberOfLabelObjects()!=1?"s":""); 
    fflush(stdout);
  }

  // compute correct histogram settings for the output image (its RGB so we need to scale the input)
  using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageType3D>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(inputImage);
  imageCalculatorFilter->Compute();
  int minGray = imageCalculatorFilter->GetMinimum();
  int maxGray = imageCalculatorFilter->GetMaximum();

  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType3D>;
  HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(inputImage);
  int histogramSize = 1024;
  histogramGenerator->SetNumberOfBins(histogramSize);
  histogramGenerator->SetHistogramMin(minGray);
  histogramGenerator->SetHistogramMax(maxGray);
  histogramGenerator->Compute();
  using HistogramType = HistogramGeneratorType::HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();
  // set in calling function
  double lowerT = overlayInfo->brightnessContrastLL; // 0.01;
  double upperT = overlayInfo->brightnessContrastUL; // 0.999;
  double t1 = -1;
  double t2 = -1;
  double sum = 0;
  double total = 0;
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    total += histogram->GetFrequency(bin, 0);
  }
  for (unsigned int bin = 0; bin < histogramSize; bin++) {
    double f = histogram->GetFrequency(bin, 0) / total;
    // fprintf(stdout, "bin %d, value is %f\n", bin, f);
    sum += f;
    if (t1 == -1 && sum > lowerT) {
      t1 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
    }
    if (t2 == -1 && sum > upperT) {
      t2 = minGray + (maxGray - minGray) * (bin / (float)histogramSize);
      break;
    }
  }

  // create empty color channels at 16bit floating point precision
  typedef float FPixelType;
  // typedef itk::Image<FPixelType, 2> FloatImageType;
  using FloatImageType = itk::Image<FPixelType, 3>;
  FloatImageType::Pointer red_channel = FloatImageType::New();
  FloatImageType::Pointer green_channel = FloatImageType::New();
  FloatImageType::Pointer blue_channel = FloatImageType::New();
  fusedRegion = inputImage->GetLargestPossibleRegion();
  red_channel->SetRegions(fusedRegion);
  red_channel->Allocate();
  red_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  red_channel->SetOrigin(inputImage->GetOrigin());
  red_channel->SetSpacing(inputImage->GetSpacing());
  red_channel->SetDirection(inputImage->GetDirection());
  green_channel->SetRegions(fusedRegion);
  green_channel->Allocate();
  green_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  green_channel->SetOrigin(inputImage->GetOrigin());
  green_channel->SetSpacing(inputImage->GetSpacing());
  green_channel->SetDirection(inputImage->GetDirection());
  blue_channel->SetRegions(fusedRegion);
  blue_channel->Allocate();
  blue_channel->FillBuffer(itk::NumericTraits<FPixelType>::Zero);
  blue_channel->SetOrigin(inputImage->GetOrigin());
  blue_channel->SetSpacing(inputImage->GetSpacing());
  blue_channel->SetDirection(inputImage->GetDirection());

  // we have two modes, either Votemap mode or color by regions of interest
  if (overlayInfo->votemapMode) {
    overlayInfo->obtained_max_votemap_value = -1;
    itk::ImageRegionIterator<MaskImageType3D> maskIterator(maskImage, fusedRegion);
    itk::ImageRegionIterator<FloatImageType> redSIterator(red_channel, fusedRegion);
    itk::ImageRegionIterator<FloatImageType> greenSIterator(green_channel, fusedRegion);
    itk::ImageRegionIterator<FloatImageType> blueSIterator(blue_channel, fusedRegion);
    maskIterator.GoToBegin();
    redSIterator.GoToBegin();
    greenSIterator.GoToBegin();
    blueSIterator.GoToBegin();
    while (!maskIterator.IsAtEnd() && !redSIterator.IsAtEnd() && !greenSIterator.IsAtEnd() && !blueSIterator.IsAtEnd()) {
      unsigned short maskValue = maskIterator.Get();
      if (maskValue > overlayInfo->obtained_max_votemap_value) {
        overlayInfo->obtained_max_votemap_value = maskValue;
      }

      std::vector<float> col = labelColorsVotemap[0];
      if (maskValue > 0) {
        if (maskValue >= (overlayInfo->votemapAgree * overlayInfo->votemapMax)) {
          col = labelColorsVotemap[2];
        } else {
          col = labelColorsVotemap[1];
        }
      }
      redSIterator.Set(col[0]/255.0f);
      greenSIterator.Set(col[1]/255.0f);
      blueSIterator.Set(col[2]/255.0f);
      ++maskIterator;
      ++redSIterator;
      ++greenSIterator;
      ++blueSIterator;
    }
  } else {
    for (unsigned int n_idx = 0; n_idx < labelMap->GetNumberOfLabelObjects(); n_idx++) {
      unsigned int n = n_idx;
      ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n); // the label number is the connected component number - not the one label as mask
      int label = labelObject->GetLabel(); // it might be that all regions have the same label and only n is a good choice for the color

      // color is
      std::vector<float> col = labelColors2[ (label % (labelColors2.size()-1)) +1 ];
      itk::Index<3U> index;
      for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
        index = labelObject->GetIndex(pixelId);
        // set this position in all three color images, values are between 0 and 1
        red_channel->SetPixel(index, col[0]/255.0); 
        green_channel->SetPixel(index, col[1]/255.0);
        blue_channel->SetPixel(index, col[2]/255.0);
      }
    }
  }

  // smooth all three channels independently from each other
  using GFilterType = itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>;
  auto gaussFilterR = GFilterType::New();
  gaussFilterR->SetInput(red_channel);
  gaussFilterR->SetVariance(1.5f);
  gaussFilterR->Update();
  FloatImageType::Pointer smoothRed = gaussFilterR->GetOutput();

  auto gaussFilterG = GFilterType::New();
  gaussFilterG->SetInput(green_channel);
  gaussFilterG->SetVariance(1.5f);
  gaussFilterG->Update();
  FloatImageType::Pointer smoothGreen = gaussFilterG->GetOutput();

  auto gaussFilterB = GFilterType::New();
  gaussFilterB->SetInput(blue_channel);
  gaussFilterB->SetVariance(1.5f);
  gaussFilterB->Update();
  FloatImageType::Pointer smoothBlue = gaussFilterB->GetOutput();

  itk::ImageRegionIterator<FloatImageType> redSIterator(smoothRed, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> greenSIterator(smoothGreen, fusedRegion);
  itk::ImageRegionIterator<FloatImageType> blueSIterator(smoothBlue, fusedRegion);
  itk::ImageRegionIterator<ImageType3D> inputIterator(inputImage, fusedRegion);
  itk::ImageRegionIterator<CImageType3D> fusedIterator(fused, fusedRegion);
  //itk::ImageRegionIterator<MaskImageType3D> maskIterator(maskImage, fusedRegion);

  inputIterator.GoToBegin();
  fusedIterator.GoToBegin();
  redSIterator.GoToBegin();
  greenSIterator.GoToBegin();
  blueSIterator.GoToBegin();
  // maskIterator.GoToBegin();
  float f = 0.6; // weight of the underlay, at 0.1 only mask is visible
  if (overlayInfo->votemapMode) {
    f = 0.3; // in votemap mode we want to see more of the underlay
  }
  float red, green, blue;
  while (!inputIterator.IsAtEnd() && !fusedIterator.IsAtEnd() && !redSIterator.IsAtEnd() && !greenSIterator.IsAtEnd() && !blueSIterator.IsAtEnd()) {
    float scaledP = ((float) inputIterator.Get() - t1) / (t2 - t1);
    CPixelType value = fusedIterator.Value();

    red = redSIterator.Get();
    green = greenSIterator.Get();
    blue = blueSIterator.Get();

    // only in votemap mode the alpha value is mask value dependent
    /*if (overlayInfo->votemapMode) {
      unsigned short maskValue = maskIterator.Get();
      if (maskValue > 0) {
        if (maskValue >= (overlayInfo->votemapAgree * overlayInfo->votemapMax)) {
          f = 0.6; // col = labelColorsVotemap[2];
        } else {
          // a yellow region of uncertainty should be more visible, use a lower alpha value for the blending of the underlay
          f = 0.3; // col = labelColorsVotemap[1];
        }
      } else {
        f = 0.6; // default
      }
    }*/

    // alpha blend
    red = f * scaledP + red * (1 - f);
    green = f * scaledP + green * (1 - f);
    blue = f * scaledP + blue * (1 - f);

    red = std::min<float>(1, std::max<float>(0, red));
    green = std::min<float>(1, std::max<float>(0, green));
    blue = std::min<float>(1, std::max<float>(0, blue));

    value.SetRed((int)(red * 255));
    value.SetGreen((int)(green * 255));
    value.SetBlue((int)(blue * 255));
    fusedIterator.Set(value);

    ++inputIterator; 
    ++fusedIterator;
    ++redSIterator;
    ++greenSIterator;
    ++blueSIterator;
  }
  return fused;
}

int main(int argc, char *argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  versionString = std::string("0.0.1.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  command.SetVersion(versionString.c_str());
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("ImageAndMask2Fused: Creates a fused image series from image and mask pair (DICOM format).");
  command.SetCategory("image conversion");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("maskdir", "Directory with input mask DICOM series (0 - background, 1 - foreground).", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for images/, labels/, fused/, and reports/ folder as DICOM. The redcap/ folder contains series folders with output.json files for REDCap imports.", MetaCommand::STRING, true);

  command.SetOption(
    "UIDFixed", "u", false,
    "If enabled identifiers are stable - will not change for a given input. This allows image series to overwrite each other - assuming that the PACS "
    "supports this overwrite mode. By default the SeriesInstanceUID and SOPInstanceUID values are generated again every time the processing is done.");
  command.SetOptionLongTag("UIDFixed", "uid-fixed");

  command.SetOption("TitleText", "t", false, "Specify the title text on the fused image.");
  command.SetOptionLongTag("TitleText", "title");
  command.AddOptionField("TitleText", "title", MetaCommand::STRING, false);

  command.SetOption("SubTitleText", "s", false, "Specify a sub-title text for version information.");
  command.SetOptionLongTag("SubTitleText", "subtitle");
  command.AddOptionField("SubTitleText", "subtitle", MetaCommand::STRING, false);

  command.SetOption("Info", "i", false, "Specify an info message that will appear in the bottom left corner. This option can be used to identify the version/container used for creating the segmentation.");
  command.SetOptionLongTag("Info", "info");
  command.AddOptionField("Info", "info", MetaCommand::STRING, false);

  command.SetOption("VotemapMax", "o", false, "Enable votemap display for overlapping regions of interest. The provided mask volume is assumed to be of a votemap type. The provided value should be the maximum agreement value possible.");
  command.SetOptionLongTag("VotemapMax", "votemapmax");
  command.AddOptionField("VotemapMax", "votemapmax", MetaCommand::FLOAT, false);

  command.SetOption("VotemapAgree", "a", false, "If vote map mode is enabled (see VotemapMax) you can specify a level for disagreement relative to the VotemapMax (in percent). A value of >= 0.75 for a VotemapMax value of 5 would color voxel in the overlay blue if the votemap is larger than 4 and yellow otherwise. Only values of 0 in the vote map are displayed with a transparent overlay.");
  command.SetOptionLongTag("VotemapAgree", "votemapagree");
  command.AddOptionField("VotemapAgree", "votemapagree", MetaCommand::FLOAT, false, "0.75");

  command.SetOption("BrightnessContrastLL", "d", false, "Set threshold for brightness / contrast based on cummulative histogram lower limit (percentage dark pixel 0.01).");
  command.SetOptionLongTag("BrightnessContrastLL", "brightness-contrast-ll");
  command.AddOptionField("BrightnessContrastLL", "value", MetaCommand::FLOAT, false);

  command.SetOption("BrightnessContrastUL", "b", false, "Set threshold for brightness / contrast based on cummulative histogram upper limit (percentage bright pixel 0.999).");
  command.SetOptionLongTag("BrightnessContrastUL", "brightness-contrast-ul");
  command.AddOptionField("BrightnessContrastUL", "value", MetaCommand::FLOAT, false);

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  if (command.GetOptionWasSet("UIDFixed"))
    stableUIDs = true;

  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  std::string input = command.GetValueAsString("indir");
  std::string mask = command.GetValueAsString("maskdir");
  std::string output = command.GetValueAsString("outdir");

  if (input.size() == 0 || output.size() == 0 || mask.size() == 0) {
    return 1;
  }

  float votemapMax = 0.0;
  bool votemapMode = false;
  if (command.GetOptionWasSet("VotemapMax")) {
    votemapMax = command.GetValueAsFloat("VotemapMax", "votemapmax");
    votemapMode = true;
  }
  float votemapAgree = 0.75;
  if (command.GetOptionWasSet("VotemapAgree")) {
    votemapAgree = command.GetValueAsFloat("VotemapAgree", "votemapagree");
    if (votemapAgree < 0 || votemapAgree > 1.0) {
      votemapAgree = std::min<float>(std::max<float>(0.0, votemapAgree), 1.0);
      fprintf(stdout, "Warning: votemap agree value not between 0 and 1. Adjusted to %f.\n", votemapAgree);
    }
  }
  float brightness_contrast_ll = 0.01;
  float brightness_contrast_ul = 0.999;
  float brightnesscontrast_ll = brightness_contrast_ll;
  float brightnesscontrast_ul = brightness_contrast_ul;
  if (command.GetOptionWasSet("BrightnessContrastLL")) {
    brightnesscontrast_ll = command.GetValueAsFloat("BrightnessContrastLL", "value");
    if (brightnesscontrast_ll < 0 || brightnesscontrast_ll > 1.0) {
      fprintf(stdout, "Warning: lower brightness values not between 0 and 1. Adjusted to 0.01.\n");
      brightnesscontrast_ll = 0.01;
    }
  }
  if (command.GetOptionWasSet("BrightnessContrastUL")) {
    brightnesscontrast_ul = command.GetValueAsFloat("BrightnessContrastUL", "value");
    if (brightnesscontrast_ul < 0 || brightnesscontrast_ul > 1.0) {
      fprintf(stdout, "Warning: upper brightness values not between 0 and 1. Adjusted to 0.999.\n");
      brightnesscontrast_ul = 0.999;
    }
  }
  if (brightnesscontrast_ul < brightnesscontrast_ll) {
    float tmp = brightnesscontrast_ll;
    brightnesscontrast_ll = brightnesscontrast_ul;
    brightnesscontrast_ul = tmp;
  }
  brightness_contrast_ll = brightnesscontrast_ll;
  brightness_contrast_ul = brightnesscontrast_ul;
  if (verbose) {
    fprintf(stdout, "create report with brightness/contrast settings %.03f %.03f\n", brightness_contrast_ll, brightness_contrast_ul);
  }

  std::string infoMessage = command.GetValueAsString("Info", "info");
  std::string titleText = command.GetValueAsString("TitleText", "title");
  std::string subTitleText = command.GetValueAsString("SubTitleText", "subtitle");

  if (infoMessage == "") {
    // set default subtitle text
    infoMessage = std::string("For Research Use Only - Not for use in diagnostic procedures.");
  }

  OverlayInfos_t overlayInfos = {titleText, subTitleText, infoMessage, brightness_contrast_ll,  brightness_contrast_ul, votemapMode, votemapMax, votemapAgree};

  MaskImageType3D::Pointer maskImage = readMaskImageSerie(mask);
  if (maskImage == nullptr) {
    std::cerr << "Could not read mask image series from: " << mask << std::endl;
    return EXIT_FAILURE;
  } 

  ImageType3D::Pointer inputImage = readImageSerie(input);
  if (inputImage == nullptr) {
    std::cerr << "Could not read image series from: " << input << std::endl;
    return EXIT_FAILURE;
  } 

  CImageType::Pointer fusedImage = computeFusedImage(inputImage, maskImage, &overlayInfos);
  if (fusedImage == nullptr) {
    std::cerr << "Could not compute fused image!" << std::endl;
    return EXIT_FAILURE;
  } 

  saveFusedImageSeries(fusedImage, output, input, &overlayInfos);

  return EXIT_SUCCESS;  

}