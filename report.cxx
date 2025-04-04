#include "report.h"
#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDefs.h"
#include "gdcmGlobal.h"
#include "gdcmImage.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmPhotometricInterpretation.h"
#include "gdcmUIDGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"


#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <gdcmFile.h>
#include <gdcmImage.h>
#include <math.h>
#include <codecvt>
#include <sstream>

#include <boost/date_time.hpp>
#include <boost/math/distributions/lognormal.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

#define WIDTH 2400
#define HEIGHT 2400

unsigned char image_buffer[HEIGHT][WIDTH];
//unsigned char image_buffer512[512][512];
//unsigned char image_buffer828[512][828];
unsigned int image_buffer_gen_size[2] = {0,0};
unsigned char **image_buffer_gen = NULL;

Report *getDefaultReportStruct() {
  Report *report = new Report;
  // create some example UIDs for the report
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->SeriesInstanceUID = fuid.Generate();
  }
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->StudyInstanceUID = fuid.Generate();
  }
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->SOPInstanceUID = fuid.Generate();
  }
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->FrameOfReferenceUID = fuid.Generate();
  }
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->StudyID = fuid.Generate();
    report->StudyID = report->StudyID.substr(0, 16);
  }
  {
    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    report->AccessionNumber = fuid.Generate();
    report->AccessionNumber = report->AccessionNumber.substr(0, 16);
  }

  report->SeriesDescription = std::string("Report");
  report->StudyDate = std::string("");
  report->StudyTime = std::string("");
  report->measures = std::vector<std::map<std::string, std::string>>();
  report->key_fact = std::string("");
  report->key_unit = std::string("");
  report->VersionString = std::string("");
  report->TitleText = std::string("MMIV.no report. Made with AI");
  return report;
}

// provide string as str.c_str()
// return if conversion was successful
int stdstring2stdwstring(std::wstring *wstr, const char* ptr) {
    wstr->clear();
    std::mbtowc(nullptr, 0, 0); // reset the conversion state
    const char* end = ptr + std::strlen(ptr);
    int ret{};
    for (wchar_t wc; (ret = std::mbtowc(&wc, ptr, end - ptr)) > 0; ptr += ret) {
        //std::wcout << wc;
        wstr->push_back(wc);
    }
    std::wcout << '\n';
    return ret;
}

// draws into image_buffer
void draw_bitmap(FT_Bitmap *bitmap, int width, int height, FT_Int x, FT_Int y) {
  FT_Int i, j, p, q;
  FT_Int x_max = x + bitmap->width;
  FT_Int y_max = y + bitmap->rows;

  for (i = x, p = 0; i < x_max; i++, p++) {
    for (j = y, q = 0; j < y_max; j++, q++) {
      if (i < 0 || j < 0 || i >= width || j >= height)
        continue;

      image_buffer[j][i] |= bitmap->buffer[q * bitmap->width + p];
    }
  }
}


void addBar(CImageType::Pointer img, int height) {
  CImageType::RegionType kregion = img->GetLargestPossibleRegion();
  using ImageSizeType = typename CImageType::SizeType;
  ImageSizeType regionSize;
  regionSize[0] = kregion.GetSize()[0];
  regionSize[1] = height;

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

  CImageType::RegionType::IndexType regionStart;
  regionStart = kregion.GetIndex();
  regionStart[1] = height;
  regionSize[0] = kregion.GetSize()[0];
  regionSize[1] = 1;

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

void addMarker(char *buffer, int posx, int posy, int fontSize) {
  unsigned char *bvals = (unsigned char *)buffer;

  // draw 4 lines to highlight the location (fontSize was 9 before)
  int line_length = fontSize/2.0; // in pixel
  int offset = 7.0/9.0*fontSize; // in pixel away from that location
  int line_width = 3.0/9.0*fontSize; // one pixel wide

  // top bar
  for (int j = -floor(line_width/2.0); j < floor(line_width/2.0); j++) {
    for (int i = 0; i < line_length; i++) {
        int idx = (posx+j) + (posy - offset - i) * image_buffer_gen_size[0];
        if (idx < 0 || idx > image_buffer_gen_size[1]*image_buffer_gen_size[0])
          continue;
        bvals[idx] = 255; //(unsigned char)std::max(0.0f, std::min(current_image_max_value, 1.0f));
    }
  }
  // bottom bar
  for (int j = -floor(line_width/2.0); j < floor(line_width/2.0); j++) {
    for (int i = 0; i < line_length; i++) {
        int idx = (posx+j) + (posy + offset + i - 1) * image_buffer_gen_size[0];
        if (idx < 0 || idx > image_buffer_gen_size[1]*image_buffer_gen_size[0])
          continue;
        bvals[idx] = 255; //(unsigned char)std::max(0.0f, std::min(current_image_max_value, 1.0f));
    }
  }

  // left bar
  for (int j = -floor(line_width/2.0); j < floor(line_width/2.0); j++) {
    for (int i = 0; i < line_length; i++) {
        int idx = (posx - offset - i) + (posy + j) * image_buffer_gen_size[0];
        if (idx < 0 || idx > image_buffer_gen_size[1]*image_buffer_gen_size[0])
          continue;
        bvals[idx] = 255; //(unsigned char)std::max(0.0f, std::min(current_image_max_value, 1.0f));
    }
  }

  // right bar
  for (int j = -floor(line_width/2.0); j < floor(line_width/2.0); j++) {
    for (int i = 0; i < line_length; i++) {
        int idx = (posx + offset + i - 1) + (posy + j) * image_buffer_gen_size[0];
        if (idx < 0 || idx > image_buffer_gen_size[1]*image_buffer_gen_size[0])
          continue;
        bvals[idx] = 255; //(unsigned char)std::max(0.0f, std::min(current_image_max_value, 1.0f));
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

// generic method using dynamic array for image_buffer_gen
void addToReportGen(char *buffer, std::string font_file, int font_size, std::string sstext, int posx, int posy, float radiants) {
  FT_Library library;


  //std::wstring stext;

  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>, wchar_t> convert;
  std::wstring stext = convert.from_bytes(sstext);

  //if (!stdstring2stdwstring(&stext, sstext.c_str())) {
  //  fprintf(stderr, "Error converting string to wstring\n");
  //  return;
  //}

  bool verbose = 1;
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
  // gdcm::ImageReader reader;

  // unsigned long len = WIDTH * HEIGHT * 8;
  //  char *buffer = new char[len];

  error = FT_Init_FreeType(&library); /* initialize library */

  if (error != 0) {
    fprintf(stderr, "\033[0;31mError\033[0m: The freetype library could not be initialized with this font.\n");
    return;
  }

  int start_px = 40;
  int start_py = 40;
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
  target_height = 512;

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




/*  const char *text = sstext.c_str();
  for (n = 0; n < num_chars; n++) {
    if (text[n] == '\n') {
      continue; // ignore newlines
    }

    FT_Set_Transform(face, &matrix, &pen);

    error = FT_Load_Char(face, stext[n], FT_LOAD_RENDER);
    if (error) {
      fprintf(stdout, "\033[0;31mError\033[0m:: [addToReportGen] could not load character: '%c'\n", text[n]);
      continue;
    }

    draw_bitmap_gen(&slot->bitmap, image_buffer_gen_size[0], image_buffer_gen_size[1], slot->bitmap_left, target_height - slot->bitmap_top);

    pen.x += slot->advance.x;
    pen.y += slot->advance.y;
  } */

  int nn = 0; 
  for (std::wstring::iterator it = stext.begin(); it != stext.end(); it++) {
    //wchar_t c = *it;
    wchar_t c;
    int ret = get_mb(&c, sstext.c_str(), nn);
    nn++;
    if (c == '\n') {
      continue; // ignore newlines
    }
    //fprintf(stdout, "char: %s\n", c);

    FT_Set_Transform(face, &matrix, &pen);

    FT_UInt glyph_index = FT_Get_Char_Index( face, *it );
    error = FT_Load_Glyph(face, glyph_index, FT_LOAD_RENDER);
    if (error) {
      fprintf(stdout, "\033[0;31mError\033[0m:: [addToReportGen] could not load character: '%ls'\n", &c);
      continue;
    }

    /*FT_Glyph glyph;
    error = FT_Get_Glyph( face->glyph, &glyph );
    if (error) {
      fprintf(stdout, "\033[0;31mError\033[0m:: [addToReportGen] could not get glyph: '%ls'\n", it);
      continue;
    } */

    // add here 
    // FT_Get_Glyph(face->glyph, &glyph);
    // FT_Glyph_To_Bitmap(&glyph, FT_RENDER_MODE_NORMAL, 0, 1);

    draw_bitmap_gen(&slot->bitmap, image_buffer_gen_size[0], image_buffer_gen_size[1], slot->bitmap_left, target_height - slot->bitmap_top);

    pen.x += slot->advance.x;
    pen.y += slot->advance.y;
  } 



  FT_Done_Face(face);

  // int px = start_px;
  // int py = start_py + 10.0 * (font_size);
  //  (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));

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

void addToReport(char *buffer, std::string font_file, int font_size, std::string stext, int posx, int posy, float radiants) {
  FT_Library library;

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
  // gdcm::ImageReader reader;

  // unsigned long len = WIDTH * HEIGHT * 8;
  //  char *buffer = new char[len];

  error = FT_Init_FreeType(&library); /* initialize library */

  if (error != 0) {
    fprintf(stderr, "\033[0;31mError\033[0m: The freetype library could not be initialized with this font.\n");
    return;
  }

  int start_px = -10;
  int start_py = -10;
  int text_lines = 1;

  float repeat_spacing = 2.0f;
  int xmax = WIDTH;
  int ymax = HEIGHT;
  num_chars = stext.size();

  int px = posx;
  // WIDTH - ((num_chars + 2) * font_size - start_px);
  // int py = start_py + (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));
  int py = posy;
  // start_py + 2.0 * (font_size);

  memset(image_buffer, 0, sizeof(unsigned char) * HEIGHT * WIDTH);

  angle = radiants;
  target_height = HEIGHT;

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
  pen.x = (num_chars * 1 * 64);
  pen.y = (target_height - 80) * 64; // the 60 here is related to the font size!
  const char *text = stext.c_str();
  for (n = 0; n < num_chars; n++) {
    if (text[n] == '\n') {
      continue; // ignore newlines
    }

    /* set transformation */
    FT_Set_Transform(face, &matrix, &pen);

    error = FT_Load_Char(face, text[n], FT_LOAD_RENDER);
    if (error) {
      fprintf(stdout, "\033[0;31mError\033[0m: Could not load character: '%c'\n", text[n]);
      continue; /* ignore errors */
    }
    /* now, draw to our target surface (convert position) draws into image_buffer */
    draw_bitmap(&slot->bitmap, WIDTH, HEIGHT, slot->bitmap_left, target_height - slot->bitmap_top);

    /* increment pen position */
    pen.x += slot->advance.x;
    pen.y += slot->advance.y;
  }

  FT_Done_Face(face);

  // int px = start_px;
  // int py = start_py + 10.0 * (font_size);
  //  (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));

  float current_image_min_value = 0.0f;
  float current_image_max_value = 255.0f;
  unsigned char *bvals = (unsigned char *)buffer;
  for (int yi = 0; yi < HEIGHT; yi++) {
    for (int xi = 0; xi < WIDTH; xi++) {
      if (image_buffer[yi][xi] == 0)
        continue;
      // I would like to copy the value from image over to
      // the buffer. At some good location...

      int newx = px + xi;
      int newy = py + yi;
      int idx = newy * xmax + newx;
      if (newx < 0 || newx >= xmax || newy < 0 || newy >= ymax)
        continue;
      if (image_buffer[yi][xi] == 0)
        continue;

      // instead of blending we need to use a fixed overlay color
      // we have image information between current_image_min_value and current_image_max_value
      // we need to scale the image_buffer by those values.
      float f = 0;
      float v = 1.0f * image_buffer[yi][xi] / 255.0; // 0 to 1 for color, could be inverted if we have a white background
      float w = 1.0f * bvals[idx] / current_image_max_value;
      float alpha_blend = (v + w * (1.0f - v));

      // fprintf(stdout, "%d %d: %d\n", xi, yi, bvals[idx]);
      bvals[idx] = (unsigned char)std::max(
          0.0f, std::min(current_image_max_value, current_image_min_value + (alpha_blend) * (current_image_max_value - current_image_min_value)));
    }
  }
}

// Create a report image
// distribution can be eiter "norm" or "lognorm" 
void saveReport(Report *report, std::string distribution, float mean_mean, float mean_stds, bool verbose) {

  std::string font_file = "Menlo.ttf";
  if (const char *env_p = std::getenv("REPORT_FONT_PATH")) {
    font_file = std::string(env_p);
  }
  if (!boost::filesystem::exists(font_file)) {
    fprintf(stderr, "\033[0;31mError\033[0m: no font provided. Set the environment variable REPORT_FONT_PATH to the location of a ttf file.\n");
    fflush(stderr);
    return;
  }

  FT_Library library;

  double angle = 0.0;
  int target_height = HEIGHT;
  int n, num_chars;

  // int font_length = 20;

  std::string font_path = font_file;
  int font_size = 11;
  int face_index = 0;

  FT_Face face;
  FT_GlyphSlot slot;
  FT_Matrix matrix; /* transformation matrix */
  FT_Vector pen;    /* untransformed origin  */
  FT_Error error;
  // gdcm::ImageReader reader;

  unsigned long len = WIDTH * HEIGHT; // in bytes
  char *buffer = new char[len];
  char *buffer_color = new char[WIDTH * HEIGHT * 3];


  error = FT_Init_FreeType(&library); /* initialize library */

  if (error != 0) {
    fprintf(stderr, "\033[0;31mError\033[0m: The freetype library could not be initialized with this font.\n");
    return;
  }

  // add the key image
  if (report->keyImage) {

    // for that keyImage    
    CImageType::RegionType kregion = report->keyImage->GetLargestPossibleRegion();
    int KWIDTH = kregion.GetSize()[0];
    int KHEIGHT = kregion.GetSize()[1];
    char *kbuffer = new char[KWIDTH * KHEIGHT];
    char *kbuffer_color = new char[KWIDTH * KHEIGHT * 3];
    memset(&kbuffer[0], 0, sizeof(char)*KWIDTH*KHEIGHT);
    memset(&kbuffer_color[0], 0, sizeof(char)*KWIDTH*KHEIGHT*3);

    // create space for image_buffer unsigned char image_buffer828[512][828];
    if (image_buffer_gen) {
      // delete its memory before we create a new array
      for (int i = 0; i < image_buffer_gen_size[1]; i++) {
        if (image_buffer_gen[i]) 
          free(image_buffer_gen[i]);
      }
      if (image_buffer_gen)
        free(image_buffer_gen);
      image_buffer_gen_size[0] = 0;
      image_buffer_gen_size[1] = 0;
    }
    // allocate new space for the image buffer
    image_buffer_gen = (unsigned char**) malloc(KHEIGHT * sizeof(unsigned char*));
    if (image_buffer_gen) { 
      for (int i = 0; i < KHEIGHT; i++) {
        image_buffer_gen[i] = (unsigned char*)malloc(KWIDTH * sizeof(unsigned char));
        if (!image_buffer_gen[i]) {
          fprintf(stderr, "Error allocating memory for report image\n");
          fflush(stderr);
          exit(-1); // give up
        }
      }
    }
    // remember that size for the next free
    image_buffer_gen_size[0] = KWIDTH;
    image_buffer_gen_size[1] = KHEIGHT;

    if (0) {
      // to compute a z-score we can use the publication:
      // W.Limthongkul et al. The Spine Journal 10 (2010) pages 153-158
      // combine male and female, combine all L1 - L5, L1 - L5
      std::vector<float> means = {38.15, 41.48, 44.21, 44.61, 42.52, 25.18, 27.37, 29.54, 30.19, 28.80 };
      std::vector<float> stds = {  9.25,  7.78, 10.14,  9.96, 10.14,  4.31,  4.53,  4.4,   3.07,  2.63 };
      float mean_mean = means[0];
      float mean_stds = stds[0];
      int sum_n = 10; // Algorithm described by Cochrane
      for (int i = 1; i < means.size(); i++) {
        // assume 10 for each group, we don't know how many participants participated in each measurement group
        float tmp_mean_mean = (sum_n * mean_mean + 10 * means[i]) / (sum_n + 10);
        // use the old mean for this computation
        // sqrt(((n1-1)*s1*s1 + (n2-1)*s2*s2 + n1 * n2 / (n1 + n2) * (m1*m1 + m2*m2 - 2 * m1 * m2)) / (n1 + n2 -1));
        mean_stds = sqrt(((sum_n-1)*mean_stds*mean_stds + (10-1)*stds[i]*stds[i] + (sum_n * 10) / 
                        (sum_n + 10) * (mean_mean*mean_mean + means[i]*means[i] - 2 * mean_mean * means[i])) / (sum_n + 10 -1));
        sum_n += 10;
        mean_mean = tmp_mean_mean;
        if (verbose)
          fprintf(stdout, "model used for z-score is: %f %f\n", mean_mean, mean_stds);
      }
      if (verbose) {
        fprintf(stdout, "model used for z-score is: %f %f\n", mean_mean, mean_stds);
        fflush(stdout);
      }
      // a better setting for mean and std is - based on some example runs in BackToBasic
      mean_mean = 23.31783;
      mean_stds = 4.539313;
    }

    // add the physical size to each label
    for (int i = 0; i < report->keyImageTexts.size(); i++) {
      for (int j = 0; j < report->measures.size(); j++) {
        if ( atoi(report->keyImageTexts[i].c_str())-1 == j ) {
          // remember what keyImageText is what label
          // orderOfRegions.insert(std::pair<int, int>(i, j));

          // TODO: measures needs to be organized like the text after a resorting
          float a = atof(report->measures[j].find("physical_size")->second.c_str()) / 1000.0;

          // If we do a log-normal distribution we need to compute the cdf (prob of being smaller than a)
          // and the survival function (compelemnt of the cdf, prob of being larger than a), we should report the smaller of the two.


          // z-score is
          float perc = 0.0f;
          float zscore = 0.0f;
          if (distribution == "norm") {
            zscore = (a - mean_mean) / mean_stds;
            //perc = 100.0f * 0.5f *  (1.0f + (boost::math::erf(zscore / sqrtf(2.0)))); // or, better behaving distribution
            //perc = 100*(/*0.5f * */ boost::math::erfc(- zscore / sqrtf(2.0)));
            perc = 100.0f *  (1.0f + (boost::math::erf(- abs(zscore) / sqrtf(2.0))));
          } else if (distribution == "lognorm") {
            zscore = (log(a) - mean_mean) / mean_stds;
            //zscore = boost::math::pdf(boost::math::lognormal(mean_mean, mean_stds), a);
            perc = 100.0 * fmin(1.0-boost::math::cdf(boost::math::lognormal(mean_mean, mean_stds), a), boost::math::cdf(boost::math::lognormal(mean_mean, mean_stds), a));
          } else {
            fprintf(stderr, "Error: unknown distribution %s specified, only norm and lognorm are supported\n", distribution.c_str());
            return;
          }

          // create the message string
          char str2[256];
          snprintf(str2, 256, "%.02f", zscore);
          //std::stringstream stream2;
          //stream2 << std::fixed << std::setprecision(2) << zscore;
          char str3[256];
          //std::stringstream stream3;
          if (perc < 1)
            snprintf(str3, 256, "<1");
          else {
            if (perc > 99) {
              snprintf(str3, 256, ">99");
            } else {
              snprintf(str3, 256, "%d", (int)perc);
            }
          }
          char str[256];
          snprintf(str, 256, "%.03f", a);
          //std::stringstream stream;
          //stream << std::fixed << std::setprecision(3) << a;
          // meaning of p is: probability of randomly drawing a volume that is further away from the mean than the z-score 
          report->keyImageTexts[i] += std::string(": ") + std::string(str) + std::string(" cm3 ") + std::string("z: ") + std::string(str2) + std::string(" p: ") + std::string(str3) + std::string(" %");
          if (verbose) {
            fprintf(stdout, "  overlay: \"%s\"\n", report->keyImageTexts[i].c_str());
          }
        }
      }
    }
    // debug orderOfRegions
    if (0) {
      //std::map<int, int>::iterator it;
      //for (it = orderOfRegions.begin(); it != orderOfRegions.end(); it++) {
      //  fprintf(stdout, "Key: %d, Value: %d\n", it->first, it->second);
      //  fflush(stdout);
      //}
    }

    // use the generic method to draw text using image_buffer_gen
    // mark the top of the image as "Generated by AI"
    // what is the size of the keyImage?
    int kw = report->keyImage->GetLargestPossibleRegion().GetSize()[0];
    int barHeight = 0.04 * kw;
    int fontSize = 0.012 * kw;
    addBar(report->keyImage, barHeight); // in dark yellow/orange
    addToReportGen(kbuffer, font_file, fontSize, report->TitleText, 10, (barHeight/2)-(fontSize), 0);
    addToReportGen(kbuffer, font_file, fontSize, std::string("[area") + (report->keyImagePositions.size()!=1?std::string("s"):std::string("")) + std::string(" of interest: ") + std::to_string(report->keyImageTexts.size()) + std::string("]"), 5, 15+barHeight, 0);  
    if (report->VersionString.size() > 0)
      addToReportGen(kbuffer, font_file, fontSize, report->VersionString, 10, (barHeight+5)*2, 0);  

    for (int k = 0; k < report->keyImagePositions.size(); k++) {
      //fprintf(stdout, "print %s at %d %d\n", report->keyImageTexts[k].c_str(), report->keyImagePositions[k][0], report->keyImagePositions[k][1]);
      //fflush(stdout);
      std::string::size_type pos = 0;
      if ( ( pos = report->keyImageTexts[k].find("z:") ) != std::string::npos) {
        std::string piece1 = report->keyImageTexts[k].substr(0, pos);
        std::string piece2 = report->keyImageTexts[k].substr(pos);
        //fprintf(stdout, "string: %s <-> %s", piece1.c_str(), piece2.c_str());
        // add a marker for the exact location, need to know how large the character is...
        addMarker(kbuffer, report->keyImagePositions[k][0], report->keyImagePositions[k][1], fontSize);

        addToReportGen(kbuffer, font_file, fontSize, piece1, report->keyImagePositions[k][0]-(5.0/9.0*fontSize) + (20/9.0*fontSize), report->keyImagePositions[k][1]-(30.0/9.0*fontSize), 0);  
        addToReportGen(kbuffer, font_file, fontSize/2, piece2, report->keyImagePositions[k][0]+(35.0/9.0*fontSize) + (20/9.0*fontSize), report->keyImagePositions[k][1]-(50.0/9.0*fontSize), 0);  

        // add TextTopRight
        std::istringstream f(report->TextTopRight);
        int c = 0;
        for (std::string line; std::getline(f, line, '\n'); c++) {
          addToReportGen(kbuffer, font_file, fontSize/2, line, KWIDTH-(fontSize*line.size()), k*(KHEIGHT / report->keyImagePositions.size())+(c * 1.5*fontSize)+(3.0*fontSize), 0);
        }

        // un-changeable text "Not for clinical use"
        addToReportGen(kbuffer, font_file, fontSize/2, "For Research Use Only – Not for use in diagnostic procedures.", 10, (k+1)*(KHEIGHT / report->keyImagePositions.size())-(2.5*fontSize), 0);
      } else {
        addToReportGen(kbuffer, font_file, fontSize/2, report->keyImageTexts[k], report->keyImagePositions[k][0]-5, report->keyImagePositions[k][1]-(30.0/9.0*fontSize), 0);  
      }
    }

    //itk::ImageRegionIterator<CImageType> kIterator(report->keyImage, kregion);
    itk::ImageRegionIteratorWithIndex<CImageType> kIterator(report->keyImage, kregion);
    kIterator.GoToBegin();
    while(!kIterator.IsAtEnd()) {
      CPixelType val = kIterator.Value();
      CImageType::IndexType kidx = kIterator.GetIndex(); 
      //fprintf(stdout, "idx %d %d: %d %d %d\n", kidx[0], kidx[1], val[0], val[1], val[2]);
      //fflush(stdout);
      int ttt = (unsigned char)kbuffer[kidx[1]*KWIDTH+kidx[0]]; // a text color other than white is 255/255,255/237,255/160 (yellowish-orange)

      int c1 = (int)val[0];
      int c2 = (int)val[1];
      int c3 = (int)val[2];

      // now use alpha blending instead of |=

      //      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+0] = (c1 |= ttt);
      //      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+1] = (c2 |= ttt);
      //      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+2] = (c3 |= ttt);
      
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

      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+0] = 255.0 * Cr;
      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+1] = 255.0 * Cg;
      kbuffer_color[3*(kidx[1]*KWIDTH+kidx[0])+2] = 255.0 * Cb;

      ++kIterator;
    }

    gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
    pixeldata.SetByteValue(kbuffer_color, KWIDTH * KHEIGHT * 3); // in bytes

    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;
    im->SetNumberOfDimensions(2);
    im->SetDimension(0, KWIDTH);
    im->SetDimension(1, KHEIGHT);
    im->SetPhotometricInterpretation(gdcm::PhotometricInterpretation::RGB); // change_image.GetPhotometricInterpretation());
    im->GetPixelFormat().SetSamplesPerPixel(3);
    im->SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);


    // gdcm::Image im = change_image;
    gdcm::File *filePtr = new gdcm::File;
    gdcm::Anonymizer anon;
    anon.SetFile(*filePtr);
    anon.Replace(gdcm::Tag(0x0008, 0x0008), "DERIVED\\SECONDARY\\OTHER"); // ImageType
    anon.Replace(gdcm::Tag(0x0028, 0x0002), "3");            // SamplesperPixel
    anon.Replace(gdcm::Tag(0x0028, 0x0004), "RGB");  // PhotometricInterpretation
    anon.Replace(gdcm::Tag(0x0028, 0x0010), std::to_string(KHEIGHT).c_str());         // Rows
    anon.Replace(gdcm::Tag(0x0028, 0x0011), std::to_string(KWIDTH).c_str());          // Columns
    anon.Replace(gdcm::Tag(0x0028, 0x0030), "1\\1"); // PixelSpacing

    anon.Replace(gdcm::Tag(0x0008, 0x0050), report->AccessionNumber.c_str());


    // anon.Replace(gdcm::Tag(0x0020, 0x0010), report->StudyID.c_str()); // we might get the wrong StudyID here from PACS, lets use the StudyInstanceUID value here
    anon.Replace(gdcm::Tag(0x0020, 0x0010), report->StudyInstanceUID.c_str());

    anon.Replace(gdcm::Tag(0x0020, 0x0052), report->FrameOfReferenceUID.c_str());
    anon.Replace(gdcm::Tag(0x0010, 0x0010), report->PatientName.c_str());
    anon.Replace(gdcm::Tag(0x0010, 0x0020), report->PatientID.c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x000d), report->StudyInstanceUID.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0090), report->ReferringPhysician.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0020), report->StudyDate.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0030), report->StudyTime.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0080), report->InstitutionName.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x1030), report->StudyDescription.c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x0011), std::to_string(1000).c_str());

    anon.Replace(gdcm::Tag(0x0020, 0x0013), std::to_string(0).c_str()); // InstanceNumber
    anon.Replace(gdcm::Tag(0x0008, 0x103e), std::string("Biomarker report (research PACS)").c_str());

    boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
    char dateOfReport[9];
    int year = timeLocal.date().year();
    int month = timeLocal.date().month();
    int day = timeLocal.date().day();
    snprintf(dateOfReport, 9, "%04d%02d%02d", year, month, day);
    std::string DateOfSecondaryCapture = std::string(dateOfReport);
    // std::to_string(timeLocal.date().year()) + std::to_string(timeLocal.date().month()) + std::to_string(timeLocal.date().day());
    char timeOfReport[7];
    snprintf(timeOfReport, 7, "%02d%02d%02d", (int)(timeLocal.time_of_day().hours()), (int)(timeLocal.time_of_day().minutes()),
            (int)(timeLocal.time_of_day().seconds()));
    std::string TimeOfSecondaryCapture = std::string(timeOfReport);
    //    std::to_string(timeLocal.time_of_day().hours()) + std::to_string(timeLocal.time_of_day().minutes()) + std::to_string(timeLocal.time_of_day().seconds());

    anon.Replace(gdcm::Tag(0x0018, 0x1012), DateOfSecondaryCapture.c_str());
    anon.Replace(gdcm::Tag(0x0018, 0x1014), TimeOfSecondaryCapture.c_str());
    anon.Replace(gdcm::Tag(0x0018, 0x1016), std::string("pr2mask").c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x4000), std::string("Region of interest shape, intensity and texture measures").c_str());

    gdcm::DataSet &ds = filePtr->GetDataSet(); // ds = reader.GetFile().GetDataSet();
    im->SetDataElement(pixeldata);
    gdcm::Attribute<0x0008, 0x18> ss;
    // adjust the SOPInstanceUID string in case we have more than one report to write
    // int size_num = std::to_string(report->summary.size()).size()+1;
    std::string marker = report->SOPInstanceUID.substr(report->SOPInstanceUID.find_last_of(".") + 1);
    std::string newMarker = marker + std::to_string(report->summary.size());
    std::string newSOPInstanceUID = report->SOPInstanceUID.substr(0,report->SOPInstanceUID.find_last_of("."));
    if (newSOPInstanceUID.size() + newMarker.size() > 62) {
      newSOPInstanceUID = newSOPInstanceUID.substr(0, newSOPInstanceUID.size()-newMarker.size()-1);
    }
    newSOPInstanceUID = newSOPInstanceUID + std::string(".") + newMarker;
    //fprintf(stdout, "%s %s\n", report->SOPInstanceUID.c_str(), newSOPInstanceUID.c_str());
    ss.SetValue(newSOPInstanceUID.c_str()); // TODO: we need a different SOPInstanceUID for each roi
    ds.Replace(ss.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000e> ss2;
    ss2.SetValue(report->SeriesInstanceUID.c_str());
    ds.Replace(ss2.GetAsDataElement());

    gdcm::Attribute<0x0010, 0x0010> ss3;
    ss3.SetValue(report->PatientName.c_str());
    ds.Replace(ss3.GetAsDataElement());

    gdcm::Attribute<0x0010, 0x0020> ss4;
    ss4.SetValue(report->PatientID.c_str());
    ds.Replace(ss4.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000d> ss5;
    ss5.SetValue(report->StudyInstanceUID.c_str());
    ds.Replace(ss5.GetAsDataElement());

    // we would like to store the numeric information in some kind of sequence here as well - step towards a structured report
    if (0) {
      //static const gdcm::Global &g = gdcm::Global::GetInstance();
      //static const gdcm::Dicts &dicts = g.GetDicts();
      //static const gdcm::Dict &pubdict = dicts.GetPublicDict();
      //const gdcm::DictEntry &dictentry = pubdict.GetDictEntry(gdcm::Tag);
      gdcm::DataElement entry;
      entry.SetTag(gdcm::Tag(0x0040, 0x0100)); // CODE VALUE
      entry.SetByteValue("44", (uint32_t)strlen("44"));
      gdcm::DataElement code;
      code.SetTag(gdcm::Tag(0x0008, 0x0104));
      code.SetByteValue("Volume", (uint32_t)strlen("Volume")); // CODE MEANING

      gdcm::SmartPointer<gdcm::SequenceOfItems> sq = new gdcm::SequenceOfItems();
      sq->SetLengthToUndefined();

      for (int measure_idx = 0; measure_idx < report->measures.size(); measure_idx++) {
        auto m = report->measures[measure_idx];
        std::string labelName = std::string("region ") + std::to_string(measure_idx);
        if (strlen(labelName.c_str()) % 2 != 0)
          labelName += std::string(" ");

        for (std::map<std::string, std::string>::iterator iter = m.begin(); iter != m.end(); ++iter) {
          gdcm::Item it;
          it.SetVLToUndefined();
          gdcm::DataSet &nds = it.GetNestedDataSet();

          // resultJSON["measures"].push_back(*iter);
          std::string code_meaning = iter->first;
          std::string code_value = iter->second;
          if (strlen(code_meaning.c_str()) % 2 != 0)
            code_meaning += std::string(" ");
          if (strlen(code_value.c_str()) % 2 != 0)
            code_value += std::string(" ");

          gdcm::DataElement entry;
          entry.SetTag(gdcm::Tag(0x0008, 0x0100)); // CODE VALUE
          entry.SetByteValue(code_value.c_str(), (uint32_t)strlen(code_value.c_str()));
          gdcm::DataElement code;
          code.SetTag(gdcm::Tag(0x0008, 0x0104));
          code.SetByteValue(code_meaning.c_str(), (uint32_t)strlen(code_meaning.c_str())); // CODE MEANING
          gdcm::DataElement region;
          region.SetTag(gdcm::Tag(0x0008, 0x0103));
          region.SetByteValue(labelName.c_str(), (uint32_t)strlen(labelName.c_str()));

          nds.Insert(entry);
          nds.Insert(code);
          nds.Insert(region);
          sq->AddItem(it);
        }
      }

      gdcm::DataElement nn;
      nn.SetTag(gdcm::Tag(0x0040, 0x0a730));
      nn.SetValue(*sq);
      nn.SetVLToUndefined();

      ds.Insert(nn);
    }

    gdcm::ImageWriter writer;
    writer.SetImage(*im);
    writer.SetFile(*filePtr);

    // file names should have a roi counter attached, see if it ends with .dcm and remove it
    std::string out_filename = std::string(report->filename);
    if (out_filename.substr(out_filename.find_last_of(".") + 1) == "dcm") {
      out_filename = out_filename.substr(0, out_filename.find_last_of("."));
    }
    writer.SetFileName((out_filename + std::string("_") + std::to_string(report->summary.size()) + std::string(".dcm")).c_str());
    if (!writer.Write()) {
      return;
    }
    delete[] kbuffer;
    delete[] kbuffer_color;
  } // end of keyImage

  for (int roi = 0; roi < report->summary.size(); roi++) {
//fprintf(stdout, "go over all roi in saveReport...%d\n", roi);
//fflush(stdout);
    //fprintf(stdout, "Start creating report page %d\n", roi+1);
    //fflush(stdout);
    // set the buffer to black (=0)
    memset(&buffer[0], 0, sizeof(char)*len);

    int start_px = 10;
    int start_py = 10;
    int text_lines = report->summary[roi].size();

    float repeat_spacing = 2.0f;
    int xmax = WIDTH;
    int ymax = HEIGHT;

    int px = start_px;
    int py = start_py + (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));
    py = start_py;

    // write one line of text
    for (int line = 0; line < report->summary[roi].size(); line++) {
      memset(image_buffer, 0, sizeof(unsigned char) * HEIGHT * WIDTH);

      // int lengths_min = placements[placement]["lengths"][0];
      // int lengths_max = placements[placement]["lengths"][1];

      int num_chars = report->summary[roi][line].size();

      error = FT_New_Face(library, font_file.c_str(), face_index, &face); /* create face object */

      if (face == NULL) {
        fprintf(stderr, "\033[0mError\033[0m: no face found, provide the filename of a ttf file...\n");
        return;
      }

      int font_size_in_pixel = font_size;
      // fprintf(stdout, "try setting size %d %d\n", font_size_in_pixel, num_chars);
      error = FT_Set_Char_Size(face, font_size_in_pixel * 64, 0, 150, 150); // font_size_in_pixel * 64, 0, 96, 0); /* set character size */
      /* error handling omitted */
      if (error != 0) {
        fprintf(stderr, "\033[0;31mError\033[0m: FT_Set_Char_Size returned error, could not set size %d.\n", font_size_in_pixel);
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
      pen.x = 1 * 64;
      pen.y = (target_height - 20) * 64;

      const char *text = report->summary[roi][line].c_str();
      for (n = 0; n < num_chars; n++) {
        if (text[n] == '\n') {
          continue; // ignore newlines
        }

        /* set transformation */
        FT_Set_Transform(face, &matrix, &pen);

        error = FT_Load_Char(face, text[n], FT_LOAD_RENDER);
        if (error)
          continue; /* ignore errors */

        /* now, draw to our target surface (convert position) draws into image_buffer */
        draw_bitmap(&slot->bitmap, WIDTH, HEIGHT, slot->bitmap_left, target_height - slot->bitmap_top);

        /* increment pen position */
        pen.x += slot->advance.x;
        pen.y += slot->advance.y;
      }

      FT_Done_Face(face);

      float current_image_min_value = 0.0f;
      float current_image_max_value = 255.0f;
      unsigned char *bvals = (unsigned char *)buffer;
      for (int yi = 0; yi < HEIGHT; yi++) {
        for (int xi = 0; xi < WIDTH; xi++) {
          if (image_buffer[yi][xi] == 0)
            continue;
          // I would like to copy the value from image over to
          // the buffer. At some good location...
          int newx = px + xi;
          int newy = py + yi;
          int idx = newy * xmax + newx;
          if (newx < 0 || newx >= xmax || newy < 0 || newy >= ymax)
            continue;

          // instead of blending we need to use a fixed overlay color
          // we have image information between current_image_min_value and current_image_max_value
          // we need to scale the image_buffer by those values.
          float f = 0;
          float v = 1.0f * image_buffer[yi][xi] / 255.0; // 0 to 1 for color, could be inverted if we have a white background
          float w = 1.0f * bvals[idx] / current_image_max_value;
          float alpha_blend = (v + w * (1.0f - v));

          // fprintf(stdout, "%d %d: %d\n", xi, yi, bvals[idx]);
          bvals[idx] = (unsigned char)std::max(
              0.0f, std::min(current_image_max_value, current_image_min_value + (alpha_blend) * (current_image_max_value - current_image_min_value)));
        }
      }
      py += (repeat_spacing * 1.0 * font_size);
    }
//fprintf(stdout, "go over all roi in saveReport...%d\n", roi);
//fflush(stdout);

    // write the key fact a little bit larger on the top right
    if (1) {
      FT_Library library2;

      double angle;
      int target_height;
      int n, num_chars;

      // int font_length = 20;

      std::string font_path = font_file;
      int font_size = 42;
      int face_index = 0;

      FT_Face face;
      FT_GlyphSlot slot;
      FT_Matrix matrix;
      FT_Vector pen;
      FT_Error error;
      // gdcm::ImageReader reader;

      unsigned long len = WIDTH * HEIGHT;
      // char *buffer = new char[len];

      error = FT_Init_FreeType(&library2);

      if (error != 0) {
        fprintf(stderr, "\033[0;31mError\033[0m: The freetype library could not be initialized with this font.\n");
        return;
      }

      int start_px = 10;
      int start_py = 10;
      int text_lines = 1;

      float repeat_spacing = 2.0f;
      int xmax = WIDTH;
      int ymax = HEIGHT;
      // current roi's size is (overwrite the overall size in this loop)
      std::stringstream stream;
      if (report->measures.size() > roi && report->measures[roi].find("physical_size") != report->measures[roi].end()) {
        stream << std::fixed << std::setprecision(2) << (atof(report->measures[roi].find("physical_size")->second.c_str())/1000.0);
      } else {
        stream << "unknown ";
      }
      report->key_fact = stream.str();
      report->key_unit = std::string("cm^3");

      num_chars = report->key_fact.size();

      int px = WIDTH - ((num_chars + 2) * font_size - start_px);
      // int py = start_py + (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));
      int py = start_py + 2.0 * (font_size);

      memset(image_buffer, 0, sizeof(char) * HEIGHT * WIDTH);

      //angle = 0;
      //target_height = HEIGHT;

      error = FT_New_Face(library2, font_file.c_str(), face_index, &face);

      if (face == NULL) {
        fprintf(stderr, "\033[0;31mError\033[0m: no face found, provide the filename of a ttf file...\n");
        return;
      }

      int font_size_in_pixel = 36;
      error = FT_Set_Char_Size(face, font_size_in_pixel * 64, 0, 150, 150); // font_size_in_pixel * 64, 0, 96, 0);
      if (error != 0) {
        fprintf(stderr, "\033[0;31mError\033[0m: FT_Set_Char_Size returned error, could not set size %d.\n", font_size_in_pixel);
        return;
      }

      slot = face->glyph;

      matrix.xx = (FT_Fixed)(cos(angle) * 0x10000L);
      matrix.xy = (FT_Fixed)(-sin(angle) * 0x10000L);
      matrix.yx = (FT_Fixed)(sin(angle) * 0x10000L);
      matrix.yy = (FT_Fixed)(cos(angle) * 0x10000L);

      pen.x = (num_chars * 1 * 64);
      pen.y = (target_height - 80) * 64; // the 60 here is related to the font size!
      const char *text = report->key_fact.c_str();
      for (n = 0; n < num_chars; n++) {
        if (text[n] == '\n') {
          continue; // ignore newlines
        }

        FT_Set_Transform(face, &matrix, &pen);

        error = FT_Load_Char(face, text[n], FT_LOAD_RENDER);
        if (error)
          continue;

        draw_bitmap(&slot->bitmap, WIDTH, HEIGHT, slot->bitmap_left, target_height - slot->bitmap_top);

        pen.x += slot->advance.x;
        pen.y += slot->advance.y;
      }

      FT_Done_Face(face);

      // int px = start_px;
      // int py = start_py + 10.0 * (font_size);
      //  (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));

      float current_image_min_value = 0.0f;
      float current_image_max_value = 255.0f;
      unsigned char *bvals = (unsigned char *)buffer;
      for (int yi = 0; yi < HEIGHT; yi++) {
        for (int xi = 0; xi < WIDTH; xi++) {
          if (image_buffer[yi][xi] == 0)
            continue;
          // I would like to copy the value from image over to
          // the buffer. At some good location...

          int newx = px + xi;
          int newy = py + yi;
          int idx = newy * xmax + newx;
          if (newx < 0 || newx >= xmax || newy < 0 || newy >= ymax)
            continue;
          if (image_buffer[yi][xi] == 0)
            continue;

          // instead of blending we need to use a fixed overlay color
          // we have image information between current_image_min_value and current_image_max_value
          // we need to scale the image_buffer by those values.
          float f = 0;
          float v = 1.0f * image_buffer[yi][xi] / 255.0; // 0 to 1 for color, could be inverted if we have a white background
          float w = 1.0f * bvals[idx] / current_image_max_value;
          float alpha_blend = (v + w * (1.0f - v));

          // fprintf(stdout, "%d %d: %d\n", xi, yi, bvals[idx]);
          bvals[idx] = (unsigned char)std::max(
              0.0f, std::min(current_image_max_value, current_image_min_value + (alpha_blend) * (current_image_max_value - current_image_min_value)));
        }
      }
      // add the units
      //    int px = WIDTH - ((num_chars + 2) * font_size - start_px);
      //    int py = start_py + 2.0 * (font_size);

      // // addToReport(buffer, font_file, 36, report->key_fact, WIDTH - ((num_chars + 2) * font_size - start_px), start_py + 2.0 * (font_size), 0);
      addToReport(buffer, font_file, 26, std::string("cm"), (WIDTH) - ((1.2) * font_size), start_py + 0.5 * (font_size), -3.1415927 / 2.0);
      addToReport(buffer, font_file, 16, std::string("3"), (WIDTH) - ((0.8) * font_size), start_py + 2.0 * (font_size), -3.1415927 / 2.0);
    }

    //
    // add a logo in the lower right corner of the report
    //
    const int bitmapWidth = 13;
    const int bitmapHeight = 10;
    unsigned char bitmap[bitmapWidth * bitmapHeight] = {
          '@', '@', '@', '@',  0 ,  0 ,  0 ,  0 ,  0 , '@', '@', '@', '@',
          '@', '@', '@', '@',  0 ,  0 ,  0 ,  0 ,  0 , '@', '@', '@', '@',
          '@',  0 ,  0 , '@',  0 , '@', '@',  0 ,  0 , '@',  0 ,  0 , '@',
          '@',  0 ,  0 , '@', '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 , '@',
          '@',  0 ,  0 ,  0 , '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 ,  0 ,
          '@',  0 ,  0 ,  0 , '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 ,  0 ,
          '@',  0 ,  0 ,  0 , '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 ,  0 ,
          '@',  0 ,  0 ,  0 , '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 ,  0 ,
          '@',  0 ,  0 ,  0 , '@',  0 ,  0 , '@',  0 , '@',  0 ,  0 ,  0 ,
          '@',  0 ,  0 ,  0 ,  0 , '@', '@',  0 ,  0 , '@',  0 ,  0 ,  0 
    };
    int startX1 = WIDTH - bitmapWidth;  // Calculate the starting X position for the bitmap
    int startY1 = HEIGHT - bitmapHeight;  // Calculate the starting Y position for the bitmap

    // Copy the bitmap into the 2D image buffer
    for (int ytt = 0; ytt < bitmapHeight; ++ytt) {
        for (int xtt = 0; xtt < bitmapWidth; ++xtt) {
            buffer[(startY1 + ytt) * WIDTH + (startX1 + xtt)] = bitmap[ytt * bitmapWidth + xtt];
        }
    }

    // we need some color values here instead of gray-scale
    for (unsigned int c = 0; c < WIDTH * HEIGHT; c++) {
        buffer_color[3*c+0] = buffer[c];
        buffer_color[3*c+1] = buffer[c];
        buffer_color[3*c+2] = buffer[c];
    }

    gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
    pixeldata.SetByteValue(buffer_color, WIDTH * HEIGHT * 3); // in bytes

    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;
    im->SetNumberOfDimensions(2);
    im->SetDimension(0, xmax);
    im->SetDimension(1, ymax);
    im->SetPhotometricInterpretation(gdcm::PhotometricInterpretation::RGB); // change_image.GetPhotometricInterpretation());
    im->GetPixelFormat().SetSamplesPerPixel(3);

//    im->GetPixelFormat().SetBitsAllocated(8); // change_image.GetPixelFormat().GetBitsAllocated());
//    im->GetPixelFormat().SetBitsStored(8);    // change_image.GetPixelFormat().GetBitsStored());
//    im->GetPixelFormat().SetHighBit(7);
//    im->GetPixelFormat().SetPixelRepresentation(gdcm::PixelFormat::UINT8);
    //im->SetSlope(1.0);
    //im->SetIntercept(0);
    im->SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);


    // gdcm::Image im = change_image;
    gdcm::File *filePtr = new gdcm::File;
    gdcm::Anonymizer anon;
    anon.SetFile(*filePtr);
    anon.Replace(gdcm::Tag(0x0008, 0x0008), "DERIVED\\SECONDARY\\OTHER"); // ImageType
    anon.Replace(gdcm::Tag(0x0028, 0x0002), "3");            // SamplesperPixel
    anon.Replace(gdcm::Tag(0x0028, 0x0004), "RGB");  // PhotometricInterpretation
    anon.Replace(gdcm::Tag(0x0028, 0x0010), std::to_string(HEIGHT).c_str());         // Rows
    anon.Replace(gdcm::Tag(0x0028, 0x0011), std::to_string(WIDTH).c_str());          // Columns
    anon.Replace(gdcm::Tag(0x0028, 0x0030), "1\\1"); // PixelSpacing

//    anon.Replace(gdcm::Tag(0x0028, 0x1050), "128"); // WindowCenter
//    anon.Replace(gdcm::Tag(0x0028, 0x1051), "255"); // WindowWidth
//    anon.Replace(gdcm::Tag(0x0028, 0x1052), "0");   // RescaleIntercept
//    anon.Replace(gdcm::Tag(0x0028, 0x1053), "1");   // RescaleSlope
                                                    //  anon.Replace(gdcm::Tag(0x0028, 0x0103), "0");   // use unsigned 0..255

    anon.Replace(gdcm::Tag(0x0008, 0x0050), report->AccessionNumber.c_str());
    // use a new StudyID tag
    //anon.Replace(gdcm::Tag(0x0020, 0x0010), report->StudyID.c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x0010), report->StudyInstanceUID.c_str());

    anon.Replace(gdcm::Tag(0x0020, 0x0052), report->FrameOfReferenceUID.c_str());
    anon.Replace(gdcm::Tag(0x0010, 0x0010), report->PatientName.c_str());
    anon.Replace(gdcm::Tag(0x0010, 0x0020), report->PatientID.c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x000d), report->StudyInstanceUID.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0090), report->ReferringPhysician.c_str());
    //  anon.Replace(gdcm::Tag(0x0008, 0x103e), report->SeriesDescription.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0020), report->StudyDate.c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0030), report->StudyTime.c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x0011), std::to_string(1000).c_str());
    anon.Replace(gdcm::Tag(0x0008, 0x0080), report->InstitutionName.c_str());


    anon.Replace(gdcm::Tag(0x0020, 0x0013), std::to_string(roi+1).c_str()); // InstanceNumber
    anon.Replace(gdcm::Tag(0x0008, 0x103e), std::string("Biomarker report (research PACS)").c_str());

    boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
    char dateOfReport[9];
    int year = timeLocal.date().year();
    int month = timeLocal.date().month();
    int day = timeLocal.date().day();
    snprintf(dateOfReport, 9, "%04d%02d%02d", year, month, day);
    std::string DateOfSecondaryCapture = std::string(dateOfReport);
    // std::to_string(timeLocal.date().year()) + std::to_string(timeLocal.date().month()) + std::to_string(timeLocal.date().day());
    char timeOfReport[7];
    snprintf(timeOfReport, 7, "%02d%02d%02d", (int)(timeLocal.time_of_day().hours()), (int)(timeLocal.time_of_day().minutes()),
            (int)(timeLocal.time_of_day().seconds()));
    std::string TimeOfSecondaryCapture = std::string(timeOfReport);
    //    std::to_string(timeLocal.time_of_day().hours()) + std::to_string(timeLocal.time_of_day().minutes()) + std::to_string(timeLocal.time_of_day().seconds());

    anon.Replace(gdcm::Tag(0x0018, 0x1012), DateOfSecondaryCapture.c_str());
    anon.Replace(gdcm::Tag(0x0018, 0x1014), TimeOfSecondaryCapture.c_str());
    anon.Replace(gdcm::Tag(0x0018, 0x1016), std::string("pr2mask").c_str());
    anon.Replace(gdcm::Tag(0x0020, 0x4000), std::string("Region of interest shape, intensity and texture measures").c_str());

    //im->GetDataElement().SetByteValue(buffer, WIDTH * HEIGHT);
    //im->GetPixelFormat().SetSamplesPerPixel(1);

    gdcm::DataSet &ds = filePtr->GetDataSet(); // ds = reader.GetFile().GetDataSet();
    im->SetDataElement(pixeldata);
    gdcm::Attribute<0x0008, 0x18> ss;
    // adjust the SOPInstanceUID string in case we have more than one report to write
    // int size_num = std::to_string(report->summary.size()).size()+1;
    std::string marker = report->SOPInstanceUID.substr(report->SOPInstanceUID.find_last_of(".") + 1);
    std::string newMarker = marker + std::to_string(roi);
    std::string newSOPInstanceUID = report->SOPInstanceUID.substr(0,report->SOPInstanceUID.find_last_of("."));
    if (newSOPInstanceUID.size() + newMarker.size() > 62) {
      newSOPInstanceUID = newSOPInstanceUID.substr(0, newSOPInstanceUID.size()-newMarker.size()-1);
    }
    newSOPInstanceUID = newSOPInstanceUID + std::string(".") + newMarker;
    //fprintf(stdout, "%s %s\n", report->SOPInstanceUID.c_str(), newSOPInstanceUID.c_str());
    ss.SetValue(newSOPInstanceUID.c_str()); // TODO: we need a different SOPInstanceUID for each roi
    ds.Replace(ss.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000e> ss2;
    ss2.SetValue(report->SeriesInstanceUID.c_str());
    ds.Replace(ss2.GetAsDataElement());

    gdcm::Attribute<0x0010, 0x0010> ss3;
    ss3.SetValue(report->PatientName.c_str());
    ds.Replace(ss3.GetAsDataElement());

    gdcm::Attribute<0x0010, 0x0020> ss4;
    ss4.SetValue(report->PatientID.c_str());
    ds.Replace(ss4.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000d> ss5;
    ss5.SetValue(report->StudyInstanceUID.c_str());
    ds.Replace(ss5.GetAsDataElement());

    gdcm::ImageWriter writer;
    writer.SetImage(*im);
    writer.SetFile(*filePtr);

    // file names should have a roi counter attached, see if it ends with .dcm and remove it
    std::string out_filename = std::string(report->filename);
    if (out_filename.substr(out_filename.find_last_of(".") + 1) == "dcm") {
      out_filename = out_filename.substr(0, out_filename.find_last_of("."));
    }
    writer.SetFileName((out_filename + std::string("_") + std::to_string(roi) + std::string(".dcm")).c_str());
    if (!writer.Write()) {
      return;
    }
  } // loop over rois   
  delete[] buffer;
  delete[] buffer_color;
  FT_Done_FreeType(library);
}

void show_image(void) {
  int i, j;

  for (i = 0; i < HEIGHT; i++) {
    for (j = 0; j < WIDTH; j++)
      putchar(image_buffer[i][j] == 0 ? ' ' : image_buffer[i][j] < 64 ? '.' : (image_buffer[i][j] < 128 ? '+' : '*'));
    putchar('\n');
  }
}

