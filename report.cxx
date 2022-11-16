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
#include <algorithm>
#include <boost/filesystem.hpp>
#include <gdcmFile.h>
#include <gdcmImage.h>
#include <math.h>

#include <ft2build.h>
#include FT_FREETYPE_H

#define WIDTH 2400
#define HEIGHT 2400

unsigned char image_buffer[HEIGHT][WIDTH];

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
  return report;
}

void draw_bitmap(FT_Bitmap *bitmap, FT_Int x, FT_Int y) {
  FT_Int i, j, p, q;
  FT_Int x_max = x + bitmap->width;
  FT_Int y_max = y + bitmap->rows;

  for (i = x, p = 0; i < x_max; i++, p++) {
    for (j = y, q = 0; j < y_max; j++, q++) {
      if (i < 0 || j < 0 || i >= WIDTH || j >= HEIGHT)
        continue;

      image_buffer[j][i] |= bitmap->buffer[q * bitmap->width + p];
    }
  }
}

void saveReport(Report *report) {

  std::string font_file = "Menlo.ttf";
  if (const char *env_p = std::getenv("REPORT_FONT_PATH")) {
    font_file = std::string(env_p);
  }
  if (!boost::filesystem::exists(font_file)) {
    fprintf(stderr, "Error: no font provided. Set the environment variable REPORT_FONT_PATH to the location of a ttf file.\n");
    return;
  }

  FT_Library library;

  double angle;
  int target_height;
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

  unsigned long len = WIDTH * HEIGHT * 8;
  char *buffer = new char[len];

  error = FT_Init_FreeType(&library); /* initialize library */

  if (error != 0) {
    fprintf(stderr, "Error: The freetype library could not be initialized with this font.\n");
    return;
  }

  int start_px = 10;
  int start_py = 10;
  int text_lines = report->summary.size();

  float repeat_spacing = 2.0f;
  int xmax = WIDTH;
  int ymax = HEIGHT;

  int px = start_px;
  int py = start_py + (text_lines * font_size + text_lines * (repeat_spacing * 0.5 * font_size));

  // write one line of text
  for (int line = 0; line < report->summary.size(); line++) {

    memset(image_buffer, 0, HEIGHT * WIDTH);

    // int lengths_min = placements[placement]["lengths"][0];
    // int lengths_max = placements[placement]["lengths"][1];

    int num_chars = report->summary[line].size();
    angle = 0;
    int target_height = HEIGHT;

    error = FT_New_Face(library, font_file.c_str(), face_index, &face); /* create face object */

    if (face == NULL) {
      fprintf(stderr, "Error: no face found, provide the filename of a ttf file...\n");
      return;
    }

    int font_size_in_pixel = font_size;
    // fprintf(stdout, "try setting size %d %d\n", font_size_in_pixel, num_chars);
    error = FT_Set_Char_Size(face, font_size_in_pixel * 64, 0, 150, 150); // font_size_in_pixel * 64, 0, 96, 0); /* set character size */
    /* error handling omitted */
    if (error != 0) {
      fprintf(stderr, "Error: FT_Set_Char_Size returned error, could not set size %d.\n", font_size_in_pixel);
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
    const char *text = report->summary[line].c_str();
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
      draw_bitmap(&slot->bitmap, slot->bitmap_left, target_height - slot->bitmap_top);

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
    py += (repeat_spacing * 1.0 * font_size);
  }

  gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
  pixeldata.SetByteValue(buffer, len);

  gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;
  im->SetNumberOfDimensions(2);
  im->SetDimension(0, xmax);
  im->SetDimension(1, ymax);
  im->SetPhotometricInterpretation(gdcm::PhotometricInterpretation::MONOCHROME2); // change_image.GetPhotometricInterpretation());
  im->GetPixelFormat().SetSamplesPerPixel(1);

  im->GetPixelFormat().SetBitsAllocated(8); // change_image.GetPixelFormat().GetBitsAllocated());
  im->GetPixelFormat().SetBitsStored(8);    // change_image.GetPixelFormat().GetBitsStored());
  im->GetPixelFormat().SetHighBit(7);
  im->GetPixelFormat().SetPixelRepresentation(gdcm::PixelFormat::INT8);
  im->SetSlope(1.0);
  im->SetIntercept(0);

  // gdcm::Image im = change_image;
  gdcm::File *filePtr = new gdcm::File;
  gdcm::Anonymizer anon;
  anon.SetFile(*filePtr);
  anon.Replace(gdcm::Tag(0x0028, 0x0002), "1");            // SamplesperPixel
  anon.Replace(gdcm::Tag(0x0028, 0x0004), "MONOCHROME2");  // PhotometricInterpretation
  anon.Replace(gdcm::Tag(0x0028, 0x0010), std::to_string(HEIGHT).c_str());         // Rows
  anon.Replace(gdcm::Tag(0x0028, 0x0011), std::to_string(WIDTH).c_str());          // Columns
  anon.Replace(gdcm::Tag(0x0028, 0x0030), "1\\1"); // PixelSpacing

  anon.Replace(gdcm::Tag(0x0028, 0x1050), "128"); // WindowCenter
  anon.Replace(gdcm::Tag(0x0028, 0x1051), "255"); // WindowWidth
  anon.Replace(gdcm::Tag(0x0028, 0x1052), "0");   // RescaleIntercept
  anon.Replace(gdcm::Tag(0x0028, 0x1053), "1");   // RescaleSlope

  anon.Replace(gdcm::Tag(0x0008, 0x0050), report->AccessionNumber.c_str());
  anon.Replace(gdcm::Tag(0x0020, 0x0010), report->StudyID.c_str());
  anon.Replace(gdcm::Tag(0x0020, 0x0052), report->FrameOfReferenceUID.c_str());
  anon.Replace(gdcm::Tag(0x0010, 0x0010), report->PatientName.c_str());
  anon.Replace(gdcm::Tag(0x0010, 0x0020), report->PatientID.c_str());
  anon.Replace(gdcm::Tag(0x0020, 0x000d), report->StudyInstanceUID.c_str());
  anon.Replace(gdcm::Tag(0x0008, 0x103e), report->SeriesDescription.c_str());

  im->GetDataElement().SetByteValue(buffer, WIDTH * HEIGHT * sizeof(int8_t));
  im->GetPixelFormat().SetSamplesPerPixel(1);

  gdcm::DataSet &ds = filePtr->GetDataSet(); // ds = reader.GetFile().GetDataSet();
  im->SetDataElement(pixeldata);
  gdcm::Attribute<0x0008, 0x18> ss;
  ss.SetValue(report->SOPInstanceUID.c_str());
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

  writer.SetFileName(report->filename.c_str());
  if (!writer.Write()) {
    return;
  }
  delete[] buffer;
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