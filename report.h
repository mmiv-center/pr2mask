#ifndef REPORT_H
#define REPORT_H

#include "itkMetaDataObject.h"
#include "itkRGBPixel.h"
#include "itkImageAdaptor.h"
#include <string>
#include <vector>
#include <map>
// a simple report builder that exports a DICOM file
// provide the information required to setup the report

  using CPixelType = itk::RGBPixel<unsigned char>;
  using CImageType = itk::Image<CPixelType, 2>;

typedef struct {
  std::string filename;
  std::string VersionString; // provided by the user on the command line, version of the AI
  std::string pr2maskVersionString; // the version of pr2mask used to generate the report images
  std::string ContainerVersionString; // the version of the docker container used on DTS (CONTAINER_VERSION env variable)
  std::string PatientName;
  std::string PatientID;
  std::string SOPInstanceUID; // 64chars max
  std::string SeriesInstanceUID;
  std::string StudyInstanceUID;
  std::string StudyDescription;
  std::string FrameOfReferenceUID;
  std::string ReferringPhysician;
  std::string StudyID; // 16 chars max
  std::string AccessionNumber;
  std::string SeriesDescription;
  std::string ReportSeriesInstanceUID;
  std::string ReportSOPInstanceUID;
  std::string StudyDate;
  std::string StudyTime;
  std::string SeriesDate;
  std::string SeriesTime;
  std::string key_fact;
  std::string key_unit;
  std::string InstitutionName;
  std::string TitleText;
  std::string TextTopRight; // try to make this multi-line
  std::map<int, std::string> TextTopRightLabels; // try to make this multi-line with labels
  float BrightnessContrastLL;
  float BrightnessContrastUL;
  std::string ReportType;
  CImageType::Pointer keyImage;
  std::vector<std::array<int, 2> > keyImagePositions;
  std::vector<std::string> keyImageTexts;
  std::vector<std::vector<std::string> > summary;
  std::vector<std::map<std::string, std::string> > measures;
} Report;

Report *getDefaultReportStruct();
void saveReport(Report *report, std::string distribution="norm", float mean_mean = 23.31783, float mean_stds = 4.539313, bool verbose = false);


#endif