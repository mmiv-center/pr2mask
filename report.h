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
  std::string key_fact;
  std::string key_unit;
  std::string InstitutionName;
  CImageType::Pointer keyImage;
  std::vector<std::array<int, 2> > keyImagePositions;
  std::vector<std::string> keyImageTexts;
  std::vector<std::vector<std::string> > summary;
  std::vector<std::map<std::string, std::string> > measures;
} Report;

Report *getDefaultReportStruct();
void saveReport(Report *report);

#endif