#ifndef REPORT_H
#define REPORT_H

#include <string>
#include <vector>
// a simple report builder that exports a DICOM file
// provide the information required to setup the report

typedef struct {
  std::string filename;
  std::string SOPInstanceUID; // 64chars max
  std::string SeriesInstanceUID;
  std::string StudyInstanceUID;
  std::string FrameOfReferenceUID;
  std::string StudyID; // 16 chars max
  std::string AccessionNumber;
  std::vector<std::string> summary;
} Report;

Report *getDefaultReportStruct();
void saveReport(Report *report);

#endif