#ifndef REPORT_H
#define REPORT_H

#include <string>
#include <vector>
#include <map>
// a simple report builder that exports a DICOM file
// provide the information required to setup the report

typedef struct {
  std::string filename;
  std::string PatientName;
  std::string PatientID;
  std::string SOPInstanceUID; // 64chars max
  std::string SeriesInstanceUID;
  std::string StudyInstanceUID;
  std::string FrameOfReferenceUID;
  std::string ReferringPhysician;
  std::string StudyID; // 16 chars max
  std::string AccessionNumber;
  std::string SeriesDescription;
  std::string StudyDate;
  std::string StudyTime;
  std::string key_fact;
  std::string key_unit;
  std::vector<std::string> summary;
  std::vector<std::map<std::string, std::string> > measures;
} Report;

Report *getDefaultReportStruct();
void saveReport(Report *report);

#endif