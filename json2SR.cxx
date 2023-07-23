#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */
#include "dcmtk/ofstd/ofconsol.h"     /* for COUT */

#include "dcmtk/dcmsr/dsrdoc.h"       /* for main interface class DSRDocument */
#include "dcmtk/dcmsr/codes/ucum.h"

#include "dcmtk/dcmdata/dctk.h"       /* for typical set of "dcmdata" headers */

#include <boost/date_time.hpp>
#include <boost/program_options.hpp>
#include "json.hpp"
#include <string>

// forward declarations
static void generate(DSRDocument *doc, OFString &studyUID_01, nlohmann::json &report);

int main(int argc, char *argv[]) {
    setlocale(LC_NUMERIC, "en_US.utf-8");

    boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();
    std::string versionString = std::string("0.0.1.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");

    namespace po = boost::program_options;
    po::options_description desc("Convert JSON created by imageAndMask2Report into DICOM-SR.\nUsage: <program> [options] <json-input-file>\n\nAllowed options");

    desc.add_options()
        ("version,v", "print the version string")
        ("help,h", "produce help message")
        ("json-input-file,j", po::value< std::vector<std::string> >(), "json-input-file")
    ;

    po::positional_options_description p;
    p.add("json-input-file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
          options(desc).positional(p).run(), vm);
    po::notify(vm);

    //po::store(po::parse_command_line(argc, argv, desc), vm);
    //po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    if (vm.count("version")) {
        std::cout << versionString.c_str() << "\n";
        return 1;
    }

    if (vm.count("json-input-file")) {
        std::vector<std::string> jsons = vm["json-input-file"].as< std::vector< std::string> >();
        for (int i = 0; i < jsons.size(); i++) {
           std::cout << "Reading input json " << jsons[i] << ".\n";
        }
    } else {
        std::cout << "No input json specified.\n";
        std::cout << desc << "\n";
        exit(-1);
    }

    std::vector<std::string> input = vm["json-input-file"].as<std::vector< std::string> >();

    // read in the config file as json
    for (int i = 0; i < input.size(); i++) {
        std::ifstream f(input[i]);
        nlohmann::json report = nlohmann::json::parse(f);

        std::string outputFilename = input[i];
        outputFilename = outputFilename.substr(0,outputFilename.find_last_of('.'))+".dcm";
 
        /* make sure data dictionary is loaded */
        if (!dcmDataDict.isDictionaryLoaded()) {
            CERR << "Warning: no data dictionary loaded, "
                << "check environment variable: "
                << DCM_DICT_ENVIRONMENT_VARIABLE << OFendl;
        }

        DSRDocument *doc = new DSRDocument();
        OFString studyUID_01;
        std::string StudyInstanceUID("");
        nlohmann::json meta;
        if (report.contains("meta-data")) {
            meta = report["meta-data"];
            if (meta.is_array() && meta.size() > 0) {
                for (int j = 0; j < meta.size(); j++) {
                    if (meta[j].contains("StudyInstanceUID")) {
                        StudyInstanceUID = meta[j]["StudyInstanceUID"];
                        StudyInstanceUID.erase(remove( StudyInstanceUID.begin(), StudyInstanceUID.end(), '\"' ),StudyInstanceUID.end());
                        studyUID_01 = OFString(StudyInstanceUID.c_str());
                        //fprintf(stderr, "Found study instance uid %s\n", value.str().c_str());
                        break;
                    }
                }
            }
        }

        generate(doc, studyUID_01, report);

        //COUT << OFString(79, '-') << OFendl;
        //COUT << "mkreport: " << outputFilename << OFendl << OFendl;
        //doc->print(COUT, DSRTypes::PF_shortenLongItemValues);
        //COUT << OFendl;

        DcmFileFormat *fileformat = new DcmFileFormat();
        DcmDataset *dataset = NULL;
        if (fileformat != NULL)
            dataset = fileformat->getDataset();
        if (dataset != NULL) {
                doc->getCodingSchemeIdentification().addPrivateDcmtkCodingScheme();
                if (doc->write(*dataset).good()) {
                    // now set the meta data for this document (make part of our study)
                    dataset->putAndInsertString(DCM_StudyInstanceUID, StudyInstanceUID.c_str());
                    std::string derivedSeriesInstanceUID = StudyInstanceUID;
                    std::string endString = ".5";
                    if (derivedSeriesInstanceUID.substr(derivedSeriesInstanceUID.size() - 2, 2) == ".5")
                        endString = ".6";
                    // change it so that we end up with a new series instance uid - always in the same way, always at most 64 characters in length
                    derivedSeriesInstanceUID = derivedSeriesInstanceUID.substr(0, 64 - 2) + endString;
                    dataset->putAndInsertString(DCM_SeriesInstanceUID, derivedSeriesInstanceUID.c_str());
                    if (meta.size() > 0 && meta[0].contains("PatientName")) {
                        std::string PatientName = meta[0]["PatientName"];
                        dataset->putAndInsertString(DCM_PatientName, PatientName.c_str());
                    }
                    if (meta.size() > 0 && meta[0].contains("PatientID")) {
                        std::string PatientID = meta[0]["PatientID"];
                        dataset->putAndInsertString(DCM_PatientID, PatientID.c_str());
                    }

                    std::cout << "Write: " << outputFilename << "..." << OFendl;
                    OFString filename(outputFilename.c_str());
                    if (fileformat->saveFile(filename, EXS_LittleEndianExplicit).bad())
                        CERR << "ERROR: could not save dataset to file '" << filename << "'" << OFendl;
                    else {
                        std::cout << "done." << OFendl;
                    }
                } else
                    CERR << "ERROR: could not write SR document into dataset" << OFendl;
        }
        delete doc;
    }    
    return 0;
}


static void generate(DSRDocument *doc, OFString &studyUID_01, nlohmann::json &report) {
    nlohmann::json meta;
    if (report.contains("meta-data")) {
        meta = report["meta-data"];
    }
    if (meta.is_array() && meta.size() > 0) {
        meta = meta[0];
    } else {
        return; // no meta
    }
    nlohmann::json measures;
    if (report.contains("measures")) {
        measures = report["measures"];
    }
    if (measures.is_array() && measures.size() > 0) {
        measures = measures[0];
    } else {
        return; // no measures
    }

    std::string ReferringPhysician("");
    if (meta.contains("ReferringPhysician"))
        ReferringPhysician = meta["ReferringPhysician"];
    std::string PatientName("");
    if (meta.contains("PatientName"))
        PatientName = meta["PatientName"];
    std::string PatientID("");
    if (meta.contains("PatientID"))
        PatientID = meta["PatientID"];
    std::string ReportSOPInstanceUID("");
    if (meta.contains("ReportSOPInstanceUID"))
        ReportSOPInstanceUID = meta["ReportSOPInstanceUID"];
    std::string StudyInstanceUID("");
    if (meta.contains("StudyInstanceUID"))
        StudyInstanceUID = meta["StudyInstanceUID"];
    std::string ReportSeriesInstanceUID;
    if (meta.contains("SeriesInstanceUID"))
        ReportSeriesInstanceUID = meta["SeriesInstanceUID"];

    doc->createNewDocument(DSRTypes::DT_BasicTextSR);
    if (!studyUID_01.empty()) {
        //fprintf(stderr, "study instance uid is not empty\n");
        // we will overwrite the series instance uid later, ignore for now... 
        doc->createNewSeriesInStudy(studyUID_01);
    }
    doc->setStudyDescription("ROR Structured report");
    doc->setSeriesDescription("ROR SR representation for report");

    doc->setPatientName(PatientName.c_str());
    doc->setPatientID(PatientID.c_str());
    doc->setPatientBirthDate("00000000");
    doc->setPatientSex("O");
    doc->setReferringPhysicianName(ReferringPhysician.c_str());

    doc->getTree().addContentItem(DSRTypes::RT_isRoot, DSRTypes::VT_Container);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("DT.01", OFFIS_CODING_SCHEME_DESIGNATOR, "MMIV Report"));

    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_PName, DSRTypes::AM_belowCurrent);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.04", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Name"));
    doc->getTree().getCurrentContentItem().setStringValue(ReferringPhysician.c_str());
    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.05", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Organization Name"));
    doc->getTree().getCurrentContentItem().setStringValue("Mohn Medical Imaging and Visualization Centre, Haukeland University Hospital, Bergen, Norway");

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Image);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_27", OFFIS_CODING_SCHEME_DESIGNATOR, "Referenced image series"));
    // doc->getTree().getCurrentContentItem().setImageReference(DSRImageReferenceValue(UID_SecondaryCaptureImageStorage, "1.2.276.0.7230010.3.1.4.123456.1.1a", OFFalse));
    doc->getTree().getCurrentContentItem().setImageReference(DSRImageReferenceValue(UID_SecondaryCaptureImageStorage, ReportSOPInstanceUID.c_str(), OFFalse)); // specify SOPInstanceUID of report
    doc->getCurrentRequestedProcedureEvidence().addItem(StudyInstanceUID.c_str(), ReportSeriesInstanceUID.c_str(), UID_SecondaryCaptureImageStorage, ReportSOPInstanceUID.c_str(), OFFalse );

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_28", OFFIS_CODING_SCHEME_DESIGNATOR, "AI Assessment"));
    doc->getTree().getCurrentContentItem().setStringValue("This report image is referenced by this SR.");

    // measures
    // walk through all the measures we find in the report
    for (auto& el : measures.items()) {
        std::string key = el.key();
        std::string value = el.value();
        doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
        doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_01", OFFIS_CODING_SCHEME_DESIGNATOR, "Description"));
        doc->getTree().getCurrentContentItem().setStringValue((std::string("Measure: ") + key).c_str());

        // doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Num);
        doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
        doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("112187", "SNM3", "Measure"));
        // doc->getTree().getCurrentContentItem().setNumericValue(DSRNumericMeasurementValue(value.c_str(), CODE_UCUM_Pixels));
        doc->getTree().getCurrentContentItem().setStringValue(value.c_str());
    }

    doc->getCodingSchemeIdentification().addItem(CODE_UCUM_CodingSchemeDesignator);
    doc->getCodingSchemeIdentification().setCodingSchemeUID(CODE_UCUM_CodingSchemeUID);
    doc->getCodingSchemeIdentification().setCodingSchemeName(CODE_UCUM_CodingSchemeDescription);
    doc->getCodingSchemeIdentification().setCodingSchemeResponsibleOrganization("Mohn Medical Imaging and Visualization Centre, Bergen, Norway");

/*

    //doc->createNewDocument(DSRTypes::DT_EnhancedSR);
    doc->createNewDocument(DSRTypes::DT_BasicTextSR);
    if (!studyUID_01.empty()) {
        //fprintf(stderr, "study instance uid is not empty\n");
        // we will overwrite the series instance uid later, ignore for now... 
        doc->createNewSeriesInStudy(studyUID_01);
    }
    doc->setStudyDescription("ROR Structured report");
    doc->setSeriesDescription("ROR SR representation for report");

    doc->setPatientName(PatientName.c_str());
    doc->setPatientID(PatientID.c_str());
    doc->setPatientBirthDate("00000000");
    doc->setPatientSex("O");
    doc->setReferringPhysicianName(ReferringPhysician.c_str());

    doc->getTree().addContentItem(DSRTypes::RT_isRoot, DSRTypes::VT_Container);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("DT.01", OFFIS_CODING_SCHEME_DESIGNATOR, "MMIV Report"));
    //doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_Text);

    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_PName, DSRTypes::AM_belowCurrent);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.04", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Name"));
    doc->getTree().getCurrentContentItem().setStringValue(ReferringPhysician.c_str());
    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.05", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Organization Name"));
    doc->getTree().getCurrentContentItem().setStringValue("Mohn Medical Imaging and Visualization Centre, Bergen, Norway");

    // reference the image (study series sop)
    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Image);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_27", OFFIS_CODING_SCHEME_DESIGNATOR, "Referenced image series"));
    doc->getTree().getCurrentContentItem().setImageReference(DSRImageReferenceValue(UID_SecondaryCaptureImageStorage, ReportSOPInstanceUID.c_str())); // specify SOPInstanceUID of report
    doc->getCurrentRequestedProcedureEvidence().addItem(StudyInstanceUID.c_str(), ReportSeriesInstanceUID.c_str(), UID_SecondaryCaptureImageStorage, ReportSOPInstanceUID.c_str());

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_28", OFFIS_CODING_SCHEME_DESIGNATOR, "AI Assessment"));
    doc->getTree().getCurrentContentItem().setStringValue("This report image is referenced by this SR.");

    doc->getTree().goUp();
    fprintf(stderr, "writing %s %s %s %s\n", StudyInstanceUID.c_str(), ReportSeriesInstanceUID.c_str(), ReportSOPInstanceUID.c_str(), UID_SecondaryCaptureImageStorage);

    //doc->getTree().addContentItem(DSRTypes::RT_isRoot, DSRTypes::VT_Container);
    //doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("DT.06", OFFIS_CODING_SCHEME_DESIGNATOR, "Consultation Report"));

    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_PName, DSRTypes::AM_belowCurrent);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.04", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Name"));
    doc->getTree().getCurrentContentItem().setStringValue("FIONA Image AI done by ror");
    doc->getTree().addContentItem(DSRTypes::RT_hasObsContext, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("IHE.05", OFFIS_CODING_SCHEME_DESIGNATOR, "Observer Organization Name"));
    doc->getTree().getCurrentContentItem().setStringValue("Mohn Medical Imaging and Visualization Centre, Bergen, Norway");

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Container);
    doc->getTree().getCurrentContentItem().setContinuityOfContent(DSRTypes::COC_Continuous);

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text, DSRTypes::AM_belowCurrent);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_01", OFFIS_CODING_SCHEME_DESIGNATOR, "Description"));
    doc->getTree().getCurrentContentItem().setStringValue("This report has been created by");

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_PName);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_04", OFFIS_CODING_SCHEME_DESIGNATOR, "Referring Physician"));
    doc->getTree().getCurrentContentItem().setStringValue("FIONA Image AI");

    // measures
    // walk through all the measures we find in the report
    for (auto& el : measures.items()) {
        std::string key = el.key();
        std::string value = el.value();
        doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
        doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_01", OFFIS_CODING_SCHEME_DESIGNATOR, "Description"));
        doc->getTree().getCurrentContentItem().setStringValue((std::string("Measure: ") + key).c_str());

        doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Num);
        doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("112187", "SNM3", "Measure"));
        doc->getTree().getCurrentContentItem().setNumericValue(DSRNumericMeasurementValue(value.c_str(), CODE_UCUM_Pixels));
    }

    doc->getTree().goUp();

    doc->getTree().addContentItem(DSRTypes::RT_contains, DSRTypes::VT_Text);
    doc->getTree().getCurrentContentItem().setConceptName(DSRCodedEntryValue("CODE_02", OFFIS_CODING_SCHEME_DESIGNATOR, "Diagnosis"));
    doc->getTree().getCurrentContentItem().setStringValue("Left to the user.");

    // add additional information on UCUM coding scheme
    doc->getCodingSchemeIdentification().addItem(CODE_UCUM_CodingSchemeDesignator);
    doc->getCodingSchemeIdentification().setCodingSchemeUID(CODE_UCUM_CodingSchemeUID);
    doc->getCodingSchemeIdentification().setCodingSchemeName(CODE_UCUM_CodingSchemeDescription);
    doc->getCodingSchemeIdentification().setCodingSchemeResponsibleOrganization("Mohn Medical Imaging and Visualization Centre, Bergen, Norway");

    */
}


