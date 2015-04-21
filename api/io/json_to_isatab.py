__author__ = 'alfie'

import os, glob, json, ntpath, string

from api.io.common_functions import CommonFunctions

class JsonToIsatabWriter():
    commonFunctions = CommonFunctions()

    def __init__(self):
        self._investigation_file_pattern = "i_*.json"
        self._study_file_pattern = "s_*.json"
        self._assay_file_pattern = "a_*.json"

        self._isatab_i_ontology_source_ref_sec = ["Term Source Name", "Term Source File",
                                                "Term Source Version", "Term Source Description"]

        self._isatab_i_investigation_sec = ["Investigation Identifier", "Investigation Title",
                                            "Investigation Description", "Investigation Submission Date",
                                            "Investigation Public Release Date", "Comment [Created With Configuration]",
                                            "Comment [Last Opened With Configuration]"]

        self._isatab_i_investigation_publications_sec = ["Investigation PubMed ID", "Investigation Publication DOI",
                                                        "Investigation Publication Author List", "Investigation Publication Title",
                                                        "Investigation Publication Status", "Investigation Publication Status Term Accession Number",
                                                        "Investigation Publication Status Term Source REF"]

        self._isatab_i_investigation_contacts_sec = ["Investigation Person Last Name", "Investigation Person First Name",
                                                    "Investigation Person Mid Initials", "Investigation Person Email",
                                                    "Investigation Person Phone", "Investigation Person Fax",
                                                    "Investigation Person Address", "Investigation Person Affiliation",
                                                    "Investigation Person Roles", "Investigation Person Roles Term Accession Number",
                                                    "Investigation Person Roles Term Source REF", "Comment [Investigation Person REF]"]

        self._isatab_i_study_sec = ["Study Identifier", "Study Title",
                                    "Study Description", "Comment [Study Grant Number]",
                                    "Comment [Study Funding Agency]", "Study Submission Date",
                                    "Study Public Release Date", "Study File Name"]

        self._isatab_i_study_design_descriptors_sec = ["Study Design Type", "Study Design Type Term Accession Number",
                                                    "Study Design Type Term Source REF"]

        self._isatab_i_study_publications_sec = ["Study PubMed ID", "Study Publication DOI",
                                                "Study Publication Author List", "Study Publication Title",
                                                "Study Publication Status", "Study Publication Status Term Accession Number",
                                                "Study Publication Status Term Source REF"]

        self._isatab_i_study_factors_sec = ["Study Factor Name", "Study Factor Type",
                                            "Study Factor Type Term Accession Number", "Study Factor Type Term Source REF"]

        self._isatab_i_study_assays_sec = ["Study Assay Measurement Type", "Study Assay Measurement Type Term Accession Number",
                                        "Study Assay Measurement Type Term Source REF", "Study Assay Technology Type",
                                        "Study Assay Technology Type Term Accession Number", "Study Assay Technology Type Term Source REF",
                                        "Study Assay Technology Platform", "Study Assay File Name"]

        self._isatab_i_study_protocols_sec = ["Study Protocol Name", "Study Protocol Type",
                                            "Study Protocol Type Term Accession Number", "Study Protocol Type Term Source REF",
                                            "Study Protocol Description", "Study Protocol URI",
                                            "Study Protocol Version", "Study Protocol Parameters Name",
                                            "Study Protocol Parameters Name Term Accession Number", "Study Protocol Parameters Name Term Source REF",
                                            "Study Protocol Components Name", "Study Protocol Components Type",
                                            "Study Protocol Components Type Term Accession Number", "Study Protocol Components Type Term Source REF"]

        self._isatab_i_study_contacts_sec = ["Study Person Last Name", "Study Person First Name",
                                            "Study Person Mid Initials", "Study Person Email",
                                            "Study Person Phone", "Study Person Fax",
                                            "Study Person Address", "Study Person Affiliation",
                                            "Study Person Roles", "Study Person Roles Term Accession Number",
                                            "Study Person Roles Term Source REF", "Comment [Study Person REF]"]

    def parsingJson(self, json_dir, output_dir):
        if os.path.isdir(json_dir):
            i_filenames = glob.glob(os.path.join(json_dir, self._investigation_file_pattern))
            s_filenames = glob.glob(os.path.join(json_dir, self._study_file_pattern))
            a_filenames = glob.glob(os.path.join(json_dir, self._assay_file_pattern))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            self.writeJsonInvestigationToIsatab(i_filenames, output_dir)
            self.writeJsonStudyAssayToIsatab(s_filenames, output_dir, "studySampleTable", "studyTableHeaders", "studyTableData")
            self.writeJsonStudyAssayToIsatab(a_filenames, output_dir, "assayTable", "assayTableHeaders", "assayTableData")

    def writeJsonInvestigationToIsatab(self, filenames, output_dir):
        my_str = ""

        assert len(filenames) == 1
        json_investigation = filenames[0]
        assert os.path.exists(json_investigation), "Did not find investigation file: %s" % json_investigation

        with open(json_investigation) as in_handle:
            jsonData = json.load(in_handle)
            # ONTOLOGY SOURCE REFERENCE
            my_str = self.writeSectionInvestigation(my_str, "ONTOLOGY SOURCE REFERENCE", jsonData["ontologySourceReference"], self._isatab_i_ontology_source_ref_sec)
            # INVESTIGATION
            my_str = my_str + "INVESTIGATION" + "\n"
            for i in self._isatab_i_investigation_sec:
                my_str = my_str + i + "\t\"" + jsonData["investigation"][self.commonFunctions.makeAttributeName(i)] + "\"" + "\n"
            # INVESTIGATION PUBLICATIONS
            my_str = self.writeSectionInvestigation(my_str, "INVESTIGATION PUBLICATIONS", jsonData["investigationPublications"], self._isatab_i_investigation_publications_sec)
            # INVESTIGATION CONTACTS
            my_str = self.writeSectionInvestigation(my_str, "INVESTIGATION CONTACTS", jsonData["investigationContacts"], self._isatab_i_investigation_contacts_sec)
            for study in jsonData["studies"]:
                # STUDY
                my_str = my_str + "STUDY" + "\n"
                for i in self._isatab_i_study_sec:
                    my_str = my_str + i + "\t\"" + study["study"][self.commonFunctions.makeAttributeName(i)] + "\"" + "\n"
                # STUDY DESIGN DESCRIPTORS
                my_str = self.writeSectionInvestigation(my_str, "STUDY DESIGN DESCRIPTORS", study["studyDesignDescriptors"], self._isatab_i_study_design_descriptors_sec)
                # STUDY PUBLICATIONS
                my_str = self.writeSectionInvestigation(my_str, "STUDY PUBLICATIONS", study["studyPublications"], self._isatab_i_study_publications_sec)
                # STUDY FACTORS
                my_str = self.writeSectionInvestigation(my_str, "STUDY FACTORS", study["studyFactors"], self._isatab_i_study_factors_sec)
                # STUDY ASSAYS
                my_str = self.writeSectionInvestigation(my_str, "STUDY ASSAYS", study["assays"], self._isatab_i_study_assays_sec)
                # STUDY PROTOCOLS
                my_str = self.writeSectionInvestigation(my_str, "STUDY PROTOCOLS", study["studyProtocols"], self._isatab_i_study_protocols_sec)
                # STUDY CONTACTS
                my_str = self.writeSectionInvestigation(my_str, "STUDY CONTACTS", study["studyContacts"], self._isatab_i_study_contacts_sec)
        # now we write out each of the study files
        with open(os.path.join(output_dir, ntpath.basename(str(json_investigation)).split(".")[0] + ".txt"), "w") as file_isatab:
            file_isatab.write(my_str)

    def writeSectionInvestigation(self, my_str, sec_header, study, sec_header_group):
        my_str = my_str + sec_header + "\n"
        for i in sec_header_group:
            my_str = my_str + i + "\t"
            for b in study:
                my_str = my_str + "\"" + b[self.commonFunctions.makeAttributeName(i)] + "\"" + "\t"
            my_str = my_str + "\n"
        return my_str

    def writeJsonStudyAssayToIsatab(self, filenames, output_dir, tableNameTitle, tableHeaderTitle, tableDataTitle):
        my_str = ""
        assert len(filenames) > 0
        for each_file in filenames:
            assert os.path.exists(each_file), "Did not find study / assay file: %s" % each_file
            with open(each_file) as in_handle:
                json_each_s = json.load(in_handle)
                # the study headers
                for node_h in json_each_s[tableNameTitle][tableHeaderTitle]:
                    my_str = my_str + "\"" + node_h["name"] + "\"" + "\t"
                    if "attributes" in node_h:
                        for n_h in node_h["attributes"]:
                            my_str = my_str + "\"" + n_h["name"] + "\"" + "\t"
                my_str = my_str + "\n"
                # now for each of the rows
                for node_r in json_each_s[tableNameTitle][tableDataTitle]:
                    for n_r in node_r:
                        my_str = my_str + "\"" + n_r + "\"" + "\t"
                    my_str = my_str + "\n"
                # now we write out each of the study files
                with open(os.path.join(output_dir, ntpath.basename(str(each_file)).split(".")[0] + ".txt"), "w") as file_isatab:
                    file_isatab.write(my_str)

########

mywriter = JsonToIsatabWriter()
folder_name = "BII-I-1"
json_dir = os.path.join("../../tests/data", folder_name + "-json")
output_dir = os.path.join("../../tests/data", folder_name + "-generatedIsatab")
mywriter.parsingJson(json_dir, output_dir)
