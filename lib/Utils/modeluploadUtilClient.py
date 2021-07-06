import os
import logging
import re
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class modeluploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_model(self):
        """
        The upload method (Incomplete documentation)
        Get the output name for the model
        params should include:
            username,
            staging_file_names,
            genes_table_ref,
            description,
            output_names
        """

        print("params: ", self.params)
        self.validate_import_model_from_staging_params()

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception("Expecting a single output name, got a different number"
                            f": {len(op_nms)}. Output Names: " + \
                            ", ".join(op_nms))
        else:
            model_name = op_nms[0]

        model_str, past_end_str = self.get_model_and_pastEnd_strs(self.params["standard_model_name"])

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()
        #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
        #        self.params['username'], str(date_time))

        # We create the data for the object
        model_data = {
            "file_type": "KBaseRBTnSeq.RBTS_Model",
            "utc_created": str(date_time),
            "standard_model_name": self.params["standard_model_name"],
            "model_string": model_str,
            "past_end_string": past_end_str,
            "description": "Manual Upload: " + self.params["description"]
        }


        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBaseRBTnSeq.RBTS_Model",
                    "data": model_data,
                    "name": model_name,
                }
            ],
        }

        logging.info("Using DFU to save the object info:")
        # save_objects returns a list of object_infos
        dfu_object_info = self.dfu.save_objects(save_object_params)[0]

        return {
            "Name": dfu_object_info[1],
            "Type": dfu_object_info[2],
            "Date": dfu_object_info[3],
        }

    def validate_import_model_from_staging_params(self):

        # check for required parameters
        for p in [
            "username",
            "description",
            "output_names",
            "standard_model_name"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        if self.params["standard_model_name"] == "Custom":
            for p in ["model_str", "past_end_str"]:
                if p not in self.params:
                    raise ValueError(f"'{p}' parameter is required when creating a custom model.")
            if self.params["model_str"] == "None" or self.params["model_str"] == "":
                raise ValueError(f" input model_str cannot be '' or 'None' when creating a custom model.")

    def get_model_and_pastEnd_strs(self, standard_model_name):
        """
        Description:
            In this function we get the two parts of the model-
            The model string, which is the part of the transposon in which the barcode sits.
            And the past end string, which is after the transposon.
        Returns:
            
        """

        if standard_model_name == "Custom":
            model_str = self.params["model_str"]
            past_end_str = self.params["past_end_str"]
            self.check_model_parts(model_str, past_end_str)
        else:
            if standard_model_name ==  "Sc_Tn5":
                model_str = "nnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCAGCTATGACATGAAGATGTGTATAAGAGACAG" 
                past_end_str = "GGAAGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATTGGGCCCAGGTACCAACTACGTCAGGTGGCACTTT"
            elif standard_model_name == "ezTn5_Tet_Bifido":
                model_str = "nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCACCTCGACAGATGTGTATAAGAGACAG" 
                past_end_str = ""
            elif standard_model_name == "ezTn5_kan1":
                model_str = "nnnnnnCTAAGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACAGATGTGTATAAGAGACAG"
                past_end_str = ""
            elif standard_model_name == "ezTn5_kanU":
                model_str = "nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACAGATGTGTATAAGAGACAG"
                past_end_str = ""
            elif standard_model_name == "magic_Tn5":
                model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCGGTTGAGATGTGTATAAGAGACAG'
                past_end_str = 'TCGACGGCTTGGTTTCATAAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGA'
            elif standard_model_name == "magic_mariner":
                model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
                past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGCTTATCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGAGCACCGCCAGGTGCGAATAAGGGACAGTGAAGAAG'
            elif standard_model_name == "magic_mariner.2":
                model_str = 'nnnnnnnnnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
                past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGCTTATCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGAGCACCGCCAGGTGCGAATAAGGGACAGTGAAGAAG'
            elif standard_model_name == "pHIMAR_kan":
                model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
                past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA' 
            elif standard_model_name == "pKMW3":
                model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
                past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA'
            elif standard_model_name == "pKMW3_universal":
                model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
                past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA'
            elif standard_model_name == "pKMW7":
                model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG' 
                past_end_str = 'TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC'
            elif standard_model_name == "pKMW7_U":
                model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG'
                past_end_str = 'TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC'


        logging.info("Model String: '{model_str}'."
                    " Past End String: '{past_end_str}'.")

        return model_str, past_end_str

    def check_model_parts(self, model_str, past_end_str):
        """
        Description:
            Making sure the model string and past end string make sense
        Args:
            model_str (str): 
            past_end_str (str):
        """
        possible_values = ["A","C","T","G","N"]

        for x in model_str:
            if x.upper() not in possible_values:
                raise Exception(f"Unknown character: '{x}'. possible values are" \
                                + " " + ", ".join(possible_values))
        for x in past_end_str:
            if x not in possible_values:
                raise Exception(f"Unknown character: '{x}', possible values are" \
                                + " " + ", ".join(possible_values))




