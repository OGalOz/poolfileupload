import os
import logging
import re
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class expsfileuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_expsfile(self):
        """
        The upload method

        We perform a number of steps:
        Get name of expsfile as it is in staging.
        Find the expsfile in /staging/expsfile_name
        Get the output name for the expsfile
        Get the column headers for the exps file for
            data and testing purposes. 
        Test if expsfile is well-formed.
        We send the file to shock using dfu.
        We get the handle and save the object with all
            the necessary information- including related genome.
        params should include:
            username,
            staging_file_names,
            genome_ref,
            description,
            output_names
        """

        print("params: ", self.params)
        self.validate_import_expsfile_from_staging_params()


        # Name of file in staging:
        stg_fs = self.params["staging_file_names"]
        if len(stg_fs) != 1:
            raise Exception("Expecting a single experiments file, got a different number"
                            f" of staging files: {len(stg_fs)}. Files: " + \
                            ", ".join(sgf_fs))
        else:
            staging_exps_fp_name = stg_fs[0]

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception("Expecting a single output name, got a different number"
                            f": {len(op_nms)}. Output Names: " + \
                            ", ".join(op_nms))
        else:
            expsfile_name = op_nms[0]


        print("expsfile_name: ", expsfile_name)
        print("top dir /:", os.listdir("/"))
        print("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet!")
        else:
            print("Succesfully recognized staging directory")

        # This is the path to the exps file
        expsfile_fp = os.path.join(self.staging_folder, staging_exps_fp_name)

        # We check correctness of exps file. Returns list and int
        column_header_list, num_rows, exps_df = self.check_exps_file(
                                                            expsfile_fp,
                                                            self.params['sep_type'])

        # Making sure every column name is a string (not int)
        column_header_list = [str(x) for x in column_header_list]


        # We copy the file from staging to scratch
        new_exps_fp = os.path.join(self.shared_folder, expsfile_name)

        # If TSV we change nothing, otherwise we convert it to TSV
        if self.params["sep_type"] == "TSV":
            shutil.copyfile(expsfile_fp, new_exps_fp)
        else:
            # Converting CSV to TSV
            exps_df.to_csv(new_exps_fp, sep="\t", index=False)

        expsfile_fp = new_exps_fp

        # We create the handle for the object:
        file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": expsfile_fp, "make_handle": True, "pack": "gzip"}
        )

        # The following var res_handle only created for simplification of code
        res_handle = file_to_shock_result["handle"]

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()
        #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
        #        self.params['username'], str(date_time))

        # We create the data for the object
        exps_data = {
            "file_type": "KBaseRBTnSeq.RBTS_ExperimentsTable",
            "expsfile": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "column_header_list": column_header_list,
            "num_cols": str(len(column_header_list)),
            "num_lines": str(num_rows),
            "related_genome_ref": self.params["genome_ref"],
            "related_organism_scientific_name": self.get_genome_organism_name(
                self.params["genome_ref"]
            ),
            "description": self.params["description"],
        }

        '''
        #DEBUGGING:
        logging.info("Printing types for data dict: ")
        for k, val in exps_data.items():
            logging.info(f"{k}: {type(val)}")
        '''

        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBaseRBTnSeq.RBTS_ExperimentsTable",
                    "data": exps_data,
                    "name": expsfile_name,
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

    def validate_import_expsfile_from_staging_params(self):

        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "genome_ref",
            "description",
            "output_names"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def check_exps_file(self, expsfile_fp, separator):
        """
        Args:
            expsfile_fp (str): Path to experiments file
            separator (str): "CSV" or "TSV", "CSV" implying comma, "TSV" tab
        """
        sep = "\t" if separator == "TSV" else ","

        required = [
            "SetName",
            "Index",
            "Description",
            "Date_pool_expt_started",
            "Group"
        ]


        exps_df = pd.read_table(expsfile_fp, sep=sep) 

        for col_name in required:
            if col_name not in exps_df.columns:
                raise Exception(f"Required column name: {col_name} missing from file headers.",
                                " File headers are: " + ", ".join(exps_df.columns))

        # Set Name and Index must be unique: (+ in pandas means vector wise)
        st_nm_ix = list(exps_df["SetName"] + exps_df["Index"])
        unq_dict = {}
        for val in st_nm_ix:
            if val in unq_dict:
                raise Exception("SetName and Index have to be unique per row."
                                f" Repeated value: {val}")
            else:
                unq_dict[val] = 1

        if "control_bool" in exps_df.columns:
            if "control_group" not in exps_df.columns:
                raise Exception("If 'control_bool' column is in Experiments file, then"
                                " 'control_group' column must also be included."
                                " File headers are currently: " + ", ".join(exps_df.columns))
        elif "control_group" in exps_df.columns:
                raise Exception("If 'control_group' column is in Experiments file, then"
                                " 'control_bool' column must also be included."
                                " File headers are currently: " + ", ".join(exps_df.columns))

        logging.info("For experiments file, num rows is: "
                    f" {exps_df.shape[0]}, num columns is {exps_df.shape[1]}")

        return [list(exps_df.columns), exps_df.shape[0], exps_df]


    def get_genome_organism_name(self, genome_ref):
        # Getting the organism name using WorkspaceClient
        ws = self.params['ws_obj'] 
        res = ws.get_objects2(
            {
                "objects": [
                    {
                        "ref": genome_ref,
                        "included": ["scientific_name"],
                    }
                ]
            }
        )
        scientific_name = res["data"][0]["data"]["scientific_name"]
        return scientific_name
    def get_genome_organism_name_from_genes_table(self, gene_table_ref):
        # Getting the organism name using WorkspaceClient
        ws = self.params['ws_obj'] 
        res = ws.get_objects2(
            {
                "objects": [
                    {
                        "ref": gene_table_ref,
                        "included": ["related_organism_scientific_name"],
                    }
                ]
            }
        )
        scientific_name = res["data"][0]["data"]["related_organism_scientific_name"]
        return scientific_name
