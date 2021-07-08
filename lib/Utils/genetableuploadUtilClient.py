
import os
import logging
import re
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class genetableuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_genes_table(self):
        """
        The upload method

        We perform a number of steps:
        Get name of genetable as it is in staging.
        Find the genetable in /staging/genetable_name
        Get the output name for the genetable
        Get the column headers for the genes file for
            data and testing purposes. Should be len 11.
        Test if genetable is well-formed.
        We send the file to shock using dfu.
        We get the handle and save the object with all
            the necessary information- including related genome.
        params should include:
            output_names,
            staging_file_names,
            ws_obj,
            workspace_id,

        """
        print("params: ", self.params)
        self.validate_import_genetable_from_staging_params()

        # Name of file in staging:
        stg_fs = self.params["staging_file_names"]
        if len(stg_fs) != 1:
            raise Exception("Expecting a single genes table file, got a different number"
                            f" of staging files: {len(stg_fs)}. Files: " + \
                            ", ".join(sgf_fs))
        else:
            staging_genes_fp_name = stg_fs[0]

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception("Expecting a single output name, got a different number"
                            f": {len(op_nms)}. Output Names: " + \
                            ", ".join(op_nms))
        else:
            genetable_name = op_nms[0]


        print("genetable_name: ", genetable_name)
        print("top dir /:", os.listdir("/"))
        print("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet! Error will be thrown")
        else:
            print("Succesfully recognized staging directory")
        # This is the path to the genes file
        genetable_fp = os.path.join(self.staging_folder, staging_genes_fp_name)

        # CHECK genes FILE:
        column_header_list, num_lines, genes_df = self.check_genes_file(genetable_fp,
                                                             self.params["sep_type"])

        if len(column_header_list) != 11:
            print(
                "WARNING: Number of columns is not 11 as expected: {}".format(
                    len(column_header_list)
                )
            )
        # We copy the file from staging to scratch
        new_genes_fp = os.path.join(self.shared_folder, genetable_name)

        if self.params["sep_type"] == "TSV":
            shutil.copyfile(genetable_fp, new_genes_fp)
        else:
            genes_df.to_csv(new_genes_fp, sep="\t", index=False)

        genetable_fp = new_genes_fp
        # We create the handle for the object:
        file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": genetable_fp, "make_handle": True, "pack": "gzip"}
        )
        # The following var res_handle only created for simplification of code
        res_handle = file_to_shock_result["handle"]

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()
        #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
        #        self.params['username'], str(date_time))

        # We create the data for the object
        genes_data = {
            "file_type": "KBaseRBTnSeq.RBTS_InputGenesTable",
            "input_genes_table": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "column_header_list": column_header_list,
            "column_headers_str": ", ".join(column_header_list),
            "num_lines": str(num_lines),
            "related_genome_ref": self.params["genome_ref"],
            "related_organism_scientific_name": self.params["organism_scientific_name"]
        }






        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBaseRBTnSeq.RBTS_InputGenesTable",
                    "data": genes_data,
                    "name": genetable_name,
                }
            ],
        }
        # save_objects returns a list of object_infos
        dfu_object_info = self.dfu.save_objects(save_object_params)[0]
        print("dfu_object_info: ")
        print(dfu_object_info)
        return {
            "Name": dfu_object_info[1],
            "Type": dfu_object_info[2],
            "Date": dfu_object_info[3],
        }

    def validate_import_genetable_from_staging_params(self):
        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "description",
            "organism_scientific_name",
            "genome_ref",
            "output_names"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def check_genes_file(self, genetable_fp, separator):
        """
        We check the genes file by initializing into dict format

        """

        
        sep = "\t" if separator == "TSV" else ","


        genes_dtypes = {
                'locusId': str,
                'sysName': str,
                'type': int,
                'scaffoldId': str,
                'begin': int,
                'end': int,
                'strand': str,
                'name': str,
                'desc': str,
                'GC': float,
                'nTA': int
                }
        genes_df = pd.read_table(genetable_fp, sep=sep, dtype=genes_dtypes)

        req_cols = ['locusId', 'sysName', 'type', 'scaffoldId', 'begin', 'end', 'strand',
                    'name', 'desc', 'GC', 'nTA']
        for x in req_cols:
            if x not in genes_df.columns:
                raise Exception(f"Required column name {x} not found in genes table.")


        # Checking for duplicates
        locusIds_dict = {}
        for locusId in genes_df["locusId"]:
            if locusId in locusIds_dict:
                raise Exception(f"Duplicate locusId: {locusId}")
            else:
                locusIds_dict[locusId] = 1
        for ix, strand in genes_df["strand"].iteritems():
            if strand not in ["+", "-"]:
                raise Exception(f"Incorrect strand value: {strand} at row {ix}")

        for ix, frac in genes_df["GC"].iteritems():
            if frac < 0 or frac > 1:
                raise Exception("GC content must be between 0 and 1"
                                f" Value at row {ix} is {frac}")


        logging.info("genetable columns are: " + ", ".join(genes_df.columns))

        return [list(genes_df.columns), genes_df.shape[0], genes_df]


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


