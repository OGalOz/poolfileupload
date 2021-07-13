import os
import logging
import re
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class poolfileuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_poolfile(self):
        """
        The upload method

        We perform a number of steps:
        Get name of poolfile as it is in staging.
        Find the poolfile in /staging/poolfile_name
        Get the output name for the poolfile
        Get the column headers for the pool file for
            data and testing purposes. Should be len 12.
        Test if poolfile is well-formed.
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
        self.validate_import_poolfile_from_staging_params()

        # Name of file in staging:
        stg_fs = self.params["staging_file_names"]
        if len(stg_fs) != 1:
            raise Exception("Expecting a single pool file, got a different number"
                            f" of staging files: {len(stg_fs)}. Files: " + \
                            ", ".join(sgf_fs))
        else:
            staging_pool_fp_name = stg_fs[0]

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception("Expecting a single output name, got a different number"
                            f": {len(op_nms)}. Output Names: " + \
                            ", ".join(op_nms))
        else:
            poolfile_name = op_nms[0]


        print("poolfile_name: ", poolfile_name)
        print("top dir /:", os.listdir("/"))
        print("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet! Error will be thrown")
        else:
            print("Succesfully recognized staging directory")
        # This is the path to the pool file
        poolfile_fp = os.path.join(self.staging_folder, staging_pool_fp_name)

        # CHECK POOL FILE:
        column_header_list, num_lines, pool_df = self.check_pool_file(poolfile_fp,
                                                             self.params["sep_type"])

        if len(column_header_list) != 12:
            print(
                "WARNING: Number of columns is not 12 as expected: {}".format(
                    len(column_header_list)
                )
            )
        # We copy the file from staging to scratch
        new_pool_fp = os.path.join(self.shared_folder, poolfile_name)

        if self.params["sep_type"] == "TSV":
            shutil.copyfile(poolfile_fp, new_pool_fp)
        else:
            pool_df.to_csv(new_pool_fp, sep="\t", index=False)

        poolfile_fp = new_pool_fp
        # We create the handle for the object:
        file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": poolfile_fp, "make_handle": True, "pack": "gzip"}
        )
        # The following var res_handle only created for simplification of code
        res_handle = file_to_shock_result["handle"]

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()
        #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
        #        self.params['username'], str(date_time))
        fastq_refs = []

        # We create the data for the object
        pool_data = {
            "file_type": "KBaseRBTnSeq.RBTS_MutantPoolFile",
            "poolfile": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "column_header_list": column_header_list,
            "column_headers_str": ", ".join(column_header_list),
            "num_lines": str(num_lines),
            "fastqs_used": fastq_refs,
            "fastqs_used_str": "NA",
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "related_genes_table_ref": self.params["genes_table_ref"],
            "model_ref": self.params["model_ref"],
            "related_organism_scientific_name": self.get_genome_organism_name(
                self.params["genes_table_ref"]
            ),
            "description": "Manual Upload: " + self.params["description"],
        }

        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBaseRBTnSeq.RBTS_MutantPoolFile",
                    "data": pool_data,
                    "name": poolfile_name,
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

    def validate_import_poolfile_from_staging_params(self):
        prms = self.params
        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "genes_table_ref",
            "description",
            "output_names"
        ]:
            if p not in prms:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        if "model_ref" not in prms or prms["model_ref"] == "" \
            or prms["model_ref"] is None:
            raise Exception("When uploading a poolfile, you must include a model ref.")

    def check_pool_file(self, poolfile_fp, separator):
        """
        We check the pool file by initializing into dict format

        """

        
        sep = "\t" if separator == "TSV" else ","

        """
        dtypes = {
                "barcode": str,
                "rcbarcode": str,
                "scaffold": str,
                "nTot": 'Int64',
                "n": 'Int64',
                "strand": str,
                "pos": 'Int64',
                "n2": 'Int64',
                "scaffold2": str,
                "strand2": str,
                "pos2": 'Int64',
                "nPastEnd": 'Int64'
        }
        """

        pool_df = pd.read_table(poolfile_fp, sep=sep)

        req_cols = ["barcode", "rcbarcode", "scaffold", 
                    "pos", "strand"] 

        for x in req_cols:
            if x not in pool_df.columns:
                raise Exception(f"Required column name {x} not found in pool file.")

        # Checking for duplicates
        barcodes_dict = {}
        for barcode in pool_df["barcode"]:
            if barcode in barcodes_dict:
                raise Exception(f"Duplicate barcode: {barcode}")
            else:
                barcodes_dict[barcode] = 1
        for ix, strand in pool_df["strand"].iteritems():
            if strand not in ["+", "-"]:
                if pool_df["scaffold"].iloc[ix] != "pastEnd":
                    raise Exception(f"Incorrect strand value: {strand} at row {ix}")
                else:
                    logging.info(f"Found pastEnd hit at row {ix}")


        for ix, pos in pool_df["pos"].iteritems():
            if pos < 0:
                raise Exception("Positions must be positive."
                                f" Value at row {ix} is {pos}")

        if "pos2" in pool_df.columns:
            for ix, pos in pool_df["pos2"].iteritems():
                if pos < 0:
                    raise Exception("Positions must be positive."
                                    f" Value at row {ix} is {pos}")


        logging.info("Poolfile columns are: " + ", ".join(pool_df.columns))

        return [list(pool_df.columns), pool_df.shape[0], pool_df]


    def get_genome_organism_name(self, gene_table_ref):
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
        logging.info("Workspace get objects 2 results:")
        logging.info(res)

        scientific_name = res["data"][0]["data"]["related_organism_scientific_name"]
        return scientific_name
