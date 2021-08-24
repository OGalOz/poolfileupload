import os
import logging
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class barcodecountfileuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_barcodecountfile(self):

        """
        The upload method

        We perform a number of steps:
        Get name of barcodecount file as it is in staging.
        Find the barcodecount file in /staging/barcodecount_name
        Get the output name for the barcodecount file
        Get the column headers for the barcode count file for
            data and testing purposes. 
        Test if barcodecount file is well-formed.
        NOTE: We use output_name as set_name - it is important that
            these are equivalent!!!!!
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
        self.validate_import_file_from_staging_params()

        # Name of files in staging (Not path):
        staging_fp_names = self.params["staging_file_names"]
        # Output name of barcodecount file:
        barcodecount_names = self.params["output_names"]

        if len(staging_fp_names) != len(barcodecount_names):
            raise Exception("The number of barcodecount staging files "
                            "must be equal to the number of output names, "
                            "as the two correspond to each other. "
                            f"Number of barcodecount staging files: {len(staging_fp_names)}."
                            f" Number of output names: {len(barcodecount_names)}.")


        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet!")
        else:
            print("Succesfully recognized staging directory")

        op_info_list = []

        for i in range(len(staging_fp_names)):
            staging_fp_name, crnt_bc_op_name = staging_fp_names[i], barcodecount_names[i]
            # This is the path to the barcodecount file in staging
            barcodecount_fp = os.path.join(self.staging_folder, staging_fp_name)
            # We check correctness of barcodecount file in staging
            column_header_list, num_lines, bc_df = self.check_barcodecount_file(barcodecount_fp, 
                                                                      self.params["sep_type"])

            # We copy the file from staging to scratch
            new_bc_fp = os.path.join(self.shared_folder, crnt_bc_op_name)

            if self.params["sep_type"] == "TSV":
                shutil.copyfile(barcodecount_fp, new_bc_fp)
            else:
                bc_df.to_csv(new_bc_fp, sep="\t", index=False)
            #barcodecount_scratch_fp is location of barcodecount file in scratch
            barcodecount_scratch_fp = new_bc_fp


            # We create the KBase handle for the object:
            file_to_shock_result = self.dfu.file_to_shock(
                {"file_path": barcodecount_scratch_fp, "make_handle": True, "pack": "gzip"}
            )
            # The following var res_handle only created for simplification of code
            res_handle = file_to_shock_result["handle"]

            # Keep track of our own datetime
            date_time = datetime.datetime.utcnow()
            #new_desc = "Uploaded by {} on (UTC) {} using Uploader. User Desc: ".format(
            #        self.params['username'], str(date_time))
            fastq_refs = []

            # We create the data for the object
            barcodecount_data = {
                "file_type": "KBaseRBTnSeq.RBTS_BarcodeCount",
                "barcodecount": res_handle["hid"],
                # below should be shock
                "handle_type": res_handle["type"],
                "shock_url": res_handle["url"],
                "shock_node_id": res_handle["id"],
                "compression_type": "gzip",
                "column_header_list": column_header_list,
                "column_headers_str": ", ".join(column_header_list),
                "mutantpool_ref": self.params["mutantpool_ref"],
                "fastqs_used": fastq_refs,
                "fastqs_used_str": "NA",
                "file_name": res_handle["file_name"],
                "utc_created": str(date_time),
                "set_name": crnt_bc_op_name, 
                "num_lines": str(num_lines),
                "protocol_used": self.params["protocol_type"],
                "related_genome_ref": self.params["genome_ref"],
                "related_organism_scientific_name": self.get_genome_organism_name(
                    self.params["genome_ref"]
                ),
                "description": "Manual Upload: " + self.params["description"],
            }
    
            # To get workspace id:
            ws_id = self.params["workspace_id"]
            save_object_params = {
                "id": ws_id,
                "objects": [
                    {
                        "type": "KBaseRBTnSeq.RBTS_BarcodeCount",
                        "data": barcodecount_data,
                        "name": crnt_bc_op_name,
                    }
                ],
            }
            # save_objects returns a list of object_infos
            dfu_object_info = self.dfu.save_objects(save_object_params)[0]
            print("dfu_object_info: ")
            print(dfu_object_info)
    
            op_info_list.append({
                "Name": dfu_object_info[1],
                "Type": dfu_object_info[2],
                "Date": dfu_object_info[3],
            })

        return op_info_list
    

    def check_barcodecount_file(self, barcodecount_fp, sep):
        """
        We check the barcodecount file by initializing into dict format
   
        Currently a weak test- should add more testing capabilities.
        """
        # Expected fields
        exp_f = "barcode rcbarcode scaffold strand pos".split(" ") 

        sep = "\t" if sep == "TSV" else ","

        barcodecount_df = pd.read_table(barcodecount_fp, sep=sep)

        for field in exp_f:
            if field not in barcodecount_df.columns:
                raise Exception(f"Expected field {field} in barcodecount but didn't get it."
                                " Current fields: " + ", ".join(barcodecount_df.columns))
        
        for ix, val in barcodecount_df["strand"].iteritems():
            if val not in ["+", "-"]:
                raise Exception("The strand column of barcodecount must be one of '+' or '-',"
                                f" current value at row number {ix + 1} is {val}.")

        for ix, val in barcodecount_df["pos"].iteritems():
            if not isinstance(val, int):
                raise Exception("The 'pos' column of barcodecount must be an integer,"
                                f" current value at row number {ix + 1} is {val}.")

    
        return [list(barcodecount_df.columns), barcodecount_df.shape[0], barcodecount_df]

    def validate_import_file_from_staging_params(self):
        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "protocol_type",
            "genome_ref",
            "mutantpool_ref",
            "description",
            "output_names",
            "ws_obj",
            "workspace_id"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

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
