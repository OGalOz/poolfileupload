import os
import logging
import re
import shutil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class poolcountfileuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_poolcountfile(self):
        raise Exception("INCOMPLETE - DO NOT USE")

        """
        The upload method

        We perform a number of steps:
        Get name of poolcount file as it is in staging.
        Find the poolcount file in /staging/poolcount_name
        Get the output name for the poolcount file
        Get the column headers for the pool count file for
            data and testing purposes. Should be len 12.
        Test if poolcount file is well-formed.
        We send the file to shock using dfu.
        We get the handle and save the object with all
            the necessary information- including related genome.

        """
        print("params: ", self.params)
        self.validate_import_poolfile_from_staging_params()

        # Name of file in staging:
        staging_pool_fp_name = self.params["staging_file_name"]

        # Output name of pool file:
        poolfile_name = self.params["poolfile_name"]

        print("poolfile_name: ", poolfile_name)
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet!")
        else:
            print("Succesfully recognized staging directory")

        # This is the path to the pool file
        poolfile_fp = os.path.join(self.staging_folder, staging_pool_fp_name)
        # We check correctness of pool file
        column_header_list = self.check_pool_file(poolfile_fp)
        if len(column_header_list) != 12:
            print(
                "WARNING: Number of columns is not 12 as expected: {}".format(
                    len(column_header_list)
                )
            )
        # We copy the file from staging to scratch
        new_pool_fp = os.path.join(self.shared_folder, poolfile_name)
        shutil.copyfile(poolfile_fp, new_pool_fp)
        poolfile_fp = new_pool_fp


        # We create the KBase handle for the object:
        file_to_shock_result = up['dfu'].file_to_shock(
            {"file_path": up['poolcount_fp'], "make_handle": True, "pack": "gzip"}
        )
        # The following var res_handle only created for simplification of code
        res_handle = file_to_shock_result["handle"]

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()
        new_desc = "Uploaded by {} on (UTC) {} using Uploader\n".format(
                self.params['username'], str(date_time))
        fastq_refs = ["Manual Upload"]

        # We create the data for the object
        pool_data = {
            "file_type": "KBasePoolTSV.PoolCount",
            "poolcount": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "column_header_list": column_header_list,
            "fastqs_used": fastq_refs,
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "set_name": up['set_name'], 
            "num_lines": str(num_lines),
            "related_genome_ref": up["genome_ref"],
            "related_organism_scientific_name": get_genome_organism_name(
                up["genome_ref"],
                up['ws_obj']
            ),
            "description": new_desc + up["poolcount_description"],
        }
    
        # To get workspace id:
        ws_id = up["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBasePoolTSV.PoolCount",
                    "data": pool_data,
                    "name": up['poolcount_name'],
                }
            ],
        }
        # save_objects returns a list of object_infos
        dfu_object_info = up['dfu'].save_objects(save_object_params)[0]
        print("dfu_object_info: ")
        print(dfu_object_info)
        return {
            "Name": dfu_object_info[1],
            "Type": dfu_object_info[2],
            "Date": dfu_object_info[3],
        }
    

    def check_poolcount_file(poolcount_fp):
        """
        We check the pool file by initializing into dict format
    
        The function "init_pool_dict" runs the tests to see if the file is
        correct.
        """
        # Expected fields
        exp_f = "barcode rcbarcode scaffold strand pos".split(" ") 
    
        with open(poolcount_fp, "r") as f:
            f_str = f.read()
        f_list = f_str.split('\n')
        num_lines = len(f_list)
        header_line = f_list[0]
    
    
        if header_line == '':
            raise Exception("File format incorrect: " + poolcount_fp)
    
        fields = header_line.split("\t")
    
        if not (len(fields) >= 6):
            raise Exception("Too few fields in " + poolcount_fp)
        for i in range(len(exp_f)):
            if not fields[i] == exp_f[i]:
                raise Exception(
                            "Expected {} but field is {}".format(exp_f[i], fields[i])
                        )
        return [fields, num_lines]

    def validate_import_poolfile_from_staging_params(self):
        # check for required parameters
        for p in [
            "username",
            "staging_file_name",
            "genome_ref",
            "description",
            "output_name",
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
