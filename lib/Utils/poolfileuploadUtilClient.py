import os
import logging
import re
import shutil
import datetime
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

        """
        print("params: ", self.params)
        self.validate_import_poolfile_from_staging_params()

        # Name of file in staging:
        staging_pool_fp_name = self.params["staging_file_name"]

        # Output name of pool file:
        poolfile_name = self.params["output_name"]

        print("poolfile_name: ", poolfile_name)
        print("top dir /:", os.listdir("/"))
        print("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet! Error will be thrown")
        else:
            print("Succesfully recognized staging directory")
        # This is the path to the pool file
        poolfile_fp = os.path.join(self.staging_folder, staging_pool_fp_name)
        # We check correctness of pool file
        column_header_list, num_lines = self.check_pool_file(poolfile_fp)
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
        fastq_refs = ["Manual Upload"]

        # We create the data for the object
        pool_data = {
            "file_type": "KBasePoolTSV.PoolFile",
            "poolfile": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "column_header_list": column_header_list,
            "num_lines": num_lines,
            "fastqs_used": fastq_refs,
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "related_genome_ref": self.params["genome_ref"],
            "related_organism_scientific_name": self.get_genome_organism_name(
                self.params["genome_ref"]
            ),
            "description": self.params["description"],
        }

        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBasePoolTSV.PoolFile",
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
        # check for required parameters
        for p in [
            "username",
            "staging_file_name",
            "genome_ref",
            "description",
            "output_name"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def check_pool_file(self, poolfile_fp):
        """
        We check the pool file by initializing into dict format

        The function "init_pool_dict" runs the tests to see if the file is
        correct.
        """
        col_header_list = []
        # Parse pool file and check for errors
        test_vars_dict = {"poolfile": poolfile_fp, "report_dict": {"warnings": []}}
        try:
            col_header_list, num_lines = self.init_pool_dict(test_vars_dict)
        except Exception:
            logging.warning(
                "Pool file seems to have errors - " + "Please check and reupload."
            )
            raise Exception
        return [col_header_list, num_lines]

    def init_pool_dict(self, vars_dict):

        # pool dict is rcbarcode to [barcode, scaffold, strand, pos]
        pool = {}
        num_lines = 0
        with open(vars_dict["poolfile"], "r") as f:
            header_str = f.readline()
            if header_str == '':
                raise Exception("Issue with pool file - first line empty")
            num_lines += 1
            column_header_list = [x.strip() for x in header_str.split("\t")]
            crnt_line = f.readline() 
            while crnt_line != '':
                num_lines += 1
                crnt_line.rstrip()
                pool = self.check_pool_line_and_add_to_pool_dict(
                    crnt_line, pool, vars_dict
                )
                crnt_line = f.readline()
        if len(pool.keys()) == 0:
            raise Exception("No entries in pool file")
        return [column_header_list, num_lines]

    def check_pool_line_and_add_to_pool_dict(self, pool_line, pool, vars_dict):
        """
        For a pool line to be correct it has to follow a few rules.

        We care about the first 7 columns of each pool line.
        The first line in the file is the headers, and the first 7 are
        barcode, rcbarcode, nTot, n, scaffold, strand, pos
        Both the barcodes and rcbarcodes must be entirely made up of
        characters from "ACTG". Position must be made up of any number
        of digits (including 0). Strand is from "+","-","".
        If the rcbarcode already exists in the pool, then there is a
        problem with the pool file. Each rcbarcode must be unique.
        """
        # We get first 7 columns of pool_line (out of 12)
        split_pool_line = pool_line.split("\t")[:7]
        # We remove spaces:
        for x in split_pool_line:
            x.replace(" ", "")
        if len(split_pool_line) >= 7:
            # We unpack
            (
                barcode,
                rcbarcode,
                undef_1,
                undef_2,
                scaffold,
                strand,
                pos,
            ) = split_pool_line
        else:
            warning_text = "pool file line with less than 7 tabs:\n{}".format(pool_line)
            vars_dict["report_dict"]["warnings"].append(warning_text)
            logging.warning(warning_text)
            barcode = "barcode"

        if barcode == "barcode":
            # Header line
            pass
        else:
            if not re.search(r"^[ACGT]+$", barcode):
                logging.debug(len(barcode))
                raise Exception("Invalid barcode: |{}|".format(barcode))
            if not re.search(r"^[ACGT]+$", rcbarcode):
                raise Exception("Invalid rcbarcode: |{}|".format(rcbarcode))
            if not (pos == "" or re.search(r"^\d+$", pos)):
                raise Exception("Invalid position: |{}|".format(pos))
            if not (strand == "+" or strand == "-" or strand == ""):
                raise Exception("Invalid strand: |{}|".format(strand))
            if rcbarcode in pool:
                raise Exception("Duplicate rcbarcode.")
            pool[rcbarcode] = [barcode, scaffold, strand, pos]
        return pool

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
