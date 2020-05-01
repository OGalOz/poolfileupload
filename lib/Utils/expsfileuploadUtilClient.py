import os
import logging
import re
import shutil
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
            staging_file_name,
            genome_ref,
            description,
            output_name
        """

        print("params: ", self.params)
        self.validate_import_expsfile_from_staging_params()

        # Name of file in staging: (file name or absolute path?)
        staging_exps_fp_name = self.params["staging_file_name"]

        # Output name of exps file:
        expsfile_name = self.params["output_name"]

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
        column_header_list, num_rows, setNames = self.check_exps_file(
                                                                    expsfile_fp)

        # We copy the file from staging to scratch
        new_exps_fp = os.path.join(self.shared_folder, expsfile_name)
        shutil.copyfile(expsfile_fp, new_exps_fp)
        expsfile_fp = new_exps_fp

        # We create the handle for the object:
        file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": expsfile_fp, "make_handle": True, "pack": "gzip"}
        )

        # The following var res_handle only created for simplification of code
        res_handle = file_to_shock_result["handle"]

        # We create the data for the object
        exps_data = {
            "file_type": "KBasePoolTSV.Experiments",
            "expsfile": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "file_name": res_handle["file_name"],
            "column_header_list": column_header_list,
            "num_lines": str(num_rows),
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
                    "type": "KBasePoolTSV.Experiments",
                    "data": exps_data,
                    "name": expsfile_name,
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

    def validate_import_expsfile_from_staging_params(self):

        # check for required parameters
        for p in [
            "staging_file_name",
            "genome_ref",
            "description",
            "output_name"
        ]:
            if p not in self.params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def check_exps_file(self, expsfile_fp):

        required = [
            "SetName",
            "Index",
            "Description",
            "Date_pool_expt_started",
        ]

        cols, num_rows, setNames = self.read_table(expsfile_fp, required)

        return [cols, num_rows, setNames]

    def read_table(self, fp, required):
        """
        Following function takes a filename and a list of required fields i
        (file is TSV)
        returns list of headers
        Does not return header line
        """
        with open(fp, "r") as f:
            file_str = f.read()
        file_list = file_str.split("\n")
        header_line = file_list[0]
        # Check for Mac Style Files
        if re.search(r"\t", header_line) and re.search(r"\r", header_line):
            raise Exception(
                (
                    "Tab-delimited input file {} is a Mac-style text file "
                    "which is not supported.\n"
                    "Use\ndos2unix -c mac {}\n to convert it to a Unix "
                    "text file.\n"
                ).format(fp, fp)
            )
        cols = header_line.split("\t")
        cols_dict = {}
        for i in range(len(cols)):
            cols_dict[cols[i]] = i
        for field in required:
            if field not in cols_dict:
                raise Exception(
                    "No field {} in {}. Must include fields".format(field, fp)
                    + "\n{}".format(" ".join(required))
                )
        rows = []
        # This is unique to Experiments
        setNames = []
        for i in range(1, len(file_list)):
            line = file_list[i]
            # if last line empty
            if len(line) == 0:
                continue
            line = re.sub(r"[\r\n]+$", "", line)
            split_line = line.split("\t")
            setNames.append(split_line[0])
            if not len(split_line) == len(cols):
                raise Exception(
                    "Wrong number of columns in:\n{}\nin {} l:{}".format(line, fp, i)
                )
            new_dict = {}
            for i in range(len(cols)):
                new_dict[cols[i]] = split_line[i]
            rows.append(new_dict)

        return [cols, len(file_list), setNames]

    def get_genome_organism_name(self, genome_ref):
        # Getting the organism name using WorkspaceClient
        ws = self.params["ws_obj"]
        res = ws.get_objects2(
            {"objects": [{"ref": genome_ref, "included": ["scientific_name"]}]}
        )
        scientific_name = res["data"][0]["data"]["scientific_name"]
        return scientific_name
