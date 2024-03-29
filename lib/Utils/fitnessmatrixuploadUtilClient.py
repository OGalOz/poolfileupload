import os
import logging
import re
import shutil
import datetime
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class fitnessmatrixuploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_fitnessmatrix(self):
        """
        The upload method

        We perform a number of steps:
        Get name of fitnessmatrix as it is in staging.
        Find the fitnessmatrix in /staging/op_datatype_name
        Get the output name for the fitnessmatrix
        Get the column headers for the pool file for
            data and testing purposes. Should be len 12.
        Test if fitnessmatrix is well-formed.
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
        self.validate_import_fitnessmatrix_from_staging_params()

        # Double checking number of files we want from staging
        strain_fit_bool = False
        stg_fs = self.params["staging_file_names"]
        if not len(stg_fs) in [2, 3]:
            raise Exception(
                "Expecting between 2/3 staging files, got a different number"
                f" of staging files: {len(stg_fs)}. Files: " + ", ".join(sgf_fs)
            )
        else:
            staging_fitness_matrix_fp_name = stg_fs[0]
            staging_t_score_matrix_fp_name = stg_fs[1]
            logging.info(
                "Using this file for the fitness matrix: "
                + staging_fitness_matrix_fp_name
                + ". "
            )
            logging.info(
                "Using this file for the t_score matrix: "
                + staging_t_score_matrix_fp_name
                + "."
            )
            if len(stg_fs) == 3:
                strain_fit_bool = True
                staging_strain_fit_table_fp_name = stg_fs[2]
                logging.info(
                    "Using this file for the strain fit matrix: "
                    + strain_fit_table_fp_name
                    + "."
                )

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception(
                "Expecting a single output name, got a different number"
                f": {len(op_nms)}. Output Names: " + ", ".join(op_nms)
            )
        else:
            op_datatype_name = op_nms[0]

        print("op_datatype_name: ", op_datatype_name)
        print("top dir /:", os.listdir("/"))
        print("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet! Cannot continue.")
        else:
            print("Succesfully recognized staging directory")
        # This is the path to the pool file
        fitnessmatrix_fp = os.path.join(
            self.staging_folder, staging_fitness_matrix_fp_name
        )
        t_scorematrix_fp = os.path.join(
            self.staging_folder, staging_t_score_matrix_fp_name
        )
        if strain_fit_bool:
            strain_fitmatrix_fp = os.path.join(
                self.staging_folder, staging_strain_fit_table_fp_name
            )

        # CHECK FILES:
        column_header_list, num_lines = self.check_matrix_files(
            fitnessmatrix_fp, t_scorematrix_fp, self.params["sep_type"]
        )
        if strain_fit_bool:
            self.check_strain_fit_table(strain_fitmatrix_fp, self.params["sep_type"])

        # We copy the files from staging to scratch
        new_fitness_matrix_fp = os.path.join(
            self.shared_folder, op_datatype_name + ".fit.tsv"
        )
        new_t_score_fp = os.path.join(
            self.shared_folder, op_datatype_name + ".t_score.tsv"
        )
        if strain_fit_bool:
            new_strain_fp = os.path.join(
                self.shared_folder, op_datatype_name + ".strain_fit.tsv"
            )

        if self.params["sep_type"] == "TSV":
            shutil.copyfile(fitnessmatrix_fp, new_fitness_matrix_fp)
            shutil.copyfile(t_scorematrix_fp, new_t_score_fp)
            if strain_fit_bool:
                shutil.copyfile(strain_fitmatrix_fp, new_strain_fp)
        else:
            # sep type is comma (CSVs)
            fit_df = pd.read_table(fitnessmatrix_fp, sep=",", keep_default_na=False)
            t_score_df = pd.read_table(t_scorematrix_fp, sep=",", keep_default_na=False)
            fit_df.to_csv(new_fitness_matrix_fp, sep="\t", index=False)
            t_score_df.to_csv(new_t_score_fp, sep="\t", index=False)
            if strain_fit_bool:
                strain_df = pd.read_table(
                    strain_fitmatrix_fp, sep=",", keep_default_na=False
                )
                strain_df.to_csv(new_strain_fp, sep="\t", index=False)

        # We create the handles for the objects:
        fitness_file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": new_fitness_matrix_fp, "make_handle": True, "pack": "gzip"}
        )
        t_score_file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": new_t_score_fp, "make_handle": True, "pack": "gzip"}
        )
        fitness_res_handle = fitness_file_to_shock_result["handle"]
        t_score_res_handle = t_score_file_to_shock_result["handle"]
        if strain_fit_bool:
            strain_fit_file_to_shock_result = self.dfu.file_to_shock(
                {"file_path": new_strain_fp, "make_handle": True, "pack": "gzip"}
            )
            strain_fit_res_handle = strain_fit_file_to_shock_result["handle"]

        # We create a better Description by adding date time and username
        date_time = datetime.datetime.utcnow()

        # We create the data for the object
        matrices_data = {
            "file_type": "KBaseRBTnSeq.RBTS_Gene_Fitness_T_Matrix",
            "fit_scores_handle": fitness_res_handle["hid"],
            "t_scores_handle": t_score_res_handle["hid"],
            # below should be shock
            "handle_type": fitness_res_handle["type"],
            "fitness_shock_url": fitness_res_handle["url"],
            "t_scores_shock_url": t_score_res_handle["url"],
            "fitness_shock_node_id": fitness_res_handle["id"],
            "t_scores_shock_node_id": t_score_res_handle["id"],
            "compression_type": "gzip",
            "fitness_file_name": fitness_res_handle["file_name"],
            "t_scores_file_name": t_score_res_handle["file_name"],
            "utc_created": str(date_time),
            "column_header_list": column_header_list,
            "num_cols": str(len(column_header_list)),
            "num_lines": str(num_lines),
            "related_genome_ref": self.params["genome_ref"],
            "poolcounts_used": [],
            "related_experiments_ref": self.params["experiments_ref"],
            "related_organism_scientific_name": self.get_genome_organism_name(
                self.params["genome_ref"]
            ),
            "description": "Manual Upload: " + self.params["description"],
        }
        if strain_fit_bool:
            matrices_data["strain_fit_handle"] = strain_fit_res_handle["hid"]
            matrices_data["strain_fit_shock_url"] = strain_fit_res_handle["url"]
            matrices_data["strain_fit_shock_node_id"] = strain_fit_res_handle["id"]
            matrices_data["strain_fit_file_name"] = strain_fit_res_handle["file_name"]

        # To get workspace id:
        ws_id = self.params["workspace_id"]
        save_object_params = {
            "id": ws_id,
            "objects": [
                {
                    "type": "KBaseRBTnSeq.RBTS_Gene_Fitness_T_Matrix",
                    "data": matrices_data,
                    "name": op_datatype_name,
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

    def validate_import_fitnessmatrix_from_staging_params(self):
        prms = self.params
        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "genome_ref",
            "experiments_ref",
            "description",
            "output_names",
        ]:
            if p not in prms:
                raise ValueError(
                    'When uploading a fitness matrix, "{}" parameter is required, but missing'.format(
                        p
                    )
                )

    def check_matrix_files(self, fitness_matrix_fp, t_score_matrix_fp, separator):
        """
        Args:
           fitness_matrix_fp (str): Path to fitness matrix file
           t_score_matrix_fp (str): Path to t score matrix file
           separator (str): "," or "\t"
        Returns:
            list<list<column_names (str)>, num_rows (int)>
        Description:
            We check the matrix files by initializing into dict format
        """

        sep = "\t" if separator == "TSV" else ","

        """
        dtypes = {
            "orgId", "locusId", "sysName", "geneName", "desc" All strings
        }
        """

        req_cols = ["locusId", "sysName", "geneName", "desc"]

        fitness_df = pd.read_table(fitness_matrix_fp, sep=sep, keep_default_na=False)
        t_score_df = pd.read_table(t_score_matrix_fp, sep=sep, keep_default_na=False)

        for x in req_cols:
            if x not in fitness_df.columns:
                raise Exception(
                    f"Required column name {x} not found in fitness file {fitness_matrix_fp}."
                )
            if x not in t_score_df.columns:
                raise Exception(
                    f"Required column name {x} not found in t score file {t_score_matrix_fp}."
                )

        for i in range(len(fitness_df.columns)):
            if fitness_df.columns[i] != t_score_df.columns[i]:
                raise Exception(
                    "Columns don't match up (fitness, t_score):"
                    f"{fitness_df.columns[i]} != {t_score_df.columns[i]} at column {i}"
                )

        # Making sure all non numerical values are the same for both files, and locusIds are unique.
        locusIds_dict = {}
        for ix, locusId in fitness_df["locusId"].iteritems():
            if locusId != t_score_df["locusId"].iloc[ix]:
                raise Exception(
                    f"locusIds not equal at index {ix} in fitness and t score files."
                    f"{str(fitness_df['locusId'])} != {str(t_score_df['locusId'])}"
                )
            if fitness_df["sysName"].iloc[ix] != t_score_df["sysName"].iloc[ix]:
                if not (
                    pd.isnull(fitness_df["sysName"].iloc[ix])
                    and pd.isnull(fitness_df["sysName"].iloc[ix])
                ):
                    raise Exception(
                        f"sysNames not equal at index {ix} in fitness and t score files."
                        f"{str(fitness_df['sysName'])} != {str(t_score_df['sysName'])}"
                    )
            if fitness_df["geneName"].iloc[ix] != t_score_df["geneName"].iloc[ix]:
                if not (
                    pd.isnull(fitness_df["geneName"].iloc[ix])
                    and pd.isnull(fitness_df["geneName"].iloc[ix])
                ):
                    raise Exception(
                        f"geneNames not equal at index {ix} in fitness and t score files."
                        f"{str(fitness_df['geneName'])} != {str(t_score_df['geneName'])}"
                    )
            if fitness_df["desc"].iloc[ix] != t_score_df["desc"].iloc[ix]:
                if not (
                    pd.isnull(fitness_df["desc"].iloc[ix])
                    and pd.isnull(fitness_df["desc"].iloc[ix])
                ):
                    raise Exception(
                        f"descriptions not equal at index {ix} in fitness and t score files."
                        f"{str(fitness_df['desc'])} != {str(t_score_df['desc'])}"
                    )
            if locusId in locusIds_dict:
                raise Exception(f"Duplicate locusIds at index {ix}")
            else:
                locusIds_dict[locusId] = 1

        logging.info("Matrices columns are: " + ", ".join(fitness_df.columns))

        return [list(fitness_df.columns), fitness_df.shape[0]]

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

