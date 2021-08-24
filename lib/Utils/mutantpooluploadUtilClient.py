import os
import logging
import shutil
import datetime
import subprocess
import pandas as pd
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from Utils.funcs import catch_NaN, DownloadGenomeToFNA


class mutantpooluploadUtil:
    def __init__(self, params):
        self.params = params
        self.callback_url = os.environ["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.data_folder = os.path.abspath("/kb/module/data/")
        # This is where files from staging area exist
        self.staging_folder = os.path.abspath("/staging/")
        self.shared_folder = params["shared_folder"]
        self.scratch_folder = os.path.join(params["shared_folder"], "scratch")

    def upload_mutantpool(self):
        """
        The upload method

        We perform a number of steps:
        Get name of mutantpool as it is in staging.
        Find the mutantpool in /staging/mutantpool_name
        Get the output name for the mutantpool
        Get the column headers for the mutant pool for
            data and testing purposes. Should be len 12.
        Test if mutantpool is well-formed.
        We send the file to shock using dfu.
        We get the handle and save the object with all
            the necessary information- including related genome.
        params should include:
            output_names,
            staging_file_names,
            ws_obj,
            workspace_id,
            gt_obj (genetable object)

        """
        print("params: ", self.params)
        self.validate_import_mutantpool_from_staging_params()

        # Name of file in staging:
        stg_fs = self.params["staging_file_names"]
        if len(stg_fs) != 1:
            raise Exception("Expecting a single mutant pool, got a different number"
                            f" of staging files: {len(stg_fs)}. Files: " + \
                            ", ".join(stg_fs))
        else:
            staging_pool_fp_name = stg_fs[0]

        op_nms = self.params["output_names"]
        if len(op_nms) != 1:
            raise Exception("Expecting a single output name, got a different number"
                            f": {len(op_nms)}. Output Names: " + \
                            ", ".join(op_nms))
        else:
            mutantpool_name = op_nms[0]


        logging.info("mutantpool_name: ", mutantpool_name)
        logging.info("top dir /:", os.listdir("/"))
        logging.info("/kb/module/:", os.listdir("/kb/module"))
        if not os.path.exists(self.staging_folder):
            raise Exception("Staging dir does not exist yet! Error will be thrown")
        else:
            logging.info("Succesfully recognized staging directory")
        # This is the path to the mutant pool
        mutantpool_fp = os.path.join(self.staging_folder, staging_pool_fp_name)

        # CHECK mutant pool:
        column_header_list, num_lines, pool_df, scf_names = self.check_mutant_pool(mutantpool_fp,
                                                             self.params["sep_type"])

        if len(column_header_list) != 12:
            print(
                "WARNING: Number of columns is not 12 as expected: {}".format(
                    len(column_header_list)
                )
            )

        self.check_scaffold_names_match(scf_names, self.params["genome_ref"])

        # We copy the file from staging to scratch
        new_pool_fp = os.path.join(self.shared_folder, mutantpool_name)

        logging.info("Creating Genes Table from Genome:")
        genes_table_fp = self.get_genes_table_from_genome_ref(self.params["genome_ref"])
        Stats_op_fp = os.path.join(self.shared_folder, "PoolStatsOutput.txt")
        logging.info(os.listdir('/kb/module'))
        is_dir = os.path.isdir('/kb/module/lib') 
        logging.info(is_dir)
        if is_dir:
            logging.info(os.listdir('/kb/module/lib'))
            logging.info(os.listdir('/kb/module/lib/Utils'))
        else:
            raise Exception("lib not a dir!")
        
        PoolStats_R_fp = '/kb/module/lib/Utils/PoolStats.R'
        is_file = os.path.isfile(PoolStats_R_fp)
        if not is_file:
            raise Exception("PoolStats R file not found!")

        self.run_poolstats_r(Stats_op_fp,
                           PoolStats_R_fp,
                           genes_table_fp,
                           mutantpool_fp,
                           str(num_lines),
                           self.shared_folder)

        gene_hit_frac = self.get_gene_hit_rate(Stats_op_fp)
        if isinstance(gene_hit_frac, float):
            gene_hit_frac = str(round(gene_hit_frac, 4))
            # otherwise it is string already


        if self.params["sep_type"] == "TSV":
            shutil.copyfile(mutantpool_fp, new_pool_fp)
        else:
            pool_df.to_csv(new_pool_fp, sep="\t", index=False)

        mutantpool_fp = new_pool_fp
        # We create the handle for the object:
        file_to_shock_result = self.dfu.file_to_shock(
            {"file_path": mutantpool_fp, "make_handle": True, "pack": "gzip"}
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
            "file_type": "KBaseRBTnSeq.RBTS_MutantPool",
            "mutantpool": res_handle["hid"],
            # below should be shock
            "handle_type": res_handle["type"],
            "shock_url": res_handle["url"],
            "shock_node_id": res_handle["id"],
            "compression_type": "gzip",
            "column_header_list": column_header_list,
            "column_headers_str": ", ".join(column_header_list),
            "num_lines": str(num_lines),
            "gene_hit_frac": gene_hit_frac,
            "fastqs_used": fastq_refs,
            "fastqs_used_str": "NA",
            "file_name": res_handle["file_name"],
            "utc_created": str(date_time),
            "related_genome_ref": self.params["genome_ref"],
            "tnseq_model_name": self.params["tnseq_model_name"],
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
                    "type": "KBaseRBTnSeq.RBTS_MutantPool",
                    "data": pool_data,
                    "name": mutantpool_name,
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
            "GenesTable_fp": genes_table_fp
        }

    def validate_import_mutantpool_from_staging_params(self):
        prms = self.params
        # check for required parameters
        for p in [
            "username",
            "staging_file_names",
            "genome_ref",
            "description",
            "output_names"
        ]:
            if p not in prms:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        if "tnseq_model_name" not in prms or prms["tnseq_model_name"] == "" \
            or prms["tnseq_model_name"] is None:
            raise Exception("When uploading a mutant pool, you must include a model name.")

    def check_mutant_pool(self, mutantpool_fp, separator):
        """
        We check the mutant pool by initializing into dict format


        Returns:
            [list(pool_df.columns), pool_df.shape[0], pool_df, unique_scaffold_names]
            which is:
                column names (str)
                number of lines
                the dataframe
                the scaffold names (list<str>)

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

        pool_df = pd.read_table(mutantpool_fp, sep=sep)

        req_cols = ["barcode", "rcbarcode", "scaffold", 
                    "pos", "strand"] 

        for x in req_cols:
            if x not in pool_df.columns:
                raise Exception(f"Required column name {x} not found in mutant pool.")

        past_end_rows = []
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
                    past_end_rows.append(ix)
        logging.info("Found pastEnd hit at rows: " + ", ".join([str(x) for x in past_end_rows]))


        for ix, pos in pool_df["pos"].iteritems():
            if not pd.isna(pos):
                if pos < 0:
                    raise Exception("Positions must be positive."
                                    f" Value at row {ix} is {pos}")

        if "pos2" in pool_df.columns:
            for ix, pos in pool_df["pos2"].iteritems():
                if not pd.isna(pos):
                    if pos < 0:
                        raise Exception("Positions must be positive."
                                        f" Value at row {ix} is {pos}")

        # Getting all scaffold names:
        unique_scaffold_names = list(pool_df["scaffold"].unique())


        logging.info("Mutant Pool columns are: " + ", ".join(pool_df.columns))

        return [list(pool_df.columns), pool_df.shape[0], pool_df, unique_scaffold_names]


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

    
    def run_poolstats_r(self, 
                      Stats_op_fp,
                      PoolStats_R_fp,
                      genes_table_fp,
                      pool_fp,
                      nMapped,
                      tmp_dir):
        """
        Description:
            We run an R script to get statistics regarding pool
            Writes to Stats_op_fp, that's the important file.
        Args:
            Stats_op_fp: (str) Path to R log
            PoolStats_R_fp (str) Path to R script 'PoolStats.R'
            pool_fp (str) (pool_fp) finished pool file
            genes_table_fp (str) genes table file path
            nMapped (str): string of an int
            tmp_dir: (str) Path to tmp_dir
        Returns:
            None
        """

        R_executable = "Rscript"

        RCmds = [R_executable, 
                PoolStats_R_fp, 
                pool_fp,
                genes_table_fp,
                nMapped]


        logging.info("Running R PoolStats")

        std_op = os.path.join(tmp_dir, "R_STD_OP.txt")
        with open(Stats_op_fp, "w") as f:
            with open(std_op, "w") as g:
                subprocess.call(RCmds, stderr=f, stdout=g)

        if os.path.exists(Stats_op_fp):
            logging.info("Succesfully ran R, log found at " + Stats_op_fp)

        return None

    def get_gene_hit_rate(self, R_log_fp):
        """
        Inputs:
            R_log_fp: (str)
        Outputs:
            gene_hit_frac (float or str ("NaN")): Fraction of non-essential genes hits
        Description:
            We take the Standard Error output of PoolStats.R and parse it
            in a crude manner in order to get fraction of non-essential genes
            hit.
        """
        with open(R_log_fp, "r") as f:
            rlog_str = f.read()
    
        res_d = {"failed": False}
    
        rlog_l = rlog_str.split("\n")
    
        logging.info("Rlog results: \n" + rlog_str)
        if rlog_str == '':
            res_d["failed"] = True
            logging.critical('PoolStats failed to create any results at ' + R_log_fp)
            return "NaN"
        elif len(rlog_l) < 11:
            res_d["failed"] = True

            logging.critical('PoolStats output does not have 11 lines as expected:\n"' + rlog_str)
            return "NaN"
    
        
    
        gene_hit_frac =  catch_NaN(rlog_l[3].split(' ')[8])
        #res_d['insertions'] = catch_NaN(rlog_l[0].split(' ')[0])
        #res_d['diff_loc'] = catch_NaN(rlog_l[0].split(' ')[6])
        #res_d['cntrl_ins'] = catch_NaN(rlog_l[1].split(' ')[1])
        #res_d['cntrl_distinct'] = catch_NaN(rlog_l[1].split(' ')[3][1:])
        #res_d['nPrtn_cntrl'] = catch_NaN(rlog_l[2].split(' ')[4])
        #res_d["essential_hit_rate"] = catch_NaN(rlog_l[3].split(' ')[6])
        #res_d["num_surp"] = catch_NaN(rlog_l[5].split(' ')[1])
        #res_d['stn_per_prtn_median'] = catch_NaN(rlog_l[7].split(' ')[5])
        #res_d['stn_per_prtn_mean'] = catch_NaN(rlog_l[7].split(' ')[-1])
        #res_d['gene_trspsn_same_prcnt'] = catch_NaN(rlog_l[8].split(' ')[-1][:-1])
        #res_d['reads_per_prtn_median'] = catch_NaN(rlog_l[9].split(' ')[5])
        #res_d['reads_per_prtn_mean'] = catch_NaN(rlog_l[9].split(' ')[7])
        #res_d['reads_per_mil_prtn_median'] = catch_NaN(rlog_l[10].split(' ')[-3])
        #res_d['reads_per_mil_prtn_mean'] = catch_NaN(rlog_l[10].split(' ')[-1])
    
        #logging.info("Results from parsing PoolStats.R output:")
        #logging.info(res_d)
    
        return gene_hit_frac 

    def get_genes_table_from_genome_ref(self, genome_ref):
        """
        Description:
            We use the installed client to get the genes table
            from the genome_ref.
        Args:
            self.params:
                gt_obj gene table object
            genome_ref (str): ref of genome 'A/B/C'
        """
        
        g2gt_results = self.params['gt_obj'].genome_to_genetable({'genome_ref': genome_ref})
        logging.info(g2gt_results)
        # This is where the genes table is downloaded to
        gene_table_fp = os.path.join(os.path.join(self.shared_folder, 'g2gt_results'), 'genes.GC')
         

        return gene_table_fp
       
    def check_scaffold_names_match(self, pool_scf_names, genome_ref):
        """Checks scaffold names from mutant pool are all in genome fna
   
        Args:
            pool_scf_names (list<str>): A list of all scaffold names
                                    from the mutant pool file.
        Returns: None
        Description:
            Every scaffold in the mutant pool has to be listed in
            the genome fna, otherwise there is a mismatch of
            scaffold names.
            Downloads genome fna, gets all scaffold names,
            checks every name in scf_names is in the scaffold
            names of the genome fna, otherwise raises an
            Exception.
        """

        #Downloading genome fna
        genome_fna_fp = DownloadGenomeToFNA(self.gfu, genome_ref, self.shared_folder)
        # scaffold name to length of scaffold dict
        scf_to_len_d = GetScaffoldLengths(genome_fna_fp)
        fna_scf_names = list(scf_to_len_d.keys())
        for scf in pool_scf_names:
            if scf not in fna_scf_names:
                raise Exception(f"Scaffold name {scf} not found in genome fna. Full list "
                                 "of scaffold names from genome fna: "
                                 ", ".join(fna_scf_names))
        


        return None





