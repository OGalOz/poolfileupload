# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import shutil
from Utils.mutantpooluploadUtilClient import mutantpooluploadUtil
from Utils.expsfileuploadUtilClient import expsfileuploadUtil
from Utils.barcodecountuploadUtilClient import barcodecountfileuploadUtil
from Utils.genetableuploadUtilClient import genetableuploadUtil
from Utils.fitnessmatrixuploadUtilClient import fitnessmatrixuploadUtil  
from Utils.modeluploadUtilClient import modeluploadUtil 
#from Utils.funcs import check_output_name
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.rbts_genome_to_genetableClient import rbts_genome_to_genetable
#END_HEADER


class poolfileupload:
    '''
    Module Name:
    poolfileupload

    Module Description:
    A KBase module: poolfileupload
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass

    def run_poolfileupload(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
            'workspace_name' (str):, 
            'workspace_id' (int): e.g. 62550, 
            'genome_ref' (str): 'A/B/C'
            'pool_file_type' (str): 'genes_table' or 'mutantpool' or 'barcodecount' or 'experiments' or 'model'
            'description' (str): Free string 
            'sep_type': 'TSV' or 'CSV'
            'protocol_type': fixed vocab
            'staging_file_names' (list<str>): list<filenames> 
            'output_names' list<str>: List<Free string> - Correlate to staging_file_names
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_poolfileupload

        params['shared_folder'] = self.shared_folder
        token = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=token)
        params['workspace_id'] =  ws.get_workspace_info({'workspace': params['workspace_name']})[0]
        params['ws_obj'] = ws
        # genetable object (converts genome to gene table)
        params['gt_obj'] = rbts_genome_to_genetable(self.callback_url)
        params['username'] = ctx['user_id']
        res_dir = os.path.join(self.shared_folder, "results")
        os.mkdir(res_dir)
        params['results_dir'] = res_dir
        #params['output_name'] = check_output_name(params['output_name'])

        # Checking basic params
        if not 'sep_type' in params:
            raise Exception("sep_type not in params.")
        elif params['sep_type'] not in ["TSV", "CSV"]:
            raise Exception(f"Did not recognize sep_type: {params['sep_type']}")


        logging.info(params)

        if 'RBTS_file_type' not in params:
            raise Exception("Did not get param RBTS_file_type")
        else:
            pft = params['RBTS_file_type']
            if pft in ['experiments', 'mutantpool', 'barcodecount', 'fitness_matrix', 'model']:
                if params['genome_ref'] == "":
                    raise Exception(f"When uploading {pft} files you must reference a genome object.")
                if "organism_scientific_name" in params and params["organism_scientific_name"] != "" and \
                        params["organism_scientific_name"] is not None and \
                        params["organism_scientific_name"] != "None":
                    logging.warning("When uploading anything besides a genes table, "
                                    "do not provide the organism's name (under Advanced Inputs)."
                                    f" Current name given: '{params['organism_scientific_name']}'."
                                    " This new scientific name will not be used.")
                if pft == 'mutantpool':
                    # requires genome_ref
                    pf_util = mutantpooluploadUtil(params)
                    result = pf_util.upload_mutantpool()
                    gene_table_fp = result["GenesTable_fp"] 
                    shutil.move(gene_table_fp, res_dir)
                elif pft == 'barcodecount':
                    if "protocol_type" not in params or params["protocol_type"] == "":
                        raise Exception("If uploading a barcodecount file, upload "
                                        "protocol type as well (under Advanced)." + \
                                        json.dumps(params))
                    if "mutantpool_ref" not in params or params["mutantpool_ref"] == "":
                        raise Exception("If uploading barcodecounts files, upload "
                                        "related mutant pool as well (under Advanced).")
                    pcf_util = barcodecountfileuploadUtil(params)
                    result = pcf_util.upload_barcodecountfile()
                elif pft == 'experiments':
                    expsf_util = expsfileuploadUtil(params)
                    result = expsf_util.upload_expsfile()
                elif pft == 'fitness_matrix':
                    num_stage = len(params['staging_file_names'])
                    if not num_stage > 2:
                        raise Exception("When uploading a fitness matrix, "
                                        " upload at least 2 files: the fitness scores "
                                        " and the T-score files. Optionally "
                                        " upload the strain fitness scores file."
                                        " The fitness score TSV file should be  "
                                        " the first one, the t-score should be  "
                                        " the second, and the strain fitness 3rd. ")
                    elif num_stage > 3:
                        raise Exception("Cannot take more than 3 files for "
                                        "this data type: Gene Fitness, T Scores"
                                        ", and Strain fitness scores.")
                    fitness_matrix_util = fitnessmatrixuploadUtil(params)
                    result = fitness_matrix_util.upload_fitnessmatrix()
                else:
                    # model
                    modelf_util  = modeluploadUtil(params)
                    result = modelf_util.upload_model()
            elif pft == "genes_table":
                if "organism_scientific_name" not in params or params["organism_scientific_name"] == "":
                    raise Exception("When uploading a genes table, you must provide the organism's scientific name (under Advanced Inputs).")
                if "genome_ref" not in params or params["genome_ref"] == "":
                    raise Exception("When uploading a genes table, you must provide a genome object reference (under Advanced).")
                gene_table_util = genetableuploadUtil(params)
                result = gene_table_util.upload_genes_table()
            else:
                raise Exception(f"Did not recognize pool_file_type {pft} for upload")

        text_message = "Finished uploading file \n"
        if pft != "barcodecount":
            text_message += "{} saved as {} on {}\n".format(result['Name'],
                        result['Type'], result['Date'])
        else:
            for pc_result in result:
                text_message += "{} saved as {} on {}\n".format(pc_result['Name'],
                        pc_result['Type'], pc_result['Date'])

        logging.info(text_message)




        # Returning file in zipped format:-------------------------------
        report_created = False
        if len(os.listdir(res_dir)) > 0:
            report_created = True
            logging.info("res dir: " + ", ".join(os.listdir(res_dir)))
            dfu = DataFileUtil(self.callback_url)
            file_zip_shock_id = dfu.file_to_shock({'file_path': res_dir,
                                                  'pack': 'zip'})['shock_id']

            dir_link = {
                    'shock_id': file_zip_shock_id, 
                   'name':  'results.zip', 
                   'label':'RBTS_UPLOAD_output_dir', 
                   'description': 'The directory of outputs from uploading' \
                    + 'RBTS table.'
            }

            report_params = {
                    'workspace_name' : params['workspace_name'],
                    'file_links':[dir_link],
                    "message": text_message
                    }

            #Returning file in zipped format:------------------------------------------------------------------
            report_util = KBaseReport(self.callback_url)
            report_info = report_util.create_extended_report(report_params)
            # ----------

        if not report_created:
            report = KBaseReport(self.callback_url)
            report_info = report.create({'report': {'objects_created':[],
                                                    'text_message': text_message},
                                                    'workspace_name': params['workspace_name']})

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_poolfileupload

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_poolfileupload return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]

