
import json
import os,logging,re
import subprocess
import h5py
import uuid
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from pprint import pprint
from shutil import copy
import subprocess


class poolfileuploadUtil:
    def __init__(self,params):
        self.params = params
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.data_folder = os.path.abspath('/kb/module/data/')
        #This is where files from staging area exist
        self.staging_folder = os.path.abspath('/staging/')
        self.shared_folder = params['shared_folder']
        self.scratch_folder = os.path.join(params['shared_folder'],"scratch")

    def upload_poolfile(self):

        print('params: ',self.params)
        self.validate_import_poolfile_from_staging_params(self, params)

        #Name of file in staging:
        staging_pool_fp_name = self.params['staging_file_subdir_path']

        #Output name of pool file:
        poolfile_name = self.params['poolfile_name']

        print('poolfile_name: ',poolfile_name)
        print("top dir /:",os.listdir('/'))
        print("/kb/module/:",os.listdir('/kb/module'))
        try:
            os.mkdir(self.staging_folder)
        except OSError:
            #We expect this error if staging_folder exists
            print ("Creation of the directory %s failed" % self.staging_folder)
        else:
            print ("Successfully created the directory %s " % self.staging_folder)


        #This is the path to the pool file
        poolfile_fp = os.path.join(self.staging_folder, staging_pool_fp_name)

        #We check correctness of pool file
        column_header_list = self.check_pool_file(poolfile_fp) 
        if len(column_header_list) != 12:
            print("WARNING: Number of columns is not 12 as expected: {}".format(
                len(column_header_list)))

        #We create the data for the object:
        file_to_shock_result = self.dfu.file_to_shock({'file_path':poolfile_fp,
            'make_handle': True,
            'pack':'gzip'})
        res_handle = file_to_shock_result['handle']

        pool_data = {
        'file_type' : 'KBaseRBTnSeq.PoolTSV',
        'handle_id' : res_handle['hid'],
        'shock_url' : res_handle['url'],
        'shock_node_id' : res_handle['id'],
        'shock_type' : res_handle['type'], #should be shock
        'compression_type' : "gzip",
        'file_name' : res_handle['file_name'],
        'run_method' : params['run_method'],
        'related_genome_ref' : params['genome_ref'],
        'related_organism_name' : get_genome_organism_name(
            params['genome_ref']),
        'description' : params['description']
        }


        #To get workspace id:
        ws_id = self.params['workspace_id']

        save_object_params = {
            'id': ws_id,
            'objects': [{
                'type': 'KBaseRBTnSeq.PoolTSV',
                'data': pool_data,
                'name': poolfile_name
            }]
        }

        # save_objects returns a list of object_infos
        dfu_object_info = self.dfu.save_objects(save_object_params)[0]

        print('dfu_object_info: ')
        print(dfu_object_info)
        return {'Name':dfu_object_info[1], 'Type': dfu_object_info[2], 
                'date' dfu_object_info[3]}


    def validate_import_poolfile_from_staging_params(self, params):
        # check for required parameters
        for p in ['staging_file_subdir_path', 'workspace_name', 'poolfile_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))


    def check_pool_file(self, poolfile_fp):
        col_header_list = []
        #Parse pool file and check for errors
        test_vars_dict = {
                "poolfile": pool_fp,
                "report_dict": {
                    "warnings": []
                    }
        }
       
        try:
            col_header_list = init_pool_dict(test_vars_dict)
        except Exception:
            logging.warning("Pool file seems to have errors - " \
                    + "Please check and reupload.")
            raise Exception
    
        return col_header_list
    
    
    def init_pool_dict(self, vars_dict):
        pool = {} # pool dict is rcbarcode to [barcode, scaffold, strand, pos]
        with open(vars_dict['poolfile'], "r") as f:
            poolfile_str = f.read()
            poolfile_lines = poolfile_str.split("\n")
            column_header_list = [x.strip() for x in poolfile_lines[0].split("\t")]
            for pool_line in poolfile_lines:
                pool_line.rstrip()
                pool = check_pool_line_and_add_to_pool_dict(pool_line, pool,
                        vars_dict)
        if len(pool.keys()) == 0:
            raise Exception("No entries in pool file")
    
        return column_header_list
    
    
    def check_pool_line_and_add_to_pool_dict(self, pool_line, pool, vars_dict):
    
        #We get first 7 columns of pool_line (out of 12)
        split_pool_line = pool_line.split("\t")[:7]
    
        #We remove spaces:
        for x in split_pool_line:
            x.replace(' ', '')
    
        if len(split_pool_line) >= 7:
            #We unpack
            barcode, rcbarcode, undef_1, undef_2, scaffold, strand, pos = split_pool_line
        else:
            warning_text = "pool file line with less than 7 tabs:\n{}".format(
                pool_line)
            vars_dict["report_dict"]["warnings"].append(warning_text)
            logging.warning(warning_text)
            barcode = "barcode"
    
        if barcode == "barcode": #Header line
            pass
        else:
            if not re.search(r"^[ACGT]+$", barcode):
                logging.debug(len(barcode))
                raise Exception("Invalid barcode: |{}|".format(barcode))
            if not re.search( r"^[ACGT]+$",rcbarcode ):
                raise Exception("Invalid rcbarcode: |{}|".format(rcbarcode))
            if not (pos == "" or re.search( r"^\d+$", pos)):
                raise Exception("Invalid position: |{}|".format(pos))
            if not (strand == "+" or strand == "-" or strand == ""):
                raise Exception("Invalid strand: |{}|".format(strand))
            if rcbarcode in pool:
                raise Exception("Duplicate rcbarcode.")
            pool[rcbarcode] = [barcode, scaffold, strand, pos]
        return pool
    
    def get_genome_organism_name(self, genome_ref):
        #Getting the organism name using WorkspaceClient 
        ws = Workspace(self.callback_url)
        res = ws.get_objects2({
                    "objects": [
                        "ref": genome_ref,
                        "included": ["scientific_name"]
                    ]
                })
        scientific_name = res['data'][0]['data']['scientific_name']
        return scientific_name 

        
    
    
    
    
    
    
    
