
import json
import os,logging,re
import subprocess
import h5py
import uuid
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
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
            print ("Creation of the directory %s failed" % self.staging_folder)
        else:
            print ("Successfully created the directory %s " % self.staging_folder)
        # reaction = self.params['input_deck_file']

        #This is the path to the pool file
        poolfile_fp = os.path.join(self.staging_folder, staging_pool_fp_name)

        #We check correctness of pool file
        check_pool_result = self.check_pool_file(poolfile_fp) 


        #We create the data for the object:
        pool_handle = self.dfu.file_to_shock({'file_path':poolfile_fp,
            'pack':'gzip'})['handle']['hid']

        db = {  "name": "poolfile", 
                "description": "Intermediary file for RBTnSeq", 
                "orig_workspace": self.params['workspace'],
                "poolfile_handle": pool_handle}


        #To get workspace id:
        ws_id = self.params['workspace_id']

        save_object_params = {
            'id': ws_id,
            'objects': [{
                'type': 'KBaseRBTnSeq.PoolTSV',
                'data': db,
                'name': poolfile_name
            }]
        }

        # save_objects returns a list of object_infos
        dfu_oi = self.dfu.save_objects(save_object_params)[0]

        print(dfu_oi)

        return {'Name':dfu_oi[1]}


    def validate_import_poolfile_from_staging_params(self, params):
        # check for required parameters
        for p in ['staging_file_subdir_path', 'workspace_name', 'poolfile_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))


    def check_pool_file(self, poolfile_fp):
    
        #Parse pool file and check for errors
        test_vars_dict = {
                "poolfile": pool_fp,
                "report_dict": {
                    "warnings": []
                    }
        }
       
        try:
            init_pool_dict(test_vars_dict)
        except Exception:
            logging.warning("Pool file seems to have errors - " \
                    + "Please check and reupload.")
            raise Exception
    
        return 0
    
    
    def init_pool_dict(self, vars_dict):
        pool = {} # pool dict is rcbarcode to [barcode, scaffold, strand, pos]
        with open(vars_dict['poolfile'], "r") as f:
            poolfile_str = f.read()
            poolfile_lines = poolfile_str.split("\n")
            for pool_line in poolfile_lines:
                pool_line.rstrip()
                pool = check_pool_line_and_add_to_pool_dict(pool_line, pool,
                        vars_dict)
        if len(pool.keys()) == 0:
            raise Exception("No entries in pool file")
    
        return pool
    
    
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
    
    
    
    
    
    
    
    
    
