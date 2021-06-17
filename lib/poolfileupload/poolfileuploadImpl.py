# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
from Utils.poolfileuploadUtilClient import poolfileuploadUtil
from Utils.expsfileuploadUtilClient import expsfileuploadUtil
from Utils.poolcountuploadUtilClient import poolcountfileuploadUtil
from Utils.funcs import check_output_name
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
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
        params['username'] = ctx['user_id']
        params['output_name'] = check_output_name(params['output_name'])

        # Checking basic params
        if not 'sep_type' in params:
            raise Exception("sep_type not in params.")
        elif params['sep_type'] not in ["TSV", "CSV"]:
            raise Exception(f"Did not recognize sep_type: {params['sep_type']}")



        if 'pool_file_type' not in params:
            raise Exception("Did not get param pool_file_type")
        else:
            pft = params['pool_file_type']
            if pft == 'poolfile':
                pfu = poolfileuploadUtil(params)
                result = pfu.upload_poolfile()
            elif pft == 'poolcount':
                pcfu = poolcountfileuploadUtil(params)
                result = pcfu.upload_poolcountfile()
            elif pft == 'experiments':
                expsfu = expsfileuploadUtil(params)
                result = expsfu.upload_expsfile()
            else:
                raise Exception("Did not recognize pool_file_type for upload")

        text_message = "Finished uploading file \n"
        text_message += "{} saved as {} on {}\n".format(result['Name'],
                        result['Type'], result['Date'])

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

