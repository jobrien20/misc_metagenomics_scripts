
import yaml
import os
import subprocess
import sys
from subprocess import Popen, PIPE, STDOUT



class workflow_manager:
    required_configurations = ['directory_of_datasets', 'single_or_multiple_datasets']


    def __init__(self, configuration_file_path):
        self.configuration_file_path = configuration_file_path
        self.config_dict = self.get_configurations()
        self.run_datasets()


    def get_configurations(self): # loads yaml file and converts it into one single dictionary with list of configurations unnested
        with open(self.configuration_file_path, "r") as file:
            yaml_dict = yaml.load(file, Loader=yaml.FullLoader)
            config_dict = {}
            list_of_dicts = [yaml_dict]

            for a_dict in list_of_dicts:
                for key,value in a_dict.items():
                    
                    if isinstance(value, dict) == True:
                        
                        list_of_dicts.append(value)
                    
                    else:

                        config_dict[key] = value
        print(config_dict)

        return config_dict

    def verify_mandatory_configurations(self):
        print("to do")

    def initialise_tools(self):
        
        workflow_tools = workflow_tools(self.configuration_dict)
        summary_tools = summary_tools(self.configuration_dict)

        return workflow_tools, summary_tools

    def run_datasets(self):
        
        directory_of_datasets = self.config_dict['directory_of_datasets']
        if self.config_dict['single_or_multiple_datasets'] == 'single':
            dataset(directory_of_datasets)

        else:
            list_of_datasets = os.listdir(directory_of_datasets)
            list_of_dataset_paths = ["%s/%s" % (directory_of_datasets, dset) for dset in list_of_datasets]
            
            for path in list_of_dataset_paths:
                dataset(path, self.config_dict)          



class dataset: # dataset object with fastq paths and attributes to be added etc.

    def __init__(self, dataset_path, configuration_dict):
        self.dataset_path = dataset_path
        self.initial_fastq_paths = self.get_fastq_paths()
        self.sample_names = self.get_sample_names()


        self.run_workflow()
    
    def run_workflow(self):

        self.run_fastqc_and_multiqc(self.initial_fastq_paths)




    def get_sample_names(self):
        sample_paths_without_ext = [".".join(fastq.split(".")[:-1]) for fastq in self.initial_fastq_paths]
        sample_names = [path.split("/")[-1] for path in sample_paths_without_ext]
        return sample_names

    def get_fastq_paths(self):
        fastq_paths = ["%s/%s" % (self.dataset_path, file) for file in os.listdir(self.dataset_path) if file[-6:] == '.fastq' or file[-3:] == '.fq']
        return fastq_paths

    def run_fastqc_and_multiqc(self, fastq_files):
        initial_fastqc_directory = f"{self.configuration_dict['output_directory']}/initial_fastqc_results"
        initial_multiqc_directory = f"{self.configuration_dict['output_directory']}initial_multiqc_results"
        
        os.mkdir(initial_fastqc_directory)
        os.mkdir(initial_fastqc_directory)

        for fastq in fastq_files:
            fastqc_args = ['fastqc', fastq, '-threads', self.configuration_dict['threads'], '-outdir', initial_fastqc_directory]
            subprocess.call(fastqc_args)
        
        multiqc_args = ['multiqc', initial_fastqc_directory, '-o', initial_multiqc_directory]

        subprocess.call(multiqc_args)

   
    def run_trimming(self, fastq_files):
        trimming_directory = f"{self.configuration_dict['output_directory']}/trimmed_fastqs"
        os.mkdir(trimming_directory)
        trim_galore_args = ['trim_galore', ]



workflow_manager(sys.argv[1])