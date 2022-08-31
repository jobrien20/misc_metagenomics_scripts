import yaml
import os
import os.path
import subprocess
import sys
from subprocess import Popen, PIPE, STDOUT
import shutil



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

                        config_dict[key] = str(value)
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
            dataset(directory_of_datasets, self.config_dict)

        else:
            list_of_datasets = os.listdir(directory_of_datasets)
            list_of_dataset_paths = ["%s/%s" % (directory_of_datasets, dset) for dset in list_of_datasets]
            
            for path in list_of_dataset_paths:
                dataset(path, self.config_dict)          


class dataset: # dataset object with fastq paths and attributes to be added etc


    def __init__(self, dataset_path, configuration_dict):
        self.dataset_path = dataset_path
        self.configuration_dict = configuration_dict
        self.initial_fastq_paths, self.fastq_ext = self.get_fastq_paths()
        self.sample_names = self.get_sample_names()
        self.metaphlan_db = self.check_and_build_dbs()

        self.run_workflow()

    def run_workflow(self):
        self.run_metaphlan_taxonomy()




    def get_sample_names(self):

        sample_paths_without_ext = [fastq.replace(f"{self.fastq_ext}", "") for fastq in self.initial_fastq_paths]

        sample_names = [path.split("/")[-1] for path in sample_paths_without_ext]
        
        if self.configuration_dict['paired_or_unpaired'] == 'paired':
            sample_names = self.get_pairs(sample_names)
       
        return sample_names

    def get_fastq_paths(self):
        if self.configuration_dict['gzip_compressed'] == 'Y':
            fastq_paths = ["%s/%s" % (self.dataset_path, file) for file in os.listdir(self.dataset_path) if file[-9:] == '.fastq.gz' or file[-6:] == '.fq.gz']
            fastq_ext = ".gz"
        else:
            fastq_paths = ["%s/%s" % (self.dataset_path, file) for file in os.listdir(self.dataset_path) if file[-6:] == '.fastq' or file[-3:] == '.fq']
            fastq_ext = ""

        if ".fastq" in fastq_paths[0]:
            fastq_ext = ".fastq" + fastq_ext
        else:
            fastq_ext = ".fq" + fastq_ext

        return fastq_paths, fastq_ext

    def get_pairs(self, list_of_fastqs):
        forward_pair = self.configuration_dict['forward_pair']
        backward_pair = self.configuration_dict['backward_pair']
        unique_fastqs_dict = {}
        for fastq in list_of_fastqs:
            if forward_pair in fastq:
                
                fastq_no_ext = fastq.replace(forward_pair, "")
            
            else:
                
                fastq_no_ext = fastq.replace(backward_pair, "")
            
            if fastq_no_ext not in unique_fastqs_dict.keys():
                
                unique_fastqs_dict[fastq_no_ext] =  fastq

            else:

                if forward_pair in fastq:

                    fastq_names = [fastq, unique_fastqs_dict[fastq_no_ext]]

                else:

                    fastq_names = [unique_fastqs_dict[fastq_no_ext], fastq]

                unique_fastqs_dict[fastq_no_ext] = fastq_names
        
        return unique_fastqs_dict

    def check_for_indices(self):
        build_file_extensions = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

        indice_count = 0
        host_indices = []
        for extent in build_file_extensions:
            for file in os.listdir(f"{self.configuration_dict['database_directory']}/metaphlan_db_directory"):
                print(file)
                if extent in file:
                    indice_count += 1
                    print("hi")
                    if indice_count == 6:
                        host_file_name = file.split("/")[-1].replace(extent, "")
                        print("found all indices")
                        return True
                    break
        return False

    def check_and_build_dbs(self):
        database_directory = self.configuration_dict['database_directory']
        try:
            os.mkdir(database_directory)
        except:
            "db dir already produced"

        metaphlan_directory = f"{database_directory}/metaphlan_db_directory"
        indice_check = self.check_for_indices()



        try:
            os.mkdir(metaphlan_directory)

        except:

            "metaphlan dir already produced"

        if indice_check == True:
            print("yep")
            return metaphlan_directory

        metaphlan_build_args = ['metaphlan', '--install', '--bowtie2db', metaphlan_directory]
        subprocess.call(metaphlan_build_args)

        indice_check = self.check_for_indices()
        if indice_check == False:
            print("metaphlan db build failed. Exitting.")
            exit()


        return metaphlan_directory

    def run_metaphlan_taxonomy(self):
        print("running metaphlan taxonomy")
        metaphlan_output_dir = f"{self.configuration_dict['output_directory']}/metaphlan_results_directory/"
        try:
            os.mkdir(metaphlan_output_dir)
        except:
            print("metaphlan output directory already produced.")
        for sample_name in self.sample_names:
            
            sample_fastq_path = f"{self.dataset_path}/{sample_name}{self.fastq_ext}"
            sample_metaphlan_path = f"{metaphlan_output_dir}/{sample_name}_metaphlan_output"
            print(sample_fastq_path)
            print(sample_metaphlan_path)

            metaphlan_args = ['metaphlan', sample_fastq_path, '--nproc', self.configuration_dict['threads'], '--input_type', 'fastq', '--bowtie2db', self.metaphlan_db, '-o', sample_metaphlan_path]
            subprocess.call(metaphlan_args)
        #metaphlan metagenome.fastq --input_type fastq -o profiled_metagenome.txt

workflow_manager(sys.argv[1])