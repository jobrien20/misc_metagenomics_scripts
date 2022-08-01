
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
        
        if self.configuration_dict['paired_or_unpaired'] == 'paired':
            sample_names = self.get_pairs(sample_names)
       
        return sample_names

    def get_fastq_paths(self):
        if self.configuration_dict['gzip_compressed'] == 'Y':
            fastq_paths = ["%s/%s" % (self.dataset_path, file) for file in os.listdir(self.dataset_path) if file[-9:] == '.fastq.gz' or file[-6:] == '.fq.gz']
        else:
            fastq_paths = ["%s/%s" % (self.dataset_path, file) for file in os.listdir(self.dataset_path) if file[-9:] == '.fastq' or file[-6:] == '.fq']
        return fastq_paths

    def get_pairs(self, list_of_fastqs):
        forward_pair = self.configuration_dict['forward_pair']
        backward_pair = self.configuration_dict['backward_pair']
        unique_fastqs_dict = []
        for fastq in self.sample_names:
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
        
        return unique_fastqs


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

   
    def run_trimming(self):
        trimming_directory = f"{self.configuration_dict['output_directory']}/trimmed_fastqs"
        os.mkdir(trimming_directory)
        if self.configuration_dict['paired_or_unpaired'] == 'Y':
            for sample_name,fwd_and_bck in self.sample_names.items():
                forward_sample = f"{self.dataset_path}/{fwd_and_bck[0]}"
                backward_sample = f"{self.dataset_path}/{fwd_and_bck[1]}"

                trim_galore_args = ['trim_galore', '-q', self.configuration_dict['trim_phred_quality'], self.configuration_dict['minimum_read_length'], '--trim-n', '--cores', self.configuration_dict['threads'], '--output_dir', trimming_directory, '--paired', forward_sample, backward_sample]
                subprocess.call(trim_galore_args)
        # trim_galore -q 20 --gzip --paired --length 50 --trim-n --output_dir --cores

    def run_merging(self):
        merging_directory = f"{self.configuration_dict['output_directory']}/merged_fastqs"
        os.mkdir(merging_directory)
        for sample_name,fwd_and_bck in self.sample_names.items():
                
                forward_sample = f"{self.dataset_path}/{fwd_and_bck[0]}"
                backward_sample = f"{self.dataset_path}/{fwd_and_bck[1]}"
                merged_fastq_path = f"{merging_directory}/merged_{sample_name}"

                merge_args = ['NGmerge', '-m', self.configuration_dict['minimum_ngmerge_overlap'], '-p', self.configuration_dict['perc_mismatches_allowed_in_overlap'], '-1', forward_sample, '-2', backward_sample, '-o', merged_fastq_path]
                subprocess.call(merge_args)

    def run_bowtie_alignment(self, fastq_paths):
        
        bowtie_directory = f"{self.configuration_dict['results_directory']}/bowtie_alignment_directory"
        os.mkdir(bowtie_directory)
        summary_directory = f"{bowtie_directory}/results_summaries"
        os.mkdir(summary_directory)
        
        for fastq in fastq_paths:
            
            fastq_name = fastq.split("/")[-1]
            unaligned_reads_path = f"{bowtie_directory}/bowtie_unaligned_{fastq_name}"
            sum_path = f"{summary_directory}/{fastq.split('.')[0]}_bowtie_sum.txt"
            bowtie_args = ['bowtie2', '-x', self.configuration_dict['bowtie_host_indice_directory'], '-U', fastq, '--very-sensitive', '-p', self.configuration_dict['threads'], '>', unaligned_reads_path, '2>', sum_path]
            subprocess.call(bowtie_args)

workflow_manager(sys.argv[1])