import sys
import os
import pandas as pd
def main():
    metaphlan_dir = sys.argv[1]
    metaphlan_files = os.listdir(metaphlan_dir)
    shared_part_of_file_name = ["bowtie_unaligned_merged_", "_metaphlan_output"] # removes shared parts of file name so can just grab the unique sample name for purpose of adding to formatting
    summary_file_path = f"{metaphlan_dir}/classification_summary.tsv"
    total_classif_info_path = f"{metaphlan_dir}/total_classifications.csv"
    prepare_summary_file(summary_file_path)
    taxa_dataframes_dict = {}
    for file in metaphlan_files:
        
        output_file = f"{metaphlan_dir}/{file}"
        sample_name = file
        for shared_part in shared_part_of_file_name:

            sample_name = sample_name.replace(shared_part, "")

        file_taxa_rank_dict = run_formatting_on_metaphlan_file(output_file)
        generate_summary(file_taxa_rank_dict, sample_name, summary_file_path)

        for rank,taxa_rank_dict in file_taxa_rank_dict.items():
            

            normalised_taxa_dict = makes_abundances_equal_one_hundred(rank, taxa_rank_dict)
            normalised_taxa_dict['Sample_name'] = sample_name # adds the name of the sample to the taxa dict, allows for easy index setting with later function
            add_to_taxa_vals_to_taxa_dataframe(rank, normalised_taxa_dict, taxa_dataframes_dict)


    taxa_dataframes_dict = set_df_indexes_to_samp_name(taxa_dataframes_dict)
    produce_final_files(taxa_dataframes_dict, metaphlan_dir)
    get_total_classifications(taxa_dataframes_dict, total_classif_info_path)
    




def prepare_summary_file(summary_file_path):
    first_row = '-\tKingdom\tPhyla\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain'
    part_for_second_row ='\tNum_taxa,%_abund_classified'
    second_row = "Sample"
    for i in range(0, len(first_row), 1):
        second_row = second_row + part_for_second_row
    with open(summary_file_path, "w") as f:
        f.write(f"{first_row}\n")
        f.write(f"{second_row}\n")

    all_classifications_dict = {'k__', ''}




def run_formatting_on_metaphlan_file(metaphlan_output_file):
    metaphlan_classification_lines = get_file_lines(metaphlan_output_file)
    taxa_rank_dict = separate_classifications(metaphlan_classification_lines)

    taxa_dict_by_rank = {}
    for taxa_rank,lines in taxa_rank_dict.items():

        taxa_dict = get_classif_and_percent_abundance(lines, taxa_rank)
        taxa_dict_by_rank[taxa_rank] = taxa_dict
    
    return taxa_dict_by_rank
        

def get_file_lines(metaphlan_output_file):
    
    with open(metaphlan_output_file) as f:
        file_lines = f.readlines()
        relevant_lines = file_lines[3:]

    return relevant_lines
def separate_classifications(metaphlan_classification_lines):

    taxa_ranks = ['t__', 's__', 'g__', 'f__', 'o__', 'c__', 'p__', 'k__']
    taxa_rank_dict = {}
    for classification_line in metaphlan_classification_lines:
        for rank in taxa_ranks:
            if rank in classification_line:
                if rank not in taxa_rank_dict:
                    taxa_rank_dict[rank] = [classification_line]
                else:
                    taxa_rank_dict[rank].append(classification_line)
                break


    return taxa_rank_dict

def get_classif_and_percent_abundance(list_of_taxa_lines_for_rank, rank):
    taxa_dict = {}

    for line in list_of_taxa_lines_for_rank:
        line_split_by_tab = line.split("\t")
        classif_part = line_split_by_tab[0]
        abundance_part = float(line_split_by_tab[-2])


        highest_classif = classif_part.split("|")[-1]

        taxa_dict[highest_classif] = abundance_part
    return taxa_dict


def generate_summary(taxa_dict_by_rank, sample_name, summary_file_path):
    taxa_ranks = ['t__', 's__', 'g__', 'f__', 'o__', 'c__', 'p__', 'k__']
    sample_taxa_abundance_and_num_info = ""
    for taxa_rank,taxa_dict in taxa_dict_by_rank.items():
        
        total_abundance = 0
        num_of_taxas = 0
        
        for taxas,abundance in taxa_dict.items():
            total_abundance += abundance
            num_of_taxas += 1
        sample_taxa_abundance_and_num_info = sample_taxa_abundance_and_num_info + f"\t{num_of_taxas},{total_abundance}"
                    
    
    with open(summary_file_path, "a") as file:
        file.write(f"{sample_name}\t{sample_taxa_abundance_and_num_info}\n")


def makes_abundances_equal_one_hundred(rank, taxa_dict):
    list_of_taxas = list(taxa_dict.keys())
    total_abund = sum(taxa_dict.values())
    multiplying_factor = 100 / total_abund
    normalised_taxa_dict = {}
    for taxa,abund in taxa_dict.items():
        new_abund = abund * multiplying_factor
        normalised_taxa_dict[taxa] = new_abund

    
    return normalised_taxa_dict

def add_to_taxa_vals_to_taxa_dataframe(rank, normalised_taxa_dict, taxa_dataframes_dict):

    if rank in taxa_dataframes_dict.keys():
        
        rank_dataframe = taxa_dataframes_dict[rank]
        current_taxas_found = rank_dataframe.columns.tolist()
        taxas_in_sample = normalised_taxa_dict.keys()
        new_taxas_found = [taxa for taxa in taxas_in_sample if taxa not in current_taxas_found]
        
        for new_taxa in new_taxas_found:
            rank_dataframe[new_taxa] = 0
        taxas_not_in_sample = [taxa for taxa in current_taxas_found if taxa not in taxas_in_sample]
        
        for not_in in taxas_not_in_sample:
            normalised_taxa_dict[not_in] = 0
        sample_df = pd.DataFrame(normalised_taxa_dict, index=[0])
        rank_dataframe = pd.concat([rank_dataframe, sample_df], ignore_index=True)
      #  rank_dataframe = rank_dataframe.append(normalised_taxa_dict, ignore_index=True)

    
    else:

        rank_dataframe = pd.DataFrame(normalised_taxa_dict, index=[0])

    
    taxa_dataframes_dict[rank] = rank_dataframe
    return taxa_dataframes_dict

def set_df_indexes_to_samp_name(taxa_dataframes_dict):

    for rank,dataframe in taxa_dataframes_dict.items():
        dataframe = dataframe.set_index('Sample_name')
        print(dataframe)
        taxa_dataframes_dict[rank] = dataframe
    
    return taxa_dataframes_dict


def produce_final_files(taxa_dataframes_dict, metaphlan_dir):
    for rank,dataframe in taxa_dataframes_dict.items():
        output_file = f"{metaphlan_dir}/{rank}_final_output_table.csv"
        dataframe.to_csv(output_file, sep=',', index=True)

def get_total_classifications(taxa_dataframes_dict, total_classif_info_path):
    taxa_rank_total_classifs_dict = {}
    taxa_rank_convert_dict = {"k__": "Kingdom", "p__" : "Phylum", "c__" : "Class", "o__" : "Order", "f__" : "family", "g__" :  "genus", "s__" : "species", "t__" : "strain"}

    for taxa_rank, dataframe in taxa_dataframes_dict.items():
        print(taxa_rank)
        print(len(dataframe.columns.tolist()))
        taxa_rank_name = taxa_rank_convert_dict[taxa_rank]
        taxa_rank_total_classifs_dict[taxa_rank_name] = len(dataframe.columns.tolist())

    with open(total_classif_info_path, "w") as file:
        file.write("Taxonomic rank,Total classifications in data")
        for taxa_name, taxa_count in taxa_rank_total_classifs_dict.items():
            file.write(f"{taxa_name},{taxa_count}\n")

            
main()