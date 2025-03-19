#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO

def run_command(command):
    """
    Run a command in the shell and return the stdout and stderr.

    Args:
        command (str): The command to be executed.

    Returns:
        tuple: A tuple containing the stdout and stderr as strings.
    """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    # Print stdout and stderr
    if stdout:
        print(f'\033[92m{stdout}\033[0m')
    if stderr:
        print(f'\033[91m{stderr}\033[0m')

    return stdout, stderr

def create_report(tabulated_data: Path, taxa_IDs: dict, output_dir: Path, aproach: str):
    """
    Create a report table containing the taxon, source, number of records and highest resolution.

    Args:
        tabulated_data (Path): The path to the tabulated data.
        taxa_IDs (dict): A dictionary containing the taxon and the corresponding sequence IDs.
        output_dir (Path): The path to the output directory.
        aproach (str): The approach used to classify the sequences.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the report table.
    """
    # Extract the data from the tabulated files
    tab_data = pd.read_csv(tabulated_data, sep='\t', skiprows=[1])

    # Crete folder for individual reports
    if not os.path.exists(f"{output_dir}/Individual_reports_{aproach}"):
        os.mkdir(f"{output_dir}/Individual_reports_{aproach}")

    # Create a report table
    taxa_levels = {'f__': 'Family', 'g__': 'Genus', 's__': 'Species'}
    report_table = pd.DataFrame(columns=['Taxon', 'Source', 'Number of records', 'Highest resolution', 'Species', 'Genus', 'Family'])

    for taxon, seqID in taxa_IDs.items():
        # Check if the taxon is in the tabulated data
        id_check = [id.replace("_NCBI", "") for id in seqID]
        if set(id_check).intersection(set(tab_data['id'].tolist())):
            # Save the matched lines in a separate file
            matched_IDs = tab_data[tab_data['id'].isin(id_check)]
            matched_IDs.to_excel(f'{output_dir}/Individual_reports_{aproach}/{taxon.replace(" ", "_")}.xlsx', index=False)

            # Get the source of the
            if any("_NCBI" in id for id in seqID):
                source = "NCBI"
            else:
                source = "Fera database"

            # Get the number of seqIds present in tab_data
            rec = tab_data[tab_data['id'].isin(id_check)].shape[0]

            # Get the highest resolution
            tab_taxa = matched_IDs['Taxon'].str.replace("; ", ";")
            tab_taxa_matches = tab_taxa[tab_taxa.apply(lambda x: any(char in x for char in taxon.split()))]

            if tab_taxa_matches.empty:
                res = 'Not found'
            else:
                top_level = taxa_levels[tab_taxa_matches.iloc[0].split(';')[-1][:3]]
                if top_level == "Species" and not any(tab_taxa_matches.apply(lambda x: all(char in x for char in taxon.split()))):
                    res = "Different species"
                else:
                    res = top_level

            # Get how many lines in tab_taxa_matches have taxa keys
            if res != 'Not found':
                prop_sp = round(tab_taxa[tab_taxa.apply(lambda x: 's__' in x)].shape[0] / tab_taxa.shape[0], 2)*100
                prop_g = round(tab_taxa[tab_taxa.apply(lambda x: 'g__' in x)].shape[0] / tab_taxa.shape[0], 2)*100
                prop_f = round(tab_taxa[tab_taxa.apply(lambda x: 'f__' in x)].shape[0] / tab_taxa.shape[0], 2)*100
            else:
                prop_sp = 'Not found'
                prop_g = 'Not found'
                prop_f = 'Not found'

            # Create the taxon row
            row = pd.DataFrame({'Taxon': taxon,
                                'Source': source,
                                'Number of records': [rec],
                                'Highest resolution': res,
                                'Species': str(prop_sp) + '%',
                                'Genus': str(prop_g) + '%',
                                'Family': str(prop_f) + '%'})
        else:
            # Create the taxon row
            row = pd.DataFrame({'Taxon': taxon,
                                'Source': 'None',
                                'Number of records': [0],
                                'Highest resolution': 'None',
                                'Species': 'None',
                                'Genus': 'None',
                                'Family': 'None'})

        # Append the row to the report table
        report_table = pd.concat([report_table, row], ignore_index=True)

        # Sort table by Source
        report_table = report_table.sort_values(by='Source', ascending=True)

    return report_table

def main():
    """
    Main function to execute the taxonomic availability and resolution check.
    """
    # Input arguments
    arg = parse_arguments()


    #################################################### Set parameters ####################################################

    if arg.config:
        df = pd.read_csv(arg.config, delimiter="\t", index_col=0)
        try:
            amplicon_row = df.loc[[arg.amplicon]]
            config_dict = amplicon_row.to_dict('records')[0]
        except KeyError:
            return None

        classifier = Path(config_dict["Classifier"])
        primer_f = config_dict["Fprimer"]
        primer_r = config_dict["Rprimer"]
        ref_seqs_path = config_dict["Database"]
        ref_taxa_path = config_dict["Taxonomy"]
    else:
        try:
            classifier = Path(arg.classifier)
        except TypeError:
            pass
        primer_f = arg.primer_f
        primer_r = arg.primer_r
        ref_seqs_path = Path(arg.ref_seqs_path)
        ref_taxa_path = Path(arg.ref_taxa_path)

    # Export taxonomy artifact and read is as pandas df
    if not os.path.exists('Reference_db/taxonomy.tsv'):
        run_command(f'qiime tools export --input-path {ref_taxa_path} --output-path Reference_db')
    else:
        print("The reference taxonomy file has already been extracted.")
    ref_taxa = pd.read_csv('Reference_db/taxonomy.tsv', sep='\t')

    # Format the taxa names
    ref_taxa['Taxon'] = ref_taxa['Taxon'].str.replace("; ?s__", " ", regex=True).str.replace("_", " ")

    # Check if arg.taxon is a file or a string
    try:
        with open(arg.taxon, 'r') as file:
            taxa = [line.strip() for line in file]
    except FileNotFoundError:
        taxa = [arg.taxon]

    # Create a directory to save the representative sequences
    if not os.path.exists('Rep_seqs'):
        os.mkdir('Rep_seqs')


    #################################################### Search taxon ####################################################

    n=0
    seq_ids_dict = {}
    for taxon in taxa:
        matched_rows = ref_taxa[ref_taxa.iloc[:, 1].str.contains(taxon)]

        if not matched_rows.empty:
            print(f"Taxon {taxon} found in the reference taxonomy.")
            # Get the Feature ID from the first column in the matched row
            seq_ids = matched_rows['Feature ID'].tolist()
            seq_ids_dict[taxon] = seq_ids

            # Save the seq_ids in a tsv file
            seq_ids_file = f"Reference_db/{taxon.replace(' ', '_')}_seq_ids.tsv"
            with open(seq_ids_file, 'w') as file:
                file.write("Feature ID\n")
                for seq_id in seq_ids:
                    file.write(f"{seq_id}\n")

            # Extract the sequences from ref_seqs and save them in fasta files
            taxon_file = f"Rep_seqs/{taxon.replace(' ', '_')}_{arg.amplicon}.qza"

            # Filter the reference sequences
            run_command(f'qiime feature-table filter-seqs \
                            --i-data {ref_seqs_path} \
                            --m-metadata-file {seq_ids_file} \
                            --o-filtered-data {taxon_file}')

            if arg.trimmed_db == 'N':
                trimmed_file = f"Rep_seqs/trimmed_{taxon.replace(' ', '_')}_{arg.amplicon}.qza"
                # Trim the sequences
                run_command(f"qiime feature-classifier extract-reads \
                                --i-sequences {taxon_file} \
                                --p-f-primer {primer_f} \
                                --p-r-primer {primer_r} \
                                --p-n-jobs {arg.thread} \
                                --p-min-length 50 \
                                --p-max-length 600 \
                                --o-reads {trimmed_file} \
                                --verbose 2>/dev/null")
                print("Sequences trimmed.")
                rep_seqs = trimmed_file
            else:
                rep_seqs = taxon_file
            n += 1

            # Merge the sequences into a single file
            if os.path.exists(rep_seqs):
                if n == 1:
                    run_command(f'cp {rep_seqs} Rep_seqs/merged-rep-seqs.qza')
                else:
                    run_command(f'qiime feature-table merge-seqs \
                                    --i-data Rep_seqs/merged-rep-seqs.qza \
                                    --i-data {rep_seqs} \
                                    --o-merged-data Rep_seqs/merged-rep-seqs.qza')

        else:
            print(f"Taxon {taxon} not found in reference taxonomy.")

            # Check if the user wants to download the sequences from NCBI
            if arg.download_all == 'N':
                download = input("Do you want to download sequences from NCBI? [Y]/N: ")
            else:
                download = arg.download_all

            if download.upper() == 'Y' or download == '':
                # Download sequences from NCBI
                if arg.amplicon == '16S':
                    amp = ' 16S'
                elif arg.amplicon == 'ITS':
                    amp = ' AND (ITS OR internal transcribed spacer)'
                elif arg.amplicon == 'COI':
                    amp = ' AND (COX1 OR CO1 OR COI)'
                elif arg.amplicon == '18S':
                    amp = ' 18S'
                else:
                    amp = arg.amplicon

                query = taxon + amp
                taxon_dir = f"{taxon.replace(' ', '_')}"
                run_command(f'qiime rescript get-ncbi-data \
                                --p-query "({query})" \
                                --p-logging-level INFO \
                                --p-n-jobs {arg.thread} \
                                --output-dir NCBIdata/{taxon_dir} \
                                --verbose 2>/dev/null')

                # Filter the sequences
                run_command(f'qiime taxa filter-seqs \
                                --i-sequences NCBIdata/{taxon_dir}/sequences.qza \
                                --i-taxonomy NCBIdata/{taxon_dir}/taxonomy.qza \
                                --p-include {taxon.split()[-1]} \
                                --o-filtered-sequences NCBIdata/{taxon_dir}/sequences_filtered.qza \
                                --verbose 2>/dev/null')

                # Trim the sequences
                rep_seqs = f"Rep_seqs/NCBI_trimmed_{taxon.replace(' ', '_')}_{arg.amplicon}.qza"
                run_command(f"qiime feature-classifier extract-reads \
                                --i-sequences NCBIdata/{taxon_dir}/sequences_filtered.qza \
                                --p-f-primer {primer_f} \
                                --p-r-primer {primer_r} \
                                --p-n-jobs {arg.thread} \
                                --p-min-length 50 \
                                --p-max-length 600 \
                                --o-reads {rep_seqs} \
                                --verbose 2>/dev/null")
                print("Sequences downloaded and trimmed.")
                n += 1

                if os.path.exists(rep_seqs):
                    # Extract the sequences from the fasta file
                    run_command(f"qiime tools export --input-path {rep_seqs} --output-path NCBIdata/{taxon_dir}")
                    ncbi_IDs = []
                    for seq_record in list(SeqIO.parse(f"NCBIdata/{taxon_dir}/dna-sequences.fasta", "fasta")):
                        # add "_NCBI" to each sequence ID
                        seq_id = seq_record.id + "_NCBI"
                        ncbi_IDs.append(seq_id)
                    seq_ids_dict[taxon] = ncbi_IDs

                    # Merge the sequences into a single file
                    if n == 1:
                        run_command(f'cp {rep_seqs} Rep_seqs/merged-rep-seqs.qza')
                    else:
                        run_command(f'qiime feature-table merge-seqs \
                                        --i-data Rep_seqs/merged-rep-seqs.qza \
                                        --i-data {rep_seqs} \
                                        --o-merged-data Rep_seqs/merged-rep-seqs.qza')

                elif not os.path.exists(rep_seqs):
                    print("The trimming resulted in no sequences.")
                    seq_ids_dict[taxon] = 'None'
                    if len(taxa) == 1:
                        return None
            else:
                if len(taxa) == 1:
                    print("Exiting the program.")
                    return None
                else:
                    continue


    #################################################### Classify sequences ####################################################

    # Classify the sequences
    if not os.path.exists('Classification'):
        os.mkdir('Classification')

    print("Classifying sequences...")
    if arg.taxa_approach == 'nb':
        run_command(f'qiime feature-classifier classify-sklearn \
                        --i-classifier {classifier} \
                        --i-reads Rep_seqs/merged-rep-seqs.qza \
                        --p-n-jobs {arg.thread} \
                        --o-classification Classification/taxonomy_nb.qza')

    elif arg.taxa_approach == 'vs':
        run_command(f'qiime feature-classifier classify-consensus-vsearch \
                        --i-query Rep_seqs/merged-rep-seqs.qza \
                        --i-reference-reads {ref_seqs_path} \
                        --i-reference-taxonomy {ref_taxa_path} \
                        --p-threads {arg.thread} \
                        --o-classification Classification/taxonomy_vs.qza \
                        --o-search-results Classification/blast')

    elif arg.taxa_approach == 'both':
        run_command(f'qiime feature-classifier classify-sklearn \
                        --i-classifier {classifier} \
                        --i-reads Rep_seqs/merged-rep-seqs.qza \
                        --p-n-jobs {arg.thread} \
                        --o-classification Classification/taxonomy_nb.qza')

        run_command(f'qiime feature-classifier classify-consensus-vsearch \
                        --i-query Rep_seqs/merged-rep-seqs.qza \
                        --i-reference-reads {ref_seqs_path} \
                        --i-reference-taxonomy {ref_taxa_path} \
                        --p-threads {arg.thread} \
                        --o-classification Classification/taxonomy_vs.qza \
                        --o-search-results Classification/blast')

    # Tabulate the data
    if arg.taxa_approach == 'both':
        approach = ['nb', 'vs']
    else:
        approach = [arg.taxa_approach]

    for app in approach:
        run_command(f'qiime metadata tabulate \
                        --m-input-file Rep_seqs/merged-rep-seqs.qza \
                        --m-input-file Classification/taxonomy_{app}.qza \
                        --o-visualization Classification/tabulated_data_{app}.qzv \
                        --verbose 2>/dev/null')

        run_command(f'qiime tools export \
                        --input-path Classification/tabulated_data_{app}.qzv \
                        --output-path Classification')

        os.rename('Classification/metadata.tsv', f'Classification/tabulated_data_{app}.tsv')

        # Remove unnecessary files
        rm_list = ['css', 'index.html', 'js', 'q2templateassets']
        for item in rm_list:
            item_path = f'Classification/{item}'
            if os.path.isdir(item_path):
                 shutil.rmtree(f'Classification/{item}')
            else:
                os.remove(f'Classification/{item}')


    #################################################### Create report table ####################################################

    # Create a report table
    if not os.path.exists('Report'):
        os.mkdir('Report')

    for app in approach:
        tabulated_data = Path(f'Classification/tabulated_data_{app}.tsv')
        report_table = create_report(tabulated_data, seq_ids_dict, Path('Report'), app)
        report_table.to_excel(f'Report/report_table_{app}.xlsx', index=False)


# Parse the command line arguments.
def parse_arguments():
    """
    Parse the command line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description='Check taxonomic availability and resolution.')
    parser.add_argument('--taxon', help='Taxon to be processed', required=True)
    parser.add_argument('--amplicon', help='Amplicon region', required=True)
    parser.add_argument('--trimmed_db', help='Has the database been trimmed?', choices=['Y', 'N'], required=True)
    parser.add_argument('--classifier', help='Classifier to be used')
    parser.add_argument('--ref_seqs_path', help='Reference sequences')
    parser.add_argument('--ref_taxa_path', help='Reference taxonomy')
    parser.add_argument('--primer_f', help='Forward primer')
    parser.add_argument('--primer_r', help='Reverse primer')
    parser.add_argument('--config', help='Configuration dictionary')
    parser.add_argument('--thread', default=round(os.cpu_count() * 0.5), help='Number of threads')
    parser.add_argument('--taxa_approach', help='Classifier to be used', default='both', choices=['nb', 'vs', 'both'])
    parser.add_argument('--download_all', help='Download all sequences from NCBI', default='N', choices=['Y', 'N'])

    return parser.parse_args()

if __name__ == '__main__':
    main()