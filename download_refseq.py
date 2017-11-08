import os
import time
import ftplib
import argparse
import multiprocessing
import pandas as pd
from accessoryFunctions import accessoryFunctions

from Bio import Entrez
from time import sleep


class RefSeqDownload(object):

    def get_download_list(self, assembly_summary_file):
        # Define the relevant fields from assembly_summary.txt
        fields = [
            'taxid',
            'ftp_path',
            'assembly_level'
        ]

        # Drop text file into pandas dataframe
        df = pd.read_csv(assembly_summary_file, delimiter='\t', low_memory=False, skiprows=1, usecols=fields)

        # Filter for requested dataset
        if self.assembly_level is not None:
            df = df[df['assembly_level'] == self.assembly_level]

        # List to store every ftp:taxonomy relationship. Used later to download + make subdirs.
        download_list = []

        # Iterates through every row of the dataframe and feeds the FTP and taxid into the get_taxonomy function.
        # ETA: ~60mins
        # Extremely slow and has no right to be, but it's good enough.
        for index, row in df.iterrows():
            try:
                to_append = self.get_taxonomy(row['ftp_path'], row['taxid'])
                download_list.append(to_append)
            except:
                print('Encountered URL error. Trying again in 1 minute.')

                # Relax for 60 seconds in the event of a timeout
                sleep(60)

                # Try again
                to_append = self.get_taxonomy(row['ftp_path'], row['taxid'])
                download_list.append(to_append)

            # NCBI doesn't like more than 3 requests per second
            sleep(0.33)

        return download_list

    # Function to fetch taxonomy by taxid via Entrez - returns a dictionary connecting the ftp link to the taxonomy
    def get_taxonomy(self, ftp, taxid):
        # efetch based on taxid
        search = Entrez.efetch(id=int(taxid),
                               db='taxonomy',
                               retmode='xml',
                               email=self.email)

        # read the search result
        search_result = Entrez.read(search)
        tax_dict = {}

        # get full taxonomy and store in dictionary
        try:
            for item in search_result[0]['LineageEx']:
                if item['Rank'] == 'no rank':
                    pass
                else:
                    tax_dict[item['Rank']] = item['ScientificName']
        except IndexError:
            pass

        tmp = {ftp: tax_dict}

        self.dump_download_list(tmp)

        return tmp

    def dump_download_list(self, download_list):
        # Drop the download list into a text file.
        download_list_file = open(os.path.join(self.output_folder, 'download_list.txt'), 'w')
        for item in download_list:
            download_list_file.write('{}\n'.format(item))

    # TODO: Implement this
    def read_download_list(self, download_list_file):
        with open(download_list_file, 'r') as f:
            download_list = [eval(line.strip()) for line in f]
        return download_list

    def download_file(self, entry):
        """
        Takes an entry from the download list and downloads it.
        """
        os.mkdir(os.path.join(self.output_folder, 'raw_files'))
        download_folder = os.path.join(self.output_folder, 'raw_files')

        # Grab FTP and taxonomy
        for ftp, taxonomy_dict in entry.items():

            # Folder structure will be created from this working directory.
            dirtree = self.output_folder

            # Simple taxonomy check
            tax_list = []
            if 'superkingdom' in taxonomy_dict:
                tax_list.append(taxonomy_dict['superkingdom'])
            if 'phylum' in taxonomy_dict:
                tax_list.append(taxonomy_dict['phylum'].replace(' ', '_'))
            if 'class' in taxonomy_dict:
                tax_list.append(taxonomy_dict['class'].replace(' ', '_'))
            if 'order' in taxonomy_dict:
                tax_list.append(taxonomy_dict['order'].replace(' ', '_'))
            if 'family' in taxonomy_dict:
                tax_list.append(taxonomy_dict['family'].replace(' ', '_'))
            if 'genus' in taxonomy_dict:
                tax_list.append(taxonomy_dict['genus'].replace(' ', '_'))
            if 'species' in taxonomy_dict:
                tax_list.append(taxonomy_dict['species'].replace(' ', '_'))
            if 'subspecies' in taxonomy_dict:
                tax_list.append(taxonomy_dict['subspecies'].replace(' ', '_'))

            # Make directory structure if it doesn't already exist
            for level in tax_list:
                dirtree = os.path.join(dirtree, level)
                try:
                    os.mkdir(dirtree)
                except OSError:
                    pass

            # Navigate to FTP directory
            self.f.cwd(ftp.replace('ftp://ftp.ncbi.nlm.nih.gov', ''))

            # Search through the directory listing for the FASTA we want
            for file in self.f.nlst():
                if file.endswith('from_genomic.fna.gz') is True:
                    pass
                elif file.endswith('_genomic.fna.gz'):
                    file_name = file

            # Get the full FTP path
            genomic_ftp = os.path.join(ftp, file_name)

            # TODO: Allow hour_start and hour_end params to be passed from cmd
            accessoryFunctions.download_file(genomic_ftp,
                                             os.path.join(download_folder, file_name),
                                             hour_start=5,  # 18
                                             hour_end=23,  # 6
                                             timeout=6000)

            # Grab filepath of fasta after download
            file_location = os.path.join(download_folder, file_name)

            # Create symlink in tax folder
            try:
                os.symlink(file_location, os.path.join(dirtree, file_name))
            except OSError:
                os.remove(os.path.join(dirtree, file_name))
                os.symlink(file_location, os.path.join(dirtree, file_name))

            # Take it easy
            sleep(0.33)
            print('Downloaded {} successfully.'.format(file_location))

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nREFSEQ DOWNLOADER' + '\033[0m')

        # Arguments
        self.args = args
        self.output_folder = args.output_folder
        self.assembly_summary = args.assembly_summary
        self.assembly_level = args.assembly_level
        self.email = args.email

        # Validation
        if self.assembly_summary is None:
            print('ERROR: Please provide the path to your assembly_summary file.')
            quit()

        if self.email is None:
            print('ERROR: Please provide an email address.')
            quit()

        # Assembly level setup
        if self.assembly_level == 'complete_genome':
            self.assembly_level = 'Complete Genome'
        elif self.assembly_level == 'all':
            self.assembly_level = None
        else:
            self.assembly_level = self.assembly_level.capitalize()

        # Retrieve a new list:
        self.download_list = self.get_download_list(self.assembly_summary)

        # FTP Setup
        url = 'ftp.ncbi.nlm.nih.gov'
        self.f = ftplib.FTP(url)
        self.f.login('anonymous')

        # Multiprocess download
        p = multiprocessing.Pool(processes=4)
        p.map(self.download_file, self.download_list)
        p.close()
        p.join()


def main():
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('output_folder',
                        help='Path to the folder you would like to download to')
    parser.add_argument('-em', '--email',
                        help='Email address required to retrieve taxonomy info from NCBI')
    parser.add_argument('-as', '--assembly_summary',
                        help='Path to assembly_summary.txt file retrieved from NCBI.'
                             'Can also pass a URL instead of a local file.')  # TODO: This isn't actually true yet
    parser.add_argument('-al', '--assembly_level',
                        help='Assembly level you would like to retrieve.'
                             'The default is complete_genome. Other options include: '
                             'assembly, contigs, chromosome, scaffold, all',
                        default='complete_genome')
    arguments = parser.parse_args()

    x = RefSeqDownload(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished RefSeqDownload functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')

if __name__ == '__main__':
    main()
