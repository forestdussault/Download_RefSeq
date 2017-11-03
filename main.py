import os
import ftplib
import multiprocessing
import pandas as pd
from accessoryFunctions import accessoryFunctions

from Bio import Entrez
from time import sleep


def get_download_list(assembly_summary_file):
    # Define the relevant fields from assembly_summary.txt
    fields = [
        'taxid',
        'ftp_path',
        'assembly_level'
    ]

    # Drop text file into pandas dataframe
    # TODO: Change hardcoded text file to a URL to NCBI's complete assembly_summary
    df = pd.read_csv(assembly_summary_file, delimiter='\t', low_memory=False, skiprows=1, usecols=fields)

    # Filter for complete genomes. Change this if you want a different dataset.
    df = df[df['assembly_level'] == 'Complete Genome']

    # List to store every ftp:taxonomy relationship. Used later to download + make subdirs.
    download_list = []

    # Iterates through every row of the dataframe and feeds the FTP and taxid into the get_taxonomy function. ~60mins.
    # Extremely slow and has no right to be, but it's good enough.
    for index, row in df.iterrows():
        try:
            to_append = get_taxonomy(row['ftp_path'], row['taxid'])
            download_list.append(to_append)
        except:
            print('Encountered URL error. Trying again in 1 minute.')

            # Relax for 60 seconds in the event of a timeout
            sleep(60)

            # Try again
            to_append = get_taxonomy(row['ftp_path'], row['taxid'])
            download_list.append(to_append)

        # NCBI doesn't like more than 3 requests per second
        sleep(0.33)

    return download_list


# Function to fetch taxonomy by taxid via Entrez - returns a dictionary connecting the ftp link to the taxonomy
def get_taxonomy(ftp, taxid):
    # efetch based on taxid
    search = Entrez.efetch(id=int(taxid),
                           db='taxonomy',
                           retmode='xml',
                           email='forest.dussault@inspection.gc.ca')

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
    return tmp


def dump_download_list(download_list):
    # Drop the download list into a text file.
    download_list_file = open('download_list.txt', 'w')
    for item in download_list:
        download_list_file.write('{}\n'.format(item))


def read_download_list(download_list_file):
    with open(download_list_file, 'r') as f:
        download_list = [eval(line.strip()) for line in f]
    return download_list


def download_file(entry):
    """
    Takes an entry from the download list and downloads it.
    """
    # Establish connection
    url = 'ftp.ncbi.nlm.nih.gov'
    f = ftplib.FTP(url)
    f.login('anonymous')

    # Point to main downloads folder
    # TODO: Remove hardcoded path
    download_folder = '/mnt/nas/Databases/RefSeq/raw_files'

    # Grab FTP and taxonomy
    for ftp, taxonomy_dict in entry.items():

        # Folder structure will be created from this working directory.
        # TODO: Remove hardcoded path
        dirtree = '/mnt/nas/Databases/RefSeq/'

        # Simple taxonomy check
        tax_list = []
        if 'superkingdom' in taxonomy_dict:
            tax_list.append(taxonomy_dict['superkingdom'])
        if 'phylum' in taxonomy_dict:
            tax_list.append(taxonomy_dict['phylum'])
        if 'class' in taxonomy_dict:
            tax_list.append(taxonomy_dict['class'])
        if 'order' in taxonomy_dict:
            tax_list.append(taxonomy_dict['order'])
        if 'family' in taxonomy_dict:
            tax_list.append(taxonomy_dict['family'])
        if 'genus' in taxonomy_dict:
            tax_list.append(taxonomy_dict['genus'])
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
        f.cwd(ftp.replace('ftp://ftp.ncbi.nlm.nih.gov', ''))

        # Search through the directory listing for the FASTA we want
        for file in f.nlst():
            if file.endswith('from_genomic.fna.gz') is True:
                pass
            elif file.endswith('_genomic.fna.gz'):
                file_name = file

        # Get the full FTP path
        genomic_ftp = os.path.join(ftp, file_name)

        accessoryFunctions.download_file(genomic_ftp,
                                         os.path.join(download_folder, file_name),
                                         hour_start=8,
                                         hour_end=12,
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


def main():
    # Retrieve a new list:
    # download_list = get_download_list('/mnt/nas/Databases/RefSeq/assembly_summary.txt')

    # Load an old list:
    # TODO: Remove hardcoded path
    download_list = read_download_list('/mnt/nas/Databases/RefSeq/download_list.txt')

    # Multiprocess download
    p = multiprocessing.Pool(processes=4)

    p.map(download_file, download_list)

    p.close()

    p.join()

    print('Job complete.')

if __name__ == '__main__':
    main()
