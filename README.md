# Download_RefSeq

### Description

This CLI script will take an assembly summary report file from NCBI and
download the contents according to user specified filtering.
The downloaded files will be available in two parent folders:
1) raw_files
2) Symbolic links to the raw_files folder with subfolder hierarchy according to taxonomy

You can retrieve an assembly summary report file from here:
_https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#asmsumfiles_

### Installation Instructions

#### Dependencies

This script has dependencies on the following external programs, which must be installed and on the
system $PATH for Program to work.

- Python (v3.5): _https://www.python.org/downloads/_

This script also depends on the python modules found in `requirements.txt`.

To install these modules, use: `pip3 install -r requirements.txt`

#### Download And Installation

To download this repository, use: `git clone https://github.com/forestdussault/download_refseq.git`

### Command Line Arguments
`output_folder` (Required): Path to the folder you would like to store your downloads.

`-em, --email` (Required): Email address required to retrieve taxonomy info from NCBI.

`-as, --assembly_summary` (Required): Path to assembly summary report file retrieved from NCBI.

`-al, --assembly_level`: Assembly level you would like to retrieve. The default is complete_genome.
Options include:
- assembly
- contigs
- chromosome
- scaffold
- all

`-t, --time`: Set a timeframe for the downloader in 24h format. Defaults to 5:00PM - 6:00AM.
The script will only proceed with downloading during the specified timeframe.
i.e. --time 17 6

`-v, --verbose`: Set this flag to receive detailed output.

### Example Usage
Command to download all scaffold level entries in assembly_summary.txt.
The downloads will only proceed between 6:00PM and 1:00AM:

`/home/username/refseq_download
--email example@example.com
--assembly_summary /home/refseq_download/assembly_summary.txt
--assembly_level scaffold
--time 18 1`