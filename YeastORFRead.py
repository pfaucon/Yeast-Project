# Return either the set of yeast genes in DNA or protein format
# If we don't have them locally then go to yeastgenome.org and pull them
# NOTE: Possibly look into NCBI if yeastgenome gives trouble
#       http://www.ncbi.nlm.nih.gov/nuccore/nc_001144.5

import urllib
import gzip
import os #for abs path
import sys #for normal print

#do a bunch of prints to ensure things are working, it is very verbose
DEBUG = False

#this isn't important data but having it cached is nice, abspath strips trailing '/'
data_directory = os.path.abspath("../../Data/Cache") + "/"

# as of now this function requires biopython to be installed on the system
# http://biopython.org
from Bio import SeqIO

def yeast_main(is_dna_seq):
    yeast_sequences = getYeastData(is_dna_seq)
    
    # This is the parsed record, the Name and Seq fields are of use
    return yeast_sequences

def getYeastData(is_dna_seq):
    
    
    #do we want to do nucleotide or amino acid analysis ?
    if is_dna_seq:
        zipped_filename = data_directory + "yeast_genome_dna.gz"
        filename = data_directory + "yeast_genome_dna.fasta"
        source_url = "http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz"
        individual_root_dir = data_directory + "yeast_dna/"
    else:
        zipped_filename = data_directory + "yeast_genome_proteins.gz"
        filename = data_directory + "yeast_genome_proteins.fasta"
        source_url = "http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz"
        individual_root_dir = data_directory + "yeast_protein/"

    #Download the data if it isn't present
    try:
        if DEBUG:
            print "trying to open yeast data file at: " + filename
        f_ptr = open(filename, "r")
        if DEBUG:
            print "File opened successfully!"
        #file_contents = f_ptr.read()
        f_ptr.close()
    except IOError:
        
        if not os.path.exists(data_directory):
            sys.stdout.write("WARNING: cache directories were not found, creating...  ")
            os.makedirs(data_directory)
            if not os.path.exists(individual_root_dir):
                os.makedirs(individual_root_dir)
            print "OK!"
        
        print "WARNING: Couldn't find the file, downloading the compressed version..."
        print "source: " + source_url
        print "destination: " + zipped_filename

        #download from previously known good loction
        urllib.urlretrieve(source_url, zipped_filename);

        if DEBUG:
            sys.stdout.write("Download completed, decompressing...")

        #unzip the results
        f = gzip.open(zipped_filename, 'rb');
        fileout = open(filename, 'wb');
        fileout.write(f.read())
        
        fileout.close();
        f.close()
        if DEBUG:
            print "file decompressed";



    # for some reason biopython wants to give fasta files to blast
    # they need to be divided for paralellization though


    yeast_sequences = extractORFNames(filename)


    # return the full contents of the fasta file
    return yeast_sequences;

# as of now this function requires biopython to be installed on the system
# http://biopython.org
def extractORFNames(filename):
    sys.stdout.write("parsing record... ")
    #Create a BioPython seq object
    records = list(SeqIO.parse(filename, "fasta"))
    print "parse complete, " + str(len(records)) + " records found!"
    if DEBUG:
        print "dumping one record of " + str(len(records)) + "..."
        print records[0]
        print "sequence: " + records[0].seq
    
    return records

# Standard boilerplate to call the main() function to begin
# the program.  This also ensures that the main isn't called when used as a module
if __name__ == '__main__':
    yeast_main(False)