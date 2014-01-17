__author__ = 'philippefaucon'
import os
import glob
import sys
import re
import sqlite3
import time

import bioservices
from Bio.Blast import NCBIXML
from Bio import SeqIO
import Bio.Blast.Applications

import BlastAggregation
import ProjectDefinitions
import AggregateToSQLite

ProjectDefinitions.make_directory_structure()


def drugbank_to_yeast_main():

    #prepare the fasta for blast searching...
    drugbank_fasta = ProjectDefinitions.data_directory + "drugbank_fasta.fasta"
    #some genes diverge, so one drugbase entry yields 2+ uniprot entries, this should be handled
    trouble_genes = ProjectDefinitions.data_directory + "unmatched_genes.txt"
    get_drugbank(drugbank_fasta,trouble_genes,override=False)

    #blast search...
    blastp_results = ProjectDefinitions.data_directory + "drugbank_blastpresults.xml"
    deltablast_results = ProjectDefinitions.data_directory + "drugbank_deltablastresults.xml"
    search_blast(fasta_file=drugbank_fasta,blastp_results=blastp_results, deltablast_results=deltablast_results)

    #reformat from BLASTXML to a tab-delimited format for easier understanding
    aggregate_blastp = ProjectDefinitions.results_directory + "aggregate_drugbank_blastpresults.txt"
    aggregate_deltablast = ProjectDefinitions.results_directory + "aggregate_drugbank_deltablastresults.txt"
    fasta_sequences = BlastAggregation.get_fasta_sequences(drugbank_fasta)
    BlastAggregation.analyzeDirectory(blastp_results,aggregate_blastp,sequences=fasta_sequences,overwrite=True)
    BlastAggregation.analyzeDirectory(deltablast_results,aggregate_deltablast,sequences=fasta_sequences,overwrite=True)

    sqlite_file = ProjectDefinitions.results_directory + "drugbank_sqlite.txt"
    AggregateToSQLite(aggregate_blastp,aggregate_deltablast, sqlite_file)


def get_drugbank(drugbank_fasta=None,trouble_genes=os.devnull, override=False):
    """
    make a fasta file from the previously generated drugbank database with only Homo Sapiens records.
    @param drugbank_fasta: primary output
    @param trouble_genes:  output file in case of difficulty with uniprot matching
    @param override: if the drugbank_fasta file exists should it be replaced?
    @return: Error code, 0 is OK
    """

    if drugbank_fasta is None or trouble_genes is None:
        print("Error: tried to call get_drugbank missing the output file")
        return -1

    u = bioservices.UniProt(verbose=False)

    if os.path.isfile(drugbank_fasta) and not override:
        print "drugbank fasta file found, not regenerating!"
    else:
        out = open(drugbank_fasta,"w")
        out_trouble = open(trouble_genes,"w")

        conn = sqlite3.connect(ProjectDefinitions.data_directory + ProjectDefinitions.drugbank_db_name)
        curs = conn.cursor()
        curs.execute("SELECT UniProt_ID FROM drugs Where Species='Homo sapiens'")
        while True:
            row = curs.fetchone()
            if row is None:
                print "done printing drug results!"
                break
            #print row[0]
            if(len(row[0])<2):
                continue
            try:
                rest = u.searchUniProtId(row[0],format='fasta')
            except Exception:
                print "ran into a problem with gene: ", row[0],". Error: ", sys.exc_info()[0]
                out_trouble.write(row[0] + "\n")
            else:
                #print "hit: ", rest
                out.write(rest)
        out.close()
        out_trouble.close()
        return 0


def search_blast(fasta_file=None,blastp_results=None,deltablast_results=None,blastpoverwrite=False,deltablastoverwrite=False):

    if(fasta_file is None):
        print("ERROR: Tried to call search_blast with no fast file...")
        return -1

    if blastp_results is not None:
        if os.path.isfile(blastp_results) and not blastpoverwrite:
            print "drugbank blastp file found, not regenerating!"
        else:
            start = time.time()
            blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(
                query=fasta_file, db="nr_yeast", outfmt=5, out=blastp_results,num_threads=16)
            print(blastp_cline)
            stdout, stderr = blastp_cline()
            print "blastp complete! total time(seconds) =" +str(time.time()-start)

    if deltablast_results is not None:
        if os.path.isfile(deltablast_results) and not deltablastoverwrite:
            print "drugbank deltablast file found, not regenerating!"
        else:
            print "running deltablast..."
            start = time.time()
            deltablast_cline = Bio.Blast.Applications.NcbideltablastCommandline(
                query=fasta_file, db="nr_yeast", outfmt=5, out=deltablast_results,num_threads=16)
            print(deltablast_cline)
            stdout, stderr = deltablast_cline()
            print "deltablast complete! total time(seconds) =" +str(time.time()-start)

    return 0


# Standard boilerplate to call the main() function to begin
# the program.  This also ensures that the main isn't called when used as a module
if __name__ == '__main__':
    drugbank_to_yeast_main()