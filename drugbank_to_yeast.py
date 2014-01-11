__author__ = 'philippefaucon'
import os
import glob
import sys
import re
import sqlite3
import ProjectDefinitions
import time

import bioservices
from Bio.Blast import NCBIXML
from Bio import SeqIO
import Bio.Blast.Applications

import BlastAggregation

u = bioservices.UniProt(verbose=False)

drugbank_fasta = ProjectDefinitions.data_directory + "drugbank_fasta.fasta"
#some genes diverge, so one drugbase entry yields 2+ uniprot entries, this should be handled
trouble_genes = ProjectDefinitions.data_directory + "unmatched_genes.txt"

if os.path.isfile(drugbank_fasta):
    print "drugbank fasta file found, not regenerating!"
else:
    out = open(drugbank_fasta,"w")
    out_trouble = open(trouble_genes,"w")

    conn = sqlite3.connect(ProjectDefinitions.data_directory + ProjectDefinitions.drugbank_db_name)
    curs = conn.cursor()
    curs.execute("SELECT UniProt_ID FROM drugs WHERE Species='Homo sapiens'")
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


blastp_results = ProjectDefinitions.data_directory + "drugbank_blastpresults.xml"
aggregate_blastp = ProjectDefinitions.results_directory + "aggregate_drugbank_blastpresults.txt"
if os.path.isfile(blastp_results):
    print "drugbank blastp file found, not regenerating!"
else:
    start = time.time()
    blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(
        query=drugbank_fasta, db="nr_yeast", outfmt=5, out=blastp_results,num_threads=16)
    print(blastp_cline)
    stdout, stderr = blastp_cline()
    print "blastp complete! total time(seconds) =" +str(time.time()-start)

deltablast_results = ProjectDefinitions.data_directory + "drugbank_deltablastresults.xml"
aggregate_deltablast = ProjectDefinitions.results_directory + "aggregate_drugbank_deltablastresults.txt"
if os.path.isfile(deltablast_results):
    print "drugbank deltablast file found, not regenerating!"
else:
    print "running deltablast..."
    start = time.time()
    blastp_cline = Bio.Blast.Applications.NcbideltablastCommandLine(
        query=drugbank_fasta, db="nr_yeast", outfmt=5, out=deltablast_results,num_threads=16)
    print(blastp_cline)
    stdout, stderr = blastp_cline()
    print "deltablast complete! total time(seconds) =" +str(time.time()-start)


fasta_sequences = BlastAggregation.get_fasta_sequences(drugbank_fasta)

BlastAggregation.analyzeDirectory(blastp_results,aggregate_blastp,sequences=fasta_sequences,overwrite=True)
BlastAggregation.analyzeDirectory(deltablast_results,aggregate_deltablast,sequences=fasta_sequences,overwrite=False)


    # Standard boilerplate to call the main() function to begin
    # the program.  This also ensures that the main isn't called when used as a module
    #if __name__ == '__main__':
    #    drugbank_to_yeast_main(False)