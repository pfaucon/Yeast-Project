# this script will take a list of biopython sequences,
# it will use the gene names and try to find human homologs using a blast search
# return a list of close matches of the yeast genes in humans

#NOTE: we only support local query (BLAST+) because DELTA-BLAST is unsupported from QBLAST(webAPI)

#http://www.biostars.org/p/73892/
#How to set up an alias to search only humans: http://www.biostars.org/p/6528/

#code requires biopython(biopython.org)
#python wrapper to make the webserver calls easier to manage
import os
import time
import sys

# as of now this function requires biopython to be installed on the system
# http://biopython.org
from Bio import SeqIO
import Bio.Blast.Applications
from Bio.Blast.Applications import NcbiblastpCommandline
#from Bio.Blast.Applications import NcbideltablastCommandline

import ProjectDefinitions

def blast_main(yeast_sequences,is_dna_seq):
    
    #ensure that our directories exist
    configureCachedirs()
    
    initialTime = time.time()
    
    # This is the parsed record, the Name and Seq fields are of use
    i = 6000
    for record in yeast_sequences:
        if i > 0:

            print "name:" + record.name
            
            filename = ProjectDefinitions.fastadir + record.name +".fasta"
            
            #write the fasta individually
            output_handle = open(filename, "w")
            SeqIO.write(record,output_handle, "fasta")
            output_handle.close()
            
            
            #do a blastp run
            resultname = ProjectDefinitions.blastpdir + record.name + "_blastpresults.xml"
            if os.path.isfile(resultname):
                print "blastp cache record found, skipping blastp!"
            else:
                start = time.time()
                blastp_cline = NcbiblastpCommandline(query=filename, db="nr_humans", outfmt=5, out=resultname,num_threads=16)
                print(blastp_cline)
                stdout, stderr = blastp_cline()
                print "blastp complete! total time(seconds) =" +str(time.time()-start)
            
            
            #do a deltablast run
            resultname = ProjectDefinitions.deltablastdir + record.name + "_deltablastresults.xml"
            if os.path.isfile(resultname):
                print "delta blast cache record found, skipping delta blast!"
            else:
                start = time.time()
                deltablast_cline = Bio.Blast.Applications.NcbideltablastCommandline(query=filename, db="nr_humans", outfmt=5, out=resultname,num_threads=16)
                print(deltablast_cline)
                stdout, stderr = deltablast_cline()
                print "delta blast complete! total time(seconds) =" +str(time.time()-start)
            
            i = i-1

    print "total time(seconds) = " + str(time.time()-initialTime) + ". Time in hours: " + str((time.time()-initialTime)/360)

def configureCachedirs():
    if not os.path.exists(ProjectDefinitions.fastadir):
        os.makedirs(ProjectDefinitions.fastadir)
    
    if not os.path.exists(ProjectDefinitions.blastpdir):
        os.makedirs(ProjectDefinitions.blastpdir)

    if not os.path.exists(ProjectDefinitions.deltablastdir):
        os.makedirs(ProjectDefinitions.deltablastdir)