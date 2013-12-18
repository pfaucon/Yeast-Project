# this script will take a list of biopython sequences,
# it will use the gene names and try to find human homologs using a blast search
# return a list of close matches of the yeast genes in humans

#NOTE: we only support local query (BLAST+) because DELTA-BLAST is unsupported from QBLAST(webAPI)

#code requires biopython(biopython.org)
#python wrapper to make the webserver calls easier to manage
import os
import time

# as of now this function requires biopython to be installed on the system
# http://biopython.org
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline


#FIXME: This isn't global between the scripts, this could be a problem in the future
#this isn't important data but having it cached is nice, abspath strips trailing '/'
data_directory = os.path.abspath("../../Data/Cache") + "/"

def blast_main(yeast_sequences,is_dna_seq):
    
    
    # This is the parsed record, the Name and Seq fields are of use
    i = 1
    for record in yeast_sequences:
        if i > 0:

            print "name:" + record.name
            #possibly try est_human database (second parm) if we can't restrict other ways
            #result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)
            #results = result_handle.read()
            #print results
            #result_handle.close()
            #debug code
            #save_file = open("blast_results.xml","w")
            #save_file.write(results)
            #save_file.close()
            
            filename = data_directory+record.name +".fasta"
            resultname = data_directory+record.name +"_results.xml"
            output_handle = open(filename, "w")
            SeqIO.write(record,output_handle, "fasta")
            output_handle.close()
  
            sys.stdout.write("blasting... ")
            start = time.time()
            blastp_cline = NcbiblastpCommandline(query=filename, db="nr", outfmt=5, out=resultname,num_threads=4,)
            print(blastp_cline)
            stdout, stderr = blastp_cline()
            print "blasted! total time =" +str(time.time()-start)
            
            i = i-1