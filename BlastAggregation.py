#given the results from blastp and deltablast do some basic analytics on the results

from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
import glob
import sys
import re


#from http://stackoverflow.com/questions/5686211/are-there-any-function-that-can-calculate-score-for-the-aligned-sequences-give-p
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from itertools import izip
blosum.update(((b,a),val) for (a,b),val in blosum.items())

import ProjectDefinitions


def blast_aggregation_main():

    seq = get_fasta_sequences()

    blastp_file_path = ProjectDefinitions.blastpdir+"*_blastpresults.xml"
    aggregate_blastp = ProjectDefinitions.data_directory + "aggregate_blastp.txt"
    analyzeDirectory(blastp_file_path,aggregate_blastp, sequences=seq, overwrite=True)

    deltablast_file_path = ProjectDefinitions.deltablastdir+"*_deltablastresults.xml"
    aggregate_deltablast = ProjectDefinitions.data_directory + "aggregate_deltablast.txt"
    analyzeDirectory(deltablast_file_path,aggregate_deltablast, sequences=seq, overwrite=True)


    return [aggregate_blastp, aggregate_deltablast]


# given a directory full of results extract all of the relevant information
# and load it into an aggregate file
def analyzeDirectory(directory,outfilename, sequences, overwrite=False):


    #if the file exists and we don't want to overwrite it then just quit out
    if(os.path.exists(outfilename) and not overwrite):
        print "INFO: file ", outfilename, " exists and overwrite is disabled, skipping..."
        return

    print "loading all sequences in ",directory, "..."
    outhandle=open(outfilename,"w")
    outhandle.write("sequence name\tsequence_length\tmatch_db\tmatch_name\tmatch_length\tpct_cover\tscore_frac\te_value\n")
    for file in glob.glob(directory):

        if(os.path.getsize(file) <= 1):
            print("File: " + file + " Appears to be  empty, skipping...")
            continue
        print("now working on:" + file)
        result_handle=open(file)
        blast_record=NCBIXML.read(result_handle)

        nalignments=1
        for alignment in blast_record.alignments:
            if nalignments<=0:
                continue
            nalignments -= 1

            nhsps=1
            for hsp in alignment.hsps:

                #print "title: ", alignment.title, " hit_id: ", alignment.hit_id, " definition:", alignment.hit_def
                if nhsps<=0:
                    continue
                nhsps -= 1
                #information regarding hsp structure present in Bio.Blast.Records.py

                query_tokens = blast_record.query.split()
                sequence_name = query_tokens[0]
                #fasta sequences seem to end with '*', just strip that character
                query_sequence = sequences[sequence_name][:-1]
                #calculate the maximum possible score for this sequence, gap cost is arbitrary
                max_score = sum(score_pairwise(query_sequence,query_sequence,blosum,-1,-1))

                #if we have less than 50% of the points it isn't relevant
                if ((hsp.score*100)/max_score) < 50:
                    continue

                match_length = hsp.query_end-hsp.query_start
                #if we have less than 80% coverage then it isn't relevant
                if ((match_length*100.0)/blast_record.query_letters)<80:
                    continue

                match_db = get_match_db(alignment.accession)

                outhandle.write(blast_record.query +"\t"+ str(blast_record.query_letters) +"\t")
                outhandle.write(match_db +"\t"+ alignment.accession +"\t"+ str(match_length) +"\t"+ str((match_length*100.0)/blast_record.query_letters) +"\t")
                outhandle.write(str((hsp.score*100)/max_score) +"\t"+ str(hsp.expect))
                outhandle.write("\n")
        result_handle.close()

        #for full results comment this line out
        #break

    outhandle.close()

def get_match_db(accession_name):
    """
    This function will find the database that an accession name (probably) came from.
    This uses a very primitive not-from-documentation technique so functionality is limited.
    @param accession_name: string from blast, it should not be preprocessed in any way
    @return: database for an accession number "refseq","genbank","pdb","uniprot","???"
    """
    #refseq seems to be NP_ or XP_ followed by strictly integers
    refseq = '^.P_[0-9]*'
    rs = re.compile(refseq)
    m = rs.match(accession_name)
    if m and m.group(0) == accession_name:
        return "refseq"

    #genbank seems to be strictly numbers followed by strictly integers
    genbank = '^[a-zA-Z]*[0-9]*'
    gb = re.compile(genbank)
    m = gb.match(accession_name)
    if m and m.group(0) == accession_name:
        return "genbank"

    #pdb is a mix of integers and characters with _A,_B,_C at the end
    pdb = '^[A-Za-z0-9]*_[A-Z]'
    pb = re.compile(pdb)
    m = pb.match(accession_name)
    if m and m.group(0) == accession_name:
        return "pdb"

    #uniprot is a shit format that is mixed integers and numbers, this will probably fail (overmatch)
    uniprot = '^[A-Za-z0-9]*'
    up = re.compile(uniprot)
    m = up.match(accession_name)
    if m and m.group(0) == accession_name:
        return "uniprot"

    return "???"



    re.compile(refseq)
    test = re.split(accession_name)

def get_fasta_sequences():
    filename = ProjectDefinitions.data_directory + "yeast_genome_proteins.fasta"
    sys.stdout.write("parsing record... ")
    #Create a BioPython seq object
    records = list(SeqIO.parse(filename, "fasta"))
    print "parse complete, " + str(len(records)) + " records found!"

    ret={}
    for record in records:
        #print record
        ret[record.name]=record.seq

    return ret


#load information from the blastp and deltablast results and then find matches


#from http://stackoverflow.com/questions/5686211/are-there-any-function-that-can-calculate-score-for-the-aligned-sequences-give-p
def score_pairwise(seq1, seq2, matrix, gap_s, gap_e, gap = True):
    try:
        for A,B in izip(seq1, seq2):
            diag = ('-'==A) or ('-'==B)
            yield (gap_e if gap else gap_s) if diag else matrix[(A,B)]
            gap = diag
    except KeyError:
        print "KeyError in one of the following sequence (typically premature stop codon):\n"
        print "seq1: ", seq1, "\n"
        print "seq2: ", seq2, "\n"
        #put a huge number in the score to make sure this result is nixed
        yield 1000000000
#seq1 = 'MSDLANSEKYYDEDPYGFEDESAPITAEDSWAVISAFFREKGLVSQQLDSFNQFVDYTLQDIICEDSTLILEQLAQHTTESDNISRKYEISFGKIYVTKPMVNESDGVTHALYPQEARLRNLTYSSGLFVDVKKRTYEAIDVPGRELKYELIAEESEDDSESGKVFIGRLPIMLRSKNCYLSEATESDLYKLKECPFDMGGYFIINGSEKVLIAQERSAGNIVQVFKKAAPSPISHVAEIRSALEKGSRFISTLQVKLYGREGSSARTIKATLPYIKQDIPIVIIFRALGIIPDGEILEHICYDVNDWQMLEMLKPCVEDGFVIQDRETALDFIGRRGTALGIKKEKRIQYAKDILQKEFLPHITQLEGFESRKAFFLGYMINRLLLCALDRKDQDDRDHFGKKRLDLAGPLLAQLFKTLFKKLTKDIFRYMQRTVEEAHDFNMKLAINAKTITSGLKYALATGNWGEQKKAMSSRAGVSQVLNRYTYSSTLSHLRRTNTPIGRDGKLAKPRQLHNTHWGLVCPAETPEGQACGLVKNLSLMSCISVGTDPMPIITFLSEWGMEPLEDYVPHQSPDATRVFVNGVWHGVHRNPARLMETLRTLRRKGDINPEVSMIRDIREKELKIFTDAGRVYRPLFIVEDDESLGHKELKVRKGHIAKLMATEYQDIEGGFEDVEEYTWSSLLNEGLVEYIDAEEEESILIAMQPEDLEPAEANEENDLDVDPAKRIRVSHHATTFTHCEIHPSMILGVAASIIPFPDHNQSPRNTYQSAMGKQAMGVFLTNYNVRMDTMANILYYPQKPLGTTRAMEYLKFRELPAGQNAIVAIACYSGYNQEDSMIMNQSSIDRGLFRSLFFRSYMDQEKKYGMSITETFEKPQRTNTLRMKHGTYDKLDDDGLIAPGVRVSGEDVIIGKTTPISPDEEELGQRTAYHSKRDASTPLRSTENGIVDQVLVTTNQDGLKFVKVRVRTTKIPQIGDKFASRHGQKGTIGITYRREDMPFTAEGIVPDLIINPHAIPSRMTVAHLIECLLSKVAALSGNEGDASPFTDITVEGISKLLREHGYQSRGFEVMYNGHTGKKLMAQIFFGPTYYQRLRHMVDDKIHARARGPMQVLTRQPVEGRSRDGGLRFGEMERDCMIAHGAASFLKERLMEASDAFRVHICGICGLMTVIAKLNHNQFECKGCDNKIDIYQIHIPYAAKLLFQELMAMNITPRLYTDRSRDF'
#seq2 = 'MSDLANSEKYYDEDPYGFEDESAPITAEDSWAVISAFFREKGLVSQQLDSFNQFVDYTLQDIICEDSTLILEQLAQHTTESDNISRKYEISFGKIYVTKPMVNESDGVTHALYPQEARLRNLTYSSGLFVDVKKRTYEAIDVPGRELKYELIAEESEDDSESGKVFIGRLPIMLRSKNCYLSEATESDLYKLKECPFDMGGYFIINGSEKVLIAQERSAGNIVQVFKKAAPSPISHVAEIRSALEKGSRFISTLQVKLYGREGSSARTIKATLPYIKQDIPIVIIFRALGIIPDGEILEHICYDVNDWQMLEMLKPCVEDGFVIQDRETALDFIGRRGTALGIKKEKRIQYAKDILQKEFLPHITQLEGFESRKAFFLGYMINRLLLCALDRKDQDDRDHFGKKRLDLAGPLLAQLFKTLFKKLTKDIFRYMQRTVEEAHDFNMKLAINAKTITSGLKYALATGNWGEQKKAMSSRAGVSQVLNRYTYSSTLSHLRRTNTPIGRDGKLAKPRQLHNTHWGLVCPAETPEGQACGLVKNLSLMSCISVGTDPMPIITFLSEWGMEPLEDYVPHQSPDATRVFVNGVWHGVHRNPARLMETLRTLRRKGDINPEVSMIRDIREKELKIFTDAGRVYRPLFIVEDDESLGHKELKVRKGHIAKLMATEYQDIEGGFEDVEEYTWSSLLNEGLVEYIDAEEEESILIAMQPEDLEPAEANEENDLDVDPAKRIRVSHHATTFTHCEIHPSMILGVAASIIPFPDHNQSPRNTYQSAMGKQAMGVFLTNYNVRMDTMANILYYPQKPLGTTRAMEYLKFRELPAGQNAIVAIACYSGYNQEDSMIMNQSSIDRGLFRSLFFRSYMDQEKKYGMSITETFEKPQRTNTLRMKHGTYDKLDDDGLIAPGVRVSGEDVIIGKTTPISPDEEELGQRTAYHSKRDASTPLRSTENGIVDQVLVTTNQDGLKFVKVRVRTTKIPQIGDKFASRHGQKGTIGITYRREDMPFTAEGIVPDLIINPHAIPSRMTVAHLIECLLSKVAALSGNEGDASPFTDITVEGISKLLREHGYQSRGFEVMYNGHTGKKLMAQIFFGPTYYQRLRHMVDDKIHARARGPMQVLTRQPVEGRSRDGGLRFGEMERDCMIAHGAASFLKERLMEASDAFRVHICGICGLMTVIAKLNHNQFECKGCDNKIDIYQIHIPYAAKLLFQELMAMNITPRLYTDRSRDF'
#print sum(score_pairwise(seq1, seq2, blosum, 0, 0))
