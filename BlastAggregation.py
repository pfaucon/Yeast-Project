#given the results from blastp and deltablast do some basic analytics on the results

import os
import glob
import sys
import re
import sqlite3

from Bio.Blast import NCBIXML
from Bio import SeqIO

#from http://stackoverflow.com/questions/5686211/are-there-any-function-that-can-calculate-score-for-the-aligned-sequences-give-p
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from itertools import izip
blosum.update(((b,a),val) for (a,b),val in blosum.items())

import ProjectDefinitions


def blast_aggregation_main():

    print("aggregating results...")


    fasta_file = ProjectDefinitions.data_directory + "yeast_genome_proteins.fasta"
    seq = get_fasta_sequences(fasta_file)

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
    outhandle = open(outfilename,"w")
    outhandle.write("sequence_name\tsequence_length\tmatch_db\tmatch_name\tmatch_length\tpct_cover\tscore_frac\tC_score\te_value\n")

    outhandle_unfound = open(outfilename[:-4] + "_unfound.txt","w")
    outhandle_unfound.write("sequence_name\tsequence_length\tmatch_db\tmatch_name\tmatch_length\tpct_cover\tscore_frac\tC_score\te_value\n")

    conn = sqlite3.connect(ProjectDefinitions.data_directory+ProjectDefinitions.drugbank_db_name)
    curs = conn.cursor()

    for file in glob.glob(directory):

        if(os.path.getsize(file) <= 1):
            print("File: " + file + " Appears to be  empty, skipping...")
            continue
        if ProjectDefinitions.DEBUG:
            print("now working on:" + file)
        result_handle = open(file)
        records = NCBIXML.parse(result_handle)

        for blast_record in records:

            nalignments=1000
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

                    sequence_name = get_gene_id_from_fasta_name(blast_record.query)
                    query_tokens = blast_record.query.split()
                    match_name = query_tokens[0]
                    #fasta sequences seem to end with '*', just strip that character
                    query_sequence = sequences[match_name][:-1]
                    #calculate the maximum possible score for this sequence, gap cost is arbitrary
                    max_score = sum(score_pairwise(query_sequence,query_sequence,blosum,-1,-1))

                    #if we have less than 50% of the points it isn't relevant
                    #if ((hsp.score*100)/max_score) < 50:
                    #    continue

                    match_length = hsp.query_end-hsp.query_start
                    #if we have less than 80% coverage then it isn't relevant
                    #if ((match_length*100.0)/blast_record.query_letters)<80:
                    #    continue

                    match_db = get_match_db(alignment.accession)
                    match_found = False
                    #FIXMEL refseq matching is currently broken, I should find a workaround
                    if match_db == "GenBank_Protein_ID" or match_db == "UniProt_ID" or match_db == "PDB_ID":

                        curs.execute("SELECT Gene_Name FROM drugs WHERE "+match_db+" LIKE '"+alignment.accession+"'")
                        #curs.execute("SELECT Gene_Name FROM drugs WHERE "+ "PDB_ID" +" LIKE '"+ "2NSI" + "'")
                        row = curs.fetchone()
                        if row != None:
                            print "Found a match for " + alignment.accession + ", it is: " + row[0]
                            match_found = True
                    pct_cover = (match_length*100)/blast_record.query_letters
                    score_frac = (hsp.score*100)/max_score
                    outstr = sequence_name +"\t"+ str(blast_record.query_letters) +"\t" + \
                             match_db +"\t"+ alignment.accession +"\t"+ str(match_length) +"\t"+ \
                             str(pct_cover) +"\t" + str(score_frac) +"\t" + str(score_frac*pct_cover/100.0) +"\t"+ \
                             str(hsp.expect) + "\n"
                    if(match_found):
                        outhandle.write(outstr)
                    else:
                        outhandle_unfound.write(outstr)

            #for full results comment this line out
            #break

        result_handle.close()



    outhandle.close()
    outhandle_unfound.close()





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
        return "GenBank_Protein_ID"

    #pdb is a mix of integers and characters with _A,_B,_C at the end
    pdb = '^[A-Za-z0-9]*_[A-Z]'
    pb = re.compile(pdb)
    m = pb.match(accession_name)
    if m and m.group(0) == accession_name:
        return "PDB_ID"

    #uniprot is a shit format that is mixed integers and numbers, this will probably fail (overmatch)
    uniprot = '^[A-Za-z0-9]*'
    up = re.compile(uniprot)
    m = up.match(accession_name)
    if m and m.group(0) == accession_name:
        return "UniProt_ID"

    return "???"



    re.compile(refseq)
    test = re.split(accession_name)

def get_fasta_sequences(fasta_file):
    sys.stdout.write("parsing record... ")
    #Create a BioPython seq object
    records = list(SeqIO.parse(fasta_file, "fasta"))
    print "fasta sequence parse complete, " + str(len(records)) + " records found!\n"

    ret={}
    for record in records:
        #print record
        ret[record.name]=record.seq

    return ret


def get_gene_id_from_fasta_name(fasta_name):
    """
    This method will extract the gene name from the fasta name, the type of ID is unknown to this method
    @param fasta_name:
    @return:
    """
    fasta_parts = fasta_name.split('|')
    return fasta_parts[1]


#FIXME: this alignment seems to have problems with the following sequences for unknown reasons, debug at some point
#MFQQFQASCLVLFFLVGFAQQTLKPQNRKVDCNKGVTGTIYEYGALTLNGEEYIQFKQFAGKHVLFVNVAAYUGLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMHWFHQAPVSTVKSDILEYLKQFNT
#MGCAEGKAVAAAAPTELQTKGKNGDGRRRSAKDHHPGKTLPENPAGFTSTATADSRALLQAYIDGHSVVIFSRSTCTRCTEVKKLFKSLCVPYFVLELDQTEDGRALEGTLSELAAETDLPVVFVKQRKIGGHGPTLKAYQEGRLQKLLKMNGPEDLPKSYDYDLIIIGGGSGGLAAAKEAAQYGKKVMVLDFVTPTPLGTRWGLGGTCVNVGCIPKKLMHQAALLGQALQDSRNYGWKVEETVKHDWDRMIEAVQNHIGSLNWGYRVALREKKVVYENAYGQFIGPHRIKATNNKGKEKIYSAERFLIATGERPRYLGIPGDKEYCISSDDLFSLPYCPGKTLVVGASYVALECAGFLAGIGLDVTVMVRSILLRGFDQDMANKIGEHMEEHGIKFIRQFVPIKVEQIEAGTPGRLRVVAQSTNSEEIIEGEYNTVMLAIGRDACTRKIGLETVGVKINEKTGKIPVTDEEQTNVPYIYAIGDILEDKVELTPVAIQAGRLLAQRLYAGSTVKCDYENVPTTVFTPLEYGACGLSEEKAVEKFGEENIEVYHSYFWPLEWTIPSRDNNKCYAKIICNTKDNERVVGFHVLGPNAGEVTQGFAAALKCGLTKKQLDSTIGIHPVCAEVFTTLSVTKRSGASILQAGCU
#MSLGRLCRLLKPALLCGALAAPGLAGTMCASRDDWRCARSMHEFSAKDIDGHMVNLDKYRGFVCIVTNVASQUGKTEVNYTQLVDLHARYAECGLRILAFPCNQFGKQEPGSNEEIKEFAAGYNVKFDMFSKICVNGDDAHPLWKWMKIQPKGKGILGNAIKWNFTKFLIDKNGCVVKRYGPMEEPLVIEKDLPHY
#MARLLQASCLLSLLLAGFVSQSRGQEKSKMDCHGGISGTIYEYGALTIDGEEYIPFKQYAGKYVLFVNVASYUGLTGQYIELNALQEELAPFGLVILGFPCNQFGKQEPGENSEILPTLKYVRPGGGFVPNFQLFEKGDVNGEKEQKFYTFLKNSCPPTSELLGTSDRLFWEPMKVHDIRWNFEKFLVGPDGIPIMRWHHRTTVSNVKMDILSYMRRQAALGVKR
#MAFIAKSFYDLSAISLDGEKVDFNTFRGRAVLIENVASLUGTTTRDFTQLNELQCRFPRRLVVLGFPCNQFGHQENCQNEEILNSLKYVRPGGGYQPTFTLVQKCEVNGQNEHPVFAYLKDKLPYPYDDPFSLMTDPKLIIWSPVRRSDVAWNFEKFLIGPEGEPFRRYSRTFPTINIEPDIKRLLKVA
#MCAARLAAAAAAAQSVYAFSARPLAGGEPVSLGSLRGKVLLIENVASLUGTTVRDYTQMNELQRRLGPRGLVVLGFPCNQFGHQENAKNEEILNSLKYVRPGGGFEPNFMLFEKCEVNGAGAHPLFAFLREALPAPSDDATALMTDPKLITWSPVCRNDVAWNFEKFLVGPDGVPLRRYSRRFQTIDIEPDIEALLSQGPSC

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
