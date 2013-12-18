import os
import sys

print "I think I'm running from '" + os.path.abspath(".")+ "', if that's not where my cache is then make sure I'm being launched properly!"

sys.stdout.write("Verifying local libraries...")

import YeastORFRead
import HumanORFFind
import DrugTargetPairRetriever

print "Libraries verified!"

is_dna_seq = False

sequences = YeastORFRead.yeast_main(is_dna_seq)
#orfs = execfile("YeastORFRead.py");
#humans = HumanORFFind.blast_main(sequences,is_dna_seq)
drugs = DrugTargetPairRetriever.drug_main()

#print len(orfs)