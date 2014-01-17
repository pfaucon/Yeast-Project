import os
import sys

import YeastORFRead
import HumanORFFind
import BlastAggregation
import BlastAnalysis
import drugbank_to_yeast

print "I think I'm running from '" + os.path.abspath(".")+ "', if that's not where my cache is then make sure I'm being launched properly!"

#dna seq or protein seq?
is_dna_seq = False
sequences = YeastORFRead.yeast_main(is_dna_seq)

#find reasonable drug/target pairs
import DrugTargetPairRetriever
drugs = DrugTargetPairRetriever.drug_main()

#data analysis mode or actual run-through mode?
lookForResults = False

if lookForResults:
    


    # Perform a blast alignment of all yeast genes against human genes
    # this takes a long time
    humans = HumanORFFind.blast_main(sequences,is_dna_seq)

    #after the alignment we'll aggregate all the results into a single CSV (for speed)
    #this file can be copied from the server, the rest of the processing is lighter lifting
    BlastAggregation.blast_aggregation_main()

else:

    #results = BlastAggregation.blast_aggregation_main()
    execfile("drugbank_to_yeast.py")

    print "fin!"




