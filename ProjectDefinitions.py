import os
import sys

#do a bunch of prints to ensure things are working, it is very verbose
DEBUG = False

data_directory = os.path.abspath("../../Data/Cache") + "/"
fastadir = data_directory + "fasta/"
blastpdir = data_directory + "blastp/"
deltablastdir = data_directory + "deltablast/"
drugbankdir = data_directory + "drugbank/"
drugbankinfo = drugbankdir + "all_target_ids_all.csv"

def makeDirsIfNecessary(Directory):
    if not os.path.exists(Directory):
        if DEBUG:
            sys.stdout.write("INFO: cache directory '" + Directory + "' was not found, creating... ")
        os.makedirs(Directory)