# Return a set of dictionaries providing drug->gene and gene->drug pairs
# If we don't have them locally then go to www.drugbank.ca/ to pull them

import urllib
import os #for abs path
import zipfile

#do a bunch of prints to ensure things are working, it is very verbose
DEBUG = False

#this isn't important data but having it cached is nice, abspath strips trailing '/'
data_directory = os.path.abspath("../../Data/Cache/DrugBank") + "/"

def drug_main():
    drug_data = getDrugData()
    
    return drug_data

def getDrugData():
    

    #all of the drug information, more than we need now
    #zipped_filename = data_directory + "drugbank.xml.zip"
    #filename = data_directory + "drugbank.xml"
    #source_url = "http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip"
    
    #all known targets of drugs
    zipped_filename = data_directory + "all_target_ids_all.csv.zip"
    filename = data_directory + "all_target_ids_all.csv"
    source_url = "http://www.drugbank.ca/system/downloads/current/all_target_ids_all.csv.zip"

    #Download the data if it isn't present
    try:
        if DEBUG:
            sys.stdout.write("trying to open yeast data file at: " + filename + " ")
        f_ptr = open(filename, "r")
        if DEBUG:
            print "File opened successfully!"
        f_ptr.close()
    except IOError:
        
        if not os.path.exists(data_directory):
            sys.stdout.write("WARNING: cache directories were not found, creating... ")
            os.makedirs(data_directory)
            print "OK!"
        
        print "WARNING: Couldn't find the file, downloading the compressed version..."
        print "source: " + source_url
        print "destination: " + zipped_filename
        
        #download from previously known good loction
        urllib.urlretrieve(source_url, zipped_filename)
        if DEBUG:
            sys.stdout.write("Download completed, decompressing...")
        
        #unzip the results
        zf_ptr = zipfile.ZipFile(zipped_filename)
        zipfile.ZipFile.extractall(zf_ptr,data_directory)
        if DEBUG:
            print "file decompressed!"
    
    # return the full contents of the fasta file
    return "test";


# Standard boilerplate to call the main() function to begin
# the program.  This also ensures that the main isn't called when used as a module
if __name__ == '__main__':
    drug_main()