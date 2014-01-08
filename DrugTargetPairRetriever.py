# Return a set of dictionaries providing drug->gene and gene->drug pairs
# If we don't have them locally then go to www.drugbank.ca/ to pull them

import ProjectDefinitions
import sys

import urllib
import zipfile
import csv,sqlite3

#all of the drug information, more than we need now
#zipped_filename = data_directory + "drugbank.xml.zip"
#filename = data_directory + "drugbank.xml"
#source_url = "http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip"

#all known targets of drugs
zipped_filename = ProjectDefinitions.drugbankdir + "all_target_ids_all.csv.zip"
filename = ProjectDefinitions.drugbankdir + "all_target_ids_all.csv"
source_url = "http://www.drugbank.ca/system/downloads/current/all_target_ids_all.csv.zip"

def drug_main():
    drug_data = getDrugData()
    
    return drug_data

def getDrugData():

    #Download the data if it isn't present
    try:
        if ProjectDefinitions.DEBUG:
            sys.stdout.write("trying to open yeast data file at: " + filename + " ")
        f_ptr = open(filename, "r")
        if ProjectDefinitions.DEBUG:
            print "File opened successfully!"
        f_ptr.close()
    except IOError:
        

        
        print "WARNING: Couldn't find the file, downloading the compressed version..."
        print "source: " + source_url
        print "destination: " + zipped_filename
        
        #download from previously known good loction
        urllib.urlretrieve(source_url, zipped_filename)
        if ProjectDefinitions.DEBUG:
            sys.stdout.write("Download completed, decompressing...")
        
        #unzip the results
        zf_ptr = zipfile.ZipFile(zipped_filename)
        zipfile.ZipFile.extractall(zf_ptr,ProjectDefinitions.drugbankdir)
        if ProjectDefinitions.DEBUG:
            print "file decompressed!"

    fileToSQLite();

    # return the full contents of the fasta file
    return "test";


#from the original drug data some of the drugs are genetic modifiers, and some are physical modifiers
#this function will separate physiological effects from genetic modifications
def fileToSQLite():
    print "trying to convert file to sqlite database..."
    import csv, sqlite3
    conn = sqlite3.connect(ProjectDefinitions.data_directory+"drugbank.sqlite")
    curs = conn.cursor()
    curs.execute("DROP TABLE IF EXISTS 'drugs'");
    curs.execute("CREATE TABLE drugs (ID INTEGER PRIMARY KEY,Name TEXT,Gene_Name TEXT ,GenBank_Protein_ID INTEGER,GenBank_Gene_ID TEXT,UniProt_ID TEXT,Uniprot_Title TEXT ,PDB_ID TEXT,GeneCard_ID TEXT,GenAtlas_ID TEXT,HGNC_ID TEXT,HPRD_ID INTEGER,Species_Category TEXT,Species TEXT,Drug_IDs TEXT);")
    #curs.execute("CREATE TABLE drugs (geneName TEXT PRIMARY KEY, type INTEGER, term TEXT, definition TEXT);")
    
    print "database created, reading in..."
    
    with open(filename, 'rU') as fin:
        dr = csv.DictReader(fin)
        to_db = [(i['ID'], i['Name'], i['Gene Name'], i['GenBank Protein ID'], i['GenBank Gene ID'], i['UniProt ID'], i['Uniprot Title'], i['PDB ID'], i['GeneCard ID'], i['GenAtlas ID'], i['HGNC ID'], i['HPRD ID'], i['Species Category'], i['Species'], i['Drug IDs']) for i in dr]
    
    curs.executemany("INSERT INTO drugs (ID,Name,Gene_Name,GenBank_Protein_ID,GenBank_Gene_ID,UniProt_ID,Uniprot_Title ,PDB_ID,GeneCard_ID,GenAtlas_ID,HGNC_ID,HPRD_ID,Species_Category,Species,Drug_IDs) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)
    conn.commit()

    curs.execute("SELECT Gene_Name FROM drugs WHERE Drug_IDs LIKE '%DB00157%'")
    while True:
        row = curs.fetchone()
        if row == None:
            print "done printing drug results!"
            break
        print row[0]


# Standard boilerplate to call the main() function to begin
# the program.  This also ensures that the main isn't called when used as a module
if __name__ == '__main__':
    drug_main()