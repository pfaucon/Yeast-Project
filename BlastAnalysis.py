"""
This module will take the results of BlastAggregation and compare them against
the drugbase db try and find matches.
"""
import sqlite3
import csv

import ProjectDefinitions

def blast_analysis_main():
    """
    Main function, take the processed blast results and try to find them
    in our drugbank database
    @return: None
    """
    conn = sqlite3.connect(ProjectDefinitions.data_directory + ProjectDefinitions.drugbank_db_name)
    curs = conn.cursor()

def search_record(db_connection,csv_line):
    """
    This function will use a connection and a csv line to search drugbase.
    If a hit is found then we will update the database with the yeast gene.
    @return:
    """
    curs = db_connection.cursor()

def parse_csv_line(csv_line):
    """
    This function returns a dictionary with the fields of the csv line parsed out
    @param csv_line:
    @return:
    """



#from the original drug data some of the drugs are genetic modifiers, and some are physical modifiers
#this function will separate physiological effects from genetic modifications
def fileToSQLite():
    print "trying to convert file to sqlite database..."
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