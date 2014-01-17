import sqlite3
import csv
import ProjectDefinitions

def aggregate_to_sqlite_main():
    sql_file = ProjectDefinitions.data_directory+ProjectDefinitions.drugbank_db_name
    aggregate_blastp = ProjectDefinitions.results_directory + "aggregate_drugbank_blastpresults_unfound.txt"
    aggregate_delatblast = ProjectDefinitions.results_directory + "aggregate_drugbank_deltablastresults_unfound.txt"

    aggregate_to_sqlite(aggregate_blastp,aggregate_delatblast,sql_file)

def aggregate_to_sqlite(aggregate_blastp_filename,aggregate_deltablast_filename,sqlite_filename):


    "sequence name\tsequence_length\tmatch_db\tmatch_name\tmatch_length\tpct_cover\tscore_frac\tC_score\te_value\n"

    conn = sqlite3.connect(sqlite_filename)
    curs = conn.cursor()
    curs.execute("DROP TABLE IF EXISTS 'aggregate'");
    curs.execute("CREATE TABLE aggregate (sequence_name TEXT NOT NULL, match_db TEXT, match_name TEXT NOT NULL, p_C_score DOUBLE, d_C_score DOUBLE);")
    curs.execute("CREATE INDEX aggregate_score_idx ON aggregate(p_C_score)")
    print "database created, reading in..."

    with open(aggregate_blastp_filename, 'rU') as fin:
        dr = csv.DictReader(fin, delimiter='\t')
        #to_db = [(i['sequence name'], i['sequence_length'], i['match_db'], i['match_name'], i['match_length'], i['pct_cover'], i['score_frac'], i['C_score'], i['e_value']) for i in dr]
        to_db = [(i['sequence_name'], i['match_db'], i['match_name'], i['C_score']) for i in dr]

    curs.executemany("INSERT INTO aggregate (sequence_name,match_db,match_name,p_C_score) VALUES (?, ?, ?, ?);", to_db)
    conn.commit()

    with open(aggregate_deltablast_filename, 'rU') as fin:
        dr = csv.DictReader(fin, delimiter='\t')
        #to_db = [(i['sequence name'], i['sequence_length'], i['match_db'], i['match_name'], i['match_length'], i['pct_cover'], i['score_frac'], i['C_score'], i['e_value']) for i in dr]
        to_db = [(i['sequence_name'], i['match_db'], i['match_name'], i['C_score']) for i in dr]

    curs.executemany("INSERT INTO aggregate (sequence_name,match_db,match_name,p_C_score) VALUES (?, ?, ?, ?);", to_db)
    conn.commit()

    #curs.execute("SELECT * FROM aggregate")
    #curs.execute("SELECT MAX(p_C_score), sequence_name, match_name FROM aggregate GROUP BY sequence_name,match_name")
    sql_query = "SELECT UniProt_ID, Drug_IDs, match_db, match_name, p_C_score FROM "
    sql_query += "drugs INNER JOIN (SELECT sequence_name, match_db, match_name, MAX(p_C_score) AS p_C_score FROM "
    sql_query += "(SELECT * FROM aggregate WHERE p_C_score > 70) GROUP BY sequence_name,match_name)aggregate"
    sql_query += " ON drugs.Uniprot_ID=aggregate.sequence_name"
    sql_query += " ORDER BY p_C_score DESC"
    print sql_query
    curs.execute(sql_query)
    count=0
    outfile = open(ProjectDefinitions.final_result_csv,"w")
    #outfile.write("Uniprot_ID,Drug_ID,match_db,match_name,p_C_score")
    writer = csv.writer(outfile)
    writer.writerow(['Uniprot_ID','Drug_ID','match_db','match_name','p_C_score'])
    writer.writerows(curs)
    """
    while True:
        row = curs.fetchone()
        if row == None:
            print "done printing drug results, ", count, " in total"
            break
        count += 1
        print row[0], row[1], row[2]
        outfile.write("")
    """

# Standard boilerplate to call the main() function to begin
# the program.  This also ensures that the main isn't called when used as a module
if __name__ == '__main__':
    aggregate_to_sqlite_main()