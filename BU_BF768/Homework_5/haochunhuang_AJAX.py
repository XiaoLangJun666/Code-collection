#!/usr/bin/env python3
import cgi 
import cgitb
import pymysql
import json
from string import Template

#establish the connection, change username and password before running the script
connection = pymysql.connect(
    host = 'bioed.bu.edu',
    user = 'haochun',
    password = 'Hhc20160406',
    db = 'miRNA',
    port = 4253
)


cgitb.enable()


def query_database(query):
    # Create the cursor object
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        return cursor.fetchall()
    except pymysql.Error as e:
        print(e)
    finally:
        cursor.close()

def main():
    print("Content-type: application/json\n")
    form = cgi.FieldStorage()

    # Check if form data is available
    if not form:
        print(json.dumps([]))  # Empty response if no form data
        return

    selector = form.getvalue('request')
    if not selector:
        print(json.dumps([]))  # Empty response if selector is not provided
        return

    # Handling different types of requests
    if selector == 'histogram':
        miRNA_name = form.getvalue('miRNA')
        if not miRNA_name:
            print(json.dumps([]))  # Empty response if miRNA name is not provided
            return
        
        # Query for histogram data
        query_1 = f"SELECT t.score FROM miRNA mr JOIN targets t USING (mid) WHERE mr.name='{miRNA_name}'"
        results_1 = query_database(query_1)
        print(json.dumps(results_1))

    elif selector == 'table':
        miRNA_seq = form.getvalue('seq')
        if not miRNA_seq:
            print(json.dumps([]))
            return
        query_2=f"select name , seq from  miRNA mr where mr.seq regexp '{miRNA_seq}' order by seq ASC"
        results_2=query_database(query_2)
        print(json.dumps(results_2))

if __name__ == "__main__":
    main()
