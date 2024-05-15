#!/usr/bin/env python3
"""
This is a python file containing three functions
read_query()
connect_database()
execute_query()

Import into another program as:

import part1_shell as tools

"""
import re #you can use this, if you want, for regular expression search to exclude comment lines in query file
import pymysql

def read_query(filename):
	#open file, read lines and concatenate into one string
	#convert newline characters to blanks
	#ignore lines that start with "--"
	#ignore ends of lines including and beyond "--"
	query=''
	with open(filename) as file :
		for line in file:
			if line.strip().startswith('--'):
				continue
			comment_start=line.find('--')
			if comment_start != -1:
				line=line[:comment_start]
			query += line.replace('\n', ' ').strip()+' '
	
	return query
	
def connect_database(database, username, password):
	#connect to database on bioed
	#create cursor object
	connection = pymysql.connect(
	host='bioed.bu.edu',
	user=username,          #use your username for MariaDB
	password=password,    #use your password for MariaDB
	db= database,            #your db is the same as your username
	port = 4253)          #required on bioed
	# Create the cursor object
	cursor = connection.cursor()
	return (connection, cursor)


def execute_query(cursor, query):
	#execute the query with try:  .. except:
	#on error, print the query and the error
	#return all rows in results
	try:
		cursor.execute(query)
	except pymysql.Error as e:
		print(e)
	results = cursor.fetchall()	

	return results

