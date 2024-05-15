#!/usr/bin/env python3

import haochun_part1  as tools
import sys


def main(database,username,password,filename):
	query=tools.read_query(filename)

	connection,cursor=tools.connect_database(database,username,password)

	results=tools.execute_query(cursor,query)

	print("Query executed: ")
	print(query)

	print("\nResults:")
	for row in results:
		print(row)

	connection.close()


if __name__ =="__main__":
	db,user,pw,fname=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

	main(db,user,pw,fname)
