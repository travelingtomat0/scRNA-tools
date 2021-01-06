import os
from os import path
from shutil import copyfile
import csv

src_base = "/home/tomap//raw/cluster/home/tomap/raw_files/"
destination_base = "/home/tomap/new_infercnv_annotations/"
destination_addr = []

direc = os.listdir()
	
tmp = []

for el in direc:
	print(el)
	if el[0] == 'M':
		tmp.append(el)
		destination_addr.append(destination_base + el[0:7] + "_output_dir")

print(tmp)

for i in range(0, len(tmp)):
	name = src_base + tmp[i]
	dst_name = destination_addr[i] + "/count_data.txt"
	new_name = destination_addr[i] + "/nexus_annotations.txt"
	print(name)
	print(dst_name)
	print(path.exists(destination_addr[i]))
	if path.exists(destination_addr[i]):
		res = []
		with open(name) as inp:
			reader = csv.reader(inp, delimiter='\t')
			a, *reader = reader
			for row in reader:
				a, b = row
				if b == "Melanoma.melanocytic":
					res.append(a)
				if b == "Melanoma.mesenchymal":
					res.append(a)
		with open(dst_name, 'w') as out:
			writer = csv.writer(out, delimiter='\t')
			for el in res:
				writer.writerow([el])
		os.rename(dst_name, new_name)
