import os
from os import path
from shutil import copyfile

src_base = "/home/tomap/raw/cluster/home/tomap/raw_files/"
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
	dst_name = destination_addr[i] + "/count_data.h5"
	print(name)
	print(dst_name)
	print(path.exists(destination_addr[i]))
	if path.exists(destination_addr[i]):
		copyfile(name, dst_name)
