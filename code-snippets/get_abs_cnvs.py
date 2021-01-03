import csv

# TODO: This could also only be used for the marker gene copy numper variations. (Would also make more sense, since only used those dimensions in umap)

tmp = []

marker_names = []

def in_names(gene):
	for el in marker_names:
		if el == gene:
			print(gene)
			return True
	return False

with open("Garnett_Markers.csv") as inp:
	reader = csv.reader(inp, delimiter=',')
	a, *reader = reader
	for row in reader:
		a, b, c, d = row
		if a != "":
			marker_names.append(a)
			print(a)
		if b != "":
			marker_names.append(b)
			print(b)
		if c != "":
			marker_names.append(c)
			print(c)
		if d != "":
			marker_names.append(d)
			print(d)

print("---")

with open("infercnv.21_denoised.observations.txt") as inp:
	reader = csv.reader(inp, delimiter=' ')
	a, *reader = reader
	for el in a:
		tmp.append([el, 0])
	hd, *tl = reader
	if in_names(hd[0]):
		for i in range(1, len(hd)):
			tmp[i-1][1] = abs(float(hd[i])-1)
	ct = 0
	for row in tl:
		if in_names(row[0]) & (ct == 0):
			ct = ct + 1
			for i in range(1, len(row)):
				tmp[i-1][1] = tmp[i-1][1] + abs(float(row[i])-1)
		elif in_names(row[0]):
			for i in range(1, len(row)):
				tmp[i-1][1] = (tmp[i-1][1] + abs(float(row[i])-1))/2.0

with open("cnv_overexpr.txt", 'w') as out:
	writer = csv.writer(out, delimiter=' ')
	for row in tmp:
		writer.writerow(row)
