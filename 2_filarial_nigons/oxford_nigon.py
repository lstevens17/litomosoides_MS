import sys

def parse_table(table_file):
	with open(table_file, 'r') as table:
		table_dict, chr_list = {}, []
		for line in table:
			if not line.startswith("#"):
				cols = line.rstrip("\n").split()
				buscoID, status = cols[0], cols[1]
				if status == 'Complete':
					chr, start, stop = cols[2].split(":")[0], int(cols[3]), int(cols[4])
					table_dict[buscoID] = [chr, str(((start+stop)/2))]
					if not chr in chr_list:
						chr_list.append(chr)
	return table_dict, sorted(chr_list)

def parse_nigon_table(nigon_table_file):
	with open(nigon_table_file, 'r') as nigon_tsv:
		busco2nigon_dict = {}
		for line in nigon_tsv:
			if not line.startswith("Orthogroup"):
				buscoID, nigon = line.rstrip("\n").split("\t")[0], line.rstrip("\n").split("\t")[1]
				busco2nigon_dict[buscoID] = nigon
	return busco2nigon_dict


def print_comparison(sp1_table_dict, sp2_table_dict, busco2nigon_dict):
	for buscoID, sp1_pos_info in sp1_table_dict.items():
		if buscoID in sp2_table_dict:
			sp2_pos_info = sp2_table_dict[buscoID]
			try:
				nigon = busco2nigon_dict[buscoID]
			except KeyError:
				nigon = "unassigned"
			print(buscoID + "\t" + nigon + "\t" + "\t".join(sp1_pos_info) + "\t" + "\t".join(sp2_pos_info))



if __name__ == "__main__":
	sp1_table_file = sys.argv[1]
	sp2_table_file = sys.argv[2]
	nigon_table_file = sys.argv[3]
	sp1_table_dict, sp1_chr_list = parse_table(sp1_table_file)
	sp2_table_dict, sp2_chr_list = parse_table(sp2_table_file)
	busco2nigon_dict = parse_nigon_table(nigon_table_file)
	print_comparison(sp1_table_dict, sp2_table_dict, busco2nigon_dict)

