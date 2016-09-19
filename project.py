# dict mapping nucleotide codons to amino acids.
CODON_TO_AA = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

COMPLEMENT_BASE = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}



def get_sequence(filename):
	'''
	Returns DNA sequence from fasta file, stripped of first line
	'''

	sequence = ''
	with open(filename, 'r') as f:
		lines = f.readlines()
	lines = lines[1 : ]
	for line in lines:
		sequence = sequence + line[0 : len(line) - 2]

	return sequence 


def reverse_complement(sequence):
	'''
	Returns reverse complement of DNA sequence
	'''
	complement = ''
	reverse = ''
	for base in sequence:
		complement = complement + COMPLEMENT_BASE[base]

	reverse = complement[::-1]
	return reverse

#[line[i:i+n] for i in range(0, len(line), n)]

def find_codons(sequence):
	codons = []

	for i in range(0, len(sequence), 3):
		codons.append(sequence[i:i+3])
	return codons

#def find_ORF(codons):


	#print codons.index("ATG")
	#return start_index

# def find_stop(sequence):
# 	return stop_index

# get_sequence('escherichia_coli_k12_mg1655.fa')


