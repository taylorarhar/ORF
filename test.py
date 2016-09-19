import project

sequence = project.get_sequence('escherichia_coli_k12_mg1655.fa')
#project.reverse_complement(sequence)
project.find_start(sequence)