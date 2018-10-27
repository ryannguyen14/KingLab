import os, sys 



file_name = "Swi6.fasta"

### read in file and store in string called amino_acids
f = open(file_name, "r")
amino_acids = ('').join([line.strip('\n') for line in f.readlines()[1:]])
f.close()

window_size = 20

partitioned_amino_acids = []
for i in range(0,len(amino_acids)-window_size+1): 
	partitioned_amino_acids.append(amino_acids[i:i+window_size])

newname = "Subset_" + str(window_size) + "_" + file_name
h = open(newname,"w")

for j in range(0,len(partitioned_amino_acids)):
	fasta_entry_name = ">Swi6 Subset " + str(j) +"\n"
	h.write(fasta_entry_name)
	
	if j != len(partitioned_amino_acids)-1:
		h.write(partitioned_amino_acids[j] + "\n")
	else:
		h.write(partitioned_amino_acids[j])

h.close()

#command = "python ipc.py " + newname
#os.system(command)