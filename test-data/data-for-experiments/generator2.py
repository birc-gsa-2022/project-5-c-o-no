import random
import sys
import names

alphabet = ["A", "C", "G", "T"]

for n in range(100, 10000, 100):

    # fasta_string = "> seq\n" + ''.join(random.choices(alphabet, k=n))
    # f = open(f"test-data/dna/fasta-{n}.fa", "w")
    # f.write(fasta_string+"\n")
    # f.close()

    fasta_string = "> seq\n" + "A"*n
    f = open(f"test-data/sameChar/fasta-{n}.fa", "w")
    f.write(fasta_string+"\n")
    f.close()

