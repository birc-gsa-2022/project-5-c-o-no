import random
import sys
import names

alphabet = ["A", "C", "G", "T"]

if (len(sys.argv) != 5):
    print("Missing argument(s):\n")
    print("- Number of sequences")
    print("- Size of sequences")
    print("- Number of patterns")
    print("- Size of patterns")
    exit(0)

num_of_seqs = int(sys.argv[1])
size_of_seqs = int(sys.argv[2])
num_of_patterns = int(sys.argv[3])
size_of_patterns = int(sys.argv[4])

patterns = []
fasta_string = ""
pattern_string = ""

for _ in range(num_of_patterns):
    patterns.append(''.join(random.choice(alphabet) for _ in range(size_of_patterns)))

for _ in range(num_of_seqs):
    spaces_prefix = ' ' * random.randint(0,10) 
    spaces_suffix = ' ' * random.randint(0,10)
    ### Generate header
    seq_header = '>' + spaces_prefix + names.get_random_name() + spaces_suffix
    fasta_string += seq_header + '\n'

    ### Generate sequence
    seq = ''.join(random.choice(alphabet) for _ in range(size_of_seqs))

    # Inject random newlines in sequences
    num_of_newlines = random.randint(0,size_of_seqs*0.1)
    indices_to_screw_up_with_newlines = [random.randint(0,size_of_seqs-1) for _ in range(num_of_newlines)]
    for i in indices_to_screw_up_with_newlines:
        seq = seq[:i] + '\n' + seq[i:]

    # Inject the patterns in sequences
    num_of_pattern_injections = num_of_patterns*5
    for i in range(0,size_of_seqs+num_of_newlines,size_of_patterns):
        probability_of_insert = num_of_pattern_injections/(size_of_seqs+num_of_newlines)
        if (random.uniform(0,1) <= probability_of_insert): 
            seq = seq[:i] + random.choice(patterns) + seq[i:]


    fasta_string += seq + '\n'

#print(fasta_string)
for i, pattern in enumerate(patterns):
    pattern_string += "@read"+str(i+1)+"\n"+pattern+"\n"

f = open("fasta-n" + str(size_of_seqs) + "-m" + str(size_of_patterns) + ".fa", "w")
f.write(fasta_string+"\n")
f.close()
    
f = open("reads-n" + str(size_of_seqs) + "-m" + str(size_of_patterns) + ".fastq", "w")
f.write(pattern_string+"\n")
f.close()








