### FUNCTIONS ###

def read_accs_and_sequences_from_fasta(infile):
    """
    Input: readfile: Fasta file.
    Outputs: accs_and_sequences: List of tuples. Containing accs and sequences, e.g. [(acc, aTHNtem..)..()].
    """
    accs = list()
    sequences = list()
    seq = ""
    read_acc = False


    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")

    infile = open(infile, "r")
    readfile = infile.readlines()
    infile.close()

    for line in readfile:
        line = line.strip()
        if line.startswith(">"):
            acc = line.split(">")[1]
            if read_acc:
                accs.append(acc)
                sequences.append(seq)
                #reset sequence string
                seq = ""
            #catch first accesion.
            else:
                accs.append(acc)
        else:
            seq += line
            read_acc = True

    #get last sequence
    sequences.append(seq)
    accs_and_sequences = tuple( zip(accs, sequences) )
    return accs_and_sequences