import sys

# amino acid mapping table
amino_acids = {
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R", # Arg
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A", #Ala
    "AAU":"N","AAC": "N", #Asn
    "GAU":"D","GAC":"D", #Asp
    "UGU":"C","UGC":"C", #Cys
    "GAA":"E","GAG":"E", #Glu
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G", #Gly
    "CAA":"Q","CAG":"Q", #Gln
    "CAU":"H","CAC":"H", #His
    "AUU":"I","AUC":"I","AUA":"I", #Ile
    "AAA":"K","AAG":"K", #Lys
    "UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L", #Leu
    "AUG":"Met", #Met
    "UUU":"F","UUC":"F", #Phe
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P", #Pro
    "UAA":"Stop","UAG":"Stop","UGA":"Stop", #Stop
    "AGU":"S","AGC":"S","UCU":"S","UCC":"S","UCA":"S","UCG":"S", #Ser
    "UGG":"W", #Trp
    "UAU":"Y","UAC":"Y", #Tyr
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T", #Thr
    "GUU":"V","GUG":"V","GUA":"V","GUG":"V" #Val
}

# DNA to mRNA mapping
DNA_to_mRNA = {
    "T":"A",
    "G":"C",
    "C":"G",
    "A":"U"
}

# mRNA complement mapping
mRNA_to_complement = {
    "A":"U",
    "U":"A",
    "C":"G",
    "G":"C"
}

def translate_sequence(frame, direction, data):
    
    local_DNA = data.upper()
    start_frame = frame-1
    
    mRNA_sequence = transcribe_to_mRNA(local_DNA)
    tRNA_sequence = transcribe_to_tRNA(mRNA_sequence)
    
    if direction is "natural":
        print "5'3' Frame "+str(frame)
    elif direction is "reverse":
        print "3'5' Frame "+str(frame)
        tRNA_sequence = get_reverse_strand(tRNA_sequence)

    map_amino_acid(tRNA_sequence[start_frame:])
    print "\n"


def get_reverse_strand(sequence):
    
    sequence = sequence[::-1]
    return ''.join([mRNA_to_complement[nucleotide] for nucleotide in sequence])
    

def map_amino_acid(sequence):
    
    sequence_length = len(sequence)
    current=0
    protein = []
    while (current < sequence_length):
        
        codon = sequence[current:current+3]
        if codon in amino_acids:
            amino_acid = amino_acids[codon]
            protein.append(amino_acid)
        current = current + 3
    
    print ' '.join(protein)
    
def transcribe_to_mRNA(DNA_sequence):
    
    return ''.join([DNA_to_mRNA[nucleotide] for nucleotide in DNA_sequence])
    
def transcribe_to_tRNA(mRNA_sequence):
    
    return ''.join([mRNA_to_complement[nucleotide] for nucleotide in mRNA_sequence])
    

def translate_rig(data):
    
    translate_sequence(1, "natural", data)
    translate_sequence(2, "natural", data)
    translate_sequence(3, "natural", data)
    
    translate_sequence(1, "reverse", data)
    translate_sequence(2, "reverse", data)
    translate_sequence(3, "reverse", data)

def main():
    
    if(len(sys.argv) <= 1):
        print "Usage: $ python openReadingFrames.py gene.txt"
        sys.exit()
    
    with open(sys.argv[1], 'r') as gene_file:
        data = gene_file.readlines()
        data = [line.replace(' ','').strip() for line in data]
        data = ''.join(data)
    
    translate_rig(data)

if __name__ == '__main__':
    main()