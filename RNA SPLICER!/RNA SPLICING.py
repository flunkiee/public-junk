#


# def fasta_parser(file_path):
#     sequences = {}  # dictionary to store sequences with headers
#     current_header = ""
#     current_sequence = ""

#     with open(file_path, 'r') as file:
#         for line in file:
#             line = line.strip()
#             if line.startswith(">"):
#                 if current_header:  # save previous sequence
#                     sequences[current_header] = current_sequence
#                     current_sequence = ""  # reset current sequence
#                 current_header = line[1:]  # extract header without '>'
#             else:
#                 current_sequence += line

#         # save the last sequence
#         if current_header:
#             sequences[current_header] = current_sequence

#         fullseq = next(iter(sequences.values()))
#         introns = list(sequences.values())[1:]
#     return fullseq, introns
# # #ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGTTCAAAGTTTGCGCCTAG
# # fullseq = "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"
# # introns = ["ATCGGTCGAA","ATCGGTCGAGCGTGT"]

# def intron_remover(fullseq, introns): #clear introns before translation
#     exons = fullseq
#     for intron in introns:
#         exons = exons.replace(intron, "")

#     print(exons)
#     return exons #be careful about return placement 

# def transcriptor(exons):
#     transseq = exons.replace("T","U")
#     print(transseq)
#     return transseq

# def translator(transseq): #why is this sequence in french.
#     residuetable = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", #table of residues
#     "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
#     "UAU":"Y", "UAC":"Y", "UAA":" STOP", "UAG":" STOP",
#     "UGU":"C", "UGC":"C", "UGA":" STOP", "UGG":"W",
#     "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
#     "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
#     "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
#     "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
#     "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
#     "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
#     "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
#     "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
#     "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
#     "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
#     "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
#     "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#     amino_acids = [] # initialize empty list
#     for i in range(0, len(transseq), 3): # for values in range of 1 - maxlength iterate in steps of 3
#         codon = transseq[i:i+3] # codon is equal to the current step of 3
#         if codon in residuetable: # if the codon is in the table
#             amino_acids.append(residuetable[codon]) #add to a list named aminoacids
#     return amino_acids

# def __main__():
#     fullseq, introns = fasta_parser("rosalind_splc.txt")
#     exons = intron_remover(fullseq, introns)
#     transseq = transcriptor(exons)
#     amino_acids = translator(transseq)
#     print(''.join(map(str, amino_acids)))

# __main__()

#ATGGTCTACATAGCTGACAAACAGCACGTAGCATCTCGAGAGGCATATGGTCACATGTTCAAAGTTTGCGCCTAG

class RNA_SPLICER:
    def __init__(self, file_path):
        self.file_path = file_path
        self.fullseq = None
        self.introns = None
        self.exons = None
        self.transseq = None
        self.current_header = ""
        self.current_sequence = ""
        self.residuetable = {
            "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
            "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
            "UAU":"Y", "UAC":"Y", "UAA":" STOP", "UAG":" STOP",
            "UGU":"C", "UGC":"C", "UGA":" STOP", "UGG":"W",
            "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
            "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
            "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
            "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        }

    def fasta_parser(self):
        sequences = {}
        with open(self.file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if self.current_header:
                        sequences[self.current_header] = self.current_sequence
                        self.current_sequence = ""
                    self.current_header = line[1:]
                else:
                    self.current_sequence += line
            if self.current_header:
                sequences[self.current_header] = self.current_sequence
        self.fullseq = next(iter(sequences.values()))
        self.introns = list(sequences.values())[1:]

    def intron_remover(self):
        self.exons = self.fullseq
        for intron in self.introns:
            self.exons = self.exons.replace(intron, "")
        return self.exons

    def transcriptor(self):
        self.transseq = self.exons.replace("T", "U")
        return self.transseq

    def translator(self):
        amino_acids = [self.residuetable[self.transseq[i:i+3]] for i in range(0, len(self.transseq), 3) if self.transseq[i:i+3] in self.residuetable]
        return amino_acids

# main
if __name__ == "__main__":
    file_path = "rosalind_splc.txt"
    splice = RNA_SPLICER(file_path)
    splice.fasta_parser()
    splice.intron_remover()
    splice.transcriptor()
    amino_acids = splice.translator()

    print("Amino Acids:", (''.join(map(str, amino_acids))))