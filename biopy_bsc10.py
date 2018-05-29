"""
This program, takes an input file from command line; runs a BLAST search (using swissprot database); stores the top ten results in hits.txt; then finds the protein sequences and runs a MAFFT alignment.
"""

#imports needed to run the program
import sys
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Entrez
Entrez.email = "bsc10@student.le.ac.uk"
from Bio.Align.Applications import MafftCommandline


#user inputs protein sequence file in commandline
user_file = sys.argv[1]
if not user_file.endswith(".fasta"): #ensures a ".fasta" file is being used
    print("Wrong file type. Please use a .fasta file")
    exit()
record = SeqIO.read(user_file, format="fasta")
result_handle = NCBIWWW.qblast("blastp", "swissprot", record.format("fasta"))
with open("my_blast.xml","w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()


#parsing the BLAST result to find top ten results, and writing output to hits.txt
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)
output_blast = open("hits.txt","w")
seq_id = []
i = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        i = i+1
        if i <11:
            output_blast.write('****Alignment '+ str(i) +'****'+"\t")
            output_blast.write('Sequence Identifier: ' + alignment.accession +"\t")
            seq_id.append(alignment.accession) # list seq_ids to find the sequence in NCBI
            output_blast.write("Score: "+ str(hsp.score)+"\t")
            output_blast.write('E-Value: '+ str(hsp.expect)+"\t")
            sequence_identity = round(((hsp.identities/alignment.length)*100),2)
            output_blast.write("Sequence Identity: " + str(sequence_identity)+ "\n")
result_handle.close()
output_blast.close()


#retrieving sequences from NCBI databases
out_fasta = open("hits.fasta","w")
for element in seq_id:
    handle = Entrez.efetch(db="protein", id=element, rettype="fasta",retmode="text")
    out_fasta.write(handle.read())
out_fasta.close()


#aligning sequences using MAFFT
out_mafft = open("msa.fasta","w")
mafft_cline = MafftCommandline(input = "hits.fasta")
stdout, stderr = mafft_cline()
out_mafft.write(stdout)
out_mafft.close()
