# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:36:45 2022

@author: neil_
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import pandas as pd
import csv
#import re
import os
from sys import argv


##Load metadata file
#metadata = pd.read_excel(argv[1])
os.chdir(r'C:\Users\neil_\Desktop\School\M. BINF\BINF Project\Project\Pipeline')
metadata = pd.read_excel('ecoli_plasmid_database_sk.xlsx')
metadata = metadata[['Run', 'epi_type']]
metadata['epi_type'] = metadata['epi_type'].map({'clinical': 'H', 'environmental/other': 'N'})


## Assess presence-absence variation from txt files
with open('present_and_absent.txt') as f:
    isolate_list = [line.rstrip().lstrip(">") for line in f]
with open('present.txt') as f:
    present_list = [line.rstrip().lstrip(">") for line in f]
       
absent_list = set(isolate_list)^set(present_list)


## Build dict where keys are isolates and values are metadata for each isolate
I_dict={}
V_dict={}
with open("aa.fasta", "w") as output_handle:
    #for rec in SeqIO.parse(argv(4), 'fasta')
    for rec in SeqIO.parse('gene_variants.fasta', 'fasta'):
        Contig = rec.id
        Isolate = Contig.split("_")[0]
        NSeq = rec.seq
        PSeq = NSeq.ungap().translate()
        #PSeq = NSeq.translate()
        Pop =  metadata.query("Run == @Isolate")['epi_type'].values[0]
        SeqIO.write(SeqRecord(seq=PSeq, id=rec.id, description="translated sequence"), output_handle, 'fasta')
        if Isolate not in I_dict.keys():
            V_dict[Contig] = [NSeq, PSeq, Pop, Isolate]
        I_dict[Isolate] = Pop


for Isolate in absent_list:
    Pop =  metadata.query("Run == @Isolate")['epi_type'].values[0]
    I_dict[Isolate] = Pop


#pad all sequences to have same length for alignment
records = SeqIO.parse('aa.fasta', 'fasta')
records = list(records)
maxlen = max(len(record.seq) for record in records)
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

#alignment should 
output_file = '{}_padded.fasta'.format(os.path.splitext('aa.fasta')[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
alignment = AlignIO.read(output_file, "fasta")

#create list of all aa sequences for gene variants
PSeq_aln = []
Variants = []
for record in alignment:
    PSeq_aln.append(str(record.seq))
    Variants.append(record.id)
    
###Create binary matrices for feht
#create matrix where each column is a gene variant sequence of translated amino acids
df = pd.DataFrame({'sequences':PSeq_aln})
SAP = df['sequences'].apply(lambda x:pd.Series(list(x)))


binarymatrix = pd.DataFrame.transpose(pd.get_dummies(SAP))
binarymatrix.set_axis(Variants, axis=1,inplace=True)
with open('sav-data.txt', 'wt', newline= '') as out_file:
    binarymatrix.to_csv('sav-data.txt', sep='\t', mode='a')

#Create matrix where each row is an allele
allelicmatrix = pd.DataFrame.transpose(pd.get_dummies(df))
allelicmatrix.set_axis(Variants, axis=1,inplace=True)
with open('allelic-data.txt', 'wt', newline= '') as out_file:
    allelicmatrix.to_csv('allelic-data.txt', sep='\t', mode='a')
    
    
#Create 1-row matrix for PAV.
PA_data = []
for isolate in isolate_list:
    if isolate in present_list:
        PA_data.append(1)
    else:
        PA_data.append(0)

PAmatrix = pd.DataFrame({'Genome':isolate_list, 'Present':PA_data})
PAmatrix = PAmatrix.set_index('Genome').T
with open('pav-data.txt', 'wt', newline= '') as out_file:
    PAmatrix.to_csv('pav-data.txt', sep='\t', mode='a')


##Write SAV/SAP data to txt file
with open('metadata-SAV.txt', 'wt', newline='') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Contig', 'Population'])
    for contig in V_dict.keys():
        tsv_writer.writerow([contig, V_dict[contig][2]])

##Write PAV data to txt file   
with open('metadata-PAV.txt', 'wt', newline='') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Isolate', 'Population'])
    for Isolate in I_dict.keys():
        tsv_writer.writerow([Isolate, I_dict[Isolate]])

