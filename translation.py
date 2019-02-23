#!usr/bin/env python

import csv
from pprint import pprint
from matplotlib import pyplot

def codon_table(filename):
    with open(filename, mode='r') as file:
        codons, aacs = [], []
        reader = csv.reader(file)
        for line in reader:
            codons.append(line[0])
            aacs.append(line[2])
        codon_table = dict(zip(codons, aacs))
        del codon_table['Codon']
    return codon_table

def get_code_tabe(genetic_codon_table):
    return codon_table(filename)

filename = 'codon_table.csv'
genetic_code = get_code_tabe(codon_table(filename))
#pprint(genetic_code)

def clean_dna(sequence):
    sequence = sequence.upper()
    alphabet = set('ACGT')
    seq = ''
    for nuc in sequence:
        if nuc in alphabet:
            seq += nuc
    return seq

def count_foreign_dna_bases(sequence):
    sequence = sequence.upper()
    alphabet = set('ACGT')
    for nuc in sequence:
        if nuc not in alphabet:
            counter = sequence.count(nuc)
            return counter

def is_rna(sequence):
    seq =  sequence.upper()
    rna_bases = set('ACGUacgu')
    return set(seq) <= rna_bases

def is_dna(sequence):
    seq =  sequence.upper()
    rna_bases = set('ACGTacgt')
    return set(seq) <= rna_bases

def get_base_count(sequence, nucleotide):
    seq = sequence.upper()
    base = seq.count(nucleotide)
    return base

def get_gc_content(sequence):
    seq = sequence.upper()
    for i in range(0, len(seq)):
        n = get_base_count(seq, 'N')
        g = get_base_count(seq, 'G')
        c = get_base_count(seq, 'C')
        total = round((g + c) / (len(seq) - n) * 100, 2)
        return total
    
def dna_reverse_complement(sequence):
    reverse_complement = ''
    complement = {'T': 'A', 'G': 'C', 'A': 'T', 'C': 'G'}
    for base in sequence[::-1]:
        reverse_complement += complement[base]
    return ''.join(reverse_complement)
        
def dna_transcribe(sequence):
    sequence = sequence.upper()
    return sequence.replace('T', 'U')

def translate_codon(triplet):
    assert is_rna(triplet)
    codon = triplet.upper()
    return genetic_code[codon]

def get_codons(sequence):
    rna = dna_transcribe(sequence.upper())
    assert is_rna(rna)
    codons = []
    for i in range(0, len(rna), 3):
         codon = rna[i:i+3]
         if len(codon) < 3:
             break
         codons.append(codon)
    return codons

def get_translation(sequence):
    rna = dna_transcribe(sequence.upper())
    stops = ['UAA', 'UAG', 'UGA']
    protein = ''
    length = len(seq) - (len(seq) % 3) - 1 
    for i in range(0, length, 3):
        codons = get_codons(rna)
    for codon in codons:
        aac = translate_codon(codon)
        if codon in stops:
            break
        else:
             protein += aac
    return protein

def get_translation_all_frames(sequence):
    rev_comp = dna_reverse_complement(sequence)
    rna = dna_transcribe(sequence.upper())
    rna_rev = dna_transcribe(rev_comp.upper())
    peptides_all_frames = {'+1': [],'+2': [],'+3': [],
                           '-1': [],'-2':[],'-3': []}
    for idx in range(0, 3):
        i_rna = rna[idx::]
        i_rna_rev = rna_rev[idx::]
        rna_translation = get_translation(i_rna)
        rna_rev_translation = get_translation(i_rna_rev)
        if idx == 0:
            peptides_all_frames['+1'] = rna_translation
            peptides_all_frames['-1'] = rna_rev_translation
        elif idx == 1:
            peptides_all_frames['+2'] = rna_translation
            peptides_all_frames['-2'] = rna_rev_translation
        elif idx == 2:
            peptides_all_frames['+3'] = rna_translation
            peptides_all_frames['-3'] = rna_rev_translation
    return peptides_all_frames

def get_peptide(sequence):
    peptides = get_translation_all_frames(sequence)
    for frame,  peptide in peptides.items():
        if peptide.startswith('M'):
            return peptide
    return max((len(val), key) for key, val in peptides.items())

def read_fasta_one_sequence(filename):
    sequence = ''
    with open(filename, 'r') as file:
        for line in file:
            if line[0] not in ['>', ';']:
                sequence += line.strip()
    return sequence

def analysis_window_sized_sequence(sequence, function, window_size=10, step=5):
    '''Return a iterator that yelds(start, end, property) tuple.
    start and end used to slice de sequence and property is the result of the
    function
    Ex.
    sequence = "attagcgcaatctaactacactactgccgcgcggcatatatttaaatata"
    for start, end, gc in sliding_window_analysis(sequence, gc_content):
        print(start, end, gc)'''
    for start in range(0, len(seq), step):
        end = start + window_size
        if end > len(seq):
            break
        yield  start, end, function(sequence[start:end])

def get_gc_content_analysis(iter_result):
    for star, end, gc in iter_result:
        return next(1), next(2), next(3)

def parse_genome(filename):
    with open(filename, 'r') as fhand:
        first_line = fhand.__next__() #scapes first line
        sequence = []
        for line in fhand:
            line_split = line.split()
            seq_str = ''.join(line_split[:-1]) #discard last elm list
            seq_lst = list(seq_str)
            sequence.extend(seq_lst)
        return ''.join(sequence)

def results_writer(analysis):
    with open('GC_content_Genome.csv', 'w') as fhand:
        header = 'start,middle,end,gc_content\n'
        fhand.write(header)
        for start, end, gc in analysis:
            print(gc)
            middle = (start + end) / 2
            print(middle)
            row = '{}, {}, {}, {}\n'.format(start, middle, end, gc)
            fhand.write(row)



#seq = read_fasta_one_sequence('human_notch.fasta')
# seq = 'ATGCCGCCGCTCCTGGCGTGATGATGACCCCTAGTAATGA'
# gc = analysis_window_sized_sequence(seq, get_gc_content, 1000, step=50)
# for start, end, gc in analysis_window_sized_sequence(seq, get_gc_content):
#     print(start, end, gc)
# pyplot.plot(gc)
# pyplot.xlabel("Window size (1000 nucleotides")
# pyplot.ylabel("GC content")
# pyplot.show()
##print(clean_dna(seq))
##print(count_foreign_dna_bases(seq))
# print(dna_reverse_complement(seq))
# print(translate_codon('Auu'))
# rna = dna_transcribe(seq)
# print(get_codons(seq))
#print(get_translation_all_frames(seq))
# print(get_peptide(seq))
#print(get_translation(seq))
# filename = '/home/paulo/Documents/IntroductionToComputationPython/Sco.dna'
# seq = parse_genome(filename)
# gc_content = analysis_window_sized_sequence(seq, get_gc_content, 100000, 50000)
# analysis = results_writer(gc_content)
