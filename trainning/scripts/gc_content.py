#!usr/bin/env python


def get_base_count(sequence, nucleotide):
    seq = sequence.upper()
    base = seq.count(nucleotide)
    return base

def get_gc_content(sequence):
    for i in range(0, len(seq)):
        n = get_base_count(seq, 'N')
        g = get_base_count(seq, 'G')
        c = get_base_count(seq, 'C')
        total = round((g + c) / (len(seq) - n) * 100, 2)
        return total

def analysis_window_sized_sequence(sequence, function, window_size=10, step=5):
    '''Return a iterator that yelds(start, end, property) tuple.
    start and end used to slice de sequence and property is the result of the
    function
    Ex.
    sequence = "attagcgcaatctaactacactactgccgcgcggcatatatttaaatata"
    for start, end, gc in sliding_window_analysis(sequence, gc_content):
        print(start, end, gc)'''
    sequence = sequence.upper()
    for start in range(0, len(sequence), step):
        end = start + window_size
        if end > len(sequence):
            break
        yield  start, end, function(sequence[start:end])

def get_gc_content_analysis(iter_result):
    for star, end, gc in iter_result:
        print(star,'\t', end, '\t', gc)

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
            middle = (start + end) / 2
            row = '{}, {}, {}, {}\n'.format(start, middle, end, gc)
            fhand.write(row)

##seq = "attagcgcaatctaactacactactgccgcgcggcatatatttaaatata"
filename = '/home/paulo/Documents/IntroductionToComputationPython/Sco.dna'
seq = parse_genome(filename)
gc_content = analysis_window_sized_sequence(seq, get_gc_content, 100000, 50000)
analysis = results_writer(gc_content)
