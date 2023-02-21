#!/usr/bin/env python

import re
import cairo
import math
import argparse
from itertools import product

def get_args():
    parser = argparse.ArgumentParser(description="This code will take a fasta file (with exons capitalized and introns lowercase), and a file of motifs and output a figure showing the motifs on an png.")
    parser.add_argument("-f",  "--fasta", help="designates absolute file path to the fasta file", required = True)
    parser.add_argument("-m", "--motifs", help = "designates absolute file path to the text file with motifs", required = True)
    return parser.parse_args()
	
args = get_args()

fasta_file = args.fasta
match = re.match(re.compile("^(.*)\..*"), args.fasta) #find the name of the file before the "."
file_name =  match.group(1) 
motif_file = args.motifs

def validate_base_seq(DNA:str, RNA_flag: bool=False):
    '''This will validate if a string that is given is a base sequence, returning True or False'''
    DNA = DNA.upper()
    thelength= len(DNA)
    return thelength == DNA.count("A") + DNA.count("G") + DNA.count("C") + DNA.count("U") + DNA.count("T")

def oneline_fasta(file):
    '''Takes a fasta file that has the sequence on multiple lines and rewrites the file with no new lines'''
    o = open("oneline", "r+")
    f = open(file, "r")
    build_seq = ""
    i = 0
    for line in f:
        if line.startswith(">") == True:
            header = line
            if i > 0:
                o.write(str(build_seq) + "\n")
            o.write(str(header))
            build_seq = ""
        if line.startswith(">") == False:
            line = line.strip('\n') 
            build_seq += line
            i += 1
    o.write(str(build_seq))
    return(o)

def draw_figure(motif,fasta, file_name):


    width,height =  int(fasta.length), 100

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width,height)
    context = cairo.Context(surface)

    context.set_line_width(1)
    ##create the intron line##
    context.move_to(0,50) #(x,y)
    context.line_to(int(fasta.length)),50)
    context.stroke()

    fasta.get_exon(fasta)
    ##creates a box precisely in the center. This will have to be edited depending on how many exons there are##
    context.rectangle(100,25,200,50)    #(x0,y0,x1,y1)
    context.fill()



    surface.write_to_png(file_name+".png")

class Motif:

    def __init__(self, the_seq,length):
        '''This is a motif'''

        ##the data##
        self.motif = the_seq
        self.length = len(the_seq)
        self.ambiguous_motifs = []


    ## the methods ##

    # def the_length(self, the_seq):
    #     self.length = len(the_seq)

    def find_motif(self, Fasta):
        pass

    def ambiguous(self,motif):
        ambiguous_bases = {"Y":["C","T","U"], "y":["c","t","u"], "g":["g"], "t":["t"],"a":["a"], "c":["c"],"u":["u"], "G":["G"], "T":["T"],"C":["C"],"U":["U"], "A":["A"]}
        self.ambiguous_motifs = list(map("".join, product(*map(ambiguous_bases.get, motif))))
        return list(map("".join, product(*map(ambiguous_bases.get, motif))))

class Fasta: 

    def __init__(self, the_seq, length, header):
        '''fasta line with introns/exons'''

        ##the data##
        self.fasta = the_seq
        self.length = len(the_seq)
        self.header = header

    ## the methods ##

    # def the_length(self, the_seq):
    #     self.length = len(the_seq)

    def find_exon(self,fasta):
        #prior = len(fasta)
        #fasta = re.split(r"(A-Z]+", fasta)
        #after = len(fasta)
        the_uppers = []
        for i in range(0,len(fasta)):
            if fasta[i] is.upper():
                the_uppers.append(i)

        return fasta

    def find_intron(self,fasta):
        pass

    def find_motif(self,fasta:list,motif:list):
        pass

def classify(fasta,motif_file):
    ''' takes fasta and motif file and puts everything into a list of classes. returns 2 lists, 1 of fastas classes and one of motif classes '''
    seqs = dict()
    fastas = []
    motifs = []
    motif_seqs = []

    motif_file = open(motif_file, "r")
    oneline = open("oneline", "r+")

    for line in oneline: #fill the classes
        line=line.strip()
        if line.startswith(">") == True:
            header = line
        else: 
            seqs[header] = line

    for header, seq in seqs.items(): #fille the list with assigning the classes
        fastas.append(Fasta(seq,len(seq), header))

    for line in motif_file: #fill the classes
        line = line.strip()
        motif_seqs.append(line)

    for motif in motif_seqs: # fill the list with assigning the classes
        motifs.append(Motif(motif,len(motif)))
    
    return fastas,motifs




oneline = oneline_fasta(fasta_file)
fastas,motifs = classify(oneline,args.motifs) #returns a list of classes. 


for motifseq in motifs: 
    ambig = motifseq.ambiguous(motifseq.motif)
    #print(motifseq.ambiguous_motifs)








# for motifseq in motifs:
#     print(motifseq.length)

# for fastaseq in fastas:
#     print(fastaseq.fasta)
#     print(fastaseq.length)
#     print(fastaseq.header)