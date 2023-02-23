#!/usr/bin/env python

import re
import cairo
import math
import argparse
from itertools import product
from random import sample

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

random_numbers = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]


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

def draw_figure(motif_list,fasta_list, color_dict, file_name):

    width,height =  1000, ((len(fasta_list)*50+100))

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width,height)
    context = cairo.Context(surface)

    ii=40
    line_value = 50
    for fastaseq in fasta_list:
        #print(fasta)
        context.set_line_width(1)
        ##create the intron line##
        context.move_to(0,(line_value)) #(x,y)
        context.line_to(fastaseq.length,(line_value))
        context.set_source_rgb(0,0,0)
        context.stroke()

        exons = fastaseq.find_exon(fastaseq.fasta)
        for start,end in exons.items():
            context.rectangle(start,ii, (end-start), 20)    #(x0,y0,x1,y1)
            context.set_source_rgb(0,0,0)
            context.fill()
        for motifseq in motif_list:
            motif_dict = motifseq.find_motif(fastaseq.fasta)
            for motif,integers in motif_dict.items():
                color1= color_dict[motif][0]
                color2= color_dict[motif][1]
                color3= color_dict[motif][2]
                context.rectangle(integers[1],ii, (integers[0]-integers[1]), 20)    #(x0,y0,x1,y1)
                context.set_source_rgb(color1,color2,color3)
                #context.rectangle()
                context.fill()
        ii+= 50
        line_value += 50

        ##make a legend
    # context.rectangle(20,ii, 960, 50)    #(x0,y0,x1,y1)
    # context.set_source_rgb(0,0,0)
    # context.fill()


    surface.write_to_png(file_name+".png")
    return 

class Motif:

    def __init__(self, the_seq,length):
        '''This is a motif'''

        ##the data##

        self.motif = the_seq
        self.length = len(the_seq)
        self.ambiguous_motif = self.ambiguous(self)


    ## the methods ##

    def find_motif(self, fasta):
        regex = self.ambiguous_motif
        motif_dict = dict()
        match = str()
        color_dict = dict()
        for match in re.finditer(rf"{regex}", str(fasta)):
            #print(match.group())
            if match != None:
                motif_dict[match.group()]=(match.end(),match.start())
                #color_dict[match.group()] = sample(random_numbers, 3)
            #if len(motif_dict)>1:
        return  motif_dict


    def ambiguous(self,motif):
        ambiguous_bases = {"Y":"[CTUctu]", "y":"[ctu]", "G":"[Gg]", "T":"[Tt]","A":"[Aa]", "C":"[Cc]","U":"[Uu]", "G":"[Gg]", "T":"[TtUu]", "W":"[AaTtUu]","S":"[CcGg]", "M":"[AaCc]", "K":"[GgTtUu]", "R": "[AaGg]", "B":"[CcGgTtUu]", "D":"[AaGgTtUu]", "H":"[AaCcTtUu]", "V":"[AaCcGg]", "N":"[AaCcGgTtUu]"}
        motif = self.motif.upper()
        finder = str()
        for letter in motif:
            finder+= ambiguous_bases[letter] 
        return finder

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
        the_start = []
        the_end = []
        exons = dict()
        ii = 0
        iii = 0
        for i in range(0,len(fasta)):
            if fasta[i].isupper()== True:
                iii = 0
                if i-ii not in the_start:
                    the_start.append(i)
                ii +=1
            if the_start != []:
                if  fasta[i].isupper()== False:
                    if i - iii not in the_end:
                        the_end.append(i)
                        ii = 0
                    iii +=1
        for end in the_end:
            for start in the_start:
                exons[start] = end

        return exons

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

for fastaseq in fastas:
    for motifseq in motifs: 
        ambig = motifseq.ambiguous(motifseq.motif)
       # print(motifseq.ambiguous_motif)

for fastaseq in fastas:
    pass
    exons = fastaseq.find_exon(fastaseq.fasta)
    #print(exons)
    #print(fastaseq.length)


color_dict=dict()
for fastaseq in fastas:
    for motifseq in motifs:
        #print(motifseq.ambiguous_motif)
        motif_dict = motifseq.find_motif(fastaseq.fasta)

        #print(motif_dict)
        #print(color_dict)

        #draw_figure(motifseq,fastaseq,file_name)



    ###if it's based off of every individual motif ###
        for motif in motif_dict:

            if motif in color_dict:
                pass
            if motif not in color_dict:
                #samps = sample(random_numbers,1)
                #print(samps)
                color_dict[motif]= sample(random_numbers, 3)


    #### if its based off of the motif file:  ####
    # for motif in motif_file:
    #     if motif in color_dict:
    #         pass
    #     if motif not in color_dict:
    #         color_dict[motif] = sample(random_numbers, 3)

draw_figure(motifs,fastas,color_dict,file_name)

# for fastaseq in fastas:
#     for motifseq in motifs:
#         #print(motifseq.ambiguous_motif)
#         motif_dict = motifseq.find_motif(fastaseq.fasta)
#         print(motif_dict)



# for motifseq in motifs:
#     print(motifseq.length)

# for fastaseq in fastas:
#     print(fastaseq.fasta)
#     print(fastaseq.length)
#     print(fastaseq.header)


        # for motifseq in motif_list:
        #     motif_dict = motifseq.find_motif(fasta)
        #     if motif_dict != None:
        #         for index, start, end in enumerate(motif_dict.items()):
        #             index = index/.1
        #             context.set_source_rgb(0.9, 0.1, 0.1)
        #             context.rectangle(start,ii,(end-start),20)
        #             context.fill()


                            # for start,end in motif_dict.items():
        #     context.set_line_width(1)
        #     ##create the intron line##
        #     context.move_to(0,(line_value)) #(x,y)
        #     context.line_to(fasta.length,(line_value))
        #     context.stroke()
        #     pass