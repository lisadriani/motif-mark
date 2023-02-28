#!/usr/bin/env python

import argparse
import re
import cairo


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
colors= [[0,0.9,1],[1,.6,1],[1,0,.4],[.5,.5,1],[.8,.1,.5]]


def oneline_fasta(file, file_name):
    '''Takes a fasta file that has the sequence on multiple lines and rewrites the file with no new lines'''
    f = open(file, "r") #open the file given
    build_seq = ""
    i = 0
    with open(file_name+"_oneline.fasta", "w+") as o: #open a new file to write the fasta file as oneline
        for line in f:
            if line.startswith(">") == True: 
                header = line
                if i > 0: #as long as it's not the first line of the file
                    o.write(str(build_seq) + "\n") #write out the sequence you've built
                o.write(str(header)) #write out the header line
                build_seq = "" #reset the sequence
            if line.startswith(">") == False: #when its not a header line
                line = line.strip('\n')  #strip the new line
                build_seq += line #add it to the sequence
                i += 1
        o.write(str(build_seq)) #write out the sequence at the end of the loop
    mylength = len(file_name) + len("_oneline.fasta")
    o = str(o)
    o = o[25:(25+mylength)]
    return(o)

def draw_figure(motifs,fastas, file_name):

    height = (len(fastas)*75+75) #make the height relative to the amount of fasta files
    width,height =  1000, height

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width,height)
    context = cairo.Context(surface)

    ii=40 #position of the motif boxes
    line_value = 50 #position of the intron/exon lines
    for fastaseq in fastas: #loop through the sequence objects

        ##create the intron line##
        context.set_line_width(1)

        context.move_to(0,(line_value)) #(x,y)
        context.line_to(fastaseq.length,(line_value))
        context.set_source_rgb(0,0,0)
        context.stroke()


        # get the start and end positions of the exons, to create the larger portions of exon
        exons = fastaseq.find_exon(fastaseq.fasta) 
        for start,end in exons.items():
            context.set_line_width(20)
            context.move_to(start+1,(line_value))
            context.line_to(end+1,line_value)
            context.stroke()
            ## write the header line as a title above the intron/exon line
            context.move_to(10,line_value-20)
            context.show_text(fastaseq.header)
            context.stroke()
        ## get the motif positions and colors. 

        for motifseq in motifs: ## loop through the motifs from the text file
            ##motif_dict[seq]:(start,end,color)
            motif_dict = motifseq.find_motif(fastaseq.fasta)
            for motif,integers in motif_dict.items():
                color1= motif_dict[motif][2][0]
                color2= motif_dict[motif][2][1]
                color3= motif_dict[motif][2][2]
                context.rectangle(integers[1]+1,ii, ((integers[0]+1)-(integers[1]+1)), 20)    #(x0,y0,x1,y1)
                context.set_source_rgb(color1,color2,color3)
                context.fill()
        ii+= 75
        line_value += 75

        ##### make a legend #####
        # write text

        context.set_source_rgb(0,0,0)
        context.move_to(35,height-25)
        context.show_text("Legend:")
        context.stroke()
        ## write in the colors and labels
        position = 100 
        for motifseq in motifs: 
            context.rectangle(position,height-33.5, 10, 10)    #(x0,y0,x1,y1)
            context.set_source_rgb(motifseq.color[0],motifseq.color[1],motifseq.color[2])
            context.fill()
            context.set_source_rgb(0,0,0)
            context.move_to(position + 15, height - 25)
            context.show_text(motifseq.motif)
            lengmo = len(motifseq.motif) + 75
            position += lengmo   ## increment by the length of the motif for pretty spacing
        ####make the box around it, line 1
        position = position + 20
        context.set_line_width(1)
        context.move_to(20, height-50)
        context.line_to(20, height-10)
        context.stroke()
        # line 2
        context.set_line_width(1)
        context.move_to(position, height-50 )
        context.line_to(position, height-10)
        context.stroke()
        # line 3
        context.set_line_width(1)
        context.move_to(20,height-50 )
        context.line_to(position, height-50)
        context.stroke()
        # line 4
        context.set_line_width(1)
        context.move_to(20,height-10 )
        context.line_to(position, height-10)
        context.stroke()
    surface.write_to_png(file_name+".png")
    return 

class Motif:

    def __init__(self, the_seq,length, color):
        '''This is a motif'''

        ##the data##

        self.motif = the_seq
        self.length = len(the_seq)
        self.ambiguous_motif = self.ambiguous(self)
        self.color = color


    ## the methods ##

    def find_motif(self, fasta): 
        '''fasta is the gene sequence, returns dictionary of the regex used: end, start, color'''
        regex = self.ambiguous_motif #get the regex to look for all the options in case of degenerate bases
        color = self.color #get the color since they're shared across the degenerate options
        motif_dict = dict()
        match = str()
        for match in re.finditer(rf"{regex}", str(fasta)): 
            if match != None:
                motif_dict[match.group()]=(match.end(),match.start(), color)

        return  motif_dict


    def ambiguous(self,motif):  
        '''this will give you a regex expression that can search for degenerate options '''
        ambiguous_bases = {"Y":"[CTUctu]", "y":"[CTUctu]", "G":"[Gg]","A":"[Aa]", "C":"[Cc]","U":"[UuTt]", "T":"[TtUu]", "W":"[AaTtUu]","S":"[CcGg]", "M":"[AaCc]", "K":"[GgTtUu]", "R": "[AaGg]", "B":"[CcGgTtUu]", "D":"[AaGgTtUu]", "H":"[AaCcTtUu]", "V":"[AaCcGg]", "N":"[AaCcGgTtUu]"}
        uppercase = self.motif.upper()
        finder = str()
        for letter in uppercase:
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

    def find_exon(self,fasta):
        the_start = []
        the_end = []
        exons = dict()
        ii = 0
        iii = 0
        for i in range(0,(len(fasta))):
            if fasta[i].isupper()== True:
                iii = 0
                if i-ii not in the_start:
                    the_start.append(i)
                ii +=1
                if i == len(fasta)-1:
                    the_end.append(i)
            if the_start != []:
                if fasta[i].isupper()== False:
                    if i - iii not in the_end:
                        the_end.append(i)
                        ii = 0
                    iii +=1
        for index, start in enumerate(the_start):
            exons[start] = the_end[index]
        return exons

def classify(oneline,motif_file):
    ''' takes fasta and motif file and puts everything into a list of classes. returns 2 lists, 1 of fastas classes and one of motif classes '''
    seqs = dict()
    fastas = []
    motifs = []
    motif_seqs = []

    motif_file = open(motif_file, "r")
    oneline = open(oneline, "r+")

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
    i = 0
    for motif in motif_seqs: # fill the list with assigning the classes
        motifs.append(Motif(motif,len(motif), colors[i]))
        i +=1
    return fastas,motifs


oneline = oneline_fasta(fasta_file, file_name) #makes the oneline file
fastas,motifs = classify(oneline,args.motifs) #returns a list of classes.
draw_figure(motifs,fastas,file_name) # draws the figure













######### FOR CHECKING #############



# for fastaseq in fastas:
#     for motifseq in motifs: 
#         ambig = motifseq.(motifseq.motif)
#        # print(motifseq.ambiguous_motif)

# for fastaseq in fastas:
#     pass
#     exons = fastaseq.find_exon(fastaseq.fasta)
#     #print(exons)
    #print(fastaseq.length)



# for fastaseq in fastas:
#     for motifseq in motifs:
#         motif_dict = motifseq.find_motif(fastaseq.fasta)





##############NOT USED ###############

# for fastaseq in fastas:
#     for motifseq in motifs:
#         #print(motifseq.ambiguous_motif)
#         motif_dict = motifseq.find_motif(fastaseq.fasta)
#         print(motif_dict)

#color_dict = assign_colors(fastas,motifs)

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



# def assign_colors(fastas,motifs):
#     for fastaseq in fastas:
#         for motifseq in motifs:
#             motif_dict, degenerate_motif = motifseq.find_motif(fastaseq.fasta)

#             for motif in degenerate_motif:
#                 if motif in color_dict:
#                     pass
#                 if motif not in color_dict:
#                     color_dict[motif]= sample(random_numbers, 3)
#     return color_dict
