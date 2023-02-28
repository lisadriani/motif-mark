# motif-mark

This code will take a: 
* text file with up to five Motifs (degenerate OK)
* fasta file 
    * exons uppercase, introns lowercase
    * genes of interest up to 1000nt in length 

 
The output will be a .png that takes the name of fasta file given to it (if "Figure_1.fasta" is given, "Figure_1.png" will be the output). 
The image is to scale with introns, exons, and motifs, with label of fasta headers and color-to-motif legend. 

The code is object oriented, creating a list with all objects of Motifs, and a second list of all objects of gene sequences (called Fasta). 
The Motif object contains information like: the sequence, a color associated with it for the motif, a regex expression associated with degenerate base motifs. 
The Motif object has methods including creating the ambiguous regex expression (ambiguous()) and finding the motif in the gene sequence (find_motif()), returning the index position of the motifs. 
The Fasta object contains information like: the sequence, the length, and the header from the original fasta file. 
The Fasta object has a single method that finds the location of the exons (uppercase letters) in the gene sequences and returns the start and end positions of the Uppercase letters to draw the image correctly (find_exon()).



Here is an example figure, attached in github, based off of the fasta file, Figure_1.fasta and Fig_1_motifs.txt. 
![Figure_1.png](https://github.com/lisadriani/motif-mark/blob/main/Figure_1.png)

