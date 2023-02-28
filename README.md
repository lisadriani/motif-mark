# motif-mark

This code will take a: 
* text file with up to five Motifs (degenerate OK)
* fasta file 
    * exons uppercase, introns lowercase
    * genes of interest up to 1000nt in length 

 
The output will be
* a .png file
    * takes the name of fasta file given to it 
        * i.e. if "Figure_1.fasta" is given, "Figure_1.png" will be the output
    * drawn to scale
    * key with color-to-motif legend 
    * fasta headers

#### How it Works: 

The code is object oriented, creating a list with all objects of Motifs, and a second list of all objects of gene sequences (called Fasta). 

* Motif Object: 
    * Information: 
        * sequence
        * color (determined by the order in the motif file)
        * regex expression (for degenerate base motifs)
    * Methods:
        * ambiguous()
            * Creates regex expression to help take care of degenerate bases and DNA to RNA changes 
        * find_motif()
            * Takes the ambiguous regex expression and the gene sequence and returns a dictionary : dict[regex] = [start,end,color]
* Fasta Object: 
    * Information: 
        * sequence
        * length of the sequence
        * the header from the fasta file 
    * Method: 
        * find_exon()
            * this will find the upppercase letters from the gene sequence and save the start and end positions to draw. 



Here is an example figure, attached in github, based off of the fasta file, Figure_1.fasta and Fig_1_motifs.txt. 
![Figure_1.png](https://github.com/lisadriani/motif-mark/blob/main/Figure_1.png)

