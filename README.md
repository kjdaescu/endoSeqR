# endoSeqR
Package to identify endo-siRNAs from small RNA 

1) endo_Indexes sub-folder includes instructions on obtaining index files needed for alignment. The index file folders for each species (elegans, musculus, rerio, etc.) may be downloaded from https://drive.google.com/drive/folders/16xoPxbzVu9Zo9_gqF16umHO_swB91kxV. Alternatively, the indexes may be built with the genome, ncRNA, and cDNA fasta species files from ensembl with the "bowtie build" commands.

2) practice sub-folder includes a practice dataset, a pdf tutorial walking through the procedure (coming soon), and the expected output files.

Note that endo_MAP and clean_up in this folder are linux executables. Upon first use, it may be necessary to modify permissions, which requires going to this the directory and type the following in command line:

chmod u+x endo_MAP

chmod u+x clean_up
