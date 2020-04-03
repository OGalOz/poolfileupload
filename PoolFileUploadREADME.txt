This uploader assists in uploading Pool Files. Pool Files are TSV files with 12 
columns. The headers are, in order from left to right: "barcode", "rcbarcode",
"nTot", "n", "scaffold", "strand", "pos", "n2", "scaffold2", "strand2", "pos2",
"nPastEnd".
The pool file is part of the RBTnSeq (Random Barcode Transposon Sequencing)
analysis program. The app "MapTnSeq" also outputs a pool file. It is used 
as a reference to find the important reads after running transposon
sequencing. It is used as an input to the app "RBTnSeq Reads to Pool Counts".
The biological goal of the programs is to quantify the abundance of each 
barcoded strain in each sample.
For more information, consider the following article:
"Rapid Quantification of Mutant Fitness in Diverse Bacteria by 
Sequencing Randomly Bar-Coded Transposons" - Kelly M. Wetmore, Morgan N. Price, 
[...], Adam P. Arkin, and Adam Deutschbauer
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4436071/
