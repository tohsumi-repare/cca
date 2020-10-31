Because of old gene name to new gene name collisions, the following genes 
keep their existing name due to name collisons:

TXNRD3NB (normally maps to TXNRD3)
C16orf47 (maps to ZFHX3)
C17orf47 (maps to SEPTIN4 - another gene, SEPT4 also maps to SEPTIN4).


Otherwise, gene names are updated to a recent set of standard names.
Once can compare the old name with the new by comparing the norm (new) files
with the non-norm in the Archive directory.


--------------------------------------------------------------------------

The below is if you need to use some library other than the included libraries.

The format of the mylibrary.fasta is very specific:
>mylibrary_gene_exon1_1
FASTA_SEQUENCE
...
where mylibrary is the name of the library and gene is the name of the gene.
The exon1 and 1 represent the exon number that the sgRNA is in and the 1
represents the first sgRNA in that gene.   These things can be changed to be
other things.   Note that each entry is separated by an underscore, so there
should be no underscores in the library name, gene name, exon (or whatever it
is), or final entry..


