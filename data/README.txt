Method:

Conservation of each position in MTS
Charge of peptide with and with our mutation
Isoelectric point
B62 score of mutated aa
Length of MTS
Presence of PTM
Run regression to predict the impact


Data:
For 10 species is ready after sorting and wrangling.


Data preprocessing:
 
Protein sequences from uniprot
get_fasta_from_big_fastafile.txt


Multiple sequence alignments with CLUSTAL)
clustalo -i /Users/aarslan/Desktop/MTSignal-Peptide-project/fasta_files/3HIDH_fasta.txt -o /Users/aarslan/Desktop/MTSignal-Peptide-project/msa/3HIDH1 --outfmt=clu
USED:
$ clustalo -i /Users/aarslan/Desktop/MTSignal-Peptide-project/fasta_files/3HIDH_fasta.txt -o /Users/aarslan/Desktop/MTSignal-Peptide-project/msa/3HIDH1 --outfmt=clu

Conservation with rate4site: (negative values, which are indicative of slowly evolving, conserved sites) are divided into 4.5 equal intervals. The same 4.5 intervals are used for the scores above the average (positive values, which are indicative of rapidly evolving, variable sites))
rate4site -s /Users/aarslan/Desktop/MTSignal-Peptide-project/msa/3HIDH -a 2 -o '3HIDH_MOUSE'

USED:
for i in /Users/aarslan/Desktop/MTSignal-Peptide-project/conservation/*; do rate4site -s $i -o ${i%_*} ;done

Precomputed isoelectric:
conservation_mouse4_precomputed

$for i in /Users/aarslan/Desktop/MTSignal-Peptide-project/conservation_mouse4/*; do sed -n '1p' $i > ${i%.*}"_iso.txt" ; done

$for i in /Users/aarslan/Desktop/MTSignal-Peptide-project/conservation_mouse4/*; do python3 /Users/aarslan/Desktop/MTSignal-Peptide-project/aaCharge.py $(sed -n '1p' $i) >> ${i%.*}"_iso.txt" ; done