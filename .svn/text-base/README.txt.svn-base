#######################################################
NA_Thermo, the Das lab in-house primer design tools

(c) R. Das 2008-2013
#######################################################

Quick README by R. Das on Apr. 18, 2011. Send questions and comments to rhiju [at] stanford.edu
  [updated 22 June, 2012]

To design primers for your sequence, just follow these easy steps:

1. Start MATLAB.  Either add this 'na_thermo' directory to your path, or stay inside it.

2. Define your sequence. For example:

 sequence = 'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC';
 tag = 'P4P6';

This sequences includes a 20-nucleotide T7 promoter sequence at the beginning, and then a construct (starting with 'GG...') encoding the P4-P6 domain of the Tetrahymena ribozyme along with flanking sequences.

3. Type:

 primers = design_primers( sequence, tag );

4. This will compute primers with minimal length, annealing temperatures above a cutoff (60 C, by default), through a recursive strategy. An additional score term helps avoid primers that share multiple 3' nucleotides with other parts of the sequence, as a heuristic to reduce mispriming. The algorithm has similarities (but was developed independently) of Thachuk & Condon (2007), BIBE, Proc. of 7th IEEE International Conf., p. 123-130.

  If you want to use our original script (in use from 2008-2011), use design_primers_OLD:. (This was slower and did not rigorously optimize length.)


6. We usually copy/paste these to a Word or Excel document for easy look-up later. If you add primer labels, you can copy/paste these to the IDT website or wherever.

 

