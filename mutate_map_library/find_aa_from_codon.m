function aa_str = find_aa_from_codon(codon_str, codon_table)
% find 3-letter AA name from input codon
aa_ind = mod(strmatch(codon_str,codon_table)+20,21)+1;
aa_str = codon_table{aa_ind,1};