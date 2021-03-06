function cod_tab = codon_table

% Genetic code table. Format in RNA but not DNA sequences. First row of
%   3-letter amino acid names, second row of 1-letter amino acid names,
%   third to eigth row of codons. Stop codons are named 'X_X' and '~'.
%
% by T47, Oct 2013.
%

cod_tab = {...
    'Ala','A','GCU', 'GCC', 'GCA', 'GCG', '', '';...
    'Arg','R','CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG';...
    'Asn','N','AAU', 'AAC', '', '', '', '';...
    'Asp','D','GAU', 'GAC', '', '', '', '';...
    'Cys','C','UGU', 'UGC', '', '', '', '';...
    'Gln','Q','CAA', 'CAG', '', '', '', '';...
    'Glu','E','GAA', 'GAG', '', '', '', '';...
    'Gly','G','GGU', 'GGC', 'GGA', 'GGG', '', '';...
    'His','H','CAU', 'CAC', '', '', '', '';...
    'Ile','I','AUU', 'AUC', 'AUA', '', '', '';...
    'Leu','L','UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG';...
    'Lys','K','AAA', 'AAG', '', '', '', '';...
    'Met','M','AUG', '', '', '', '', '';...
    'Phe','F','UUU', 'UUC', '', '', '', '';...
    'Pro','P','CCU', 'CCC', 'CCA', 'CCG', '', '';...
    'Ser','S','UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC';...
    'Thr','T','ACU', 'ACC', 'ACA', 'ACG', '', '';...
    'Trp','W','UGG', '', '', '', '', '';...
    'Tyr','Y','UAU', 'UAC', '', '', '', '';...
    'Val','V','GUU', 'GUC', 'GUA', 'GUG', '', '';...
    'X_X','~','UAA', 'UGA', 'UAG', '', '', ''};

