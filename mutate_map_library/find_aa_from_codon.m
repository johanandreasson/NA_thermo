function aa_str = find_aa_from_codon(codon_str, cod_tab)

% aa_str = FIND_AA_FROM_CODON (codon_str, cod_tab);
%
% Find amino acid of the input codon.
%
% Input
% =====
%   codon_str               The input codon, in 3-letter string.
%   cod_tab                 Genetic code table.
%
% Output
% ======
%   aa_str                  Amino acid, in 3-letter string.
%
%
% by T47, Oct 2013.
%

if ~exist('cod_tab','var'); cod_tab = codon_table; end;
codon_str = strrep(codon_str, 'T', 'U');

aa_ind = mod(strmatch(codon_str,cod_tab)+20,21)+1;
aa_str = cod_tab{aa_ind,1};