function codon_str = find_codon_from_codon(codon_org, cod_tab)

% codon_str = FIND_CODON_FROM_CODON (codon_org, cod_tab);
%
% Find all synonymous codons of the input codon.
%
% Input
% =====
%   codon_org               The input codon, in 3-letter string.
%   cod_tab                 Genetic code table.
%
% Output
% ======
%   codon_str               All synonymous codons, excluding the input one.
%                               Format in cell string.
%
%
% by T47, Oct 2013.
%

if ~exist('cod_tab','var'); cod_tab = codon_table; end;
codon_org = strrep(codon_org, 'T', 'U');

aa_ind = mod(strmatch(codon_org,cod_tab)+20,21)+1;
codon_str = {};
count = 0;
for i = 3:8
    if strcmp(cod_tab(aa_ind,i),'');
        return;
    else
        if ~strcmp(cod_tab(aa_ind,i), codon_org);
            count = count + 1;
            codon_str{count} = cod_tab(aa_ind,i);
        end;
    end;
end;