function codon_str = find_codon_from_codon(codon_org, cod_tab)
% find all synonymous codons of the input codon
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