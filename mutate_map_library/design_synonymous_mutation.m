function [syn_names, syn_seqs] = design_synonymous_mutation (sequence, seq_offset, cod_offset, std_dir)

if ~exist('seq_offset','var') || isempty(seq_offset); seq_offset = 0; end;
if ~exist('cod_offset','var') || isempty(cod_offset); cod_offset = 0; end;
if ~exist('std_dir','var') || isempty(std_dir) || (std_dir ~= 1 && std_dir ~= -1); std_dir = 1; end;

cod_tab = codon_table;
[~, codons] = show_ORF(sequence, cod_offset, std_dir);

syn_names = {};
syn_seqs = {};
syn_codons = {};
mut_count = 0;
for i = 1:length(codons)
    codon_mutations = find_codon_from_codon(codons{i}, cod_tab);
    for j = 1:length(codon_mutations)
        name_count = 0;
        mut = {};
        for k = 1:3
            if ~strcmp(codon_mutations{j}{1}(k),codons{i}(k));
                name_count = name_count + 1;
                if std_dir == 1
                    mut{name_count} = [codons{i}(k),...
                        num2str(seq_offset + cod_offset + (i-1)*3 + k),...
                        codon_mutations{j}{1}(k)];
                else
                    mut{name_count} = [strrep(complement(codons{i}(k)),'T','U'),...
                        num2str(seq_offset + length(sequence) - cod_offset - ((i-1)*3 + k -1)),...
                        strrep(complement(codon_mutations{j}{1}(k)),'T','U')];
                end;
            end;
        end;
        mut_count = mut_count + 1;
        syn_names(mut_count) = {mut};
        syn_seqs(mut_count) = design_mutants({}, sequence, syn_names(mut_count), -seq_offset, '','',0);
        for k = 1:length(codons)
            if k == i;
                syn_codons(mut_count,k) = codon_mutations{j};
            else
                syn_codons(mut_count,k) = codons(k);
            end;
        end;
    end;
end;

fprintf('Output\n======\n');
syn_seqs_temp = [sequence, syn_seqs];
syn_codons_temp = cell(size(syn_codons,1)+1,length(sequence(1+cod_offset: length(sequence)- mod((length(sequence)- cod_offset),3)))/3);
syn_codons = [rot90(codons); syn_codons];
for i = 1:size(syn_codons,1)
    for j = 1:size(syn_codons,2)
        syn_codons_temp(i,j) = syn_codons(i,j);
    end;
end;
output_triplet_seq(syn_names, syn_codons_temp, syn_seqs_temp, cod_offset, 0, std_dir);
fprintf(['Total = ', num2str(length(syn_names)),'; excluding WT.\n']);
