function [syn_names, syn_seqs] = check_synonymous_mutation (wt_sequence, seq_offset, cod_offset, std_dir, all_names)

if ~exist('seq_offset','var') || isempty(seq_offset); seq_offset = 0; end;
if ~exist('cod_offset','var') || isempty(cod_offset); cod_offset = 0; end;
if ~exist('std_dir','var') || isempty(std_dir) || (std_dir ~= 1 && std_dir ~= -1); std_dir = 1; end;
if ~exist('all_names','var') || isempty(all_names); return; end;

cod_tab = codon_table;
show_ORF(wt_sequence, cod_offset, std_dir);

% prepare sequence and ORF
sequences = cell(length(all_names) + 1 ,1);
sequences{1} = wt_sequence;
cod_seqs = cell(length(all_names) + 1 ,1);
for i = 1:length(all_names)
    sequences{i+1} = design_mutants({}, wt_sequence, all_names(i), -seq_offset,'','', 0);
end;
for i = 1:length(sequences)
    if std_dir == -1;
        sequences{i} = reverse_complement(sequences{i});
    end;
    sequences{i} = strrep(sequences{i}, 'T', 'U');
    if iscell(sequences{i}); sequences{i} = sequences{i}{1}; end;
    cod_seqs{i} = sequences{i}(1+cod_offset: length(wt_sequence)- mod((length(wt_sequence)- cod_offset),3));
end;

codons = cell(length(cod_seqs),length(cod_seqs{1})/3);
for j = 1:size(codons,1)
    for i =1:3:length(cod_seqs{1})
        codons{j,(i-1)/3+1} = cod_seqs{j}(i:i+2);
        if strcmp(find_aa_from_codon(codons{j,(i-1)/3+1}, cod_tab),'X_X');
            break;
        end;
    end;
end;

% display all sequences in triplets
fprintf('Input\n=====\n');
for i = 1:length(sequences)
    if std_dir == -1;
        sequences{i} = reverse_complement(sequences{i});
    end;
end;
output_triplet_seq(all_names, codons, sequences, cod_offset, 1, std_dir);
fprintf(['Total = ', num2str(length(all_names)),'; excluding WT.\n']);

%%%%%%%%%%%%
syn_names = {};
syn_seqs = {};
syn_codons = {};
mut_count = 0;
mut_flag = 1;
for i = 2:size(codons,1)
    for j = 1:size(codons,2)
        if ~isempty(codons{i,j}) && ~isempty(codons{1,j});
            if ~strcmp(find_aa_from_codon(codons{i,j}, cod_tab), ...
                    find_aa_from_codon(codons{1,j}, cod_tab));
                mut_flag = 0;
                break;
            end;
        elseif ~(isempty(codons{i,j}) && isempty(codons{1,j}));
            mut_flag = 0;
            break;
        end;
    end;
    if mut_flag == 1;
        mut_count = mut_count + 1;
        syn_names(mut_count) = all_names(i-1);
        syn_seqs(mut_count) = sequences(i-1);
        for j = 1:size(codons,2)
            syn_codons(mut_count,j) = codons(i,j);
        end;
    end;
    mut_flag = 1;
end;

% display all sequences in triplets
fprintf('\nOutput\n======\n');
syn_seqs_temp = [wt_sequence, syn_seqs];
syn_codons_temp = [codons(1,:); syn_codons];
output_triplet_seq(syn_names, syn_codons_temp, syn_seqs_temp, cod_offset, 0, std_dir);
fprintf(['Total = ', num2str(length(syn_names)),'; excluding WT.\n']);

