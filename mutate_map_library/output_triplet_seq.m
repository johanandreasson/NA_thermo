function output_triplet_seq(all_names, codons, sequences, cod_offset, wt_flag, std_dir)
% display all sequences in triplets

cod_tab = codon_table;

if std_dir == -1;
    for i = 1: length(sequences)
        sequences(i) = strrep(reverse(sequences(i)),'T','U');
    end;
    for i = 1: size(codons,1)
        for j = 1:size(codons,2)
            codons{i,j} = strrep(complement(codons{i,j}),'T','U');
        end;
    end;
end;

mut_labels = cell(size(codons,1),1);
mut_labels{1} = 'WT';
max_label = 0;
for j = 2:size(codons,1)
    mut_labels{j} = combine_mutant_label(all_names{j-1});
    if length(mut_labels{j}) > max_label; max_label = length(mut_labels{j}); end;
end;
for j = 1+(1-wt_flag):size(codons,1)
    fprintf(['<strong>',mut_labels{j}, repmat(' ', 1,abs(length(mut_labels{j})-max_label)+4),'</strong>']);
    fprintf([sequences{j}(1:cod_offset),' ']);
    for i = 1:size(codons,2)
        if ~isempty(codons{j,i})
            if (strcmp(find_aa_from_codon(codons{j,i}, cod_tab),'X_X') && std_dir == 1) || ...
                    (strcmp(find_aa_from_codon(strrep(complement(codons{j,i}),'T','U'), cod_tab),'X_X') && std_dir == -1);
                clr = 2;
            else
                clr = 1;
            end;
            fprintf(clr,[codons{j,i},' ']);
        else
            fprintf(sequences{j}((cod_offset+3*(i-1)+1):...
                (length(sequences{j})- mod((length(sequences{j})- cod_offset),3))));
            break;
        end;
    end;
    fprintf([sequences{j}((length(sequences{j})- mod((length(sequences{j})- cod_offset),3)+1):end), '\n']);
end;


function mut_str = combine_mutant_label(mut_label)
% convert mutation label to string
mut_str = '';
for i = 1:length(mut_label)
    mut_str = [mut_str, mut_label{i},';'];
end;
if isempty(mut_str);
    mut_str = 'WT';
else
    mut_str = mut_str(1:end-1);
end;