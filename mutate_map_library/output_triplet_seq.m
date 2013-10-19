function output_triplet_seq(all_names, codons, sequences, cod_offset, wt_flag, std_dir)

%
% OUTPUT_TRIPLET_SEQ (all_names, codons, sequences, cod_offset, wt_flag, std_dir);
%
% Displays all sequences in triplets to the command line window. Stop
%   codon will be colored red. Input strand will be displayed, in RNA, not
%   DNA sequence. Wild-type should always be first entry in all_names and
%   sequences.
%
% Input
% =====
%   all_names               Set of mutatants. Format as cell. Each construct 
%                               is specified by string annotation, e.g. 
%                               'A200C'.  Multiple mutations within one 
%                               construct is recognied by a cell of strings,
%                               e.g. {'A200C', 'G201T'}. Numbering includes 
%                               offset. Use {{'WT'}} for wild-type.
%   codons                  Codons of mutants in all_names. Each row is one
%                               one mutant, each column is one amino acid
%                               codon. Will be converted to RNA sequence.
%   sequences               Sequences of mutants in all_names. Format as
%                               cell string. Will be converted to RNA
%                               sequence.
%   cod_offset              Offset for ORF. It could be any non-negative
%                               integers within sequence length. Always
%                               count from 5' of input sequence regardless
%                               of std_dir. For ORF starting from 1st nt,
%                               cod_offset = 0;
%   wt_flag                 Flag for whether wild-type will not be 
%                               displayed. 0 for no, 1 for yes.
%   std_dir                 Strand direction. 1 as input sequence is sense,
%                               -1 as anti-sense.
%
%
% by T47, Oct 2013.
%

cod_tab = codon_table;

for i = 1: length(sequences)
    if std_dir == -1;
        sequences(i) = strrep(reverse(sequences(i)),'T','U');
    else
        sequences(i) = strrep(sequences(i),'T','U');
    end;
end;
for i = 1: size(codons,1)
    for j = 1:size(codons,2)
        if std_dir == -1
            codons{i,j} = strrep(complement(codons{i,j}),'T','U');
        else
            codons{i,j} = strrep(complement(complement(codons{i,j})),'T','U');
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


