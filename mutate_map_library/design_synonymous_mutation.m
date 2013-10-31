function [syn_names, syn_seqs] = design_synonymous_mutation (sequence, seq_offset, cod_offset, std_dir, mode_flag, cod_mut, total_mut)

%
% [syn_names, syn_seqs] = DESIGN_SYNONYMOUS_MUTATION (sequence, ...
%                           [seq_offset], [cod_offset], [std_dir], ...
%                           [mode_flag], [cod_mut], [total_mut]);
%
% Enumerates possible synonymous mutations for a given RNA sequence.
%
% Input
% =====
%   sequence                The wild-type sequence, format in string.
%   [seq_offset]            Offset for sequence numbering. seq_offset is 
%                               natural numbering minums final numbering. 
%                               Format in double, default is 0. e.g., A1 in 
%                               sequence, final numbering is 50, offset is 
%                               -49.
%   [cod_offset]            Offset for codon triplet numbering. cod_offset
%                               is the natural numbering of the 1st nucelotide 
%                               in 1st codon from the 5' of the sense-strand 
%                               (regardless of STD_DIR). Format in double,
%                               default is 0.
%   [std_dir]               Flag of strand direction. 1 as SEQUENCE is 
%                               sense-strand for ORF, -1 as antisense-strand
%                               for ORF. Format in double, default is 1.
%   [mode_flag]             Flag for mode choices, format in cell of strings, 
%                               default is {'123','only'}.
%                           '123' denotes allowed numbers of nucleotide 
%                               mutations within each codon. e.g., use '1' 
%                               if only single-nucleotide change is allowed,
%                               use '23' if only double and triple
%                               nucleotide changes are allowed.
%                           'only/full' denotes multi-codon mode. Use
%                               'only' if only single-codon mutants are 
%                               allowed, use 'full' for a full scan of 
%                               multi-codon mutants list.
%   [cod_mut]               Collection for allowed numbers of codons to
%                               mutate. Format in double array, default is 
%                               [1:4], only valid under 'full' mode. e.g.,
%                               use [2,4:6] if 2 or 4 or 5 or 6 codons are
%                               allowed to mutate at the same time.
%   [total_mut]             Collection for allowed numbers of total point
%                               mutations. Format in double array, default
%                               is 0, only valid under 'full' mode. 0 means
%                               unlimited. e.g. use [1:9] if any numbers 
%                               under 10 of nucleotide mutations are 
%                               allowed.
%
% Output
% ======
%   syn_names               Mutant list of labels of all mutants. Format in
%                               cell.
%   syn_seqs                Mutated sequence list of SYN_NAMES.
%
% WARNING: 'full' mode can be extremely slow. Please proceed with caution.
%
%
% by T47, Oct 2013.
%


if nargin == 0; help( mfilename ); return; end;

if ~exist('seq_offset','var') || isempty(seq_offset); seq_offset = 0; end;
seq_offset = -seq_offset;
if ~exist('cod_offset','var') || isempty(cod_offset); cod_offset = 0; end;
if ~exist('std_dir','var') || isempty(std_dir) || (std_dir ~= 1 && std_dir ~= -1); std_dir = 1; end;

if ~exist('mode_flag','var') || isempty(mode_flag); mode_flag = {'123','only'}; end;
% read flag
mode_str = lower([mode_flag{1:end}]);
multi_flag = [];
if strfind(mode_str, '1'); multi_flag = [multi_flag, 1]; end;
if strfind(mode_str, '2'); multi_flag = [multi_flag, 2]; end;
if strfind(mode_str, '3'); multi_flag = [multi_flag, 3]; end;
if strfind(mode_str, 'full');
    full_flag = 1;
else
    full_flag = 0;
end;

if ~exist('cod_mut','var') || isempty(cod_mut); cod_mut = 1:4; end;
if ~isinteger(cod_mut); cod_mut = round(cod_mut); end;
if sign(cod_mut) == -1; cod_mut = 0; end;
if ~exist('total_mut','var') || isempty(total_mut); total_mut = 0; end;
if ~isinteger(total_mut); total_mut = round(total_mut); end;
if sign(total_mut) == -1; total_mut = 0; end;

% parse all codons
cod_tab = codon_table;
[~, codons] = show_ORF(sequence, cod_offset, std_dir);

syn_names = {};
syn_seqs = {};
syn_codons = {};
group_codons = {};
group_temp = {};
% index of output mutations cell
mut_count = 0;
for i = 1:length(codons)
    
    % get list of target mutation codons
    codon_mutations = find_codon_from_codon(codons{i}, cod_tab);
    for j = 1:length(codon_mutations)
        name_count = 0;
        mut = {};
        % make the mutation set cell of names
        for k = 1:3
            if ~strcmp(codon_mutations{j}{1}(k),codons{i}(k));
                name_count = name_count + 1;
                % numbering for strand direction
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
        
        % screen only the mutants specified by mode_flag
        if ismember(length(mut), multi_flag);
            if full_flag;
                group_temp = [group_temp, {mut}];
            else
                mut_count = mut_count + 1;
                syn_names(mut_count) = {mut};
                syn_seqs{mut_count} = get_mutant_sequence(sequence, seq_offset, mut);
                for k = 1:length(codons)
                    if k == i;
                        syn_codons(mut_count,k) = codon_mutations{j};
                    else
                        syn_codons(mut_count,k) = codons(k);
                    end;
                end;
            end;
        end;
    end;
    
    % group choices of codons
    if full_flag;
        group_codons = [group_codons, {group_temp}];
        group_temp = {};
    end;
end;

% if full_flag;
%     group_codons = group_codons(~cellfun('isempty',group_codons));
%     % choices of each codon
%     size_all = [];
%     for i = 1:length(group_codons)
%         size_all(i) = length(group_codons{i}) + 1;
%     end;
%     tag_all = zeros(size_2,length(size_all));
%     for i = 1:size_2
%         for j = 1:length(size_all)
%             tag_all(i,j) = mod(floor(i/prod(size_all(j+1:end))), size_all(j));
%         end;
%     end;
%     % tag_all is all posible combinations of codons
%     tag_all = tag_all(1:end-1,:)
%
%     % go over each combination
%     mut_count = 0;
%     for i = 1:size_2-1
%         tag_sub = tag_all(i,:);
%         if length(tag_sub(tag_sub ~= 0)) > max_cod_mut; continue; end;
%         mut = {};
%         for j = 1:length(tag_sub)
%             if tag_sub(j) ~= 0;
%                 mut = [mut, group_codons{j}{tag_sub(j)}];
%             end;
%         end;
%         if max_total_mut == 0 || length(mut) <= max_total_mut;
%             mut_count = mut_count + 1;
%             syn_names(mut_count) = {mut};
%             syn_seqs{mut_count} = get_mutant_sequence(sequence, seq_offset, mut);
%         end;
%     end;
% end;

if full_flag;
    tic;
    
    % create codon choice combination list
    group_codons = group_codons(~cellfun('isempty',group_codons));
    cod_list = [];
    for i = cod_mut
        cod_list_sub = combnk(1:length(group_codons), i);
        cod_list_sub(:, size(cod_list_sub,2)+1:max(cod_mut)) = 0;
        cod_list = [cod_list; cod_list_sub];
    end;
    % calculate total of mutants to be examined
    count = 0;
    for i = 1:size(cod_list,1)
        tag_size = [];
        for j = 1:size(cod_list,2)
            if cod_list(i,j) ~= 0;
                tag_size = [tag_size, length(group_codons{cod_list(i,j)})];
            end;
        end;
        count = count + prod(tag_size);
    end;
    toc; fprintf('\n');
    fprintf(['Number of point mutations for each codon: ', num2str(multi_flag),'\n']);
    fprintf(['Number of codon mutations: ', num2str(cod_mut),'\n']);
    fprintf(['Number of total oint mutations: ', num2str(total_mut),'\n']);
    fprintf('A total of '); fprintf(2, [num2str(count),' ']); fprintf('mutants will be calculated.\n');
    fprintf('Press any key to continue...\n');
    pause;
    
    % examine each combiantion
    syn_names = cell(1, count);
    syn_seqs = cell(1, count);
    ind = 0; count = 0; reverseStr = repmat(' ',1,29); tic;
    for i = 1:size(cod_list,1)
        % the number of possible choices for each chosen codon
        tag_size = [];
        for j = 1:size(cod_list,2)
            if cod_list(i,j) ~= 0;
                tag_size = [tag_size, length(group_codons{cod_list(i,j)})];
            end;
        end;
        % the combination of choices for chosen codon set
        tag_all = zeros(prod(tag_size),length(tag_size));
        for j = 1: size(tag_all,1)
            for k = 1: size(tag_all,2)
                tag_all(j,k) = mod(floor(j/prod(tag_size(k+1:end))), tag_size(k)) + 1;
            end;
        end;
        tag_all = tag_all([end,1:end-1],:);
        
        % list all
        for j = 1: size(tag_all,1)
            % ind for writing index, count for examination numbering
            ind = ind + 1;
            count = count + 1;
            if mod(count, ceil(size(syn_names,2)/ 10^4)) == 0 || count == size(syn_names,2);
                reverseStr = lprintf(reverseStr,['Calculating ', num2str(count), ' of ', num2str(size(syn_names,2)), ' ...\n']);
            end;
            
            mut = {};
            for k = 1:size(tag_all,2)
                mut = [mut, group_codons{cod_list(i,k)}{tag_all(j,k)}];
            end;
            if ismember(length(mut), total_mut) | total_mut == 0;
                syn_names(ind) = {mut};
                syn_seqs{ind} = get_mutant_sequence(sequence, seq_offset, mut);
            else
                ind = ind - 1;
            end;
        end;
    end;
    
    fprintf('\n'); toc; fprintf('\n');
    % delete empty and flip
    syn_names = syn_names(~cellfun('isempty',syn_names))';
    syn_seqs = syn_seqs(~cellfun('isempty',syn_seqs));
end;


% output all mutants
fprintf('Output\n======\n');
% not display if multiple codons
if full_flag;
    fprintf('Too many entries for listing, not displayed here.\n');
else
    % add WT to first row because output_triplet_seq requires that
    syn_seqs_temp = [sequence, syn_seqs];
    syn_codons_temp = cell(size(syn_codons,1)+1,length(sequence(1+cod_offset: length(sequence)- mod((length(sequence)- cod_offset),3)))/3);
    syn_codons = [rot90(codons); syn_codons];
    for i = 1:size(syn_codons,1)
        for j = 1:size(syn_codons,2)
            syn_codons_temp(i,j) = syn_codons(i,j);
        end;
    end;
    output_triplet_seq(syn_names, syn_codons_temp, syn_seqs_temp, cod_offset, 0, std_dir);    
end;
fprintf(['Total = ', num2str(length(syn_names)),'; excluding WT.\n']);

