function [mr_names, mr_seqs, non_names, non_seqs, mis_names, mis_seqs] = check_rescue_mutation (wt_sequence, wt_structure, seq_offset, all_names)

%
% [mr_names, mr_seqs, non_names, non_seqs, mis_names, mis_seqs] = ...
%   CHECK_RESCUE_MUTATION (wt_sequence, wt_structure, [seq_offset], all_names);
%
% Checks for mutate/rescue mutations from a given mutation list. Result
%   will be 3 groups:
%       mutate/rescue: mutation(s) present in helices and secondary structure
%                      rescued.
%       non-mutate: mutation(s) present only in loops and secondary
%                   structure not pertubed.
%       mutate/non-rescue: mutation(s) present in helices and secondry 
%                          structure disrupted.
%
% Input
% =====
%   wt_sequence             The wild-type sequence, format in string.
%   wt_structure            The dot-bracket annotation of secondary
%                               structure, format in string.
%   [seq_offset]            Offset for sequence numbering. seq_offset is 
%                               natural numbering minums final numbering. 
%                               Format in double, default is 0. e.g., A1 in 
%                               sequence, final numbering is 50, offset is 
%                               -49.
%   all_names               Mutation list, format in cell. 
%
% Output
% ======
%   mr_names               Mutant list of labels of all mutate/rescue 
%                               mutants. Format in cell.
%   mr_seqs                Mutated sequence list of MR_NAMES.
%   non_names               Mutant list of labels of all non-mutate 
%                               mutants. Format in cell.
%   non_seqs                Mutated sequence list of NON_NAMES.
%   mis_names               Mutant list of labels of all mutate/non-rescue 
%                               mutants. Format in cell.
%   mis_seqs                Mutated sequence list of MIS_NAMES.
%
% WARNING: Depends on the size of all_names, the run can be extremely slow. 
%          Please proceed with caution.
%
% NOTE: G-U pair is accepted as mutate/rescue.
%
%
% by T47, Oct 2013.
%

if nargin == 0; help( mfilename ); return; end;

if length(wt_sequence) ~= length(wt_structure);
    fprintf('WARNING: length of sequence and structure does not match!\n');
end;

if ~exist('seq_offset','var') || isempty(seq_offset); seq_offset = 0; end;
seq_offset = -seq_offset;
if ~exist('all_names','var') || isempty(all_names); return; end;

tic;
wt_sequence = strrep(wt_sequence, 'T','U');
fprintf('\n'); reverseStr = (' ');
for i = 1:length(all_names)
    if mod(i, ceil(length(all_names)/ 10^4)) == 0 || i == length(all_names);
        reverseStr = lprintf(reverseStr,['Checking mutation labels ', num2str(i), ' of ', num2str(length(all_names)), ' ...\n']);
    end;
    for j = 1:length(all_names{i})
        all_names{i}{j} = strrep(strrep(all_names{i}{j}, 'T', 'U'),'WU', 'WT');
    end;
end;
toc;

% get base-pair map
tic;
bps_all = convert_structure_to_bps( wt_structure ) + seq_offset;
base_pair = {'AU','GC','GU'};
fprintf('\n'); fprintf('\n'); reverseStr = (' ');

mr_names = {}; mr_seqs = {};
non_names = {}; non_seqs = {};
mis_names = {}; mis_seqs = {};

for i = 1:length(all_names);
    
    if mod(i, ceil(length(all_names)/ 10^4)) == 0 || i == length(all_names);
        reverseStr = lprintf(reverseStr,['Calculating ', num2str(i), ' of ', num2str(length(all_names)), ' ...\n']);
    end;
    
    is_pair_mut = 0;
    for j = 1:length(all_names{i})
        mut_char = all_names{i}{j};
        if find(str2num(mut_char(2:end-1)) == bps_all);
            is_pair_mut = 1;
            break;
        end;
    end;
    if is_pair_mut ==1;
        mut_seq = get_mutant_sequence(wt_sequence, seq_offset, all_names{i});
        is_rescue_mut = 1;
        for j = 1:size(bps_all,1)
            pair_char_1 = mut_seq(bps_all(j,1) -seq_offset);
            pair_char_2 = mut_seq(bps_all(j,2) -seq_offset);
            pair_char = [pair_char_1, pair_char_2];
            res1 = strfind(base_pair,pair_char);
            res2 = strfind(base_pair,fliplr(pair_char));
            res = [res1{1},res1{2},res1{3},res2{1},res2{2},res2{3}];
            if isempty(res);
                mis_names = [mis_names, {all_names{i}}];
                is_rescue_mut = 0;
                break;
            end;
        end;
        if is_rescue_mut; mr_names = [mr_names, {all_names{i}}]; end;
    else
        non_names = [non_names, {all_names{i}}];
    end;
end;

fprintf(2,'\n mutate/rescue mutants:\n');
if length(mr_names) <= 100;
    for i = 1:length(mr_names)
        mr_seqs{i} = get_mutant_sequence(wt_sequence, seq_offset, mr_names{i});
        fprintf([combine_mutant_label(mr_names{i}),'\n']);
    end;
else
    fprintf('Total of ');
    fprintf(2,[num2str(length(mr_names)), ' ']);
    fprintf('mutants.\n');
end;

fprintf(2,'\n non-mutate mutants:\n');
if length(non_names) <= 100;
    for i = 1:length(non_names)
        non_seqs{i} = get_mutant_sequence(wt_sequence, seq_offset, non_names{i});
        fprintf([combine_mutant_label(non_names{i}),'\n']);
    end;
else
    fprintf('Total of ');
    fprintf(2,[num2str(length(non_names)), ' ']);
    fprintf('mutants.\n');
end;

fprintf(2,'\n mutate/non-rescue mutants:\n');
if length(mis_names) <= 100;
    for i = 1:length(mis_names)
        mis_seqs{i} = get_mutant_sequence(wt_sequence, seq_offset, mis_names{i});
        fprintf([combine_mutant_label(mis_names{i}),'\n']);
    end;
else
    fprintf('Total of ');
    fprintf(2,[num2str(length(mis_names)), ' ']);
    fprintf('mutants.\n');
end;

fprintf('\n'); toc; fprintf('\n');
