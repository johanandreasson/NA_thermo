function all_names = combine_mutant_list (sub_names, if_wt)

%
% all_names = COMBINE_MUTANT_LIST (sub_names, [if_wt])
%
% Combines multiple mutant lists into one non-repeat list. All WT entries
%   will be removed and WT can be add back if asked by user.
%
% Input
% =====
%   sub_names           Mutant lists. Format in cell, each element is one
%                           mutations list. e.g. to combine a(1x5 mutants),
%                           b(10x1 mutants), and c(1x15 mutants), use
%                           sub_names = {a, b, c}; a, b, c could be output
%                           from design_synonymous_mutation or
%                           design_rescue.
%   [if_wt]             Flag of whether include WT as 1st in output. Format
%                           in double, 0 is no, 1 is yes. Default is 1. 
%
% Output
% ======
%   all_names           Mutant list of all non-repeat mutants from
%                           sub_names. Format in cell.
%
%
% by T47, Oct 2013.
%

if ~exist('if_wt','var') || isempty(if_wt); if_wt = 1; end;

% combine all entries in one cell
all_names = {};
count = 1;
for i = 1:length(sub_names)
    all_names(count:(count + length(sub_names{i}) -1)) = sub_names{i}(1:end);
    count = count + length(sub_names{i});
end;

% sort mutation labels in ascending
labels = cell(1,length(all_names));
for i = 1:length(all_names)
    for j = 1:length(all_names{i})
        all_names{i}{j} = all_names{i}{j}([2:end-1,1,end]);
    end;
    all_names{i} = sort(all_names{i});
    for j = 1:length(all_names{i})
        all_names{i}{j} = all_names{i}{j}([end-1,1:end-2,end]);
    end;
    labels{i} = combine_mutant_label(all_names{i});
end;

% erase repeated elements
for i = 1:length(labels)
    for j = i+1:length(labels)
        if strcmp(labels(i), labels(j)) || strcmp(labels(i), 'WT'); all_names{i} = ''; end;
    end;
end;
all_names = all_names(~cellfun('isempty',all_names));
all_names = all_names';
if if_wt; all_names = [{{'WT'}}; all_names]; end;