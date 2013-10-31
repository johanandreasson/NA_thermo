function construct_names = design_rescue(sequence, structure, offset, mode_flag, seqpos_sub)

% construct_names = DESIGN_RESCUE (sequence, structure, [offset], ...
%                                   [mode_flag], [seqpos_sub]);
%
% Designs mutation/rescue double mutants on a construct based on structure
%   provided. Insertion or deletion not supported. Modes including swapping
%   /crossing base-pairs, single/double/triple base-pair mutations, include/not
%   single mutations for quartet. Output mutant names in cell.
% Wild-type sequence will always be in the first well. For G-U pairs,
%   swapping mode will mutate to C-G pairs instead of mismatch.
%
% =Input=
%   sequence            Sequence of WT, format in string.
%   structure           Dot-bracket annotation of secondary structure to be
%                           tested. Same length as sequence, format in
%                           string.
%   [offset]            Offset for sequence numbering. Offset is natural 
%                           numbering minums final numbering. Format in 
%                           double, default is 0. e.g., A1 in sequence,
%                           final numbering is 50, offset is -49.
%   [mode_flag]         Mode choices, format in cell of strings, default is
%                           {'single','swap','quartet'}.
%                       'single/double/triple' denotes base-pair mutants, 
%                           whether single or double or triple base-pair.
%                       'swap/cross' denotes nucleotide identity, whether
%                           swapping (A>T,T>A) or crossing (A>C,T>G).
%                       'only/quartet' denotes mutants library, whether
%                           only the double mutant, or quartet including
%                           single mutants as well.
%   [seqpos_sub]        Mutation region. Format in double array, default is
%                           [1 : length(sequence)]. Value includes offset
%                           already. e.g. [1:59, 69:80].
% =Output=
%   construct_names     Label of all mutants. Format in cell.
%
%
% by T47, Oct 2013.
%

if nargin == 0; help( mfilename ); return; end;

if length(sequence) ~= length(structure);
    fprintf('WARNING: length of sequence and structure does not match!\n');
end;

if ~exist('offset','var') || isempty(offset); offset = 0; end;
if ~exist('mode_flag','var') || isempty(mode_flag); mode_flag = {'single','swap','quartet'}; end;
if ~exist('seqpos_sub','var') || isempty(seqpos_sub) || length(seqpos_sub) < 2;
    seqpos_sub = [1:length(sequence)];
else
    seqpos_temp = [];
    count = 1;
    for i = 1:length(seqpos_sub)
        if seqpos_sub(i) + offset <= length(sequence) && seqpos_sub(i) + offset >= 1;
            seqpos_temp(count) = seqpos_sub(i) + offset;
            count = count + 1;
        end;
    end;
    seqpos_sub = seqpos_temp;
end;
sequence = strrep(sequence,'U','T');

% read flag
mode_str = lower([mode_flag{1:end}]);
if strfind(mode_str, 'cross'); 
    mut_flag = 0; 
else
    mut_flag = 1; 
end;
if strfind(mode_str, 'triple'); 
    multi_flag = 2; 
elseif strfind(mode_str, 'double');
    multi_flag = 1;
else
    multi_flag = 0; 
end;
if strfind(mode_str, 'quartet'); 
    qrt_flag = 1; 
else
    qrt_flag = 0; 
end;

% get base-pair map
bps_all = convert_structure_to_bps( structure );
count = 0;
% trim with seqpos_sub
for i = 1:size(bps_all,1)
    if ~isempty(find(seqpos_sub == bps_all(i,1),1)) && ~isempty(find(seqpos_sub == bps_all(i,2),1));
        count = count + 1;
        bps(count,1) = bps_all(i,1);
        bps(count,2) = bps_all(i,2);
    end;
end;
% if empty set then quit
if count == 0;
    fprintf('WARNING: no base-pairs found in the seqpos_sub region.\n');
    construct_names = {};
    return;
end;
constructs = cell((size(bps,1)-multi_flag)*(qrt_flag*2+1)+1, 1);
construct_names = cell((size(bps,1)-multi_flag)*(qrt_flag*2+1)+1, 1);

% fill in double mutations, no.1 is always WT
tic;
constructs{1,1} = {'WT'};
construct_names{1} = constructs{1,1};
for i = 1:(size(bps,1) - multi_flag)
    % store mutation format to cell, e.g. 'A150T'
    constructs{i+1,1} = [sequence(bps(i,1)),num2str(bps(i,1) - offset),sequence_mutate_table(sequence(bps(i,1)),mut_flag)];
    if multi_flag == 2;
        constructs{i+1,2} = [sequence(bps(i+1,1)),num2str(bps(i+1,1) - offset),sequence_mutate_table(sequence(bps(i+1,1)),mut_flag)];
        constructs{i+1,3} = [sequence(bps(i+2,1)),num2str(bps(i+2,1) - offset),sequence_mutate_table(sequence(bps(i+2,1)),mut_flag)];
        constructs{i+1,4} = [sequence(bps(i,2)),num2str(bps(i,2) - offset),sequence_mutate_table(sequence(bps(i,2)),mut_flag)];
        constructs{i+1,5} = [sequence(bps(i+1,2)),num2str(bps(i+1,2) - offset),sequence_mutate_table(sequence(bps(i+1,2)),mut_flag)];
        constructs{i+1,6} = [sequence(bps(i+2,2)),num2str(bps(i+2,2) - offset),sequence_mutate_table(sequence(bps(i+2,2)),mut_flag)];
    elseif multi_flag == 1;
        constructs{i+1,2} = [sequence(bps(i+1,1)),num2str(bps(i+1,1) - offset),sequence_mutate_table(sequence(bps(i+1,1)),mut_flag)];
        constructs{i+1,3} = [sequence(bps(i,2)),num2str(bps(i,2) - offset),sequence_mutate_table(sequence(bps(i,2)),mut_flag)];
        constructs{i+1,4} = [sequence(bps(i+1,2)),num2str(bps(i+1,2) - offset),sequence_mutate_table(sequence(bps(i+1,2)),mut_flag)];
    else
        constructs{i+1,2} = [sequence(bps(i,2)),num2str(bps(i,2) - offset),sequence_mutate_table(sequence(bps(i,2)),mut_flag)];
    end;
    % correct G-U pair mutations, i.e. G-U to C-G
    temp = check_GU_pair({constructs{i+1,1:end}});
    for j = 1: length(temp); constructs{i+1,j} = temp{j}; end;
    construct_names{i+1} = {constructs{i+1,1:(end)}};
end;

% fill in single mutations if asked
start_i = i+1;
if qrt_flag;
    for i = (start_i+1):2:(start_i + (size(bps,1) - multi_flag)*2-1)
        % num is index of constructs array
        num = (i-start_i+1)/2;
        if multi_flag == 2;
            constructs{i,1} = constructs{num+1,1};
            constructs{i,2} = constructs{num+1,2};
            constructs{i,3} = constructs{num+1,3};
            constructs{i+1,1} = constructs{num+1,4};
            constructs{i+1,2} = constructs{num+1,5};
            constructs{i+1,3} = constructs{num+1,6};
        elseif multi_flag == 1;
            constructs{i,1} = constructs{num+1,1};
            constructs{i,2} = constructs{num+1,2};
            constructs{i+1,1} = constructs{num+1,3};
            constructs{i+1,2} = constructs{num+1,4};
        else
            constructs{i,1} = constructs{num+1,1};
            constructs{i+1,1} = constructs{num+1,2};
        end;
        construct_names{i} = {constructs{i,1:(end-multi_flag-1)}};
        construct_names{i+1} = {constructs{i+1,1:(end-multi_flag-1)}};
    end;
    
    % reshuffle to quartet, i.e. 1=WT, 3n-1:3n+1=quartet(s,s,d)
    construct_temp = construct_names;
    quart_num = (length(construct_names)-1)/3;
    for i = 2: quart_num + 1
        construct_names((i-1)*3-1) = construct_temp((i-1)*2+quart_num);
        construct_names((i-1)*3) = construct_temp((i-1)*2+1+quart_num);
        construct_names((i-1)*3+1) = construct_temp(i);
    end;
end;

% show construct_names as output
for i = 1:size(construct_names,1);
    name_temp = [construct_names{i,:}];
    name = '';
    for j = 1:length(name_temp)
        name = strcat(name,';',name_temp{j});
    end;
    fprintf([name(2:end),'\n']);
end;

fprintf('\n'); toc; fprint('\n');


function mut_char = sequence_mutate_table(wt_char, mut_flag)

% get mutated nucleotide
% mut_flag: 1 means swapping A>T,T>A; 0 means cross A>C,T>G
num_char = strfind('ATCG', wt_char);
if mut_flag; 
    ref_char = 'TAGC';
else
    ref_char = 'CGAT'; 
end;
mut_char = ref_char(num_char);

function mut_out = check_GU_pair(mut_in)

% find 'GU', change to 'CG'
mut_out = mut_in;
interval = length(mut_in)/2;
for i = 1:interval
    res_old = strrep([mut_in{i}(1), mut_in{i+interval}(1)],'T','U');
    res_new = strrep([mut_in{i}(end), mut_in{i+interval}(end)],'T','U');
    if strcmp(res_old,'GU') && strcmp(res_new,'CA');
        mut_out{i+interval} = [mut_out{i+interval}(1:end-1), 'G'];
    elseif strcmp(res_old,'UG') && strcmp(res_new,'AC');
        mut_out{i} = [mut_out{i}(1:end-1), 'G'];
    end;
end;



