function mut_label = parse_mutant_label( mut_str, delim)

%
% mut_label = PARSE_MUTANT_LABEL (mut_str, delim)
%
% Parse mutation label to a cell of each mutation. Wild-type ('WT') will be
%   output as empty {}.
%
% Input
% =====
%   mut_str             Mutation name in one string.
%   delim               Delimiter of mutations. Format in string, default 
%                           is ';'.
%
% Output
% ======
%   mut_label           Mutation sets input. Format as cell string.
%
%
% by T47, Oct 2013.
%

mut_label = {};
if strcmp(mut_str, 'WT');
    mut_label = {};
else
    mut_label = strread(mut_str, '%s', 'delimiter', delim)';
end;
