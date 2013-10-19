function mut_str = combine_mutant_label(mut_label)

%
% mut_str = COMBINE_MUTANT_LABEL (mut_label)
%
% Convert mutation labels to one string. ';' as delimiter. Empty input ('')
%   is interpreted as 'WT'.
%
% Input
% =====
%   mut_label           Mutation sets input. Format as cell string.
%
% Output
% ======
%   mut_str             Mutation name in one string.
%
%
% by T47, Oct 2013.
%

mut_str = '';
for i = 1:length(mut_label)
    mut_str = [mut_str, mut_label{i},';'];
end;
if isempty(mut_str);
    mut_str = 'WT';
else
    mut_str = mut_str(1:end-1);
end;