function mut_seq = get_mutant_sequence (wt_seq, offset, mutations)

%
% mut_seq = GET_MUTANT_SEQUENCE (wt_seq, offset, mutations)
%
% Returns mutated sequence.
%
% Input
% =====
%   wt_seq              The wild-type sequence, format in string.
%   offset              Offset for sequence numbering. Offset is natural 
%                          numbering minums final numbering. Format in 
%                          double, default is 0. e.g., A1 in sequence,
%                          final numbering is 50, offset is -49.
%   mutations           Set of mutatants. Format as cell. Each construct is
%                          specified by string annotation, e.g. 'A200C'. 
%                          Multiple mutations within one construct is
%                          recognied by a cell of strings, e.g. {'A200C',
%                          'G201T'}. Numbering includes offset.
%
% Output
% ======
%   mut_seq             The mutant sequence, format in string.
%
%
% by T47, Oct 2013.
%

mut_seq = wt_seq;
for i = 1:length(mutations)
    ind = str2num(mutations{i}(2:end-1));
    org_char = mutations{i}(1);
    mut_char = mutations{i}(end);
    if strcmp(mut_seq(ind-offset), org_char)
        mut_seq(ind-offset) = mut_char;
    else
        fprintf(['PROBLEM with ', mutations{i}, '.\n']);
    end;
end;