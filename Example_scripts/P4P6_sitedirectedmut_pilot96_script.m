sequence_P4P6 = ...
    'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC';

start_primer = 'TTCTAATACGACTCACTATAGGCCAAAACAACG';
end_primer = 'TTGTTGTTGTTGTTTCTTTGGTTTGGTTTTGAA';

% In future version of this script need to be very careful to also
% flag when new start/end primers need to be ordered.

% To get the numbers to work out right, order wt + 159 mutants.
% Skip the very first 'G' in the standard P4-P6 construct.
which_muts = [103:261]; % Tetrahymena ribozyme numbering.
offset = 33 - 102;  % to skip over T7 promoter, buffers on construct.
which_libraries = [ 1 , 2, 3 ];

% Quick test.
%which_muts = [103:105]; % Tetrahymena ribozyme numbering.
%offset = 33 - 102;  % to skip over T7 promoter, buffers on construct.
%which_libraries = [ 1  ];

construct_name = 'P4P6_Tail2';
design_mismatch_primers_library( sequence_P4P6, start_primer, end_primer, ...
				 offset, which_muts, ...
				 which_libraries, construct_name );


