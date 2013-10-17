sequence_P4P6 = ...
    'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC';

% In this tiling, the overall cost for the single mutants
% is less (by a few hundred dollars), but the assembly is more
% complicated -- 10 primers instead of 6. 
% Also, there are more "leftover" primer sets -- less than 24 would
% end up on a plate, and IDT requires a minimum of 24 to get the
% cheap price! 
% [Ann has checked: final assembled DNA/RNA is indistinguishable in 
% quality from assembly A, below.]
construct_name = 'P4P6-AssB';
primers = {'TTCTAATACGACTCACTATAGGCCAAAACAACGGAAT','GGCTGTTGACCCCTTTCCCGCAATTCCGTTGTTTTGGCCT','GAAAGGGGTCAACAGCCGTTCAGTACCAAGTCT','GGCCATCTCAAAGTTTCCCCTGAGACTTGGTACTGAACGGC','GGGAAACTTTGAGATGGCCTTGCAAAGGGTAT','AGGACCATGTCCGTCAGCTTATTACCATACCCTTTGCAAGGCCA','CTGACGGACATGGTCCTAACCACGC','GTTGACTTAGGACTTGGCTGCGTGGTTAGGACCATG','CAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTT','GTTGTTGTTGTTGTTTCTTTGGTTTGGTTTTGAACTGCATCCATATCAACAGAAG'};


construct_name = 'P4P6-AssA';
primers = {'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCG','GGCCATCTCAAAGTTTCCCCTGAGACTTGGTACTGAACGGCTGTTGACCCCTTTCCCG','GGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGAC','GTTGACTTAGGACTTGGCTGCGTGGTTAGGACCATGTCCGTCAGCTTATTACCATAC','CAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTT','GTTGTTGTTGTTGTTTCTTTGGTTTGGTTTTGAACTGCATCCATATCAACAGAAG'};
 



how_many_variants = zeros(1,length(sequence_P4P6));
how_many_variants(33:192) = 3;

% To get the numbers to work out right, order wt + 159 mutants.
% Skip the very first 'G' in the standard P4-P6 construct.
which_muts = [103:261]; % Tetrahymena ribozyme numbering.
offset = 33 - 102;  % to skip over T7 promoter, buffers on construct.
which_libraries = [ 1 , 2, 3 ];
[sequences_to_order, construct_names] = single_mutant_library( primers, sequence_P4P6, offset, which_muts, which_libraries,construct_name );

output_sequences_to_order_96well_diagram( sequences_to_order, primers, construct_name );

