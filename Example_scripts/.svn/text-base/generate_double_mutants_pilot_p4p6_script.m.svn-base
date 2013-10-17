primers = {'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCG','GGCCATCTCAAAGTTTCCCCTGAGACTTGGTACTGAACGGCTGTTGACCCCTTTCCCG','GGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGAC','GTTGACTTAGGACTTGGCTGCGTGGTTAGGACCATGTCCGTCAGCTTATTACCATAC','CAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTT','GTTGTTGTTGTTGTTTCTTTGGTTTGGTTTTGAACTGCATCCATATCAACAGAAG'};

 
sequence = 'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC';

mutations = { {},{'C109G'},{'G212C'},{'C109G','G212C'},...
	         {'G110C'},{'C211G'},{'G110C','C211G'},...
	         {'C143G'},{'G160C'},{'C143G','G160C'},...
	         {'C109G','G160C'},... % control
	         {'G175C'},{'C165G'},{'C166G'},...
	         {'G175C','C165G'},{'G175C','C166G'} };


offset =  -89 +20;

mutfile = 'p4p6_dblmuts.txt';
fid = fopen(  mutfile, 'w' );

mut_label_start = 'P4P6-';
for i = 1:length( mutations )

  sequence_mut = sequence;
  mutation_set = mutations{i};
  mut_label = mut_label_start;
  
  for a = 1:length( mutation_set )
    mutation = mutation_set{a};
    seqpos = str2num( mutation(2:4) ) + offset;
    startchar = mutation(1);
    endchar = mutation(end);
    if ( sequence( seqpos ) ~= startchar ) fprintf( 'PROBLEM %s\n', mutation ); end;
    sequence_mut( seqpos ) = endchar;

    if ( a> 1 ); mut_label = [mut_label,';']; end;
    mut_label  = [mut_label, mutation];
  
  end
  
  if strcmp( mut_label, mut_label_start ); mut_label = [mut_label, 'WT']; end;
  
  fprintf( 1, '%s %s\n', sequence_mut, mut_label );
  fprintf( fid, '%s %s\n', sequence_mut, mut_label );
end

fclose( fid );

mutate_primers( primers, mutfile )
