function mutate_primers( primers, muts );

if ~iscell( muts )
  [muts, sequence_names] = read_sequences( muts );
end

% squeezes out gaps:
sequence = get_primer_seq( muts{1}, 1, length( muts{1}), 1 );

primer_pos = order_primers_along_template( primers, sequence);

seq_map = get_mapping( muts{1}, sequence );


output_primers = {};
primer_names = {};

for k = 1:length( muts )
  fprintf( 2, '%s\n', sequence_names{k} );
  seq_mut = muts{k};
  for i = 1:length(primers)
    start_pos = seq_map( primer_pos(1,i) );  
    end_pos = seq_map( primer_pos(2,i) );
    dir = primer_pos(3,i);
    primer_seq = get_primer_seq( seq_mut, start_pos, end_pos, dir );

    if ~( strcmp( primer_seq, primers{i}) )
      primer_name = sprintf( '%s-%d',sequence_names{k}, i );
      fprintf( 1,'%s\t%s', primer_name,primer_seq);

      % signal if this primer is already in the order list...
      match_previous = '';
      for m = 1:length( output_primers )
	if strcmp( output_primers{m}, primer_seq )
	  match_previous = primer_names{m};
	  break;
	end
      end
      if length( match_previous ) > 0 
	fprintf( '-- matches %s', match_previous );
      else
	primer_names = [primer_names, primer_name ];
	output_primers = [output_primers, primer_seq ];
      end

      fprintf( '\n' );
    
    end
  end
  fprintf( 1, '\n' );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_map = get_mapping( seq1, sequence )

count = 0;
for i = 1:length( seq1)
  if ~(seq1( i ) == ' ' | seq1(i) == '-')
    count = count + 1;
    if ~( sequence( count ) == seq1( i ) )
      fprintf(1,['First mutant sequence should match overall sequence!\n']);
    end
    seq_map( count ) = i;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function primer_seq = get_primer_seq( seq_mut, start_pos, end_pos, dir );
primer_seq = '';
for i = start_pos:end_pos
  if ~(seq_mut(i)==' ' | seq_mut(i)=='-')
    primer_seq = [primer_seq, seq_mut(i) ]; 
  end
end
  
if (dir == -1)
  primer_seq = reverse_complement( primer_seq );
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [ sequences, sequence_names ] = read_sequences( filename );

fid = fopen( filename );
k=0;
while ~feof(fid)
  k=k+1;
  sequence_name = fgetl(fid);
  %fprintf(1,'%s\n',sequence_name );
  [sequence, remain ] = strtok(sequence_name);
  sequences{k} = sequence;
  
  if length(remain)<1
    sequence_names{k} = ['Sequence ',num2str(k)];
  else
    sequence_names{k} = strtok(remain);
  end
end

fclose(fid);
