function primer_sequences = output_RT_primers( all_primer_pos, sequence, primer_length );

count = 0;
fprintf(1,'\n');

num_primers = length( all_primer_pos );

blank_line = '';
COLWIDTH = 100;
for k = 1:max(length(sequence),COLWIDTH)
  blank_line(k) = ' ';
end

seq_line = blank_line;
bp_line  = blank_line;

for j = 1:num_primers
  
  seg_start = all_primer_pos( j );
  seg_end =  all_primer_pos( j ) + primer_length - 1;
  
  for k = seg_start:seg_end;
    seq_line(k) = strrep( reverse_complement( sequence(k) ), 'U', 'T' );
    bp_line(k) = '|';
  end

  primer_sequences{j} = reverse_complement(sequence( seg_start:seg_end )); 
  
  if ( seg_start - 1 >= 1 )
    seq_line( seg_start-1 ) = '-';
  end
  if ( seg_start - 2 >= 1 )
    seq_line( seg_start-2 ) = '<';
  end
  text_out = num2str( j ) ;
  if ( seg_start - 2 - length(text_out)  >= 1 )
    seq_line( seg_start - 3 - length(text_out) + [1:length(text_out)] ) = text_out;
  end
  
end


blank_line = blank_line( 1:COLWIDTH);
for n = 1: floor( (length(sequence)-1)/ COLWIDTH) + 1
  start_pos = COLWIDTH*(n-1) + 1; end_pos   = min( COLWIDTH*(n-1) + COLWIDTH, length(sequence));
  out_line = [sequence(start_pos:end_pos)];
  fprintf( '%s', out_line);
  fprintf( ' %d\n', end_pos );
  
  bp_line_out  = bp_line (start_pos:end_pos );
  seq_line_out = seq_line(start_pos:end_pos );
  fprintf( '%s\n%s\n',bp_line_out, seq_line_out );
  fprintf( '\n\n');
end


char(primer_sequences);

