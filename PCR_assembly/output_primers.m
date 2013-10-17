function primer_sequences = output_primers( primers_all, sequence, tag );

if ~exist( 'tag' ) tag = 'primer'; end;

% Assume primers are alternating in directionality.

count = 0;
fprintf(1,'\n');

num_primers = size( primers_all,2);

blank_line = '';
COLWIDTH = 150;
for k = 1:max(length(sequence),COLWIDTH)
  blank_line(k) = ' ';
end

seq_line_prev = blank_line;

for j = 1:num_primers
  
  primers = primers_all( :, j );
  
  seq_line = blank_line;
  
  seg_start = primers(1);
  seg_end = primers(2);
  direction = primers(3);
  
  if ( direction == 1 )
    for k = seg_start:seg_end;
      seq_line(k) = sequence(k);
    end
    primer_sequences{j} = sequence( seg_start:seg_end );   
    if ( seg_end + 1 <= length(sequence ))
      seq_line( seg_end+1 ) = '-';
    end
    if ( seg_end + 2 <= length(sequence ))
      seq_line( seg_end+2 ) = '>';
    end
    text_out = num2str( j ) ;
    if ( seg_end + 2 + length(text_out)  <= length( sequence) )
      seq_line( seg_end + 2 + [1:length(text_out)] ) = text_out;
    end
  else 
    for k = seg_start:seg_end;
      seq_line(k) = reverse_complement( sequence(k) );
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
  
   
  bp_line = blank_line;
  overlap_seq = '';
  last_bp_pos = 1;
  for k = 1:length(sequence)
    if ( isempty(strfind( 'ACGT',seq_line_prev(k))) | ...
	 isempty(strfind( 'ACGT',seq_line(k))) )
      bp_line(k) = ' ';
    else
      bp_line(k) = '|';
      last_bp_pos = k;
      overlap_seq = [ overlap_seq, sequence(k) ];
    end
  end

  if (last_bp_pos > 1)
    Tm = calc_Tm( overlap_seq, 0.2e-6,0.1,0.0015);
    Tm_text = num2str( Tm, '%6.1f' );
      bp_line( last_bp_pos + [1:length(Tm_text)+1] ) = [' ',Tm_text];
  end

  bp_lines{j} = bp_line;
  seq_lines{j} = seq_line;

  %fprintf( '%s\n%s\n', bp_line, seq_line );

  seq_line_prev = seq_line;
end

OUTPUT_STAGGER = 1;
if ~OUTPUT_STAGGER

  fprintf('%s\n',sequence);
  for k = 1:length(seq_lines)
    fprintf( '%s\n%s\n',bp_lines{k}, seq_lines{k} );  
  end
else  

  blank_line = blank_line( 1:COLWIDTH);
  for n = 1: floor( (length(sequence)-1)/ COLWIDTH)+1
    start_pos = COLWIDTH*(n-1) + 1;
    end_pos   = min( COLWIDTH*(n-1) + COLWIDTH, length(sequence));
    out_line = [sequence(start_pos:end_pos)];
    fprintf( '%s\n', out_line);
    for k = 1:length(seq_lines)
      bp_line  = bp_lines{k} (start_pos:end_pos );
      seq_line = seq_lines{k}(start_pos:end_pos );
      fprintf( '%s\n%s\n',bp_line, seq_line );
    end
    fprintf( '\n\n');
  end
end

for j = 1:length( primer_sequences )
  fprintf( '%s-%d\t%s\n', tag, j, primer_sequences{j} );
end

