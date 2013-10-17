function primer_pos = order_primers_along_template( primers, template )

num_primers = length( primers );

for n = 1:num_primers
  i = strfind( template, primers{n} );

  if ~isempty( i ) 
    pos1 = i(1);
    pos2 = i(1)+length( primers{n} )-1;    
    dir = 1;
  else
    i = strfind( template, reverse_complement(primers{n}) );
    if isempty(i)
      fprintf( 1, 'Problem with primer %d: %s',n,primers{n} );
    end
    pos1 = i(1);
    pos2 = i(1)+length( primers{n} )-1;    
    dir = -1;
  end
  primer_pos(:,n) = [pos1 pos2 dir ];
end

output_primers( primer_pos, template );
