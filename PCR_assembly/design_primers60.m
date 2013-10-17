function primers = design_primers60( sequence, tag );

sequence = strrep( sequence, 'U', 'T' );

primers{1} = sequence(1:60);

s = reverse_complement( sequence );
primers{2} = s(1:60);

fprintf( 1, '%s-F\t%s\n', tag, primers{1} )
fprintf( 1, '%s-R\t%s\n', tag, primers{2} )

