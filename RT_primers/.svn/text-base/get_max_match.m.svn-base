function max_match = get_max_match( sequence, template_sequence );
% For each residue ,find the closest exact match site elsewhere in the sequence,
% and return number of matches.
% There's probably a way to make this really fast in MATLAB-ese.
if ~exist( 'template_sequence' ); template_sequence = sequence; end;
nres = length( sequence );
nres_template = length( template_sequence );
max_match = zeros( 1, nres );

for i = 1:nres
  for j = 1:nres_template

    if (i ~= j )
      matching = 1; offset = 0;
      while matching & (i+offset <= nres) & j+offset <= nres_template
	if ( sequence(i+offset) ~= template_sequence(j+offset) ); 
	  matching = 0;
	end
	offset = offset + 1;
      end
      num_match = offset - 1;
      max_match(i) = max( max_match(i), num_match );
    end
  
  end
end

