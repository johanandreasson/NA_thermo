function [ allow_forward_endpoint, allow_reverse_endpoint ] = ...
    check_for_mispriming( sequence, LENGTH_CHECK );

if ~exist( 'LENGTH_CHECK' ) LENGTH_CHECK = 4; end;

% In the future, can do this based on deltaG's at temperature of interest!
allow_forward_endpoint = check_it( sequence, LENGTH_CHECK );

allow_reverse_endpoint = check_it( reverse_complement(sequence), LENGTH_CHECK );
allow_reverse_endpoint = allow_reverse_endpoint( end:-1:1 );

N_BP = length( sequence );


COL_SIZE = 150;

for k=1:(floor(N_BP/COL_SIZE)+1)

  allow_forward_line = '';
  allow_reverse_line = '';
  sequence_line = '';
  
  for n = 1:COL_SIZE
    pos = (k-1) * COL_SIZE+n;
    
    if (pos < N_BP)
      sequence_line(n) = sequence(pos);
      if allow_forward_endpoint(pos)
	allow_forward_line(n) = ' ';
      else
	allow_forward_line(n) = 'X';
      end
      
      if allow_reverse_endpoint(pos)
	allow_reverse_line(n) = ' ';
      else
	allow_reverse_line(n) = 'X';
      end
    end
  end
        
  fprintf('%s\n%s\n%s\n\n',allow_forward_line,sequence_line, allow_reverse_line);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function allow_endpoint = check_it( sequence, LENGTH_CHECK );

N_BP = length( sequence );

allow_endpoint = ones(1,N_BP);

for n = LENGTH_CHECK:N_BP
  for k = LENGTH_CHECK:N_BP

    if ( n~=k & ...
	 sequence( n-LENGTH_CHECK+[1:LENGTH_CHECK] ) == ... 
	 sequence( k-LENGTH_CHECK+[1:LENGTH_CHECK] ) )
      allow_endpoint( n ) = 0;
      break;
    end
  
  end
end

