function [ num_match_forward, num_match_reverse, best_match_forward, best_match_reverse ] = ...
    check_for_mispriming_complete( sequence, MATCH_REVERSE_COMPLEMENT );

if ~exist( 'MATCH_REVERSE_COMPLEMENT' ); MATCH_REVERSE_COMPLEMENT = 1; end;

N = length( sequence );

[num_match_forward1, best_match_forward1] = get_num_match( sequence, sequence );
[num_match_forward2, best_match_forward2] = get_num_match( sequence, reverse_complement(sequence) );
[num_match_forward, best_match_forward] = get_max( num_match_forward1, num_match_forward2, best_match_forward1, N+1-best_match_forward2 );

[num_match_reverse1, best_match_reverse1] = get_num_match( reverse_complement(sequence), sequence );
[num_match_reverse2, best_match_reverse2] = get_num_match( reverse_complement(sequence), reverse_complement(sequence) );
[num_match_reverse, best_match_reverse] = get_max( num_match_reverse1, num_match_reverse2, best_match_reverse1, N+1-best_match_reverse2 );

if ~MATCH_REVERSE_COMPLEMENT
  num_match_forward = num_match_forward1;   best_match_forward = best_match_forward1;
  num_match_reverse = num_match_reverse2;   best_match_reverse = best_match_reverse2;
end

num_match_reverse = num_match_reverse( end:-1:1 );

output_num_match( sequence, num_match_forward, num_match_reverse );

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [max_match, best_match] = get_num_match( sequence, template_sequence );

nres = length( sequence );
nres_template = length( template_sequence );
max_match = zeros( 1, nres );
best_match = max_match;

for i = 1:nres
  for j = 1:nres_template

    if ( i ~= j )
      matching = 1; offset = 0;
      while matching & (i-offset >= 1) & (j-offset >= 1)
	if ( sequence(i-offset) ~= template_sequence(j-offset) ); 
	  matching = 0; break;
	end
	offset = offset + 1;
      end
      num_match = offset;
      if ( num_match > max_match(i) )
	max_match(i) = num_match;
	best_match(i) = j;
      end
    end
  
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num_match, best_match] = get_max( num_match1, num_match2, best_match1, best_match2 );

num_match  = 0 * num_match1;
best_match = 0 * best_match1;

for i = 1:length( num_match1 )
  if ( num_match1(i) > num_match2(i) )
    num_match(i) = num_match1(i);
    best_match(i) = best_match1(i);
  else
    num_match(i) = num_match2(i);
    best_match(i) = best_match2(i);
  end
end