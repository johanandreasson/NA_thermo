function [ misprime_score_forward, misprime_score_reverse, best_match_forward, best_match_reverse ] = ...
    get_misprime_score( sequence, MATCH_REVERSE_COMPLEMENT );

if ~exist( 'MATCH_REVERSE_COMPLEMENT' ); MATCH_REVERSE_COMPLEMENT = 1; end;

N = length( sequence );

[misprime_score_forward1, best_match_forward1] = get_score( sequence, sequence );
[misprime_score_forward2, best_match_forward2] = get_score( sequence, reverse_complement(sequence) );
[misprime_score_forward, best_match_forward] = get_max( misprime_score_forward1, misprime_score_forward2, best_match_forward1, N+1-best_match_forward2 );

[misprime_score_reverse1, best_match_reverse1] = get_score( reverse_complement(sequence), sequence );
[misprime_score_reverse2, best_match_reverse2] = get_score( reverse_complement(sequence), reverse_complement(sequence) );
[misprime_score_reverse, best_match_reverse] = get_max( misprime_score_reverse1, misprime_score_reverse2, best_match_reverse1, N+1-best_match_reverse2 );


if ~MATCH_REVERSE_COMPLEMENT
  misprime_score_forward = misprime_score_forward1;   best_match_forward = best_match_forward1;
  misprime_score_reverse = misprime_score_reverse2;   best_match_reverse = best_match_reverse2;
end


misprime_score_reverse = misprime_score_reverse( end:-1:1 );

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
      allow_forward_line(n) = num2str(min(round(misprime_score_forward(pos)),9));
      allow_reverse_line(n) = num2str(min(round(misprime_score_reverse(pos)),9));
    end
  end
        
  fprintf('%s\n%s\n%s\n\n',allow_forward_line,sequence_line, allow_reverse_line);
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [max_misprime_score, best_match] = get_score( sequence, template_sequence );

nres = length( sequence );
nres_template = length( template_sequence );
max_misprime_score = zeros( 1, nres );
best_match = zeros(1, nres );

for i = 1:nres
  for j = 1:nres_template

    if ( i == j ); continue; end;
    
    matching = 1; offset = 0; misprime_score = 0;
    while matching & (i-offset >= 1) & (j-offset >= 1)
      if ( sequence(i-offset) ~= template_sequence(j-offset) ); 
	matching = 0; break;
      end
      
      misprime_score = misprime_score + 1;
      % C-G is worse than A-U
      if ( sequence(i-offset) == 'C' | sequence(i-offset) == 'G') misprime_score = misprime_score + 0.25; end;
      
      offset = offset + 1;
    end
    
    if ( misprime_score > max_misprime_score(i) )
      max_misprime_score(i) = misprime_score;
      best_match(i) = j;
    end
  
  end
end

% return quadratic score -- really penalize each additional match!
%max_misprime_score  = max_misprime_score.^2;


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