function [ num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse ] = ...
    check_for_mispriming_FAST( sequence, MATCH_REVERSE_COMPLEMENT );
% [ num_match_forward, num_match_reverse, ...
%   best_match_forward, best_match_reverse, ...
%   misprime_score_forward, misprime_score_reverse ] = ...
%    check_for_mispriming_FAST( sequence, MATCH_REVERSE_COMPLEMENT );
%
% Fast algorithm that takes advantage of N log N sort to 
% figure out closest match of each subsequence to other subsequences
% in the DNA or in the reverse complement.
%
% The 'misprime_score' sums up 1 for A/T matches and 1.25 for G/C matches.
%

N = length( sequence );
m = 20; % length of string subsets
subset_strings = {};
% match to sequence
for i = 1:N 
  startpos = max(i - m,1);
  endpos(i) = i;
  dir(i) = 1;
  subset_strings{i} = sequence( i:-1:startpos ); 
end

% match to reverse_complement of sequence
sequence_rc = reverse_complement( sequence );
for i = 1:N % length of sequence
  endpos(N+i) = i;
  dir(N+i) = -1;
  
  endpos_rc = N + 1 - i;
  startpos_rc = max( endpos_rc - m, 1 );
  subset_strings{N+i} = sequence_rc( endpos_rc:-1:startpos_rc ); 
end

[strings_sorted, sortidx] = sort( subset_strings );

% how close is match to neighbor?
for i = 1:(2*N-1)

  count = 0; misprime_score = 0;

  string1 = strings_sorted{i};
  string2 = strings_sorted{i+1};

  while ( count < length( string1 ) & count < length( string2 ) )
    if (string1(count+1) ~= string2(count+1) ) break; end;
    count = count + 1; 

    if ( string1(count) == 'G' | string1(count) == 'C' ) 
      misprime_score = misprime_score + 1.25;
    else
      misprime_score = misprime_score + 1;
    end
  
  end;

  match_to_next(i) = count;
  misprime_score_to_next(i) = misprime_score;

end

% compare both neighbors
match_max(1)   = match_to_next(1);
best_match(1) = 2;
misprime_score_max(1) = misprime_score_to_next(1);

match_max(2*N) = match_to_next(2*N-1);
best_match(2*N) = 2*N-1;
misprime_score_max(2*N) = misprime_score_to_next(2*N-1);

for i = 2:2*N-1;  
  if ( match_to_next(i-1) > match_to_next(i) )
    best_match(i) = i-1;
    match_max(i) = match_to_next(i-1);
    misprime_score_max(i) = misprime_score_to_next(i-1);
  else
    best_match(i) = i+1;
    match_max(i) = match_to_next(i); 
    misprime_score_max(i) = misprime_score_to_next(i);
  end
end


for i = 1:(2*N)
  if ( sortidx(i) <= N )
    num_match_forward( sortidx(i) ) = match_max(i);
    misprime_score_forward( sortidx(i) ) = misprime_score_max(i);
    best_match_forward( sortidx(i) ) = mod( sortidx( best_match(i) )-1, N) + 1;
  else
    num_match_reverse( sortidx(i)-N ) = match_max(i);
    misprime_score_reverse( sortidx(i)-N ) = misprime_score_max(i);
    best_match_reverse( sortidx(i)-N ) = mod( sortidx( best_match(i) )-1, N) + 1;
  end
end

output_num_match( sequence, num_match_forward, num_match_reverse );

return;

