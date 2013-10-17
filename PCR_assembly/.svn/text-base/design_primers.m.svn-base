function [primer_sequences ] = design_primers( sequence, min_Tm, NUM_PRIMERS, misprime_mode, tag );
% DESIGN_PRIMERS
%
% [primer_sequences ] = design_primers( sequence, tag );
%
%   or
%
% [primer_sequences ] = design_primers( sequence, min_Tm, NUM_PRIMERS, misprime_mode, tag );
%
%
% New dynamic programming (recursion) to figure out
%  DNA sequences to order for PCR assembly -- trying to
%  minimize:
%
%    Total length of primers   +  10 *  misprime_score( primer endpoint )
%
% The misprime_score is the maximum number of complementary nucleotide matches of the primer 3' end with
%   any other sequence segment in the desired sequence or its complement. G-C pairs are worth 
%    1.25 points, A-T are worth 1 point.
%
% Order (N) algorithm.  [But I need to work on speeding up pre-calculation of mispriming score, which is O(N^2) .]
%
% primer_sequences =  design_primers_NEW( sequence, min_Tm, NUM_PRIMERS, misprime_mode );
%
%  inputs:
%   sequence = sequence (use A,G,C,T). [required]
%
%  optional inputs:
%   min_Tm   = minimum melting temperature for any overlap (60 C, by
%                   default)
%   NUM_PRIMERS = desired number of primers (must be an even number). Leave blank or give 0,
%                    and script will figure this out.
%   misprime_mode = [default 0]
%                   0: use new mispriming score as described above.
%                  -1: no mispriming considerations
%                  >1: disallow any primer endpoint to be complementary to anywhere 
%                        else in the sequence at misprime_mode or more locations
%                  <0: use old-style primer 
%   output:
%    primer_sequences = the designed primers 
%                           [cell of strings -- type char(primer_sequences) for output. ]
%
% Note: oligos are assumed to be lengths between 15 and 60 nts.
%
% (C) Rhiju Das, 2012
if nargin == 0; help( mfilename ); return; end;

primer_sequences = {};
N_BP = length( sequence );

%Convert Us to Ts and lowercase to uppercase
sequence=RNA2DNA(sequence);
if exist( 'min_Tm' ) & ~isnumeric( min_Tm )
  tag = min_Tm; clear min_Tm;
end
if ~exist( 'min_Tm' )
  min_Tm = 60;
end
if exist( 'MAX_LENGTH' ) & ~isnumeric( MAX_LENGTH )
  tag = MAX_LENGTH; clear MAX_LENGTH;
end
if ~exist('MAX_LENGTH'); MAX_LENGTH = 60; end
if ~exist('MIN_LENGTH'); MIN_LENGTH = 15; end
if  ~exist( 'NUM_PRIMERS' ) | isempty( NUM_PRIMERS);  NUM_PRIMERS = 0; end;
if ~exist( 'misprime_mode' ); misprime_mode = 0; end;

DP_scores = zeros( N_BP, N_BP );
choice = zeros( N_BP, N_BP );
split_choice = cell( N_BP, N_BP );

tic
fprintf( 1, 'Precalculating Tm matrix ...' );
Tm_precalculated = precalculate_Tm( sequence );
fprintf( 1, ' done.\n' );
toc

MAX_NUM_PRIMERS = ceil( N_BP /  MIN_LENGTH );

allow_forward_endpoint = ones( 1, N_BP );
allow_reverse_endpoint = ones( 1, N_BP );

tic
%[ num_match_forward, num_match_reverse, best_match_forward, best_match_reverse ] = ...
%    check_for_mispriming_complete( sequence );
fprintf( 'Precalculating misprime score ...\n' );
[ num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse ] = ...
    check_for_mispriming_FAST( sequence );
toc

NEW_SCORE = 0; misprime_score_weight = 10.0;
if misprime_mode < -1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This was the old-style check for mispriming. Implemented here just to check
  % consistency with old code (design_primers_Nsquared.m).
  LENGTH_CHECK = abs(misprime_mode);
  [ allow_forward_endpoint, allow_reverse_endpoint ] = ...
      check_for_mispriming( sequence, LENGTH_CHECK );
elseif misprime_mode >= 0
  % new style:
  % If misprime_mode == 0, all endpoints are OK, but there is a 
  %    penalty quadratic in the number of misprime matches.
  %
  % If misprime_mode > 0, it specifies the number of misprime matches that are OK.
  %
  if ( misprime_mode > 0 )
    allow_forward_endpoint = (num_match_forward < misprime_mode ); 
    allow_reverse_endpoint = (num_match_reverse < misprime_mode );
  else
    % this is the new scoring scheme, hopefully default.
    % penalty is based on number of mismatches, with stronger weights to g/c pairings.
    NEW_SCORE = 1;
    %tic
    %fprintf( '\n precalculating misprime score -- this could be sped up\n' );
    %[ misprime_score_forward, misprime_score_reverse ] = get_misprime_score( sequence ); 
    %toc
  end
end

tic
fprintf( 'Doing dynamics programming calculation ...' );
% First set is special.
%  |                     p
%  ---------------------->
%                   ||||||
%                   <-----...
%                   q
%
num_primer_sets = (NUM_PRIMERS/2); % could be zero, meaning user does not know.
num_primer_sets_max = ceil( N_BP / MIN_LENGTH );

MAX_SCORE = N_BP * 2 + 1;  % maximum total length of primers...
if (NEW_SCORE) 
  MAX_SCORE = MAX_SCORE +  misprime_score_weight * max( [misprime_score_forward, misprime_score_reverse] ) * 2 * num_primer_sets_max;
end;

scores_start = MAX_SCORE * ones( N_BP, N_BP, num_primer_sets_max );
scores_stop  = MAX_SCORE * ones( N_BP, N_BP, num_primer_sets_max );
scores_final = MAX_SCORE * ones( N_BP, N_BP, num_primer_sets_max );

% used for backtracking:
choice_start_p = zeros( N_BP, N_BP, num_primer_sets_max );
choice_start_q = zeros( N_BP, N_BP, num_primer_sets_max );
choice_stop_i = zeros( N_BP, N_BP, num_primer_sets_max );
choice_stop_j = zeros( N_BP, N_BP, num_primer_sets_max );


% basic setup -- first primer.
for p = MIN_LENGTH:MAX_LENGTH  % STOP[forward](1)
  if ~allow_forward_endpoint(p); continue;end;
  
  % STOP[reverse](1)
  q_min = max( 1, p - MAX_LENGTH + 1 );
  q_max = p;

  for q = q_min:q_max
    if ~allow_reverse_endpoint(q); continue;end; 
  
    if ( Tm_precalculated( q, p ) ) > min_Tm
      scores_stop( p,q,1) =  (q-1) +  2 * (p - q + 1 );

      if (NEW_SCORE) scores_stop( p,q,1) = scores_stop( p,q,1) + misprime_score_weight * ( misprime_score_forward(p) +  misprime_score_reverse(q) ); end;
    
    end
    
  end
end

%colormap( jet( MAX_SCORE ) );
%image( scores_stop(:,:,1) ); pause

best_min_score = MAX_SCORE;
n = 1;
while ( n <= num_primer_sets_max  )
  
  % final scoring -- let's see if we can 'close' at the end of the sequence.
  % 
  %                 p
  %  --------------->
  %            ||||||
  %            <---------------------
  %            q                    N_BP
  %
  %
  for p = 1:N_BP % STOP[forward]
    if ~allow_forward_endpoint(p); continue;end;
    
    q_min = max( 1, p - MAX_LENGTH + 1 );  
    q_max = p;
    for q = q_min:q_max % STOP[reverse]
      if ~allow_reverse_endpoint(q); continue;end;
      
      if ( scores_stop(p,q,n) < MAX_SCORE ) % previous primer ends had overlap with good Tm and were scored
	i = N_BP + 1;
	j = N_BP;
	last_primer_length = j - q + 1;
	if last_primer_length <= MAX_LENGTH &  last_primer_length >= MIN_LENGTH
	  scores_final( p, q, n ) = scores_stop(p,q,n) + (i - p - 1);
	
	  if (NEW_SCORE) scores_final( p,q,n) = scores_final( p,q,n) + misprime_score_weight * ( misprime_score_forward(p) +  misprime_score_reverse(q) ); end;

	end
      end    
      
    end
  end

  min_score = min( min ( scores_final(:,:,n) ) );

  if  (min_score < best_min_score | n == 1 )
    best_min_score = min_score;
    best_n = n;
  end
  
  if ( n >= num_primer_sets_max ) break; end;
  if ( num_primer_sets > 0 & n == num_primer_sets ) break; end;

  % considering another primer set
  n = n + 1;

  % 
  %        p              i
  %  ------>              ------ ... ->
  %    |||||              ||||||
  %    <------------------------
  %    q                       j
  %
  for p = 1:N_BP % STOP[forward](1)

    q_min = max(1, p - MAX_LENGTH + 1); 
    q_max = p;
    for q = q_min:q_max % STOP[reverse](1)  

      if ( scores_stop(p,q,n-1) < MAX_SCORE ) % previous primer ends had overlap with good Tm and were scored
      
	% START[reverse](1)
	min_j = max( (p+1), q + MIN_LENGTH - 1 );
	max_j = min( N_BP,  q + MAX_LENGTH - 1 ) ;
	for j = min_j:max_j
	  %if (misprime_mode>0) & ~allow_forward_endpoint(j); continue;end; % at some PCR stage this will be an endpoint!	  
	  %START[reverse](2)
	  min_i = max( (p+1), j - MAX_LENGTH + 1 );
	  max_i = j;
	  for i = min_i:max_i
	    %if (misprime_mode>0) & ~allow_reverse_endpoint(i); continue;end; % at some PCR stage this will be an endpoint!

	    if ( Tm_precalculated( i, j ) ) > min_Tm 
	      
	      potential_score = scores_stop(p,q,n-1) + (i - p - 1) + 2 * (j - i + 1 );

	      if ( potential_score < scores_start(i,j,n-1)  )
		scores_start(i,j,n-1) = potential_score;
		choice_start_p(i,j,n-1) = p;
		choice_start_q(i,j,n-1) = q;
	      end
	      
	    end
	    
	  end
	end
      end
    end
  end

  %image( scores_start(:,:,1) ); pause

  % 
  %             i                     p
  %             ---------------------->
  %             ||||||           ||||||
  %  <----------------           <----- ...
  %                  j           q
  %
  
  %START[reverse](1)
  for j = 1:N_BP 
    
    %START[reverse](2)
    min_i = max( 1, j - MAX_LENGTH + 1 );
    max_i = j;
    for i = min_i:max_i;  % could also just make this 1:N_BP, but that would waste a little time.
      
      if ( scores_start(i,j,n-1) < MAX_SCORE ) % reasonable starting point
	
	% STOP[reverse](1)
	min_p = max( (j+1), i + MIN_LENGTH - 1 );
	max_p = min(  N_BP, i + MAX_LENGTH - 1 );
	for p = min_p:max_p 
	  if ~allow_forward_endpoint(p); continue;end; 
	  
	  %STOP[reverse](2)
	  min_q = max( (j+1), p - MAX_LENGTH + 1 );
	  max_q = p;
	  for q = min_q:max_q
	    if ~allow_reverse_endpoint(q); continue;end; 

	    if ( Tm_precalculated( q, p ) ) > min_Tm 
	      
	      potential_score = scores_start(i,j,n-1) + (q - j - 1) + 2 * (p - q + 1 );

	      if (NEW_SCORE) potential_score = potential_score + misprime_score_weight * ( misprime_score_forward(p) +  misprime_score_reverse(q) ); end;

	      if ( potential_score < scores_stop(p,q,n)  )
		scores_stop(p,q,n) = potential_score;
		choice_stop_i(p,q,n) = i;
		choice_stop_j(p,q,n) = j;
	      end
	      
	    end
	  end
	end
      end
    end
  end

    
end


%image( scores_stop(:,:,2) ); pause
fprintf( 1, ' done.\n' );
toc

if ( num_primer_sets > 0 )
  n = num_primer_sets;
else
  n = best_n;
end

% backtrack
[y, idx ] = min( scores_final(:,:,n) );
[min_score, q ] = min( y );
p = idx( q );

if ( min_score == MAX_SCORE )
  primer_sequence = {};
  fprintf( '\n\nNo solution found\n' );
  return;
end

primers( :, 2*n ) =  [q, N_BP, -1];

for m = n: -1 :2
  i = choice_stop_i(p,q,m);
  j = choice_stop_j(p,q,m);
  primers( :, 2*m-1) =   [i, p, 1 ];

  p = choice_start_p(i,j,m-1);
  q = choice_start_q(i,j,m-1);
  primers( :, 2*m-2 ) =  [q, j, -1];
end

primers( :, 1) =   [1, p, 1 ];

if ~exist( 'tag' ); tag = 'primer'; end;
primer_sequences = output_primers( primers, sequence, tag );


% mispriming "report"
CUTOFF = 3;
for m = 1:size(primers,2)
  
  if ( primers(3,m) == 1 ) % forward dir
    endpos =  primers(2,m);

    if ( num_match_forward( endpos ) > CUTOFF )
      problem_primers = find_primers_affected( primers, best_match_forward( endpos ) );
      fprintf( 'WARNING: Primer %d can misprime with %d-residue overlap to position %d, which is covered by primers: ', m, num_match_forward(endpos), best_match_forward( endpos )  );
      for q = 1:length( problem_primers); fprintf( ' %d', problem_primers(q) ); end;
      fprintf( '\n' );
    end
    
  else
    
    endpos =  primers(1,m);

    if ( num_match_reverse( endpos ) > CUTOFF )
      problem_primers = find_primers_affected( primers, best_match_reverse( endpos ) );
      fprintf( 'WARNING: Primer %d can misprime with %d-residue overlap to position %d, which is covered by primers: ', m, num_match_reverse(endpos), best_match_reverse( endpos )  );
      for q = 1:length( problem_primers); fprintf( ' %d', problem_primers(q) ); end;
      fprintf( '\n' );
    end
  
  end
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function primer_list = find_primers_affected( primers, pos );
primer_list = [];
for n = 1:size( primers, 2 )
  if ( pos >= primers(1,n) ) & ( pos <= primers(2,n) )
    primer_list = [ primer_list, n ];
  end
end


 