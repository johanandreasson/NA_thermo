function [primer_sequences, DP_scores, Tm_precalculated, choice, ...
	  split_choice] = design_primers_OLD( sequence, min_Tm, MAX_LENGTH, ...
					  how_many_variants,check_mispriming)
% Dynamic programming (recursion) to figure out
%  DNA sequences to order for PCR assembly -- trying to
%  minimize cost (proportional to length).
%
% primer_sequences = design_primers_OLD( sequence, min_Tm, MAX_LENGTH, ...
%					  how_many_variants,check_mispriming)
%
%  inputs:
%   sequence = sequence (use A,G,C,T). [required]
%
%  optional inputs:
%   min_Tm   = minimum melting temperature for any overlap (60 C, by
%                   default)
%   MAX_LENGTH = maximum oligo length (60 by default)
%   how_many_variants = vector containing number of mutation positions if
%                        the eventual goal is to make every single mutant 
%                        (zeros by default)
%   check_mispriming = disallow the primer to end on a pentamer that is
%                       repeated elsewhere in the strand.
%
%   output:
%    primer_sequences = the designed primers 
%                           [cell of strings -- type char(primer_sequences) for output. ]
%
% Rhiju Das --May 31, 2009

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
if ( ~exist('MAX_LENGTH') )
  MAX_LENGTH = 60;
end
if ( ~exist('how_many_variants') | isempty( how_many_variants) )
  how_many_variants = zeros( 1, N_BP );
end

DP_scores = zeros( N_BP, N_BP );
choice = zeros( N_BP, N_BP );
split_choice = cell( N_BP, N_BP );

tic
fprintf( 1, 'Precalculating Tm matrix ...' );
Tm_precalculated = precalculate_Tm( sequence );
fprintf( 1, ' done.\n' );
toc

allow_forward_endpoint = ones(1,N_BP);
allow_reverse_endpoint = ones(1,N_BP);
if ( ~exist('check_mispriming') | check_mispriming > 0 )
  LENGTH_CHECK = 4; % could be an input option
  [ allow_forward_endpoint, allow_reverse_endpoint ] = ...
      check_for_mispriming( sequence, LENGTH_CHECK );
end


tic
% Fill in matrix starting from main-diagonal and moving up.
for L = 1:N_BP % Length of window

  fprintf( 1, 'Calculating score for window lengths %d\n',L);

  for i = 1:(N_BP - L + 1) % start of window
    
    j = i + L - 1; %end of window

    %fprintf( 1, 'Calculating score for [ %3d, %3d ]\n',i,j);

    % Can we synthesize this DNA from a reverse and forward primer?
    [ primer1, primer2, cost_of_design, segment1_end, segment2_begin ] = ...
	design_two_primers( sequence(i:j), min_Tm, MAX_LENGTH , ...
			    how_many_variants(i:j), ...
			    allow_forward_endpoint(i:j),...
			    allow_reverse_endpoint(i:j),...
			    Tm_precalculated([i:j],[i:j]) );   
 
    % Shift to whole-sequence numbering
    segment1_end = segment1_end + i - 1;
    segment2_begin = segment2_begin + i - 1;
    
    % Alternatively try to synthesize this DNA from two smaller
    % pieces, whose costs have been computed previously.
    count = 0;

    % This needs to be sped up....
    %[count, costs, split_boundaries1 ] = ...
    %	get_split_costs_slow( Tm_precalculated, DP_scores,i,j, min_Tm); 

    [count, min_cost, split_bounds ] = ...
    	get_split_costs( Tm_precalculated, DP_scores,i,j, min_Tm,...
				      allow_forward_endpoint,...
				      allow_reverse_endpoint); 

    %[count1, count, max(max(Tm_precalculated([i:j],[i:j])))]
    
    if ( ( count == 0 ) & cost_of_design < 0  )
       %Failure!
       DP_scores( i, j ) = -1; % Marks that we have a problem
       choice( i,j )  = 0;
       split_choice{ i, j } = []; 
     else 
       if ( (cost_of_design > 0 ) & ( count < 1 | min_cost > cost_of_design ))
	 % Synthesize this piece!
	 DP_scores( i, j ) = cost_of_design;
	 choice( i, j ) = 1;
	 split_choice{ i, j} = [ segment1_end, segment2_begin ];
       else
	 % Best decision is to split, yo.
	 DP_scores( i, j ) = min_cost;
	 choice( i, j ) = 2;
	 split_choice{ i, j } = split_bounds;
       end
     end

   end
 end
 toc

 primer_sequences = {};

 % Backtrack to figure out what the best solution was.
 if ( DP_scores( 1, N_BP )  < 0 )
   fprintf( 1, 'Crap, could not find a solution!\n' );
   return;
 end

% return;

clf;
image( DP_scores );
colormap(jet)
hold on
plot( [1 N_BP],[1 N_BP],'color',[0.5 0.5 0.5],'linewidth',2);

% Record begin, end, and direction (-1/1)
primers = []; 
primers = track_choice( primers, choice, split_choice, 1, N_BP );
hold off;
set(gca,'xtick',1:N_BP,'xticklabel',sequence',...
	'ytick',1:N_BP,'yticklabel',sequence',...
	'tickdir','out');

primer_sequences = output_primers( primers, sequence );

 [ allow_forward_endpoint, allow_reverse_endpoint ] = ...
    check_for_mispriming( sequence, min_Tm );

fprintf( 'The designed primers are:\n' )
if exist( 'tag' )
  for i = 1:length( primer_sequences )
    fprintf( '%s-%d\t%s\n', tag,i,primer_sequences{i} );
  end
else
  for i = 1:length( primer_sequences )
    fprintf( '%s\n', primer_sequences{i} );
  end
end

return


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [count, costs, split_boundaries ] = get_split_costs_slow( Tm_precalculated, DP_scores,i,j, ...
							   min_Tm); 

count = 0;
costs = [];

MAX_OVERLAP = 1000000;

split_boundaries = [];
for m = (i+1):j % end of first segment
  if ( DP_scores( i, m ) > 0)
    for n = i:(m-1) % beginning of second segment
      if ( DP_scores( n, j ) > 0 & ( m - n ) < MAX_OVERLAP)
	%Tm = calc_Tm( sequence(n:m) ); % This should be precalculated for speed
	Tm = Tm_precalculated(n,m); 
	if ( Tm > min_Tm ) 
	       cost_of_split = DP_scores( i, m) + DP_scores( n, j );
	       
	       count = count+1;
	       costs( count ) = cost_of_split;	       
	       split_boundaries(:, count ) = [m,n];
	       
	end
      end
    end
  end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [count, min_cost, split_bounds ] = ...
    get_split_costs( Tm_precalculated, DP_scores,i,j, min_Tm,...
				  allow_forward_endpoint, allow_reverse_endpoint); 

count = 0;
min_cost = 0;
split_bounds = [0 0];

Tm_chunk = Tm_precalculated( [(i+1):j], [(i+1):j] );

[grid_m,grid_n] = meshgrid( DP_scores( i,(i+1):j), DP_scores((i+1):j,j) );

mask =  (Tm_chunk > min_Tm) .* ...
	(grid_m > 0) .* ...
	(grid_n > 0) ;

[grid_m_allow, grid_n_allow] = ...
    meshgrid( allow_reverse_endpoint( (i+1):j ), ...
	      allow_forward_endpoint( (i+1):j ) );

mask = mask .* grid_m_allow .* grid_n_allow;

count = sum( sum( mask ));

if ( count>0 )

  cost_matrix = (grid_m + grid_n);
  cost_matrix = cost_matrix + max(max(cost_matrix))*(1-mask);

  [min_cost, index1 ] = min( min( cost_matrix ));
  [min_cost, index2 ] = min( cost_matrix(:,index1) );
  split_bounds = [index1, index2] + i;

end  

