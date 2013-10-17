function [primers,Tm_precalculated] = design_mismatch_primers( sequence, mutpos, ...
						  mutseq, min_Tm, ...
						  allow_forward_endpoint, allow_reverse_endpoint );

primers = {};

WINDOWSIZE = 40;
minpos = max( mutpos - WINDOWSIZE, 1);
maxpos = min( mutpos + WINDOWSIZE, length( sequence ) );

mutpos_shift = mutpos - minpos + 1;


% Look in a limited window.
subsequence = sequence( minpos:maxpos );

subsequence_mutant = subsequence;
subsequence_mutant( mutpos_shift ) = mutseq;

% Forward primer with mutation, on reverse complement of sequence.
subsequence_comp = reverse_complement( subsequence );

required_length_forward = ...
    get_length_of_primer( subsequence_mutant, ...
			  subsequence_comp, ...
			  mutpos_shift, ...
			  min_Tm );


% Reverse primer with mutation, on reverse complement of sequence.
num_pos = length( subsequence_mutant );
subsequence_mutant = reverse_complement( subsequence_mutant );
subsequence_comp   = reverse_complement( subsequence_comp );
mutpos_shift = num_pos - mutpos_shift + 1;

required_length_reverse = ...
    get_length_of_primer( subsequence_mutant, ...
			  subsequence_comp, ...
			  mutpos_shift, ...
			  min_Tm );

required_length_reverse = required_length_reverse(end:-1:1);

%clf;
%plot( required_length_forward );
%hold on;
%plot( required_length_reverse,'r' );
%hold off;

[gridx,gridy] = meshgrid( required_length_reverse, required_length_forward );
cost_matrix = gridx+gridy;

% Need Tm of overlap.
fprintf('Calculating overlap Tm matrix...\n');
Tm_precalculated = precalculate_Tm( subsequence );


% Make a plot?
subplot(1,2,1)
image( cost_matrix )
subplot(1,2,2)
image( Tm_precalculated );

% Look around in total cost, but ensure the overlap exceeds the
% desired Tm.
mask = ( Tm_precalculated > min_Tm );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Figuring out where to avoid misprimes...\n');
N_BP   = length( sequence );

if ~exist( 'allow_forward_endpoint' ) | ~exist('allow_reverse_endpoint')
  allow_forward_endpoint = ones(1,N_BP);
  allow_reverse_endpoint = ones(1,N_BP);
  [ allow_forward_endpoint, allow_reverse_endpoint ] = ...
      check_for_mispriming( sequence );
end

[grid_m_allow, grid_n_allow] = ...
    meshgrid( allow_reverse_endpoint( minpos:maxpos ), ...
	      allow_forward_endpoint( minpos:maxpos ) );

mask = mask .* grid_m_allow .* grid_n_allow;


cost_matrix = cost_matrix + max(max(cost_matrix))*(1-mask);

[min_cost, index1 ] = min( min( cost_matrix ));
[min_cost, index2 ] = min( cost_matrix(:,index1) );

subsequence_mutant = reverse_complement( subsequence_mutant );
reverse_primer = subsequence_mutant( ...
    index1 + 1 + [-required_length_reverse(index1):-1]  );
reverse_primer = reverse_complement( reverse_primer );

forward_primer = subsequence_mutant( (index2 - 1) + [1:required_length_forward(index2)]);

primers = {forward_primer, reverse_primer};
char( primers );

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function length_cost = get_length_of_primer( subsequence_mutant, subsequence_comp, ...
					    mutpos_shift, min_Tm);

numpos = length( subsequence_mutant );

length_cost = numpos * ones(1,numpos);

% To allow primer to actually prime.
MIN_MATCH = 5; 

%Typical Phusion conditions
DNA_concentration = 0.2e-6;
monovalent_concentration = 0.1;
divalent_concentration = 0.0015;

[delH_NN, delS_NN, delG_NN, ...
 delH_AT_closing_penalty, delS_AT_closing_penalty, ...
 delG_AT_closing_penalty,...
 delH_mismatch, delS_mismatch, delG_mismatch, ...
 delH_init, delS_init ] = get_NN_parameters();

delS_DNAconc = 1.987 * log( DNA_concentration/2 ); 
delS_init =  delS_init + delS_DNAconc;



numerical_sequence  = convert_sequence( subsequence_mutant );
numerical_sequence_comp  = convert_sequence( subsequence_comp(end:-1:1) );


% potential starting locations
for i = 1 : mutpos_shift

  N_GC = (numerical_sequence(i) == 2 || numerical_sequence(i) == 3  ) ;
  delH_sum = delH_init + delH_AT_closing_penalty( numerical_sequence(i));
  delS_sum = delS_init + delS_AT_closing_penalty( numerical_sequence(i));
  L = 1;
  
  % potential ending locations of forward primer
  for j = (i+1): numpos
    if ( ~is_RC( numerical_sequence(j), numerical_sequence_comp(j) ) )
	delH_sum = delH_sum + delH_mismatch( numerical_sequence(j), ...
					     numerical_sequence_comp(j), ...
					     numerical_sequence(j-1));
	delS_sum = delS_sum + delS_mismatch( numerical_sequence(j), ...
					     numerical_sequence_comp(j), ...
					     numerical_sequence(j-1));
    elseif  ( ~is_RC( numerical_sequence(j-1), numerical_sequence_comp(j-1) ) )
      delH_sum = delH_sum + delH_mismatch( numerical_sequence_comp(j-1), ...
					   numerical_sequence(j-1), ...
					   numerical_sequence_comp(j));
      delS_sum = delS_sum + delS_mismatch( numerical_sequence_comp(j-1), ...
					   numerical_sequence(j-1), ...
					   numerical_sequence_comp(j));
    elseif ( is_RC( numerical_sequence(j), numerical_sequence_comp(j) )...
	     & ...
	     is_RC( numerical_sequence(j), numerical_sequence_comp(j) ))
      delH_sum = delH_sum + delH_NN( numerical_sequence(j-1), numerical_sequence(j));
      delS_sum = delS_sum + delS_NN( numerical_sequence(j-1), numerical_sequence(j));
    else
      fprintf( 'PROBLEM!!!!!!!!!!\n\n' );
      return;
    end

    L = L+1;
    N_GC = N_GC + (numerical_sequence(i) == 2 || numerical_sequence(i) == 3  ) ;
    
    Tm = 1000 * (delH_sum + delH_AT_closing_penalty( numerical_sequence(j))) ./ ...
	 (delS_sum + delS_AT_closing_penalty( numerical_sequence(j)));
    
    Tm = ionic_strength_correction( Tm, monovalent_concentration, ...
				    divalent_concentration, (N_GC/L), L );

    Tm = Tm - 273.15;
    if ( j > mutpos_shift+MIN_MATCH & Tm > min_Tm )
      break;
    end
  end
  
  length_cost(i) = ( j - i + 1);
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = is_RC( c1, c2 )

x = ( (c1+c2) == 5 ) ;

return
