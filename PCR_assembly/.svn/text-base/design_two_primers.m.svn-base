function [primer1, primer2, cost_of_design, segment1_end, segment2_begin ] = ...
    design_two_primers( sequence, min_Tm, MAX_LENGTH, ...
			how_many_variants, ...
			allow_forward_endpoint, ...
			allow_reverse_endpoint, ...
			Tm_precalculated );

DNA_concentration = 0.2e-6; % Typical for Phusion reaction
%Proprietary -- this is the value needed to
%reproduce Finnzyme Tm calculator results
monovalent_concentration = 0.1; 
% Nominal Mg2+ in phusion literation.
% Do we need to correct for [dNTP] sequestering of Mg2+?
divalent_concentration = 0.0015; 

N_BP = length( sequence );

primer1 = '';
primer2 = '';
cost_of_design = -1;
segment1_end = 0;
segment2_begin = 0;


if ~exist( 'min_Tm' )
  min_Tm = 65;
end
if ( ~exist('how_many_variants') )
  how_many_variants = zeros( 1, N_BP );
end
if ( ~exist('MAX_LENGTH') )
  MAX_LENGTH = 60;
end
if ( ~exist('MIN_LENGTH') ) % IDT plates require >= 15 nts.
  MIN_LENGTH = 15;
end
if ( ~exist('allow_forward_endpoint') )
  allow_forward_endpoint = ones( 1, N_BP );
end
if ( ~exist('allow_reverse_endpoint') )
  allow_reverse_endpoint = ones( 1, N_BP );
end



OUTPUT_DESIGN = 0;
if ~exist( 'Tm_precalculated' )
  Tm_precalculated = zeros( N_BP, N_BP );
  for i = 2:MAX_LENGTH
    if ( i <= N_BP )
      for j = (N_BP-MAX_LENGTH):N_BP
	if ( j > 0  & j < i)
	  Tm_precalculated(j,i) = calc_Tm( sequence(j:i), ...
					   DNA_concentration, ...
					   monovalent_concentration, ...
					   divalent_concentration );	
	end
      end
    end
  end
  OUTPUT_DESIGN = 1;
end


count = 0;

% Apply temperature cutoff..
mask = ( Tm_precalculated > min_Tm );

% Apply length cutoffs.
[xgrid,ygrid] = meshgrid(1:N_BP,1:N_BP);
L1 = xgrid;
L2 = (N_BP-ygrid+1);
mask = mask .* (L1 <= MAX_LENGTH) .* (L2 <= MAX_LENGTH );
mask = mask .* (L1 >= MIN_LENGTH) .* (L2 >= MIN_LENGTH );
mask = mask .* (xgrid > ygrid );

[allow_gridx,allow_gridy] = meshgrid( allow_forward_endpoint, allow_reverse_endpoint );
mask = mask .* allow_gridx;
mask = mask .* allow_gridy;

count = sum(sum(mask));

if (count > 0 )
  
  % Following is just 1 if we're just making one sequence.
  % However if we're trying to make every single mutant, 
  %  this effectively says that the cost is the sum of the
  %  squares of each length.
  N_variants_forward = 1 + cumsum( how_many_variants );
  N_variants_reverse = 1 + cumsum( how_many_variants(end:-1:1) );
  N_variants_reverse = N_variants_reverse( end:-1:1 );
  [ngrid_forward, ngrid_reverse ] = meshgrid( N_variants_forward, N_variants_reverse );

  cost = L1.*ngrid_forward +  + L2.*ngrid_reverse;

  % Mask out lengths that correspond to good Tm's:
  A_BIG_PENALTY = 9999999999;
  cost = cost + A_BIG_PENALTY * (1-mask);

  min_cost = min(min(cost));
  cost_of_design = min_cost;

  % If there are several solutions, pick the one with the highest Tm.
  Tms = ( cost == min_cost) .* Tm_precalculated .* mask;
  [Tm, index1 ] = max( max( Tms ) );
  [dummy, index2] = max( Tms( :,index1) );
  
  segment1_end = index1;
  segment2_begin = index2;
  
  if OUTPUT_DESIGN
    output_design( segment1_end, segment2_begin,sequence);
    primer1 = sequence(1:segment1_end);
    primer2 = reverse_complement( sequence(segment2_begin:end) );
  end
  
  
end


%cost_of_design

return;


function output_design( i, j, sequence )

N_BP = length( sequence );
for k = 1:N_BP
  if (k <= i & k >= j) 
    line1(k)=sequence(k);
    line2(k)='|';
    line3(k)=reverse_complement(sequence(k));
  else
    if ( k > i )
	    line1(k)=' ';
	    line2(k)=' ';
	    line3(k)=reverse_complement(sequence(k));
    else
      line1(k)=sequence(k);
      line2(k)=' ';
      line3(k)=' ';
    end
  end
end
fprintf(1,'%s\n%s\n%s\n',line1,line2,line3);

return;

