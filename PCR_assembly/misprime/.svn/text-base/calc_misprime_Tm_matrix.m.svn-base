function delG_matrix = calc_misprime_Tm_matrix( sequence1, sequence2, PRIMER_SIZE,  ...
				  DNA_concentration, ...
				  monovalent_concentration, ...
				  divalent_concentration, T...
				  );

if ~exist('DNA_concentration')
  DNA_concentration = 1e-5;
end
if ~exist( 'monovalent_concentration' )
  monovalent_concentration = 0.1;
end
if ~exist( 'divalent_concentration' )
  divalent_concentration = 0.0015;
end
if ~exist( 'T' )
  T = 65;
end

if ~exist( 'PRIMER_SIZE' )
  PRIMER_SIZE = 20;
end

num1 = length( sequence1 );
num2 = length( sequence2 );

delG_matrix = 99.9 * ones( num1, num2 );

MIN_ANNEAL_SEGMENT = 3;
MAX_ANNEAL_SEGMENT = 40;

need_calcs = [];
for i = 20:num1
  for j = MIN_ANNEAL_SEGMENT:num2
    if (sequence1((i-MIN_ANNEAL_SEGMENT+1):i) == ...
	reverse_complement( sequence2(  (j-MIN_ANNEAL_SEGMENT+1):j ) ) ...
	)
      need_calcs = [ need_calcs; i j];
    end
  end
end

num_calcs = size( need_calcs, 1 );

if matlabpool( 'size' ) == 0
  matlabpool( 2 );
end

calc_range = [1:num_calcs];

% Parallelized...
delG_compile = [];
parfor k = calc_range
  tic
  i = need_calcs(k,1);
  j = need_calcs(k,2);

  fprintf( 'Doing calculation %d of %d...\n',k,num_calcs);
  delG_compile( k ) = 0.0;
  %delG_matrix(i,j) = 0.0;
  
  primer = sequence1( (i-PRIMER_SIZE+1):i);
  
  landing_pad_start = j - MIN_ANNEAL_SEGMENT+1;
  landing_pad_end = ...
      min((landing_pad_start + MAX_ANNEAL_SEGMENT - 1), num2);
  landing_pad = sequence2( landing_pad_start: landing_pad_end);
  
  %primer
  %landing_pad
  
  delG_compile(k) = calc_misprime_delG( ...
      landing_pad, primer,...
      DNA_concentration, T, ...
      monovalent_concentration, divalent_concentration );	
  toc
end


for k = calc_range
  i = need_calcs(k,1)
  j = need_calcs(k,2)
  delG_matrix(i,j) = delG_compile(k);
end
