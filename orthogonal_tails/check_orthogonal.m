function [delG, Tm, free_end_F, free_end_R, bps ] = check_orthogonal( sequence_ref , sequences );

num_sequences = length( sequences );

T = 37;
DNA_concentration = 1e-6;
monovalent_concentration = 0.1;
divalent_concentration = 0.0015;

bps = {};

which_sequences = 1:num_sequences;

if matlabpool( 'size' ) == 0
  matlabpool( 2 );
end

parfor i = which_sequences
  fprintf('Doing sequence %d out of %d\n', i, num_sequences );
  [delG_out, Tm_out, bps_out] = calc_misprime_delG( reverse_complement( sequence_ref ), ...
					   sequences{i}, DNA_concentration, ...
					   T,monovalent_concentration,...
					   divalent_concentration );
  fprintf('delG( %3.1f C ) = %4.1f kcal/mol  ; Tm = %3.1f \n\n', T, delG_out, Tm_out );
					   
  if ( ~isempty( bps_out ) & ~isempty( find( bps_out(:,1) == length( sequence_ref) ) ) )
    free_end_F( i ) = 0;
  else
    free_end_F( i ) = 1;
  end

  if ( ~isempty( bps_out ) &  ~isempty( find( bps_out(:,2) == 1 )) )
    free_end_R( i ) = 0;
  else
    free_end_R( i ) = 1;
  end

    
	   
  
  delG(i) = delG_out;
  Tm(i) = Tm_out;
  bps{i} = bps_out;
  
end