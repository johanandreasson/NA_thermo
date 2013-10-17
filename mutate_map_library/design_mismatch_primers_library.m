function design_mismatch_primers_library( sequence, start_primer, end_primer, ...
				 offset, which_muts, ...
				 which_libraries, construct_name );

primer_pos = order_primers_along_template( {start_primer, end_primer}, sequence);

NUM_LIBRARIES = 3;

num_constructs = ( 1+ length( which_muts ) ) * NUM_LIBRARIES;
num_plates = floor( (num_constructs-1) / 96) + 1;

WELLS = 'ABCDEFGH';

for i = 1:num_plates;    
  sequences_to_order{ i }= {};
end

fprintf( 'Filling out sequences...\n');
tic

min_Tm = 65;
N_BP = length( sequence );
allow_forward_endpoint = ones(1,N_BP);
allow_reverse_endpoint = ones(1,N_BP);
[ allow_forward_endpoint, allow_reverse_endpoint ] = ...
    check_for_mispriming( sequence );


which_plates = [ 1 ]; %For pilot, just check plate 1.

for l_pos = 1:length( which_libraries)
  
  % lib should be a number: 1, 2,or 3 for the different possible mutations.
  lib = which_libraries( l_pos );
  
  for m_pos = 0:length(which_muts)
    
    fprintf( 'Figuring out mismatch position %d of %d, in library %d of %d...\n ',m_pos, ...
	     length( which_muts ), l_pos, length( which_libraries) );
    % which construct is this?
    n = (1 + length( which_muts )) * (lib-1)  + ( m_pos ) + 1;
    
    
    plate_num = floor( (n-1) / 96 ) + 1;
    plate_pos = mod( n-1, 96 ) + 1;
    
    plate_col = floor( (plate_pos-1) / 8 ) + 1;      
    plate_row = mod( plate_pos-1, 8) + 1;
    
    well_tag = [WELLS(plate_row),num2str(plate_col)];
    
    %%%%%%%%%%%%%%%%%
    % HACK!
    if ( ~isempty(find(plate_num == which_plates)) )
    
      % m is actual position along sequence.
      if (m_pos == 0 ) % wild type
	m = 0; 
      else 
	m = offset + which_muts( m_pos );
      end
      
      primer_F = start_primer;
      primer_R = end_primer;
      
      if (m==0) 
	well_name = 'WT';
      else
	mutseq = get_mutation( sequence( m ), lib );
	
	primers_mismatch =  design_mismatch_primers( sequence, m, ...
						     mutseq, min_Tm,...
						     allow_forward_endpoint, ...
						     allow_reverse_endpoint ...
						     );
	primer_F = primers_mismatch{1};
	primer_R = primers_mismatch{2};
	
	% Name, e.g., "C75A".
	well_name = ['Lib',num2str(lib),'-SDM-',sequence(m),...
		     num2str( which_muts( m_pos) ),...
		     mutseq];
	
    end
    
    count = length( sequences_to_order{ plate_num } ) + 1;
    sequences_to_order{ plate_num }{count} =  { well_tag, well_name, primer_F, primer_R };
    
    
  end    
  
  end
  
end

toc

plate_files = {};
count = 0;
for k = 1:1%num_plates

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Forward primers.
  primer_sequences = sequences_to_order{k};
  fprintf( '\n\n%s-SDM-Plate%d-F\n', construct_name, k );

  count = count + 1;
  plate_files{count} = [construct_name,'-SDM-Plate',num2str(k),'-F.xls'];
  fid = fopen( plate_files{count},'w');
  fprintf( fid, '%s\t%s\t%s\t%s\n', 'WellPosition','Name','Sequence','Notes');
  
  for i = 1:length( primer_sequences );
    well_tag = primer_sequences{i}{1};
    well_name = primer_sequences{i}{2};
    mut_primer = primer_sequences{i}{3};
    fprintf( '%s;%s-F;%s;\n', well_tag, well_name, mut_primer);	
    fprintf( fid, '%s\t%s-F\t%s\t\n', well_tag, well_name, mut_primer);
  end
  fclose( fid );

  %%%%%%%%%%%%%%%%%%%%%%%
  % Reverse primers.
  primer_sequences = sequences_to_order{k};
  fprintf( '\n\n%s-SDM-Plate%d-R\n', construct_name, k );

  count = count + 1;
  plate_files{count} = [construct_name,'-SDM-Plate',num2str(k),'-R.xls'];
  fid = fopen( plate_files{count},'w');
  fprintf( fid, '%s\t%s\t%s\t%s\n', 'WellPosition','Name','Sequence','Notes');
  
  for i = 1:length( primer_sequences );
    well_tag = primer_sequences{i}{1};
    well_name = primer_sequences{i}{2};
    mut_primer = primer_sequences{i}{4};
    fprintf( '%s;%s-R;%s;\n', well_tag, well_name, mut_primer);	
    fprintf( fid, '%s\t%s-R\t%s\t\n', well_tag, well_name, mut_primer);
  end
  fclose( fid );

end

fprintf('\n\n');
for k = 1:length( plate_files )
  fprintf( 'Created plate file: %s\n', plate_files{k} );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_mut = get_mutation( c, lib );
switch lib
 case 1
  switch c
   case 'A'
    c_mut = 'T';
   case 'T'
    c_mut = 'A';
   case 'C'
    c_mut = 'G';
   case 'G'
    c_mut = 'C';
  end
 case 2
  switch c
   case 'A'
    c_mut = 'C';
   case 'T'
    c_mut = 'C';
   case 'C'
    c_mut = 'A';
   case 'G'
    c_mut = 'A';
  end
 case 3
  switch c
   case 'A'
    c_mut = 'G';
   case 'T'
    c_mut = 'G';
   case 'C'
    c_mut = 'T';
   case 'G'
    c_mut = 'T';
  end
end
return
