function [ sequences_to_order, construct_names ] = single_mutant_library( primers, sequence, offset, which_muts, ...
    which_libraries, construct_name, shift );

primer_pos = order_primers_along_template( primers, sequence);

num_primers = length( primers );

NUM_LIBRARIES = 3;

num_constructs = ( 1+ length( which_muts ) ) * NUM_LIBRARIES;
num_plates = floor( (num_constructs-1) / 96) + 1;

if ~exist( 'shift' )
    shift = 0;
end

for p = 1:num_primers
    sequences_to_order{ p } = {};
    for k = 1:num_plates
        sequences_to_order{ p }{ k }  = {};
    end
end

WELLS = 'ABCDEFGH';

construct_names = {};

fprintf( 'Filling out sequences...\n');
tic
for p = 1:num_primers
    
    for l_pos = 1:length( which_libraries)
        
        % lib should be a number: 1, 2,or 3 for the different possible mutations.
        lib = which_libraries( l_pos );
        
        for m_pos = 0:length(which_muts)
            
            % which construct is this?
            n = (1 + length( which_muts )) * (lib-1)  + ( m_pos ) + 1;
            
            
            plate_num = floor( (n-1) / 96 ) + 1;
            plate_pos = mod( n-1, 96 ) + 1;
            
            plate_col = floor( (plate_pos-1) / 8 ) + 1;
            plate_row = mod( plate_pos-1, 8) + 1;
            
            well_tag = [WELLS(plate_row),num2str(plate_col)];
            
            
            % m is actual position along sequence.
            if (m_pos == 0 ) % wild type
                m = 0;
            else
                m = offset + which_muts( m_pos );
            end
            
            if ( ( m >= primer_pos( 1, p ) & ...
                    m <= primer_pos( 2, p ) ) | ...
                    m == 0 )
                count = length( sequences_to_order{ p }{ plate_num } );
                count = count + 1;
                
                wt_primer = primers{p};
                mut_primer = wt_primer;
                if (m==0)
                    well_name = 'WT';
                else
                    if (primer_pos(3,p) == -1 ) % reverse primer
                        wt_primer = reverse_complement( wt_primer );
                        mut_primer = reverse_complement( mut_primer );
                    end
                    
                    m_shift = ( m - primer_pos( 1,p ) + 1 );
                    mut_primer( m_shift ) = get_mutation( wt_primer( m_shift ), lib );
                    
                    % Name, e.g., "C75A".
                    well_name = ['Lib',num2str(lib),'-',wt_primer(m_shift),...
                        num2str(which_muts(m_pos)+shift),...
                        mut_primer(m_shift)];
                    
                    if (primer_pos(3,p) == -1 ) % reverse primer
                        wt_primer = reverse_complement( wt_primer );
                        mut_primer = reverse_complement( mut_primer );
                    end
                    
                end
                
                construct_names{n} = well_name;
                sequences_to_order{ p }{ plate_num }{ count } =  { well_tag,  well_name, mut_primer };
                
            end
            
        end
        
    end
    
end
toc


output_sequences_to_order_excel_file( sequences_to_order, primers, construct_name );

output_construct_names( construct_names, construct_name );

return

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


function output_construct_names( construct_names, construct_name );
keys_file = [construct_name,'_keys.txt'];
fprintf( 'Creating Keys file: %s\n',keys_file);
fid = fopen( keys_file,'w');

for i = 1:length( construct_names )
    fprintf( fid, '%s\n', construct_names{i} );
end

fclose( fid );