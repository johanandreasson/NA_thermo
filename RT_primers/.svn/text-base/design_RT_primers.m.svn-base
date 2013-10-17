function primers = design_RT_primers( sequence, N, primer_length, tag, add_A, bpp );

% base pair exposure
if ~exist( 'bpp' ) | isempty( bpp ); bpp = partition( sequence ); end;
if ~exist( 'primer_length' ); primer_length = 20; end;
if ~exist( 'tag' ); tag = 'RTprim'; end;
if ~exist( 'add_A' ); add_A = 1; end;

% num of potential cross-priming sites
max_match = get_max_match( sequence );

all_primer_pos = design_RT_primers_by_DP( N, primer_length, bpp, max_match );

primers = output_RT_primers( all_primer_pos, sequence, primer_length );

if add_A; primers = add_A_to_primers( primers ); end;
primers = add_FAM_to_primers( primers );

output_with_tags( primers, tag );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function primers = add_A_to_primers( primers );
NUM_A = 20;
A_string = '';
for n = 1:NUM_A; A_string = [A_string, 'A' ]; end;

for j = 1:length( primers )
  primers{j} = [A_string, primers{j} ];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function primers = add_FAM_to_primers( primers );
FAM_string = '/56-FAM/';
for j = 1:length( primers )
  primers{j} = [FAM_string, primers{j} ];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_with_tags( primers, tag )

for j = 1:length( primers )
  fprintf( 1, '%s-%d\t%s\n', tag, j, primers{j} );
end