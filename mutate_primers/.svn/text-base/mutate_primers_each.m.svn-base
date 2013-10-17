function [sequences_to_order,sequence_names] = mutate_primers_each( primers, muts, WT_mode )

if ~iscell( muts ); [muts, sequence_names] = read_sequences( muts ); end;
if ~exist('WT_mode','var') || isempty(WT_mode); WT_mode = 0; end;

% squeezes out gaps:
sequence = get_primer_seq( muts{1}, 1, length( muts{1}), 1 );
primer_pos = order_primers_along_template( primers, sequence);
seq_map = get_mapping( muts{1}, sequence );

WELLS = 'ABCDEFGH';
output_primers = {};
primer_names = {};

for p = 1:length(primers)
    sequences_to_order{ p }{1} = {};
end
count = zeros(1,length(primers));

for k = 1:length( muts )
    plate_col = floor( (k-1) / 8 ) + 1;
    plate_row = mod( k-1, 8) + 1;
    well_tag = [WELLS(plate_row),num2str(plate_col)];
    seq_mut = muts{k};
    for i = 1:length(primers)
        start_pos = seq_map( primer_pos(1,i) );
        end_pos = seq_map( primer_pos(2,i) );
        dir = primer_pos(3,i);
        primer_seq = get_primer_seq( seq_mut, start_pos, end_pos, dir );
        primer_names = sprintf( '%s-%d',sequence_names{k}, i );
        if ~strcmp(primer_seq, primers{i}) || k == 1 || WT_mode;
            count(i) = count(i) + 1;
            sequences_to_order{ i }{ 1 }{ count(i) } =  { well_tag,  primer_names, primer_seq };
        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_map = get_mapping( seq1, sequence )

count = 0;
for i = 1:length( seq1)
    if ~(seq1( i ) == ' ' || seq1(i) == '-')
        count = count + 1;
        if ~( sequence( count ) == seq1( i ) )
            fprintf(1,['First mutant sequence should match overall sequence!\n']);
        end
        seq_map( count ) = i;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function primer_seq = get_primer_seq( seq_mut, start_pos, end_pos, dir )
primer_seq = '';
for i = start_pos:end_pos
    if ~(seq_mut(i)==' ' || seq_mut(i)=='-')
        primer_seq = [primer_seq, seq_mut(i) ];
    end
end

if (dir == -1)
    primer_seq = reverse_complement( primer_seq );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sequences, sequence_names ] = read_sequences( filename )

fid = fopen( filename );
k=0;
while ~feof(fid)
    k=k+1;
    sequence_name = fgetl(fid);
    %fprintf(1,'%s\n',sequence_name );
    [sequence, remain ] = strtok(sequence_name);
    sequences{k} = sequence;
    
    if length(remain)<1
        sequence_names{k} = ['Sequence ',num2str(k)];
    else
        sequence_names{k} = strtok(remain);
    end
end

fclose(fid);
