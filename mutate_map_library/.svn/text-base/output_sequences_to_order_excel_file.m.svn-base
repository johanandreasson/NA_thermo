function output_sequences_to_order_excel_file( sequences_to_order, primers, construct_name )

num_primers =  length( primers );

plate_files = {};
count = 0;
num_plates = length( sequences_to_order{1});

for k = 1:num_plates
    for p = 1:num_primers
        primer_sequences = sequences_to_order{p}{k};
        num_primers_on_plate = length( primer_sequences );
        if (num_primers_on_plate > 0 )
            fprintf( '\n\n%s-Plate%d-Primer%d\n', construct_name, k, p );
            
            count = count+1;
            plate_files{count} = [construct_name,'-Plate',num2str(k),'-Primer',num2str(p),'.xls'];
            fid = fopen( plate_files{count},'w');
            fprintf( fid, '%s\t%s\t%s\t%s\n', 'WellPosition','Name','Sequence','Notes');
            
            for i = 1:num_primers_on_plate
                well_tag = primer_sequences{i}{1};
                well_name = primer_sequences{i}{2};
                mut_primer = primer_sequences{i}{3};
                fprintf( '%s;%s;%s;\n', well_tag, well_name, mut_primer);
                fprintf( fid, '%s\t%s\t%s\t\n', well_tag, well_name, mut_primer);
            end
            fclose( fid );
        end
    end
end

for k = 1:count
    fprintf( 'Created plate file: %s\n', plate_files{k} );
end
