function output_sequences_to_order_96well_diagram( sequences_to_order, primers, construct_name )

num_primers =  length( primers );
num_plates = length( sequences_to_order{1});

for k = 1:num_plates
    figure(); clf;
    set(gcf, 'PaperPositionMode', 'Manual', 'Color', 'White');
    set(gcf, 'Position', [(k-1)*50 0 600 800]);
    set(gcf, 'PaperOrientation', 'Portrait', ...
        'PaperSize', [8.5 11], 'PaperPosition', [0 1 8.5 10.5]);
    
    for p = 1:num_primers
        
        primer_WT = [];
        primer_sequences = sequences_to_order{p}{k};
        num_primers_on_plate = length( primer_sequences );
        if (num_primers_on_plate > 0 )
            
            subplot( (num_primers/2), 2, p );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make 96 well grid.
            cla;
            for i = 1:12;
                for j = 1:8;
                    rectangle( 'Position',[i-0.2 j-0.2 0.9 0.9],'Curvature',[1 1],'Linewidth',1);
                    hold on;
                end;
            end;
            for i = 1:12;
                h = text( i,0.2,num2str(i) );
                set(h,'fontsize',10);
            end;
            letters = 'ABCDEFGH';
            for j = 1:8;
                h = text( 0.2,j+0.2,letters(j) );
                set(h,'fontsize',10);
            end;
            set(gca,'xtick',[],'ytick',[],'ydir','reverse');
            axis equal;
            axis([-0.5 14 -0.5 9]);
            axis off;
            
            h = title( ['\bf',construct_name,'\rm-Plate\bf{\color{magenta}',num2str(k),'}\rm-Primer\bf{\color{blue}', num2str(p),'}'] );
            set(h,'fontsize',12);
            
            for i = 1:num_primers_on_plate
                well_tag = primer_sequences{i}{1};
                well_name = primer_sequences{i}{2};
                mut_primer = primer_sequences{i}{3};
                
                [row, col ] = get_row_col(well_tag);
                % T47, color WT primers in yellow
                if strcmp(mut_primer, primers{p});
                    if row == 1 && col == 1;
                        rectangle( 'Position',[col-0.2 row-0.2 0.9 0.9],'Curvature',[1 1],'Linewidth',1,...
                            'FaceColor','r');
                    else
                        rectangle( 'Position',[col-0.2 row-0.2 0.9 0.9],'Curvature',[1 1],'Linewidth',1,...
                            'FaceColor','y');
                    end;
                    primer_WT = [primer_WT, (col-1)*8+row];
                else
                    rectangle( 'Position',[col-0.2 row-0.2 0.9 0.9],'Curvature',[1 1],'Linewidth',1,...
                        'FaceColor','g');
                end;
                %text( col, row, well_name );
            end;

            primer_gray = primer_WT(diff(primer_WT) == 1)+1;
            for i =1:length(primer_gray)
                well_tag = primer_sequences{primer_gray(i)}{1};
                [row, col ] = get_row_col(well_tag);
                rectangle( 'Position',[col-0.2 row-0.2 0.9 0.9],'Curvature',[1 1],'Linewidth',1,...
                    'FaceColor',[0.8,0.8,0.8]);
            end;
        end;
    end;
    
    print_save_figure(gcf, [construct_name,'-Plate',num2str(k)], '/');
end;

% show legend of colors
fprintf('\n'); fprintf('Legend\n'); fprintf('======\n');
fprintf('RED:    Wild-type primer\n');
fprintf('GREEN:  Mutant primers\n');
fprintf('YELLOW: Wild-type primer inserts\n');
fprintf('GRAY:   Wilt-type primer deletes\n');
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [row,col] = get_row_col( well_tag )


row = strfind( 'ABCDEFGH', well_tag(1));
col = str2num( well_tag(2:end) );

return