sequence_P4P6 = ...
    'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC'
sequence_P4P6_RC = reverse_complement( sequence_P4P6);
numres = length( sequence_P4P6); % 222

% Takes a loooong time.
% Assumed [DNA] = 10 uM. Temperature = 65 C;
delG_matrix = calc_misprime_Tm_matrix( sequence_P4P6, ...
				       sequence_P4P6_RC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
image(10*min(delG_matrix,10)) % cutoff of 10 kcal/mol
colormap( 1 - gray(100))
h = title('P4P6 -- \Delta G for mispriming 20nt chunks with reverse complement. [DNA] = 10 uM. T = 65 C')
set(h,'fontsize',10);

set(gca,'ytick',1:numres, 'yticklabel',sequence_P4P6','fontsize',4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Block out super-close mispriming events.
figure(1)
delG_test = delG_matrix;
for i= 1:numres; 
  for j = -5:5;
    if ( numres-i+j >= 1 & numres-i+j <= numres )
      delG_test(i, numres-i+j) = 100;
    end
  end
end
    
delG_misprime = min( delG_test' );
plot( delG_misprime )
ylim([-10 10])
set(gca,'xtick',1:numres, 'xticklabel',sequence_P4P6','fontsize',4 )