function calc_misprime_Tm( sequence1, sequence2, DNA_concentration, ...
			   monovalent_concentration, divalent_concentration, ...
			   T );
% Assumes first residue of sequence1 base pairs with last residue
% of sequence 2. This assumption could/should be relaxed, although 
% I'm mostly interested in the case in which several base pairs are formed.
%  It should actually be possible to calculate all possible
%  mishybridizing interactions in a single shot, with much savings
%  in computational cost, but I haven't figured this out yet.

num1 = length( sequence1);
num2 = length( sequence2);

% Need two matrices for the dynamic programming:
% Q_1 is the best delG  
% 
Q_1 = zeros( num1, num2 );
Q_2 = zeros( num1, num2 );




