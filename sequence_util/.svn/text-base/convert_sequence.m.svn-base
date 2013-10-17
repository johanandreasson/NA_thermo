function numerical_sequence = convert_sequence( sequence );

numerical_sequence = double( sequence );
% Conversion of A,C,G,T to 1,2,3,4 ... following is much faster
% than switch/case, though harder to understand!
numerical_sequence( numerical_sequence == 65 ) = 1; % A
numerical_sequence( numerical_sequence == 67 ) = 2; % C
numerical_sequence( numerical_sequence == 71 ) = 3; % G
numerical_sequence( numerical_sequence == 84 ) = 4; % T
numerical_sequence( numerical_sequence == 85 ) = 4; % U
return;

