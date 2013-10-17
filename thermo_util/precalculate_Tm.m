function [Tm_precalculated,Tm_precalculated_old] = precalculate_Tm( sequence );
% This could be sped up further by doing some matlab matrixology.

DNA_concentration = 0.2e-6; % Typical for Phusion reaction
%Proprietary -- this is the value needed to
%reproduce Finnzyme Tm calculator results
monovalent_concentration = 0.1; 

% Nominal Mg2+ in phusion literation.
% Do we need to correct for [dNTP] sequestering of Mg2+?
divalent_concentration = 0.0015; 

%This could be sped up significantly, since many of the sums of
%delH, delG are shared between calculations.
N_BP = length(sequence );

 if ~exist( 'DNA_concentration' )
  DNA_concentration = 0.2e-6;
end

if ~exist( 'monovalent_concentration' )
  monovalent_concentration = 0.1;
end

if ~exist( 'divalent_concentration' )
  divalent_concentration = 0.0015;
end
% Easier to keep track of integers in Matlab
% A,C,G,T --> 1,2,3,4.
numerical_sequence  = convert_sequence( sequence );

[delH_NN, delS_NN, delG_NN, ...
 delH_AT_closing_penalty, delS_AT_closing_penalty, ...
 delG_AT_closing_penalty,...
 delH_mismatch, delS_mismatch, delG_mismatch, ...
 delH_init, delS_init ] = get_NN_parameters();


delS_DNAconc = 1.987 * log( DNA_concentration/2 ); 
delS_init =  delS_init + delS_DNAconc;

N_BP = length( numerical_sequence );

% Start at diagonal and work our way to the top right.
delH_matrix = delH_init * ones( N_BP, N_BP );
delS_matrix = delS_init * ones( N_BP, N_BP );
f_GC = zeros( N_BP, N_BP );
len_BP = ones( N_BP, N_BP );

for i = 1:N_BP
  if (numerical_sequence(i) == 2 || numerical_sequence(i) == 3 ) 
    f_GC(i,i) = 1;
  end
end

fprintf(1, 'Filling delH, delS matrix...\n' );
for i = 1:N_BP
  for j = (i+1):N_BP
    delH_matrix(i,j) = delH_matrix(i,j-1) + delH_NN( numerical_sequence(j-1), numerical_sequence(j));
    delS_matrix(i,j) = delS_matrix(i,j-1) + delS_NN( numerical_sequence(j-1), numerical_sequence(j));      
    len_BP(i,j) = len_BP(i,j-1) + 1;

    f_GC(i,j) = f_GC(i,j-1);
    if (numerical_sequence(j) == 2 || numerical_sequence(j) == 3 ) 
      f_GC(i,j)   = f_GC(i,j) + 1;
    end
  
  end
end


ADD_TERMINAL_PENALTY = 1;
if ADD_TERMINAL_PENALTY
  fprintf(1, 'Terminal penalties... \n');
  for i = 1:N_BP
    for j = (i+1):N_BP
      delH_matrix(i,j) = delH_matrix(i,j) + delH_AT_closing_penalty( numerical_sequence(i));
      delH_matrix(i,j) = delH_matrix(i,j) + delH_AT_closing_penalty( numerical_sequence(j));    
      
      delS_matrix(i,j) = delS_matrix(i,j) + delS_AT_closing_penalty( numerical_sequence(i));
      delS_matrix(i,j) = delS_matrix(i,j) + delS_AT_closing_penalty( numerical_sequence(j));    
    end
  end
end

Tm = 1000 * (delH_matrix ./ delS_matrix);
f_GC = f_GC ./ len_BP;

IONIC_STRENGTH_CORRECTION = 1;
if IONIC_STRENGTH_CORRECTION
  fprintf(1, 'Ionic strength corrections... \n');
  for i = 1:N_BP
    for j = i:N_BP;
      Tm(i,j) = ionic_strength_correction( Tm(i,j), monovalent_concentration, ...
					   divalent_concentration, f_GC(i,j), len_BP(i,j) );
    end
  end
end

% Convert to Celsius
Tm_precalculated = Tm - 273.15; 

Tm_precalculated_old = [];

return;

% Skip this by default -- this is the old, slow calculation.
tic
for i = 1:N_BP
  for j = i:N_BP
    Tm_precalculated_old(i,j) = calc_Tm( sequence(i:j),...
					 DNA_concentration, ...
					 monovalent_concentration, ...
					 divalent_concentration );
   end
end
toc

return
