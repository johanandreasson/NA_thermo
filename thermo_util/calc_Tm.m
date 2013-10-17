function Tm = calc_Tm( sequence, DNA_concentration, ...
		       monovalent_concentration, divalent_concentration);

if ~exist( 'DNA_concentration' )
  DNA_concentration = 1e-5;
end

if ~exist( 'monovalent_concentration' )
  monovalent_concentration = 1.0;
end

if ~exist( 'divalent_concentration' )
  divalent_concentration = 0.00;
end
% Easier to keep track of integers in Matlab
% A,C,G,T --> 1,2,3,4.
numerical_sequence  = convert_sequence( sequence );
[delH_NN, delS_NN, delG_NN] = get_NN_parameters();

delH_AT_closing_penalty = [ 2.2 0 0 2.2];
delS_AT_closing_penalty = [ 6.9 0 0 6.9];

%T = 273.15 + 25;
%delG = delH - (T * delS)/1000;

delH_sum = 0.2;
delS_init = -5.7;

DNA_concentration_ref = 1e-6;
delS_DNAconc = 1.987 * log( DNA_concentration/2 ); 
delS_sum =  delS_init + delS_DNAconc;

N_BP = length( numerical_sequence );
for k = 1:(N_BP-1)
  delH_sum = delH_sum + delH_NN( numerical_sequence(k), numerical_sequence(k+1));
  delS_sum = delS_sum + delS_NN( numerical_sequence(k), numerical_sequence(k+1));
end

ADD_TERMINAL_PENALTY = 1;
if ADD_TERMINAL_PENALTY
  delH_sum = delH_sum + delH_AT_closing_penalty( numerical_sequence(1));
  delH_sum = delH_sum + delH_AT_closing_penalty( numerical_sequence(end));
  
  delS_sum = delS_sum + delS_AT_closing_penalty( numerical_sequence(1));
  delS_sum = delS_sum + delS_AT_closing_penalty( numerical_sequence(end));
end

Tm = 1000 * (delH_sum / delS_sum);

f_GC = sum( numerical_sequence == 2 | numerical_sequence == 3) / N_BP;

IONIC_STRENGTH_CORRECTION = 1;
if IONIC_STRENGTH_CORRECTION
  Tm = ionic_strength_correction( Tm, monovalent_concentration, ...
				divalent_concentration, f_GC, N_BP );
end
%delH_sum
%delS_sum

delG_NN = delH_NN - (273.15+37)*delS_NN/1000;

delG_sum = delH_sum - (273.15+37) * delS_sum/1000;
% Convert to Celsius
Tm = Tm - 273.15; 

return

