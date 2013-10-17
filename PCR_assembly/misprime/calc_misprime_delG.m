function [delG, Tm, bps, Q1, Q2, choice1, choice2 ] = ...
    calc_misprime_delG( sequence1_input, sequence2_input, DNA_concentration, ...
			T, monovalent_concentration, divalent_concentration  );
% Hybridization of two DNA primers.

if ~exist('T')
  T = 37;
end
T = 273.15 + T;

[delH_NN, delS_NN, delG_NN, ...
	  delH_AT_closing_penalty, delS_AT_closing_penalty, ...
	  delG_AT_closing_penalty,...
	  delH_mismatch, delS_mismatch, delG_mismatch, ...
	  delH_init, delS_init ...
	 ] = get_NN_parameters( T );

if ~exist('DNA_concentration')
  DNA_concentration = 1e-5;
end
if ~exist( 'monovalent_concentration' )
  monovalent_concentration = 1.0;
end
if ~exist( 'divalent_concentration' )
  divalent_concentration = 0.0;
end

% Initialization of helix.
delS_DNA_conc = 1.987 * log( DNA_concentration/2 ); 
delG_init = delH_init - T * ( delS_init + delS_DNA_conc )/1000; 

% delG_mismatch( X,Y,A/C/G/U) -->  delG_mismatch( 2,4,1 ) corresponds to AC/TT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For convenience, flip sequence2, and 
% use numerical symbols 1/2/3/4 for A/C/G/T.
sequence1 = convert_sequence( sequence1_input );
sequence2 = convert_sequence( sequence2_input( end:-1:1 ) );

num1 = length( sequence1);
num2 = length( sequence2);

% Need two matrices for the dynamic programming:
% Q_2(i,j) is the best delG given that the "last" base pair
%  is at (i,j) and is immediately preceded by one base pair.
% Q_1(i,j) is the best delG given that the "last" base pair does
%  *not* have an immediately preceding base pair (these solutions
%  are not acceptable in the end, but we need to keep track of them
%  during the calculation.)
A_BIG_NUMBER = 99999999;
Q1 = A_BIG_NUMBER * ones( num1, num2 );
Q2 = A_BIG_NUMBER * ones( num1, num2 );

choice1 = cell(num1,num2);
choice2 = zeros( num1, num2 );

delG = 0;
Tm = 0;


%Start at left-hand side of num1, right-hand side of num2.
%
%for L = 3: (num1+num2)
%  for i = 3:L
%    j = L + 3 - i;
for i = 1:num1
  for j = 1:num2

    % This could perhaps be sped up dramatically if we kept track of (i,j) actually permit a single base pair and a double base pair!
    if (i <= num1 & j <= num2 & ... % Better be in the sequence
	RC( sequence1(i), sequence2(j ) )... % Better be complementary.
	)

      % For Q2, can add on to a single base pair or double base pair.
      if ( i > 1 & j > 1 & RC( sequence1(i-1), sequence2(j-1) ) )
	[best_score, choice ] = ...
	    min( [Q1( i-1, j-1 ), Q2( i-1, j-1)] );
	choice2(i,j) = choice; 
	Q2(i,j) = best_score + delG_NN( sequence1(i-1), sequence1(i) );
      end
      
      %For Q1, can add on to anything except a base pair.
      count = 0;
      scores = [];

      % This could be the first base pair...
      count = count + 1;
      index1( count ) = 0;
      index2( count ) = 0;
      scores( count ) = delG_init;
      scores( count ) = scores( count ) + delG_AT_closing_penalty( sequence1(i) );

      % Add in terminal mismatch
      if ( i > 1 & j > 1 ) 
	if (~RC( sequence1(i-1), sequence2(j-1) ) )
	  score( count ) = scores( count ) + delG_mismatch( sequence2(j-1),sequence1(i-1),sequence2(j) ); % Ending base pair.      
	end
      end

      for h = 2:(i-1)
	for l = 2:(j-1)
	  if ~( (h == i-1) & (l == j-1) )
	    count = count + 1;
	    index1( count ) = h;
	    index2( count ) = l;
	    loop_p = loop_penalty( sequence1(h:i),...
				   sequence2(l:j),...
				   T, delG_NN, ...
				   delG_mismatch, ...
				   delG_AT_closing_penalty);
	    scores( count ) = Q2(h,l) + loop_p;
	  end
	end
      end

      [score_min, count_min ] = min( scores );
      choice1{i,j} = [index1(count_min), index2( count_min )] ;
      Q1(i,j) = score_min;

      
    end
  end
end

% Add in any A_T closing penalties. 
%Q2 = Q2 + delG_AT_closing_penalty( sequence1(1) ); % Starting base pair.

for i = 2:num1
  for j = 2:num2
    if ( RC( sequence1(i), sequence2(j ) ) ... % Better be complementary.
	 & RC( sequence1(i-1), sequence2( j- 1 ) ) )
      Q2(i,j) = Q2(i,j) + delG_AT_closing_penalty( sequence1(i) ); % Ending base pair.
      
      % Add in any terminal mismatch or dangling ends.
      % How about dangles?!!
      if ( i < num1 & j < num2 )
	if (~RC( sequence1(i), sequence2(j) ) )
	  Q2(i,j) = Q2(i,j) + delG_mismatch( sequence1(i+1),sequence2(j+1),sequence1(i) ); % Ending base pair.      
	end
      end
    end
  end
end

% Figure out best closing base pair ( base pair doublet, actually).
[ best_score, j ]  = min( min( Q2 ) );
[ best_score, i ]  = min( Q2(:,j) ); 

pos1 = i;
pos2 = j;

if (pos1 > 0 | pos2 > 0 )
  bps = [];
  bps2 = [];
  in_doublet = 2;
  while( pos1 > 0 & pos2 > 0 )
    bps = [ bps; pos1, pos2];
    
    if (in_doublet==2) 
      bps2 = [ bps2; pos1, pos2];
      in_doublet = choice2( pos1, pos2 ); 
      % Did we get here from a previous double-base-pair or from a
      % lone base pair?
      pos1 = pos1 - 1;
      pos2 = pos2 - 1;
    else 
      prev_bp = choice1{ pos1, pos2 };
      in_doublet = 2;
      if length( prev_bp ) < 2
	fprintf(1,'Problem in trackback\n');
	return;
      end
      pos1 = prev_bp( 1 );
      pos2 = prev_bp( 2 );
    end
  end
  %bps = [ bps; pos1, pos2];
  
  output_bp( bps, sequence1_input, sequence2_input );
  
  % Tm estimate for this optimal solution.
  % Currently assuming that loops just contribute to
  % deltaS -- typical assumption.
  %
  % Also, however, assuming that mismatches contribute only
  % to deltaS. NEED TO FIX BY FINDING deltaS VALUE IN LITERATURE!
  %
  %
  delH_tot = delH_init;
  for i = 1:length( bps2 )
    bp_doublet_end = bps2(i);
    delH_tot =  delH_tot + delH_NN( sequence1( bp_doublet_end-1), ...
				    sequence1( bp_doublet_end) );
  end
  delH = delH_tot;
  
  delS = ( delH_tot - best_score)*1000/T;
  
  %delS = delS + delS_init + delS_DNA_conc;
  %delH = delH_tot + delH_init;
  
  Tm = 1000 * (delH / delS ) - 273.15;
  
  delG = delH - (T * delS/1000);
  
  N_BP = 0;
  f_GC = 0.0;
  for i = 1:length( bps )
    N_BP = N_BP + 1;
    if (sequence1( bps(i,1) ) == 2 | sequence1( bps(i,1) ) == 3 ) 
      f_GC = f_GC + 1;
    end
  end
  
  if (N_BP > 0 )
    f_GC = f_GC/ N_BP;
    Tm = ionic_strength_correction( Tm, monovalent_concentration, ...
				    divalent_concentration, f_GC, N_BP );
  end
else
  delG = 1000;
  Tm = -999;
  bps = [];
end
  
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = RC( i1, i2 );
x = ( (i1+i2) == 5 );
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = loop_penalty( sequence1, sequence2, ...
			   T, delG_NN, ...
			   delG_mismatch, delG_AT_closing_penalty);
%Useful for tests:
x = 0.0;

num1 = length( sequence1 );
num2 = length( sequence2 );

R = 1.9978 / 1000; % Gas constant 
n_JS = 2.44; % Jacobsen/Stockmayer coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( num1 == 2 & num2 == 2 )
  % No loop?
  fprintf( 'Hey, should not be calculating a loop penalty for this chunk!\n');
  x = 999999;
  return;
elseif ( num1 == 2 | num2 == 2 ) 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Bulge.
  
  n_bulge = max( num1, num2 ) - 2;
  if (n_bulge == 1 ) % Assume ends stack.
    x = x + delG_NN( sequence1(1), sequence1(end) );
  end
  
  % From SantaLucia, Ann Rev 2004.
  bulge_penalty = [4.0,2.9,3.1,3.2,3.3,3.5,3.7];
  n_bulge_known = length( bulge_penalty );
  if ( n_bulge <= n_bulge_known  )
    x = x + bulge_penalty( n_bulge );
  else 
    x = x + bulge_penalty( n_bulge_known ) + n_JS * R * T * log( n_bulge/n_bulge_known );
  end
  
  x = x + delG_AT_closing_penalty( sequence1(1) ) + delG_AT_closing_penalty( sequence1(2) );
  return;
elseif ( num1 == 3 & num2 == 3 )
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Mismatch  
  x = x + delG_mismatch( sequence1(2), sequence2(2), sequence1(1) );
  x = x + delG_mismatch( sequence2(2), sequence1(2), sequence2(3) );
  if ( RC( sequence1(2), sequence2(2)) )
    x = 999.999; %This case should not be considered a mismatch.
  end
  return;
else
  % From SantaLucia, Ann Rev 2004.
  internal_loop_penalty = [999,999,3.2,3.6,4.0,4.4,4.6,4.8,4.9,4.9];
  n_internal_loop_known = length( internal_loop_penalty );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Internal loops
  N = (num1+num2-4);
  if N <= n_internal_loop_known
    x = x + internal_loop_penalty( N );
  else 
    x = x + internal_loop_penalty( n_internal_loop_known ) + n_JS * R * T * log( N/n_internal_loop_known );
  end

  x = x + abs( num1 - num2 ) * 0.3;
  
  x = x + delG_mismatch( sequence1(2), sequence2(2), sequence1(1) );
  x = x + delG_mismatch( sequence2(end-1), sequence1(end-1), sequence2(end) );
  return;
end


return ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_bp( bps, sequence1_input, sequence2_input );

num1 = length( sequence1_input );
num2 = length( sequence2_input );

%Print out second sequence in reverse order.
sequence1 = sequence1_input;
sequence2 = sequence2_input(end:-1:1);
% Reverse order because first base pairs show up last in backtracking.
bps = bps( end:-1:1, : );

prev_pos1 = 0;
prev_pos2 = 0;

line_forward = '';
line_reverse = '';
line_bp = '';

for i = 1:size(bps,1)
  pos1 = bps(i,1);
  pos2 = bps(i,2);
  overhang = (pos1 - prev_pos1) - (pos2 - prev_pos2);
  if overhang > 0
    for k= 1:overhang
      if ( prev_pos1 > 0 )
	line_reverse =[ line_reverse, '-'];
      else 
	line_reverse =[ line_reverse, ' '];
      end
    end
  else    
    for k= 1:abs(overhang)
      if ( prev_pos1 > 0 )
	line_forward =[ line_forward, '-'];
      else
	line_forward =[ line_forward, ' '];      
      end      
    end    
  end

  line_forward = [line_forward, sequence1((prev_pos1+1):pos1) ];
  line_reverse = [line_reverse, sequence2((prev_pos2+1):pos2) ];

  for k= 1:(max( pos1-prev_pos1, pos2-prev_pos2)-1)
    line_bp = [ line_bp, ' '];   
  end
  line_bp = [line_bp,'|'];

  prev_pos1 = pos1;
  prev_pos2 = pos2;
end

line_forward =[ line_forward, sequence1( (prev_pos1+1):end)];
line_reverse =[ line_reverse, sequence2( (prev_pos2+1):end)];

fprintf(1,'%s\n%s\n%s\n\n',line_forward, line_bp, line_reverse );
