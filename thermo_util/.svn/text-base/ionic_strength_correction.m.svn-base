function Tm_corrected = ionic_strength_correction( Tm, monovalent_concentration, ...
						  divalent_concentration, ...
						  f_GC, N_BP);
% From Owczarzy et al., Biochemistry, 2008.

R = sqrt( divalent_concentration ) / monovalent_concentration;

Tm_corrected = Tm;
if ( R < 0.22 )
  %Monovalent dominated
  x = log( monovalent_concentration );
  Tm_corrected = 1 / ...
      ( (1/Tm) + (4.29*f_GC - 3.95)*1e-5 * x + 9.40e-6*x*x);
else 
  %Divalent dominated
  a = 3.92e-5;
  b = -9.11e-6;
  c = 6.26e-5;
  d = 1.42e-5;
  e = -4.82e-4;
  f = 5.25e-4;
  g = 8.31e-5;

  if ( R < 6.0)
    % Some competition from monovalent
    y = monovalent_concentration;
    a = 3.92e-5 * ( 0.843 - 0.352*sqrt(y)*log(y));
    d = 1.42e-5 * ( 1.279 - 4.03e-3 * log(y) - ...
		    8.03e-3 * (log( y ))^2 );
    g = 8.31e-5 * ( 0.486 - 0.258 * log(y) + ...
		    5.25e-3 * (log( y ))^3 );
  end

  x = log( divalent_concentration );
  Tm_corrected = 1 / ...
      ( (1/Tm) + a + b *x  + f_GC*(c + d*x) + ...
	(1/(2*(N_BP-1)))*(e + f*x+g*x*x));
end