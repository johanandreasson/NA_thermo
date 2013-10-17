tail2 = 'AAAGAAACAACAACAACAAC';
sequences = get_possible_tails();

[delG, Tm , free_end, bps] = check_orthogonal( tail2 , sequences );

plot( Tm );
xlabel( 'sequence number')
ylabel( 'Tm at 1 uM, 0.1 M NaCl, 1.5 mM MgCl_2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tm_cut = 25;
good_seqs = find( Tm < Tm_cut & free_end );

% Need to check again -- on reverse_complements! 
[delG_RC, Tm_RC , free_end_RC, bps_RC] = check_orthogonal( reverse_complement(tail2) , ...
			reverse_complement( {sequences{good_seqs}} ));

plot( Tm( good_seqs), 'b' );
hold on
plot( Tm_RC, 'r');
good_seqs_final = good_seqs( find( Tm_RC < Tm_cut & free_end_RC ) );


[y,i ] = min( Tm( good_seqs_final ) );
free_end( good_seqs_final( i ) )  % better be 1
tail3 = sequences{ good_seqs_final( i ) }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sequences2 = { sequences{ good_seqs_final } };
[delG2, Tm2 , free_end2, bps2] = check_orthogonal( tail3 , sequences2 );
[delG2_RC, Tm2_RC , free_end2_RC, bps2_RC] = check_orthogonal( reverse_complement(tail3) , ...
						  reverse_complement( ...
						      sequences2 ) );
good_seqs2 = find( Tm2 < Tm_cut & free_end2 & Tm2_RC < Tm_cut & ...
		   free_end2_RC );
for i = good_seqs2
  calc_misprime_delG( reverse_complement( tail3 ), sequences2{i} );
end
length( good_seqs2 )
[y,i ] = min( Tm2( good_seqs2) );
tail4 = sequences2{ good_seqs2( i ) };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sequences3 = { sequences2{ good_seqs2 } };
[delG3, Tm3 , free_end3, bps3] = check_orthogonal( tail4 , sequences3 );
[delG3_RC, Tm3_RC , free_end3_RC, bps3_RC] = ...
    check_orthogonal( reverse_complement(tail4) ,  reverse_complement( ...
						      sequences3 ) );
good_seqs3 = find( Tm3 < Tm_cut & free_end3 & Tm3_RC < Tm_cut & ...
		   free_end3_RC );
for i = good_seqs3
  calc_misprime_delG( reverse_complement( tail3 ), sequences3{i} );
end
length( good_seqs3 )

% Problem with GAAC sequence!
%good_seqs3_filter = [];
%for k = good_seqs3;
%  if isempty(  strfind( sequences3{ k }, 'GAAC' ) ) 
%    good_seqs3_filter = [ good_seqs3_filter k ];
%  end
%end

[y,i ] = min( Tm3( good_seqs3) );
tail5 = sequences3{ good_seqs3( i ) };

tail5b = sequences3{ good_seqs3(3) };
tail5c = sequences3{ good_seqs3(5) };
tail5d = sequences3{ good_seqs3(6) };

tail5e = sequences3{ good_seqs3(1) };
tail5f = sequences3{ good_seqs3(2) };
tail5g = sequences3{ good_seqs3(7) };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updated scripts -- now check for mispriming in F and R more
% selectively -- just end of *primer*
sequences4 = { sequences3{ good_seqs3 } };
sequences4 = {sequences4{2}};
[delG4, Tm4 , free_end4_F, free_end4_R, bps4] = check_orthogonal( tail5 , sequences4 );
[delG4_RC, Tm4_RC , free_end4_F_RC, free_end4_R_RC, bps4_RC] = ...
    check_orthogonal( reverse_complement(tail5) ,  reverse_complement( ...
						      sequences4 ) );
good_seqs4 = find( Tm4 < Tm_cut & free_end4_F & Tm4_RC < Tm_cut & free_end4_R_RC );

for i = good_seqs4
  calc_misprime_delG( reverse_complement( tail4 ), sequences4{i} );
end
length( good_seqs4 )
[y,i ] = min( Tm4( good_seqs4) );
tail6 = sequences4{ good_seqs4( i ) };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for more orthogonal sequences, this time with more G's,
% fewer A's.
% New pool of sequences.
% Check against tail2, tail3, tail4, tail5, tail6, one-by-one.
sequences_new = get_possible_tails2();
length( sequences_new )

[delG_new, Tm_new , free_end_new_F, free_end_new_R] = ...
    check_orthogonal( tail2 , sequences_new );
[delG_new_RC, Tm_new_RC , free_end_new_F_RC, free_end_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail2) ,  reverse_complement( ...
	sequences_new ) );

good_seqs_new = find( Tm_new < Tm_cut & free_end_new_F & Tm_new_RC < Tm_cut & free_end_new_R_RC );
sequences2_new = { sequences_new{ good_seqs_new }  };

[delG2_new, Tm2_new , free_end2_new_F, free_end2_new_R] = ...
    check_orthogonal( tail3 , sequences2_new );
[delG2_new_RC, Tm2_new_RC , free_end2_new_F_RC, free_end2_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail3) ,  reverse_complement( ...
	sequences2_new ) );

good_seqs2_new = find( Tm2_new < Tm_cut & free_end2_new_F & Tm2_new_RC < Tm_cut & free_end2_new_R_RC );
sequences3_new = { sequences2_new{ good_seqs2_new }  };

[delG3_new, Tm3_new , free_end3_new_F, free_end3_new_R] = ...
    check_orthogonal( tail4 , sequences3_new );
[delG3_new_RC, Tm3_new_RC , free_end3_new_F_RC, free_end3_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail4) ,  reverse_complement( ...
	sequences3_new ) );

good_seqs3_new = find( Tm3_new < Tm_cut & free_end3_new_F & Tm3_new_RC < Tm_cut & free_end3_new_R_RC );
sequences4_new = { sequences3_new{ good_seqs3_new }  };


[delG4_new, Tm4_new , free_end4_new_F, free_end4_new_R] = ...
    check_orthogonal( tail5 , sequences4_new );
[delG4_new_RC, Tm4_new_RC , free_end4_new_F_RC, free_end4_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail5) ,  reverse_complement( ...
	sequences4_new ) );

good_seqs4_new = find( Tm4_new < Tm_cut & free_end4_new_F & Tm4_new_RC < Tm_cut & free_end4_new_R_RC );
sequences5_new = { sequences4_new{ good_seqs4_new }  };


[delG5_new, Tm5_new , free_end5_new_F, free_end5_new_R] = ...
    check_orthogonal( tail6 , sequences5_new );
[delG5_new_RC, Tm5_new_RC , free_end5_new_F_RC, free_end5_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail6) ,  reverse_complement( ...
	sequences5_new ) );

good_seqs5_new = find( Tm5_new < Tm_cut+1 & free_end5_new_F & Tm5_new_RC < Tm_cut+1 & free_end5_new_R_RC );
sequences6_new = { sequences5_new{ good_seqs5_new }  };

[y,i ] = min( Tm5_new( good_seqs5_new ) );
tail7 = sequences5_new{ good_seqs5_new( i ) }

%%%%%%%%%%%%%%%%

[delG6_new, Tm6_new , free_end6_new_F, free_end6_new_R] = ...
    check_orthogonal( tail7 , sequences6_new );
[delG6_new_RC, Tm6_new_RC , free_end6_new_F_RC, free_end6_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail7) ,  reverse_complement( ...
	sequences6_new ) );

good_seqs6_new = find( Tm6_new < Tm_cut+1 & free_end6_new_F & Tm6_new_RC < Tm_cut+1 & free_end6_new_R_RC );
sequences7_new = { sequences6_new{ good_seqs6_new }  };

[y,i ] = min( Tm6_new( good_seqs6_new ) );
tail8 = sequences6_new{ good_seqs6_new( i ) }

%%%%%%%%%%%%%%%%%
% Nothing found in following filter.
[delG7_new, Tm7_new , free_end7_new_F, free_end7_new_R] = ...
    check_orthogonal( tail8 , sequences7_new );
[delG7_new_RC, Tm7_new_RC , free_end7_new_F_RC, free_end7_new_R_RC] = ...
    check_orthogonal( reverse_complement(tail8) ,  reverse_complement( ...
	sequences7_new ) );

good_seqs7_new = find( Tm7_new < Tm_cut+1 & free_end7_new_F & Tm7_new_RC < Tm_cut+1 & free_end7_new_R_RC );
sequences8_new = { sequences7_new{ good_seqs7_new }  };

[y,i ] = min( Tm7_new( good_seqs7_new ) );
tail9 = sequences7_new{ good_seqs7_new( i ) }


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Had another good idea -- shouldn't we end the primers with C and
% G and U -- not A?
sequences_new3 = get_possible_tails3();

[delG_new3, Tm_new3 , free_end_new3_F, free_end_new3_R] = ...
    check_orthogonal( tail2 , sequences_new3 );
[delG_new3_RC, Tm_new3_RC , free_end_new3_F_RC, free_end_new3_R_RC] = ...
    check_orthogonal( reverse_complement(tail2) ,  reverse_complement( ...
	sequences_new3 ) );

good_seqs_new3 = find( Tm_new3 < Tm_cut & free_end_new3_F & Tm_new3_RC < Tm_cut & free_end_new3_R_RC );
sequences2_new3 = { sequences_new3{ good_seqs_new3 }  };

[delG2_new3, Tm2_new3 , free_end2_new3_F, free_end2_new3_R] = ...
    check_orthogonal( tail3 , sequences2_new3 );
[delG2_new3_RC, Tm2_new3_RC , free_end2_new3_F_RC, free_end2_new3_R_RC] = ...
    check_orthogonal( reverse_complement(tail3) ,  reverse_complement( ...
	sequences2_new3 ) );

good_seqs2_new3 = find( Tm2_new3 < Tm_cut & free_end2_new3_F & Tm2_new3_RC < Tm_cut & free_end2_new3_R_RC );
sequences3_new3 = { sequences2_new3{ good_seqs2_new3 }  };

[delG3_new3, Tm3_new3 , free_end3_new3_F, free_end3_new3_R] = ...
    check_orthogonal( tail4 , sequences3_new3 );
[delG3_new3_RC, Tm3_new3_RC , free_end3_new3_F_RC, free_end3_new3_R_RC] = ...
    check_orthogonal( reverse_complement(tail4) ,  reverse_complement( ...
	sequences3_new3 ) );

good_seqs3_new3 = find( Tm3_new3 < Tm_cut & free_end3_new3_F & Tm3_new3_RC < Tm_cut & free_end3_new3_R_RC );
sequences4_new3 = { sequences3_new3{ good_seqs3_new3 }  };

[y,i ] = min( Tm3_new( good_seqs3_new3 ) );
tail10 = sequences3_new3{ good_seqs3_new3( i ) }


[delG4_new3, Tm4_new3 , free_end4_new3_F, free_end4_new3_R] = ...
    check_orthogonal( tail10 , sequences4_new3 );
[delG4_new3_RC, Tm4_new3_RC , free_end4_new3_F_RC, free_end4_new3_R_RC] = ...
    check_orthogonal( reverse_complement(tail10) ,  reverse_complement( ...
	sequences4_new3 ) );

good_seqs4_new3 = find( Tm4_new3 < Tm_cut+1 & free_end4_new3_F & Tm4_new3_RC < Tm_cut+1 & free_end4_new3_R_RC );
sequences5_new3 = { sequences4_new3{ good_seqs4_new3 }  };

[y,i ] = min( Tm4_new3( good_seqs4_new3 ) );
tail11 = sequences4_new3{ good_seqs4_new3( i ) }
