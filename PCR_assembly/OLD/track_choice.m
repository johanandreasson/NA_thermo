function primers = track_choice( primers, choice, split_choice, pos1, pos2 );

if ( choice( pos1, pos2 ) == 1 ) 
  count = size( primers, 2 );
  segment1_end = split_choice{pos1,pos2}(1);
  segment2_begin = split_choice{pos1,pos2}(2);
  primers(:, count+1 ) = [pos1,segment1_end,1];
  primers(:, count+2 ) = [segment2_begin,pos2,-1];

  plot( [pos1, pos2],[pos1 pos1],'k','linewidth',2);
  plot( [segment1_end, pos2],[segment1_end, pos1],'k','linewidth',2);
  plot( [pos2, pos2],[pos2 pos1],'w','linewidth',2);
  plot( [segment2_begin, pos2],[segment2_begin, pos1],'w','linewidth',2);
  
  plot( pos2,pos1,'ko','markerfacecolor','k');
  plot( pos1,pos1,'ko','markerfacecolor','k');
  plot( segment1_end, segment1_end,'ko','markerfacecolor','k');
  plot( segment2_begin, segment2_begin,'ko','markerfacecolor','w');
  plot( pos2,pos2,'ko','markerfacecolor','w');
  
else
  if ( choice( pos1, pos2 ) < 2 )
    fprintf( 1, 'Problem with code!\n');
  end
  segment1_end = split_choice{pos1,pos2}(1);
  segment2_begin = split_choice{pos1,pos2}(2);
  plot( [pos2, segment1_end],[pos1, pos1],'linewidth',2, ...
	'color',[0.5 0.5 0.5]);
  plot( [pos2, pos2],[pos1, segment2_begin],'linewidth',2, ...
	'color',[0.5 0.5 0.5]);

  plot( pos2,pos1,'ko','markerfacecolor','k');
  plot( segment1_end,pos1,'ko','markerfacecolor','k');
  plot( pos2, segment2_begin,'ko','markerfacecolor','k');

  primers = track_choice( primers, choice, split_choice, pos1, segment1_end );
  primers = track_choice( primers, choice, split_choice, segment2_begin, pos2 );
end

return;
