function sequences = get_possible_tails();


%words = {'CA','CAA','CAAA','CAAAA','CCA','CCAA','CCAAA','CCAAAA', ...
%	 'GAA','GAAA','GAAAA'};
%words = {'CA','CCAAAA' };
words = {'AC','AAC','AAAC','AACC','AAACC','AG','AAG','AAAG','AGG',};
start_words = {'AC','AAC','AAAC','AAAAC'};


MAX_LENGTH = 20;

sequences = {};

for i = 1:length( start_words )
  sequences = get_words( start_words{i}, words, sequences, MAX_LENGTH );
end


% Yo recursion
function sequences = get_words( current_word, words, sequences, MAX_LENGTH )

if ( length( current_word ) >  MAX_LENGTH ) 
  return;
end
if ( length( current_word ) ==  MAX_LENGTH ) 
  if ( length( strfind( current_word,'A') ) > 13 ) 
    sequences{ length( sequences ) + 1 } = current_word;
  end
  return;
end

for i = 1: length( words )
  sequences = get_words( [words{i},current_word], ...
			 words, sequences, MAX_LENGTH  );
end
