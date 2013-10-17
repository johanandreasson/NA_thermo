function output_num_match( sequence, num_match_forward, num_match_reverse );

COL_SIZE = 150;
N_BP = length( sequence );

for k=1:(floor(N_BP/COL_SIZE)+1)

  allow_forward_line = '';
  allow_reverse_line = '';
  sequence_line = '';
  
  for n = 1:COL_SIZE
    pos = (k-1) * COL_SIZE+n;
    
    if (pos < N_BP)
      sequence_line(n) = sequence(pos);
      allow_forward_line(n) = num2str(min(num_match_forward(pos),9));
      allow_reverse_line(n) = num2str(min(num_match_reverse(pos),9));
    end
  end
        
  fprintf('%s\n%s\n%s\n\n',allow_forward_line,sequence_line, allow_reverse_line);
end
