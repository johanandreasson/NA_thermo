function [aas, codons] = show_ORF(sequence, cod_offset, std_dir)

cod_tab = codon_table;

% prepare sequence and ORF
if std_dir == -1; sequence = reverse_complement(sequence); end;
sequence = strrep(sequence, 'T', 'U');
cod_seq = sequence(1+cod_offset: length(sequence)- mod((length(sequence)- cod_offset),3));
codons = cell(length(cod_seq)/3,1);

% print sense sequence in triplets
fprintf(2,'\n(sense strand)\n');
fprintf(repmat(' ',cod_offset+1,1));
for i =1:3:length(cod_seq)
    codons{(i-1)/3+1} = cod_seq(i:i+2);
    fprintf([codons{(i-1)/3+1},' ']);
    if strcmp(find_aa_from_codon(cod_seq(i:i+2),cod_tab),'X_X');
        codons = codons(1:(i-1)/3+1);
        break;
    end;
end;
fprintf('\n');

% print 3-letter peptide chain
aas = cell(length(cod_seq)/3,1);
fprintf([repmat(' ',1,length(1:cod_offset)),' ']);
for i =1:3:length(cod_seq)
    aas{(i-1)/3+1} = find_aa_from_codon(cod_seq(i:i+2),cod_tab);
    if strcmp(aas{(i-1)/3+1},'X_X');
        aas = aas(1:(i-1)/3+1);
        fprintf(2, 'X_X ');
        break;
    end;
    fprintf([aas{(i-1)/3+1},' ']);
    
end;
fprintf('\n');
% print 1-letter peptide chain
fprintf([repmat(' ',1,length(1:cod_offset)),' ']);
for i =1:length(aas)
    fprintf([' <strong>', aa_name_convert(aas{i}, cod_tab),'</strong>  ']);
end;
fprintf('\n');

% print input sequence in triplets
if std_dir == -1;
    fprintf([strrep(complement(sequence(1:cod_offset)),'T','U'),' ']);
    for i =1:3:length(cod_seq)
        fprintf([strrep(complement(cod_seq(i:i+2)),'T','U'),' ']);
    end;
    fprintf([strrep(complement(sequence((length(sequence)- mod((length(sequence)- cod_offset),3)+1):end)),'T','U'), '\n']);
else
    fprintf([sequence(1:cod_offset),' ']);
    for i =1:3:length(cod_seq)
        fprintf([cod_seq(i:i+2),' ']);
    end;
    fprintf([sequence((length(sequence)- mod((length(sequence)- cod_offset),3)+1):end), '\n']);
end;
% print input label
fprintf(2,'(input strand,');
if std_dir == -1;
    fprintf(2,' <strong>REVERSED</strong>)\n\n');
else
    fprintf(2,' <strong>ORIGINAL</strong>)\n\n');
end;


function aa_short = aa_name_convert(aa_long, codon_table)
% convert 3-letter AA name to 1-letter
aa_ind = mod(strmatch(aa_long,codon_table)+20,21)+1;
aa_short = codon_table{aa_ind, 2};

