function [mut_seq, mut_lab] = design_mutants(primers, sequence, mutations, offset, mutfile, mut_label_start, is_print)
% [mut_seq, mut_lab] = DESIGN_MUTANTS (primers, sequence, mutations, [offset],
%                                       [mutfile], [mut_label_start])
%
% Designs constructs with unlimited specified mutations. Insertion or
%   deletion not supported. Write to sequences and contruct names to text
%   file.
%
% =Input=
%   primers             Primers to assemble the targeted sequence, generated
%                          from design_primers_NEW. Format as cell.
%   sequence            Sequence of WT. Will be converted to DNA sequence.
%   mutations           Set of mutatants. Format as cell. Each construct is
%                          specified by string annotation, e.g. 'A200C'. 
%                          Multiple mutations within one construct is
%                          recognied by a cell of strings, e.g. {'A200C',
%                          'G201T'}. Numbering includes offset.
%   [offset]            The offset of sequence. Format as double, default 0.
%   [mutfile]           Filename of text output file. Format as string, default
%                          'Mutants'. All sequences will be write in this
%                          file, '.txt' extension automatically append.
%   [mut_label_start]   Prefix to all construct names, e.g. '16S-'. Format
%                          as string, default none.
%
% =Output=
%   [mut_seq]           Sequence of all mutants. Format in cell.
%   [mut_lab]           Label of all mutants. Format in cell.
%
% by T47, Aug 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('mut_label_start','var');  mut_label_start = ''; end;
if ~exist('mutfile','var') || isempty(mutfile);  mutfile = 'Mutants'; end;
mutfile = [mutfile, '.txt'];
if ~exist('offset','var');   offset = 0; end;
if ~exist('is_print','var') || isempty(is_print); is_print = 1; end;

sequence = strrep(upper(sequence),'U','T');
for i = 1:length(mutations)
    for j = 1:length(mutations{i})
        mutations{i}{j} = strrep(mutations{i}{j}, 'U', 'T');
    end;
end;

fid = fopen(  mutfile, 'w' );

for i = 1:length( mutations )

  sequence_mut = sequence;
  mutation_set = mutations{i};
  mut_label = mut_label_start;
  
  for a = 1:length( mutation_set )
    mutation = mutation_set{a};
    seqpos = str2num( mutation(2:(end-1)) ) + offset;
    startchar = mutation(1);
    endchar = mutation(end);
    if ( sequence( seqpos ) ~= startchar ) fprintf( 'PROBLEM %s\n', mutation ); end;
    sequence_mut( seqpos ) = endchar;

    if ( a> 1 ); mut_label = [mut_label,';']; end;
    mut_label  = [mut_label, mutation];

    seq_mut{i}=sequence_mut;
    lab_mut{i}=mut_label;
  end;
 
  if strcmp( mut_label, mut_label_start ); mut_label = [mut_label, 'WT']; end;
  if isempty(mutation_set)
      seq_mut{i}=sequence;
      lab_mut{i}=mut_label;
  end;
  
  if is_print; 
      fprintf( 1, '%s %s\n', sequence_mut, mut_label ); 
      fprintf( fid, '%s %s\n', sequence_mut, mut_label );
  end;
end;

fclose( fid );

if ~isempty(primers); mutate_primers( primers, mutfile ); end;

mut_seq = rot90(seq_mut,3);
mut_lab = rot90(lab_mut,3);
