function design_rescue_plate (sequence, structure, offset, mode_flag, seqpos_sub, Tm, mut_file, mut_label_start)

% construct_names = DESIGN_RESCUE_PLATE (sequence, structure, [offset], ...
%                                       [mode_flag], [seqpos_sub], [Tm], ...
%                                       [mut_file], [mut_label_start]);
%
% Designs mutation/rescue single/double/triple mutants in 96-well format.
%
% =Input=
%   sequence            Sequence of WT, format in string.
%   structure           Dot-bracket annotation of secondary structure to be
%                           tested. Same length as sequence, format in
%                           string.
%   [offset]            Offset for sequence numbering. Offset is natural 
%                           numbering minums final numbering. Format in 
%                           double, default is 0. e.g., A1 in sequence,
%                           final numbering is 50, offset is -49.
%   [mode_flag]         Mode choices, format in cell of strings, default is
%                           {'single','swap','quartet'}.
%                       'single/double/triple' denotes base-pair mutants, 
%                           whether single or double or triple base-pair.
%                       'swap/cross' denotes nucleotide identity, whether
%                           swapping (A>T,T>A) or crossing (A>C,T>G).
%                       'only/quartet' denotes mutants library, whether
%                           only the double mutant, or quartet including
%                           single mutants as well.
%   [seqpos_sub]        Mutation region. Format in double array, default is
%                           [1 : length(sequence)]. Value includes offset
%                           already. e.g. [1:59, 69:80].
%   [Tm]                Minimum melting temeperature for primer assembly.
%                           Format in double, defaut is 60.
%   [mut_file]          Filename of text output file. Format as string, default
%                          'Mutants'. All sequences will be write in this
%                          file, '.txt' extension automatically append.
%   [mut_label_start]   Prefix to all construct names, e.g. '16S-'. Format
%                          as string, default none.
%
%
% by T47, Oct 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('offset','var'); offset = []; end;
if ~exist('mode_flag','var'); mode_flag = {}; end;
if ~exist('seqpos_sub','var'); seqpos_sub = []; end;
if ~exist('Tm','var'); Tm = 60; end;
if ~exist('mut_file','var'); mut_file = ''; end;
if ~exist('mut_label_start','var'); mut_label_start = ''; end;

% find double mutants
mutations = design_rescue (sequence, structure, offset, mode_flag, seqpos_sub);
% find primer assembly
primers = design_primers(sequence, Tm);
% write mutants list
design_mutants(primers, sequence, mutations, offset, mut_file, mut_label_start, 1);
% find mutant assembly
[sequences_to_order,sequence_names] = mutate_primers_each( primers, [mut_file,'.txt'] );
% write layout files
output_sequences_to_order_excel_file( sequences_to_order, primers, mut_label_start(1:end-1) );
output_sequences_to_order_96well_diagram( sequences_to_order, primers, mut_label_start(1:end-1) );
output_construct_names( sequence_names, mut_label_start(1:end-1) );

