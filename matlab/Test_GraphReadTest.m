%
%%

addpath('./util');

INPUT_PATH = 'D:\usr\bghos\code\ts3_common_re\ts3_hicpp\graph.tsv';
SAVE_PATH = 'D:\usr\bghos\code\ts3_common_re\matlab\matgraph.mat';
PTR_SIZE_BYTES = 4;

%symTable = ProgramGraph.readProgramGraph(INPUT_PATH, true);
refTable = ProgramGraph.symTable2RefTable(symTable, true, 8);
%[symTable, refTable, funcTables] = ProgramGraph.processRODataTables(symTable, refTable, PTR_SIZE_BYTES, true);

ProgramGraph.saveTables(SAVE_PATH, symTable, refTable, funcTables);