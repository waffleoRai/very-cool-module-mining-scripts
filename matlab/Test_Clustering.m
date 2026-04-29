%
%%

addpath('./util');
SAVE_PATH = 'D:\usr\bghos\code\ts3_common_re\matlab\matgraph.mat';

heatmapRegion = ProgramGraph.genRegionSpecStruct();
heatmapRegion.section = '.text';
heatmapRegion.startAddress = uint64(0x107728B0);
heatmapRegion.endAddress = uint64(0x1099A6C0);

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(SAVE_PATH);

cAnalyzer = ClusterAnalyzer;
[symList, boundScores] = cAnalyzer.evaluate(refTable, heatmapRegion, true);