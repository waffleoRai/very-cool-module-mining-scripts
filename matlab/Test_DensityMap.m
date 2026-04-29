%
%%

addpath('./util');
SAVE_PATH = 'D:\usr\bghos\code\ts3_common_re\matlab\matgraph.mat';

BIN_SIZE = 32;

heatmapRegion = ProgramGraph.genRegionSpecStruct();
heatmapRegion.section = '.text';
heatmapRegion.startAddress = uint64(0x107728B0);
heatmapRegion.endAddress = uint64(0x1099A6C0);

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(SAVE_PATH);
refDensityMap = ReferenceHeatmap.genHeatmap(refTable, BIN_SIZE, 0, false, heatmapRegion);

refDensityMap.renderRawDensityMap(1);
refDensityMap = refDensityMap.prepForClusterAnalysis();
refDensityMap.renderPearsonDistanceMap(2);