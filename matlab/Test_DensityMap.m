%
%%

addpath('./util');
SAVE_DIR = 'D:\usr\bghos\code\ts3_common_re\matlab';
SAVE_PATH = [SAVE_DIR '\matgraph.mat'];

BIN_SIZE = 16;

heatmapRegion = ProgramGraph.genRegionSpecStruct();
heatmapRegion.section = '.text';
heatmapRegion.startAddress = uint64(0x10575C90);
heatmapRegion.endAddress = uint64(0x107D85C0);

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(SAVE_PATH);
refDensityMap = ReferenceHeatmap.genHeatmap(refTable, BIN_SIZE, 0, false, heatmapRegion);

refDensityMap.renderRawDensityMap(100);
refDensityMap = refDensityMap.prepForClusterAnalysis();
refDensityMap.renderPearsonDistanceMap(200);