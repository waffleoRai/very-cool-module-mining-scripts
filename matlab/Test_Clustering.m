%
%%

addpath('./util');
SAVE_DIR = 'D:\usr\bghos\code\ts3_common_re\matlab';
SAVE_PATH = [SAVE_DIR '\matgraph.mat'];

heatmapRegion = ProgramGraph.genRegionSpecStruct();
heatmapRegion.section = '.text';
heatmapRegion.startAddress = uint64(0x10575C90);
heatmapRegion.endAddress = uint64(0x107D85C0);

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(SAVE_PATH);

cAnalyzer = ClusterAnalyzer;
[symList, boundScores] = cAnalyzer.evaluate(refTable, heatmapRegion, true);

figure(2);
clf;
plot(symList, boundScores);

symListShort = symList;
clusterDataSavePath = [SAVE_DIR filesep sprintf('clusterscores_%08x_%08x.mat', heatmapRegion.startAddress, heatmapRegion.endAddress)];
save(clusterDataSavePath, 'symListShort', 'boundScores');

%Print to text file
clusterDataOutputPath = [SAVE_DIR filesep sprintf('clusterscores_%08x_%08x.csv', heatmapRegion.startAddress, heatmapRegion.endAddress)];
symCount = size(symList, 2);
outputTable = table('VariableNames', {'AddressString' 'BoundScore'}, 'VariableTypes', {'string' 'double'}, 'Size', [symCount 2]);
[~, symTableLookup] = ismember(symList, symTable{:, 'Index'});
outputTable{:, 'AddressString'} = symTable{symTableLookup, 'AddressString'};
outputTable{:, 'BoundScore'} = boundScores';
writetable(outputTable, clusterDataOutputPath);