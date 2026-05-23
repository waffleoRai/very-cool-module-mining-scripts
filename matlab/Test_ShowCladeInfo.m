%
%%
addpath('./util');
addpath('./gui');

SAVE_DIR = 'D:\usr\bghos\code\ts3_common_re\matlab';
PROG_INFO_PATH = [SAVE_DIR '\matgraph.mat'];
LINK_TABLE_SAVE_PATH = [SAVE_DIR 'ftlinks.mat'];

CLADE_ID = 146;

emptyCallAddress = uint64(0x10a908dc);

%------------------ Load Precalculated Data ------------------

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(PROG_INFO_PATH);

classAnalyzer = ClassFamilyAnalyzer;
classAnalyzer.emptyCallAddress = emptyCallAddress;

load(LINK_TABLE_SAVE_PATH, 'nodeLinkTable', 'linkTable');
classAnalyzer.ftNodeTable = nodeLinkTable;
classAnalyzer.ftLinkTable = linkTable;
classAnalyzer.srcSymTable = symTable;
classAnalyzer.srcFTTable = struct2table(funcList);
clear linkTable nodeLinkTable

%------------------ Run Updates ------------------

classAnalyzer = classAnalyzer.updateCommonFuncs();

linkTable = classAnalyzer.ftLinkTable;
nodeLinkTable = classAnalyzer.ftNodeTable;

%------------------ Analysis and Visualization ------------------

myNodes = nodeLinkTable(nodeLinkTable{:, 'clade'} == CLADE_ID, :);

% Graph Render
classAnalyzer.visualizeClade(CLADE_ID, 1);

% Boolean Child (Row) to Parent (Column) Map
myNodesStr = table2struct(myNodes)';
myNodesStr = ClassFamilyAnalyzer.flagPossibleLinkDirections_NodeArray(myNodesStr);

nCount = size(myNodes, 1);
xlbl = myNodes{:, 'addressString'};
ylbl = xlbl;
pcBoolMtx = false(nCount, nCount);
mCountMtx = zeros(nCount, nCount);
%TODO Vectorize
for cidx = 1:nCount
    ltbl = myNodesStr(cidx).linkTable;
    if isempty(ltbl)
        continue;
    end
    lCount = size(ltbl, 1);
    for lidx = 1:lCount
        psym = ltbl{lidx, 'symbolId'};
        pidx = find((myNodes{:, 'symbolId'} == psym), 1);
        if(ltbl{lidx, 'cpFlag'})
            pcBoolMtx(cidx, pidx) = true;
        end
        mCountMtx(cidx, pidx) = ltbl{lidx, 'fullMatches'};
    end
end
clear ltbl cidx pidx lidx lCount psym

figHandle = figure(2);
clf;
hm = heatmap(xlbl, ylbl, uint8(pcBoolMtx));
hm.Colormap = [1 0 0; 0 1 0];
hm.ColorLimits = [0 1];
hm.CellLabelColor = 'none';
hm.GridVisible = 'off';
title('Possible Parent/Child Flags');
xlabel('Parent');
ylabel('Child');

% Full Match Heatmap

[xx, yy] = meshgrid(1:nCount, 1:nCount);
diagMask = (xx == yy);
mCountMtx(diagMask) = NaN;
clear xx yy

figHandle = figure(3);
clf;
hm = heatmap(xlbl, ylbl, mCountMtx);
%hm.Colormap = turbo;
%hm.ColorLimits = [0 1];
%hm.CellLabelColor = 'none';
hm.GridVisible = 'off';
title('Function Matches');
xlabel('Parent');
ylabel('Child');
clear hm figHandle xlbl ylbl mCountMtx pcBoolMtx

% Function Match Table (With matches to consensus colored)
[classAnalyzer, funcSymIdx, entryFreq] = classAnalyzer.getCladeConsensusFunctionTable(CLADE_ID);
if ~isempty(funcSymIdx)
    isNullEmpty = or((funcSymIdx == 0), (funcSymIdx == classAnalyzer.emptyCallSymbolId));
    fprintf('--------- Consensus Table ---------\n');
    conTableSize = size(funcSymIdx, 2);
    for i = 1:conTableSize
        fprintf('\t');
        myVal = funcSymIdx(i);
        myFreq = entryFreq(i);
        freqPercent = (myFreq ./ nCount) * 100;
        if ~isNullEmpty(i)
            tidx = find((symTable{:, 'Index'} == myVal), 1);
            symaddrstr = symTable{tidx, 'AddressString'};
            fprintf('%s (Sym %d) [Frequency: %d/%d (%.2f%%)]', symaddrstr, myVal, myFreq, nCount, freqPercent);
        else
            if myVal == 0
                fprintf('<NULL>');
            else
                fprintf('<EMPTY> [Frequency: %d/%d (%.2f%%)]', myFreq, nCount, freqPercent);
            end
        end
        fprintf('\n');
    end
    clear myVal conTableSize tidx symaddrstr 
else
    fprintf('WARNING: No Function Consensus Table could be derived\n');
end

figHandle = VisualizeCladeFuncTable(myNodes, funcSymIdx, classAnalyzer.emptyCallSymbolId);