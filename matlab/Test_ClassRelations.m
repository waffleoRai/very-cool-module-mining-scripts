%
%%

addpath('./util');
SAVE_DIR = 'D:\usr\bghos\code\ts3_common_re\matlab';
PROG_INFO_PATH = [SAVE_DIR '\matgraph.mat'];
LINK_TABLE_SAVE_PATH = [SAVE_DIR 'ftlinks.mat'];

emptyCallAddress = uint64(0x10a908dc);

[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(PROG_INFO_PATH);

% funcList = fixFuncList(funcList);
% ProgramGraph.saveTables(PROG_INFO_PATH, symTable, refTable, funcList);

classAnalyzer = ClassFamilyAnalyzer;
classAnalyzer.emptyCallAddress = emptyCallAddress;

if ~isfile(LINK_TABLE_SAVE_PATH)
    classAnalyzer = classAnalyzer.initialize(symTable, funcList, true);
    nodeLinkTable = classAnalyzer.ftNodeTable;
    linkTable = classAnalyzer.ftLinkTable;
    save(LINK_TABLE_SAVE_PATH, 'nodeLinkTable', 'linkTable');
else
    load(LINK_TABLE_SAVE_PATH, 'nodeLinkTable', 'linkTable');

    % nodeLinkTable = fixNodeTable(nodeLinkTable);
    % save(LINK_TABLE_SAVE_PATH, 'nodeLinkTable', 'linkTable');

    classAnalyzer.ftNodeTable = nodeLinkTable;
    classAnalyzer.ftLinkTable = linkTable;
    classAnalyzer.srcSymTable = symTable;
    classAnalyzer.srcFTTable = struct2table(funcList);
end
clear linkTable nodeLinkTable

% testAddr = uint64(0x10b52d34);
% ftIdx = find(obj.ftNodeTable{:, 'address'} == testAddr, 1);

%classAnalyzer.renderGraph(500);
%classAnalyzer.visualizeClades(501, 5);

classAnalyzer = classAnalyzer.updateCommonFuncs();

% classAnalyzer = classAnalyzer.refineSmallClade(244);
% classAnalyzer = classAnalyzer.updateCladeMemberFunctionTerritories(244, 1);
%classAnalyzer = classAnalyzer.analyzeClade(244);
%classAnalyzer = classAnalyzer.updateCladeMemberFunctionTerritories(244, 2);

% cCount = size(classAnalyzer.cladeTable, 1);
% for c = 1:cCount
%     classAnalyzer = classAnalyzer.analyzeClade(c);
% end

nodeLinkTable = classAnalyzer.ftNodeTable;
linkTable = classAnalyzer.ftLinkTable;
cladeTable = classAnalyzer.cladeTable;

%printCladeReports(SAVE_DIR, classAnalyzer);

%---------------------------------------------------------------------

function printCladeReports(outputDir, classAnalyzer)
%Clade table
    tblPath = [outputDir filesep 'cladeTable.csv'];
    writetable(classAnalyzer.cladeTable, tblPath);
    clear tblPath

%Table for tiny clades (<5)
    isRealClade = classAnalyzer.cladeTable{:, 'MemberCount'} > 1;
    tblPath = [outputDir filesep 'tinyClades.tsv'];
    isTinyClade = and(isRealClade, classAnalyzer.cladeTable{:, 'MemberCount'} < 5);
    fh = fopen(tblPath, 'w');
    fprintf(fh, 'Name\tIndex\tSize\tMembers\n');
    tcList = find(isTinyClade)';
    for i = 1:size(tcList, 2)
        c = tcList(i);
        cn = classAnalyzer.cladeTable{c, 'Autoname'};
        memBool = classAnalyzer.ftNodeTable{:, 'clade'} == c;
        memList = classAnalyzer.ftNodeTable{memBool, 'addressString'};
        memCount = size(memList, 1);
        memPrint = [];
        for m = 1:memCount
            if isempty(memPrint)
                memPrint = memList{m};
            else
                memPrint = [memPrint ';' memList{m}];
            end
        end
        fprintf(fh, '%s\t%d\t%d\t%s\n', cn, c, memCount, memPrint);
        clear cn c memPrint memCount memBool m
    end
    fclose(fh);
    clear tblPath fh tcList i isTinyClade

%Table for medium clades (5-100)
    tblPath = [outputDir filesep 'mediumClades.tsv'];
    isMediumClade = and(isRealClade, classAnalyzer.cladeTable{:, 'MemberCount'} >= 5);
    isMediumClade = and(isMediumClade, classAnalyzer.cladeTable{:, 'MemberCount'} <= 100);
    fh = fopen(tblPath, 'w');
    fprintf(fh, 'Name\tIndex\tSize\tMembers\n');
    mcList = find(isMediumClade)';
    for i = 1:size(mcList, 2)
        c = mcList(i);
        cn = classAnalyzer.cladeTable{c, 'Autoname'};
        memBool = classAnalyzer.ftNodeTable{:, 'clade'} == c;
        memList = classAnalyzer.ftNodeTable{memBool, 'addressString'};
        memCount = size(memList, 1);
        memPrint = [];
        for m = 1:memCount
            if isempty(memPrint)
                memPrint = memList{m};
            else
                memPrint = [memPrint ';' memList{m}];
            end
        end
        fprintf(fh, '%s\t%d\t%d\t%s\n', cn, c, memCount, memPrint);
        clear cn c memPrint memCount memBool m
    end
    fclose(fh);
    clear isMediumClade tblPath fh mcList i
%TODO

%Info about mega clades (100+)
%TODO
end

% function funcList = fixFuncList(funcList)
%     funcTable = struct2table(funcList);
% 
%     funcTable{:, 'len'} = 0;
%     ftCount = size(funcTable, 1);
%     for i = 1:ftCount
%         tt = funcTable{i, 'table'};
%         ttt = tt{1};
%         funcTable{i, 'len'} = size(ttt, 1);
%     end
% 
%     funcList = table2struct(funcTable);
% end

% function nodeTable = fixNodeTable(nodeTable)
%     nodeTable{:, 'len'} = 0;
%     ftCount = size(nodeTable, 1);
%     for i = 1:ftCount
%         tt = nodeTable{i, 'table'};
%         ttt = tt{1};
%         nodeTable{i, 'len'} = size(ttt, 1);
%     end
% end

