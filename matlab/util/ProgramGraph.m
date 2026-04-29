%
%%

%%
classdef ProgramGraph

    %%
    methods (Static)

        %%
        function regInfoStruct = genRegionSpecStruct()
            regInfoStruct = struct();
            regInfoStruct.section = '.text';
            regInfoStruct.startAddress = uint64(0);
            regInfoStruct.endAddress = uint64(0);
        end

        %%
        function funcTableInfo = genFuncTableStruct(tableAllocRows)
            funcTableInfo = struct();
            funcTableInfo.address = uint64(0);
            funcTableInfo.symbolId = int32(0);
            funcTableInfo.addressString ="NULL";
            funcTableInfo.table = ProgramGraph.genFuncTableContents(tableAllocRows);
        end

        %%
        function funcTable = genFuncTableContents(allocRows)
            varNames = {'Address' 'SymbolId' 'AddressString'};
            varTypes = {'uint64' 'int32' 'string'};
            TableSize = [allocRows size(varNames, 2)];
            funcTable = table(Size=TableSize, VariableTypes=varTypes, VariableNames=varNames);
        end

        %%
        function [varNames, varTypes] = getRefTableColumns()
            varNames = {'FromAddressString' 'ToAddressString' 'FromSymbol' 'ToSymbol' 'FromAddress' 'ToAddress' 'FromSection' 'ToSection'};
            varTypes = {'string' 'string' 'int32' 'int32' 'uint64' 'uint64' 'string' 'string'};
        end

        %%
        function [symTable, refTable, funcTables] = processRODataTables(symTable, refTable, ptrSizeBytes, verbose)
            symTable = sortrows(symTable, 'Address');
            symCount = size(symTable, 1);

            % rtAppendPos = 1;
            % [varNames, varTypes] = ProgramGraph.getRefTableColumns();
            % TableSize = [2048 size(varNames, 2)];
            % rtAppend = table(Size=TableSize, VariableTypes=varTypes, VariableNames=varNames);
            addrFmtChar = 8;
            if ptrSizeBytes == 2
                addrFmtChar = 4;
            elseif ptrSizeBytes == 8
                addrFmtChar = 16;
            end
            addrFmtStr = ['%0' num2str(addrFmtChar) 'X'];

            %Filter for references that are from rdata/rodata to text
            isFromRData = or(strcmp(refTable{:, 'FromSection'}, '.rodata'), strcmp(refTable{:, 'FromSection'}, '.rdata'));
            isToText = strcmp(refTable{:, 'ToSection'}, '.text');
            filterBool = and(isFromRData, isToText);
            ftRefs = refTable(filterBool, :);

            ftSyms = unique(ftRefs{:, 'FromSymbol'})';
            ftCount = size(ftSyms,2);
            funcTables(ftCount) = ProgramGraph.genFuncTableStruct(1);
            if verbose
                fprintf('%d function table candidates found!\n', ftCount);
            end

            for i = 1:ftCount
                if verbose & (mod(i, 10) == 0)
                    percamt = (i ./ ftCount) .* 100;
                    fprintf('\tProcessing candidate %d of %d (%.2f%% complete)...\n', i, ftCount, percamt);
                end

                myId = ftSyms(i);
                mySymIdx = find(symTable{:, 'Index'} == myId, 1);
                myRefs = ftRefs((ftRefs{:, 'FromSymbol'} == myId), :);

                ftStruct = ProgramGraph.genFuncTableStruct(1);
                ftStruct.address = symTable{mySymIdx, 'Address'};
                ftStruct.symbolId = myId;
                ftStruct.addressString = sprintf(addrFmtStr, ftStruct.address);

                myReferees = refTable((refTable{:, 'ToSymbol'} == myId), :);
                myReferees = myReferees(strcmp(myReferees{:, 'FromSection'}, '.text'), :);
                refereeList = myReferees{:, 'FromSymbol'}';
                trefereeCount = size(myReferees, 1);

                %Get address of subsequent symbol so can estimate table end
                tblEndAddr = uint64(0);
                nextSymIdx = mySymIdx + 1;
                if nextSymIdx <= symCount
                    tblEndAddr = symTable{nextSymIdx, 'Address'};
                end
                entryCount = floor((tblEndAddr - ftStruct.address) ./ ptrSizeBytes);
                ftStruct.table = ProgramGraph.genFuncTableContents(entryCount);

                %Isolate referenced symbols
                [~, refSymIdxs] = ismember(myRefs{:, 'ToSymbol'}, symTable{:, 'Index'});
                refdSymTable = symTable(refSymIdxs, :);
                rSymCount = size(refdSymTable, 1);
                for j = 1:rSymCount
                    srefListCell = refdSymTable{j, 'RefereesAddr'};
                    srefList = srefListCell{1};
                    keepBool = srefList >= ftStruct.address;
                    keepBool = and(keepBool, (srefList < tblEndAddr));
                    srefList = srefList(keepBool);
                    srefList_i = floor((srefList - ftStruct.address) ./ ptrSizeBytes) + 1;
                    ftStruct.table{srefList_i, 'Address'} = refdSymTable{j, 'Address'};
                    ftStruct.table{srefList_i, 'SymbolId'} = refdSymTable{j, 'Index'};

                    %Add any indirect referees from other parts of .text...
                    if trefereeCount > 0
                        bigTableIdx = refSymIdxs(j);
                        srefListCell = symTable{bigTableIdx, 'Referees'};
                        srefList = srefListCell{1};
                        srefList = unique([srefList refereeList]);
                        srefListCell{1} = srefList;
                        symTable{bigTableIdx, 'Referees'} = srefListCell;
                    end
                end

                ftStruct.table{:, 'AddressString'} = arrayfun(@(x) sprintf(addrFmtStr, x), ftStruct.table{:, 'Address'}, 'UniformOutput', false);
                funcTables(i) = ftStruct;
            end

            %Regenerate ref table from updated sym table
            if verbose
                fprintf('Updating reference table...\n');
            end
            refTable = ProgramGraph.symTable2RefTable(symTable, verbose, addrFmtChar);
        end

        %%
        function refTable = symTable2RefTable(symTable, verbose, fmtChars)
            if (nargin < 3); fmtChars = 8; end
            [varNames, varTypes] = ProgramGraph.getRefTableColumns();

            if verbose
                fprintf('Counting references...\n');
            end
            totalReferences = 0;
            symCount = size(symTable, 1);
            for s=1:symCount
                srefListCell = symTable{s, 'Referees'};
                srefList = srefListCell{1};
                srefCount = size(srefList, 2);
                totalReferences = totalReferences + srefCount;
            end

            TableSize = [totalReferences size(varNames, 2)];
            refTable = table(Size=TableSize, VariableTypes=varTypes, VariableNames=varNames);
            cpos = 1;
            for s=1:symCount
                if verbose & (mod(s, 1000) == 0)
                    percamt = (s ./ symCount) .* 100;
                    fprintf('\tProcessing symbol %d of %d (%.2f%% complete)...\n', s, symCount, percamt);
                end
                srefListCell = symTable{s, 'Referees'};
                srefList = srefListCell{1};
                [matchBool, matchIdx] = ismember(srefList, symTable{:, 'Index'});
                matchCount = nnz(matchBool);
                if matchCount > 0
                    epos = cpos + matchCount - 1;
                    refTable{cpos:epos, 'ToSymbol'} = symTable{s, 'Index'};
                    refTable{cpos:epos, 'ToAddress'} = symTable{s, 'Address'};
                    refTable{cpos:epos, 'ToSection'} = symTable{s, 'Section'};

                    goodMatch = matchIdx(matchBool);
                    refTable{cpos:epos, 'FromSymbol'} = symTable{goodMatch, 'Index'};
                    refTable{cpos:epos, 'FromAddress'} = symTable{goodMatch, 'Address'};
                    refTable{cpos:epos, 'FromSection'} = symTable{goodMatch, 'Section'};
                    cpos = epos + 1;
                end
            end

            %Trim any unused rows
            cpos = cpos - 1;
            if cpos < totalReferences
                refTable = refTable(1:cpos,:);
            end

            addrFmtStr = ['%0' num2str(fmtChars) 'X'];
            refTable{:, 'ToAddressString'} = arrayfun(@(x) sprintf(addrFmtStr, x), refTable{:, 'ToAddress'}, 'UniformOutput', false);
            refTable{:, 'FromAddressString'} = arrayfun(@(x) sprintf(addrFmtStr, x), refTable{:, 'FromAddress'}, 'UniformOutput', false);
        end

        %%
        function filteredRefTable = filterFROMReferences(refTable, regionSpec)
            filteredRefTable = refTable;
            if isempty(refTable); return; end
            if isempty(regionSpec); return; end

            keepBool = strcmp(refTable{:, 'FromSection'}, regionSpec.section);
            if(regionSpec.startAddress > 0)
                keepBool = and(keepBool, refTable{:, 'FromAddress'} >= regionSpec.startAddress);
            end
            if(regionSpec.endAddress > 0)
                keepBool = and(keepBool, refTable{:, 'FromAddress'} < regionSpec.endAddress);
            end

            filteredRefTable = refTable(keepBool, :);
        end

        %%
        function filteredRefTable = filterTOReferences(refTable, regionSpec)
            filteredRefTable = refTable;
            if isempty(refTable); return; end
            if isempty(regionSpec); return; end

            keepBool = strcmp(refTable{:, 'ToSection'}, regionSpec.section);
            if(regionSpec.startAddress > 0)
                keepBool = and(keepBool, refTable{:, 'ToAddress'} >= regionSpec.startAddress);
            end
            if(regionSpec.endAddress > 0)
                keepBool = and(keepBool, refTable{:, 'ToAddress'} < regionSpec.endAddress);
            end

            filteredRefTable = refTable(keepBool, :);
        end
        
        %%
        function [symTable, refTable] = filterToSection(symTable, refTable, sectionName, startAddress, endAddress, inclRefOutside)
            if nargin < 4; startAddress = 0; end
            if nargin < 5; endAddress = 0; end
            if nargin < 6; inclRefOutside = false; end

            %Symbol table
            if ~isempty(symTable)
                keepBool = strcmp(symTable{:, 'Section'}, sectionName);
                if (startAddress > 0)
                    keepBool = and(keepBool, symTable{:, 'Address'} >= startAddress);
                end
                if (endAddress > 0)
                    keepBool = and(keepBool, symTable{:, 'Address'} < endAddress);
                end

                symTable = symTable(keepBool, :);
            end

            if ~isempty(refTable)
                fromBool = strcmp(refTable{:, 'FromSection'}, sectionName);
                if (startAddress > 0)
                    fromBool = and(fromBool, refTable{:, 'FromAddress'} >= startAddress);
                end
                if (endAddress > 0)
                    fromBool = and(fromBool, refTable{:, 'FromAddress'} < endAddress);
                end

                toBool = strcmp(refTable{:, 'ToSection'}, sectionName);
                if (startAddress > 0)
                    toBool = and(toBool, refTable{:, 'ToAddress'} >= startAddress);
                end
                if (endAddress > 0)
                    toBool = and(toBool, refTable{:, 'ToAddress'} < endAddress);
                end

                if(inclRefOutside)
                    keepBool = or(fromBool, toBool);
                else
                    keepBool = and(fromBool, toBool);
                end
                refTable = refTable(keepBool, :);
            end

        end

        %%
        function addrNum = addrString2Number(addrString)
            %myString = "";
            if iscell(addrString)
                myString = addrString{1};
            else
                myString = addrString;
            end

            addrNum = uint64(str2num(myString));
        end

        %%
        function symbolId = findContainingSymbol(myAddress, symTable, symTableSorted)
            if (nargin < 3); symTableSorted = false; end
            symbolId = 0;
            if isempty(symTable); return; end

            if ~isnumeric(myAddress)
                if iscell(myAddress)
                    myString = myAddress{1};
                else
                    myString = myAddress;
                end

                myAddress = uint64(str2num(myString));
            end

            filtrRow = (symTable{:, 'Address'} <= myAddress);
            partTable = symTable(filtrRow,:);
            if symTableSorted
                symbolId = partTable{size(partTable,1), 'Index'};
            else
                [~, midx] = max(partTable{:, 'Address'}, [], 'all', 'omitnan');
                symbolId = partTable{midx, 'Index'};
            end
        end

        %%
        function [refAddrList, refSymList] = resolveRawRefList(myString, symTable, symTableSorted)
            if (nargin < 3); symTableSorted = false; end
            %if (nargin < 4); verboseIndex = 0; end
            refSymList = [];
            refAddrList = [];
            if isempty(symTable); return; end
            if isempty(myString); return; end
            if strcmp(myString, ""); return; end

            % if verboseIndex > 0 & (mod(verboseIndex, 100) == 0)
            %     fprintf('Working on symbol %d...\n', verboseIndex);
            % end

            splStr = split(myString, ';')';
            refAddrList = arrayfun(@(x) ProgramGraph.addrString2Number(x), splStr);
            refSymList = arrayfun(@(x) ProgramGraph.findContainingSymbol(x, symTable, symTableSorted), refAddrList);
            %refSymList = arrayfun(@(x) ProgramGraph.findContainingSymbol(x, symTable, symTableSorted), splStr);
            refSymList = unique(refSymList);
        end

        %%
        function symTable = readProgramGraph(inputPath, verbose)
            if nargin < 2; verbose = false; end
            %fmtString = '%s0x%x%d%d%d';

            if verbose
                fprintf('Reading file...\n');
            end
            if endsWith(inputPath, '.tsv')
                rawTable = readtable(inputPath, 'FileType', 'delimitedtext', 'Delimiter', '\t');
            else
                rawTable = readtable(inputPath);
            end

            if verbose
                fprintf('Processing data...\n');
            end
            symCount = size(rawTable, 1);
            varNames = {'Name', 'Address', 'AddressString', 'Index', 'Section', 'Type', 'Referees', 'RefereesAddr', 'RefereesRaw'};
            varTypes = {'string', 'uint64', 'string', 'int32', 'string', 'string', 'cell', 'cell', 'string'};
            TableSize = [symCount size(varNames, 2)];
            symTable = table(Size=TableSize, VariableTypes=varTypes, VariableNames=varNames);

            vnames = rawTable.Properties.VariableNames;
            colName = 'SYMBOL';
            if ~ismember(colName, vnames)
                colName = 'x_SYMBOL';
            end
            symTable{:, 'Name'} = rawTable{:, colName};

            colName = 'ADDRESS';
            if ~ismember(colName, vnames)
                colName = 'x_ADDRESS';
            end
            symTable{:, 'Address'} = uint64(rawTable{:, colName});
            symTable{:, 'AddressString'} = arrayfun(@(x) sprintf('%08X', x), symTable{:, 'Address'}, 'UniformOutput', false);

            colName = 'SECTION';
            if ~ismember(colName, vnames)
                colName = 'x_SECTION';
            end
            symTable{:, 'Section'} = rawTable{:, colName};

            colName = 'REFEREES';
            if ~ismember(colName, vnames)
                colName = 'x_REFEREES';
            end
            symTable{:, 'RefereesRaw'} = rawTable{:, colName};

            colName = 'TYPE';
            if ~ismember(colName, vnames)
                colName = 'x_TYPE';
            end
            symTable{:, 'Type'} = rawTable{:, colName};

            symTable = sortrows(symTable, 'Address');
            symTable{:, 'Index'} = [1:symCount]';
            clear rawTable colName vnames TableSize varTypes varNames

            %Handle referees
            if verbose
                fprintf('Resolving references...\n');
            end
            %refSymListList = arrayfun(@(x,y) ProgramGraph.resolveRawRefList(x, symTable, true,y), symTable{:, 'RefereesRaw'}, symTable{:, 'Index'}, 'UniformOutput', false);

            for s = 1:symCount
                if verbose & (mod(s, 100) == 0)
                    perccomp = (s ./ symCount) .* 100;
                    fprintf('\tSymbol %d of %d (%.2f%% complete)...\n', s, symCount, perccomp);
                end
                rawrefString = symTable{s, 'RefereesRaw'};
                [refAddrList, refList] = ProgramGraph.resolveRawRefList(rawrefString, symTable, true);
                myCell = cell(1,1);
                myCell{1} = refList;
                symTable{s, 'Referees'} = myCell;

                myCell = cell(1,1);
                myCell{1} = refAddrList;
                symTable{s, 'RefereesAddr'} = myCell;
            end
        end

        %%
        function [symTable, refTable, funcList] = loadSavedTables(filePath)
            load(filePath, 'fileVersion', 'symbolTable', 'referenceGraph');
            symTable = symbolTable; clear symbolTable;
            refTable = referenceGraph; clear referenceGraph;
            if fileVersion >= 2
                load(filePath, 'functionTableList');
                funcList = functionTableList;
            else
                funcList = [];
            end
        end

        %%
        function saveTables(filePath, symTable, refTable, funcList)
            symbolTable = symTable;
            referenceGraph = refTable;
            functionTableList = funcList;
            fileVersion = 2;
            save(filePath, 'fileVersion', 'symbolTable', 'referenceGraph', 'functionTableList');
        end

    end

end