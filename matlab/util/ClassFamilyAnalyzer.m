%
%%
classdef ClassFamilyAnalyzer

    properties
        addrLenFmt = 8;

        emptyCallAddress = uint64(0);
        emptyCallSymbolId = 0;
        minTableSize = 2;
        ftLenCutoff = 100;
        supercladeMinSize = 100;

        srcSymTable;
        srcFTTable;

        ftNodeTable;
        ftLinkTable;
        ftGraph;
        cladeTable;

        commonFuncs = [];
    end

    methods

        %%
        function obj = initialize(obj, symTable, funcList, verbose)
            addrFmtStr = ['%0' num2str(obj.addrLenFmt) 'X'];
            obj.srcSymTable = symTable;
            obj.srcFTTable = struct2table(funcList);

            if obj.emptyCallAddress > 0
                stIdx = find((symTable{:, 'Address'} == obj.emptyCallAddress), 1);
                if ~isempty(stIdx) & (stIdx > 0)
                    obj.emptyCallSymbolId = symTable{stIdx, 'Index'};
                end
                clear stIdx;
            end
            if verbose
                pstr = ['[ClassFamilyAnalyzer.initialize] Empty call set to 0x' addrFmtStr ' (SymbolId: %d)\n'];
                fprintf(pstr, obj.emptyCallAddress, obj.emptyCallSymbolId);
            end

            %Remove any tables that are super long (probably not a manually
            %written class, probably exception vector or some shit)
            %ftLens = arrayfun(@(x) size(x, 1), obj.srcFTTable{:, 'table'});
            ftCount = size(obj.srcFTTable, 1);
            ftLens = zeros(1, ftCount);
            if verbose
                fprintf('[ClassFamilyAnalyzer.initialize] %d function table candidates input...\n', ftCount);
            end
            for i = 1:ftCount
                val = obj.srcFTTable{i, 'table'};
                val = val{1};
                ftLens(i) = size(val,1);
            end
            clear val ftCount i
            
            okayBool = ftLens <= obj.ftLenCutoff;
            okayBool = and(okayBool, (ftLens >= obj.minTableSize));
            if nnz(okayBool) < 1
                return;
            end
            obj.srcFTTable = obj.srcFTTable(okayBool, :);
            clear okayBool ftLens

            obj.srcFTTable = sortrows(obj.srcFTTable, 'symbolId');
            ftCount = size(obj.srcFTTable, 1);
            ftLinkList(ftCount) = ClassFamilyAnalyzer.genFuncTableNodeStruct();
            if verbose
                fprintf('[ClassFamilyAnalyzer.initialize] %d candidates retained. Now comparing...\n', ftCount);
            end
            for i = 1:ftCount
                ftNode = ftLinkList(ftCount);
                ftNode.fttIndex = i;
                ftNode.symbolId = obj.srcFTTable{i, 'symbolId'};
                ftNode.address = obj.srcFTTable{i, 'address'};
                val = obj.srcFTTable{i, 'addressString'};
                ftNode.addressString = val{1};
                val = obj.srcFTTable{i, 'table'};
                ftNode.table = val{1};
                ftNode.len = size(ftNode.table, 1);

                if verbose
                    fprintf('[ClassFamilyAnalyzer.initialize] Working on node %d of %d (table size: %d)...\n', i, ftCount, size(ftNode.table, 1));
                end
                for j = 1:(i-1)
                    %Scan matches for all in behind...
                    otherNode = ftLinkList(j);
                    [ftNode, otherNode] = obj.compareNodes(ftNode, otherNode);
                    ftLinkList(j) = otherNode;
                end
                %Save ftNode links to master link list
                if ~isempty(ftNode.linkTable)
                    nlt = ftNode.linkTable;
                    nlt{:, 'pfttIndex'} = i;
                    nlt{:, 'psymbolId'} = ftNode.symbolId;
                    nlt{:, 'paddressString'} = ftNode.addressString;
                    nlt{:, 'pfullLen'} = size(ftNode.table, 1);

                    if isempty(obj.ftLinkTable)
                        obj.ftLinkTable = nlt;
                    else
                        obj.ftLinkTable = [obj.ftLinkTable; nlt];
                    end
                    clear nlt
                end

                ftLinkList(i) = ftNode;
            end

            obj.ftNodeTable = struct2table(ftLinkList);

            obj.ftLinkTable = sortrows(obj.ftLinkTable, 'fullMatches', 'descend');
        end

        %%
        function [nodeA, nodeB] = compareNodes(obj, nodeA, nodeB)
            sizeA = size(nodeA.table, 1);
            sizeB = size(nodeB.table, 1);

            %Must be at least one full match
            commonLen = min(sizeA, sizeB);
            fullMatchCount = 0;
            nullMatchCount = 0;
            for i = 1:commonLen
                entryA = nodeA.table{i, 'SymbolId'};
                entryB = nodeB.table{i, 'SymbolId'};
                if (entryA == entryB)
                    if (entryA == 0) | (entryA == obj.emptyCallSymbolId)
                        nullMatchCount = nullMatchCount + 1;
                    else
                        fullMatchCount = fullMatchCount + 1;
                    end
                end
            end
            clear entryA entryB i

            if fullMatchCount < 1
                return;
            end

            linkAB = ClassFamilyAnalyzer.genLinkStructFromReferencedNode(nodeB);
            linkBA = ClassFamilyAnalyzer.genLinkStructFromReferencedNode(nodeA);
            linkAB.commonLen = commonLen; linkBA.commonLen = commonLen;
            linkAB.fullMatches = fullMatchCount; linkBA.fullMatches = fullMatchCount;
            linkAB.nullMatches = nullMatchCount; linkBA.nullMatches = nullMatchCount;

            %"Null override" matches
            apOverrideCount = 0; %Treating A as parent
            bpOverrideCount = 0; %Treating B as parent
            for i = 1:commonLen
                entryA = nodeA.table{i, 'SymbolId'};
                entryB = nodeB.table{i, 'SymbolId'};
                if (entryA == 0) | (entryA == obj.emptyCallSymbolId)
                    if (entryB ~= 0) & (entryB ~= obj.emptyCallSymbolId)
                        apOverrideCount = apOverrideCount + 1;
                    end
                end

                if (entryB == 0) | (entryB == obj.emptyCallSymbolId)
                    if (entryA ~= 0) & (entryA ~= obj.emptyCallSymbolId)
                        bpOverrideCount = bpOverrideCount + 1;
                    end
                end
            end
            clear entryA entryB i

            linkAB.pcNullOverrideMatches = apOverrideCount; linkBA.cpNullOverrideMatches = apOverrideCount;
            linkAB.cpNullOverrideMatches = bpOverrideCount; linkBA.pcNullOverrideMatches = bpOverrideCount;

            %Add to node link tables
            nodeA.linkTable = ClassFamilyAnalyzer.addEntryToLinkTable(nodeA.linkTable, linkAB);
            nodeB.linkTable = ClassFamilyAnalyzer.addEntryToLinkTable(nodeB.linkTable, linkBA);
        end

        %%
        function [pcPossible, cpPossible] = flagPossibleLinkDirections(obj)
            % pSmaller = obj.ftLinkTable{:, 'pfullLen'} < obj.ftLinkTable{:, 'fullLen'};
            % cSmaller = obj.ftLinkTable{:, 'fullLen'} < obj.ftLinkTable{:, 'pfullLen'};
            % eqSize = obj.ftLinkTable{:, 'fullLen'} == obj.ftLinkTable{:, 'pfullLen'};
            % hasPCOvr = obj.ftLinkTable{:, 'pcNullOverrideMatches'} > 0;
            % hasCPOvr = obj.ftLinkTable{:, 'cpNullOverrideMatches'} > 0;
            % 
            % pcPossible = or(pSmaller, eqSize);
            % pcPossible = and(pcPossible, ~cSmaller);
            % pcPossible = and(pcPossible, ~hasCPOvr);
            % 
            % cpPossible = or(cSmaller, eqSize);
            % cpPossible = and(cpPossible, ~pSmaller);
            % cpPossible = and(cpPossible, ~hasPCOvr);

            [pcPossible, cpPossible] = ClassFamilyAnalyzer.flagPossibleLinkDirectionsStatic(obj.ftLinkTable);
        end

        %%
        function [obj, myGraph] = generateGraph(obj)
            [pcPossible, cpPossible] = obj.flagPossibleLinkDirections();
            
            pcSrcNodes = obj.ftLinkTable{pcPossible, 'pfttIndex'}';
            pcTrgNodes = obj.ftLinkTable{pcPossible, 'fttIndex'}';
            pcWeights = obj.ftLinkTable{pcPossible, 'fullMatches'}';

            cpSrcNodes = obj.ftLinkTable{cpPossible, 'fttIndex'}';
            cpTrgNodes = obj.ftLinkTable{cpPossible, 'pfttIndex'}';
            cpWeights = obj.ftLinkTable{cpPossible, 'fullMatches'}';

            nodeNames = obj.ftNodeTable{:, 'addressString'};
            clear pcPossible cpPossible

            srcNodes = [pcSrcNodes cpSrcNodes];
            trgNodes = [pcTrgNodes cpTrgNodes];
            eWeights = [pcWeights cpWeights];
            %eWeights(:) = 1;

            % testAddr = uint64(0x10b52d34);
            % ftIdx = find(obj.ftNodeTable{:, 'address'} == testAddr, 1);
            % insrc = find(srcNodes == ftIdx);
            % intrg = find(trgNodes == ftIdx);
            % srcPrt = trgNodes(insrc);
            % trgPrt = srcNodes(intrg);

            myGraph = digraph(srcNodes, trgNodes, eWeights, nodeNames);
            clear pcSrcNodes pcTrgNodes pcWeights
            clear cpSrcNodes cpTrgNodes cpWeights

            obj.ftGraph = myGraph;
        end

        %%
        function figHandle = renderGraph(obj, figno)
            if isempty(obj.ftGraph)
                [obj, ~] = obj.generateGraph();
            end

            maxMatch = max(obj.ftLinkTable{:, 'fullMatches'}, [], 'omitnan');
            maxDrawWeight = 5.0;
            maxNodeSize = 5.0;

            nodeSizes = obj.ftNodeTable{:, 'len'} ./ maxNodeSize;
            drawWeights = (obj.ftGraph.Edges{:, 'Weight'} ./ maxMatch) .* maxDrawWeight;

            figHandle = figure(figno);
            plot(obj.ftGraph, 'LineWidth', drawWeights, 'NodeColor', 'r', 'EdgeColor', 'b');
        end

        %%
        function obj = assignClades(obj)
            if isempty(obj.ftNodeTable); return; end
            obj.ftNodeTable{:, 'clade'} = 0;
            obj.ftNodeTable{:, 'subclade'} = 0;
            if isempty(obj.ftGraph)
                [obj, ~] = obj.generateGraph();
            end
            cclades = conncomp(obj.ftGraph, 'Type', 'weak')';
            obj.ftNodeTable{:, 'clade'} = cclades;

            cladeCount = max(cclades, [], 'all', 'omitnan');
            %Clade table
            [varNames, varTypes] = ClassFamilyAnalyzer.getCladeTableColumns();
            tblSz = [cladeCount size(varNames, 2)];
            obj.cladeTable = table('Size', tblSz, 'VariableTypes', varTypes, 'VariableNames', varNames);

            %Table and Break up large clades
            for c = 1:cladeCount
                obj.cladeTable{c, 'CladeId'} = c;

                isMember = (obj.ftNodeTable{:, 'clade'} == c);
                memberCount = nnz(isMember);
                obj.cladeTable{c, 'MemberCount'} = memberCount;
                if memberCount >= obj.supercladeMinSize
                    sg = subgraph(obj.ftGraph, find(isMember));
                    scclades = conncomp(sg)';
                    obj.ftNodeTable{isMember, 'subclade'} = scclades;
                    obj.cladeTable{c, 'SubcladeCount'} = max(scclades, [], 'all', 'omitnan');
                else
                    obj.cladeTable{c, 'SubcladeCount'} = 0;
                end

                if memberCount > 1
                    obj.cladeTable{c, 'Autoname'} = "";
                else
                    memberName = obj.ftNodeTable{isMember, 'addressString'};
                    if iscell(memberName)
                        memberName = memberName{1};
                    end
                    obj.cladeTable{c, 'Autoname'} = string(['C' memberName]);
                end
            end

            %Name clades
            an = 1;
            gl = ClassFamilyAnalyzer.getGreekLetterNames();
            obj.cladeTable = sortrows(obj.cladeTable, 'MemberCount', 'descend');
            for c = 1:cladeCount
                cladeId = obj.cladeTable{c, 'CladeId'};
                isMember = (obj.ftNodeTable{:, 'clade'} == cladeId);
                memberCount = obj.cladeTable{c, 'MemberCount'};
                if memberCount > 1
                    usestr = ClassFamilyAnalyzer.glNameGen(gl, an);
                    obj.cladeTable{c, 'Autoname'} = string(['C' usestr]);
                    an = an + 1;
                else
                    memberName = obj.ftNodeTable{isMember, 'addressString'};
                    if iscell(memberName)
                        memberName = memberName{1};
                    end
                    obj.cladeTable{c, 'Autoname'} = string(['C' memberName]);
                end
            end
            obj.cladeTable = sortrows(obj.cladeTable, 'CladeId', 'ascend');
        end

        %%
        function [obj, figHandle] = visualizeClade(obj, cladeId, figno)
            figHandle = figure(figno);
            clf;

            if isempty(obj.ftNodeTable); return; end
            if isempty(obj.ftGraph)
                [obj, ~] = obj.generateGraph();
            end

            cmembers = find(obj.ftNodeTable{:, 'clade'} == cladeId);
            if isempty(cmembers); return; end
            mCount = nnz(cmembers);
            if mCount < 1; return; end

            %Isolate graph...
            sg = subgraph(obj.ftGraph, cmembers);
            lineWidths = sg.Edges{:, 'Weight'};
            lineWidths(lineWidths > 10) = 10;

            %Draw figure...
            plot(sg, 'LineWidth', lineWidths, 'NodeColor', 'r', 'EdgeColor', 'b');
            title(['Clade ' num2str(cladeId) ' (Size = ' num2str(mCount) ')']);
        end

        %%
        function figHandles = visualizeClades(obj, initFigNo, minCladeSize)
            if isempty(obj.ftNodeTable); return; end
            if ~tableHasField(obj.ftNodeTable, 'clade')
                obj = obj.assignClades();
            end

            if minCladeSize < 2; minCladeSize = 2; end

            figHandles = [];
            figNo = initFigNo;
            cladeCount = max(obj.ftNodeTable{:, 'clade'}, [], 'all', 'omitnan');
            for c = 1:cladeCount
                cmembers = find(obj.ftNodeTable{:, 'clade'} == c);
                if isempty(cmembers); continue; end
                mCount = nnz(cmembers);
                if mCount < minCladeSize; continue; end

                %Isolate graph...
                sg = subgraph(obj.ftGraph, cmembers);
                lineWidths = sg.Edges{:, 'Weight'};
                lineWidths(lineWidths > 10) = 10;

                %Draw figure...
                fh = figure(figNo);
                clf;
                plot(sg, 'LineWidth', lineWidths, 'NodeColor', 'r', 'EdgeColor', 'b');
                title(['Clade ' num2str(c) ' (Size = ' num2str(mCount) ')']);
                clear sg cmembers mCount

                figHandles = [figHandles fh];
                figNo = figNo + 1;
            end
        end

        %%
        function obj = updateCommonFuncs(obj)
            %Functions foun in more than one clade.
            if ~tableHasField(obj.ftNodeTable, 'clade')
                obj = obj.assignClades();
            end

            symCount = size(obj.srcSymTable, 1);
            foundInClades = zeros(1, symCount);
            cladeCount = max(obj.ftNodeTable{:, 'clade'}, [], 'all', 'omitnan');
            for c = 1:cladeCount
                fpool = [];
                isMember = (obj.ftNodeTable{:, 'clade'} == c);
                memberTable = obj.ftNodeTable(isMember, :);
                memberCount = size(memberTable, 1);
                for m = 1:memberCount
                    tt = memberTable{m, 'table'};
                    if iscell(tt)
                        tt = tt{1};
                    end
                    if isempty(fpool)
                        fpool = tt{:, 'SymbolId'};
                    else
                        fpool = [fpool; tt{:, 'SymbolId'}];
                    end
                end
                fpool = unique(fpool);
                %Remove 0
                fpool = fpool(fpool > 0);

                foundInClades(fpool) = foundInClades(fpool) + 1;
            end

            obj.commonFuncs = find(foundInClades > 1);
        end

        %%
        function obj = analyzeClade(obj, c)
            %TODO
            isMember = (obj.ftNodeTable{:, 'clade'} == c);
            mCount = nnz(isMember);
            if (mCount > 1) & (mCount < obj.supercladeMinSize)
                obj = obj.refineSmallClade(c);
            end

        end
        
        %%
        function obj = refineSmallClade(obj, c)
            if ~tableHasField(obj.ftNodeTable, 'clade')
                obj = obj.assignClades();
            end
            obj = obj.updateCommonFuncs();

            isMember = (obj.ftNodeTable{:, 'clade'} == c);
            %memberTable = obj.ftNodeTable(isMember, :);
            %mCount = size(memberTable, 1);

            %Cull weak edges (low full match, common functions)
            %And group these into subclades
            cmembers = find(obj.ftNodeTable{:, 'clade'} == c);
            sg = subgraph(obj.ftGraph, cmembers);
            eWeights = sg.Edges{:, 'Weight'};
            weightstd = std(eWeights, 0, 'all', 'omitnan');
            weightmean = mean(eWeights, 'all', 'omitnan');
            minWeight = weightmean - (2 * weightstd);
            weakEdges = eWeights < minWeight;
            %sg.Edges = sg.Edges(~weakEdges, :);
            sg = rmedge(sg, find(weakEdges));
            sclades = conncomp(sg, 'Type', 'weak')';

            [~, ftidx] = ismember(sg.Nodes{:, 'Name'}, obj.ftNodeTable{:, 'addressString'});
            obj.ftNodeTable{ftidx, 'subclade'} = sclades;

            obj.cladeTable{c, 'SubcladeCount'} = max(obj.ftNodeTable{isMember, 'subclade'}, [], 'all', 'omitnan');
        end

        %%
        function obj = updateCladeMemberFunctionTerritories(obj, c, subc)
            %TODO!!!!!
            if ~tableHasField(obj.ftNodeTable, 'clade')
                obj = obj.assignClades();
            end
            if ~tableHasField(obj.ftNodeTable, 'territoryStart')
                obj.ftNodeTable{:, 'territoryStart'} = int32(0);
                obj.ftNodeTable{:, 'territoryEnd'} = int32(0);
            end
            obj = obj.updateCommonFuncs();

            isMember = (obj.ftNodeTable{:, 'clade'} == c);
            if subc > 0
                isMember = and(isMember, obj.ftNodeTable{:, 'subclade'} == subc);
            end
            memberTable = obj.ftNodeTable(isMember, :);
            mCount = size(memberTable, 1);
            %Pool member functions.
            fpool = [];
            for m =  1:mCount
                tt = memberTable{m, 'table'};
                if iscell(tt)
                    tt = tt{1};
                end
                fpool = [fpool tt{:, 'SymbolId'}'];
            end
            fpool = unique(fpool);
            clear m tt

            %Remove null and common functions
            isNull = (fpool == 0);
            isCommon = ismember(fpool, obj.commonFuncs);
            isOkay = and(~isNull, ~isCommon);
            if nnz(isOkay) < 1
                %Nothing to work with
                return;
            end
            fpool = fpool(isOkay);
            clear isNull isCommon isOkay

            if (mCount > 1)
                %Build a check matrix
                poolSize = size(fpool, 2);
                checkMtx = false(mCount, poolSize);
                for m =  1:mCount
                    tt = memberTable{m, 'table'};
                    if iscell(tt)
                        tt = tt{1};
                    end
                    [~, poolIdx] = ismember(tt{:, 'SymbolId'}, fpool);
                    poolIdx = poolIdx(poolIdx > 0);
                    checkMtx(m, poolIdx) = true;
                end
                fprintf('DEBUG HOLD');
            else
                obj.ftNodeTable{isMember, 'territoryStart'} = int32(min(fpool, [], 'all', 'omitnan'));
                obj.ftNodeTable{isMember, 'territoryEnd'} = int32(max(fpool, [], 'all', 'omitnan'));
            end

        end

        %%
        function [obj, funcSymIdx, symFreq] = getCladeConsensusFunctionTable(obj, c, subc)
            if (nargin < 3); subc = 0; end
            if ~tableHasField(obj.ftNodeTable, 'clade')
                obj = obj.assignClades();
            end

            inClade = (obj.ftNodeTable{:, 'clade'} == c);
            if subc > 0
                inClade = and(inClade, (obj.ftNodeTable{:, 'subclade'} == subc));
            end
            mCount = nnz(inClade);

            if mCount < 1
                funcSymIdx = [];
                return;
            end

            cladeNodes = obj.ftNodeTable(inClade, :);
            fMaxCount = max(cladeNodes{:, 'len'});
            %funcSymIdx = zeros(1, fMaxCount);
            allmtx = zeros(mCount, fMaxCount);
            for m = 1:mCount
                mtbl = cladeNodes{m, 'table'};
                if isempty(mtbl)
                    continue;
                end
                if iscell(mtbl)
                    mtbl = mtbl{1};
                    if isempty(mtbl)
                        continue;
                    end
                end
                ml = size(mtbl, 1);
                allmtx(m, 1:ml) = mtbl{:, 'SymbolId'};
            end
            clear m mtbl ml

            %Find mode of each column and number of times it appears
            [modeM, modeF] = mode(allmtx, 1);
            modeM(modeF < 2) = 0;
            modeF(modeF < 2) = 0;

            %Trim output table, if needed
            passBool = (modeM > 0);
            lastIdx = find(passBool, 1, 'last');
            funcSymIdx = modeM(1:lastIdx);
            symFreq = modeF(1:lastIdx);
        end
        
    end

    methods (Static)

        %%
        function greekLetterNames = getGreekLetterNames()
            greekLetterNames = {'Alpha' 'Beta' 'Gamma' 'Delta' ...
                'Epsilon' 'Zeta' 'Eta' 'Theta' ...
                'Iota' 'Kappa' 'Lambda' 'Mu' ...
                'Nu' 'Xi' 'Omicron' 'Pi' ...
                'Sigma' 'Tau' 'Upsilon' 'Phi' ...
                'Chi' 'Psi' 'Omega'};
        end

        %%
        function name = glNameGen(gl, i)
            glCount = size(gl, 2);
            if (i <= glCount)
                name = gl{i};
                return;
            end

            ii = i - 1;
            digits = zeros(1, 16);
            digitsUsed = 0;
            while ii > 0
                d = mod(ii, glCount);
                di = digitsUsed + 1;
                digits(di) = d+1;
                digitsUsed = di;
                ii = floor(ii ./ glCount);
            end

            name = [];
            for j = digitsUsed:-1:1
                name = [name gl{digits(j)}];
            end

        end

        %%
        function [varNames, varTypes] = getCladeTableColumns()
            varNames = {'CladeId' 'SubcladeCount' 'MemberCount' 'Autoname'};
            varTypes = {'int32' 'int32' 'int32' 'string'};
        end

        %%
        function ftNodeStruct = genFuncTableNodeStruct()
            ftNodeStruct = struct();
            ftNodeStruct.fttIndex = 0;
            ftNodeStruct.symbolId = 0;
            ftNodeStruct.address = uint64(0);
            ftNodeStruct.addressString = "";
            ftNodeStruct.len = 0;
            ftNodeStruct.table = table.empty();
            ftNodeStruct.linkTable = table.empty();
        end

        %%
        function linkStruct = genLinkStruct()
            linkStruct = struct();
            linkStruct.fttIndex = 0;
            linkStruct.symbolId = 0;
            linkStruct.addressString = ""; %For easy debug inspection
            linkStruct.commonLen = 0;
            linkStruct.fullLen = 0;
            linkStruct.fullMatches = 0;
            linkStruct.nullMatches = 0;
            linkStruct.pcNullOverrideMatches = 0; %"this" as parent, other as child
            linkStruct.cpNullOverrideMatches = 0; %"this" as child, other as parent
        end

        %%
        function linkStruct = initLinkStruct(linkStruct, referenceNode)
            linkStruct.fttIndex = referenceNode.fttIndex;
            linkStruct.symbolId = referenceNode.symbolId;
            linkStruct.addressString = referenceNode.addressString;
            linkStruct.fullLen = size(referenceNode.table, 1);
        end

        %%
        function linkStruct = genLinkStructFromReferencedNode(referenceNode)
            linkStruct = ClassFamilyAnalyzer.genLinkStruct();
            linkStruct = ClassFamilyAnalyzer.initLinkStruct(linkStruct, referenceNode);
        end

        %%
        function linkTable = addEntryToLinkTable(linkTable, linkNode)
            appTable = struct2table(linkNode);
            if isempty(linkTable)
                linkTable = appTable;
            else
                linkTable = [linkTable; appTable];
            end
        end

        %%
        function [pcPossible, cpPossible] = flagPossibleLinkDirectionsStatic(linkTable)
            pSmaller = linkTable{:, 'pfullLen'} < linkTable{:, 'fullLen'};
            cSmaller = linkTable{:, 'fullLen'} < linkTable{:, 'pfullLen'};
            eqSize = linkTable{:, 'fullLen'} == linkTable{:, 'pfullLen'};
            hasPCOvr = linkTable{:, 'pcNullOverrideMatches'} > 0;
            hasCPOvr = linkTable{:, 'cpNullOverrideMatches'} > 0;

            pcPossible = or(pSmaller, eqSize);
            pcPossible = and(pcPossible, ~cSmaller);
            pcPossible = and(pcPossible, ~hasCPOvr);

            cpPossible = or(cSmaller, eqSize);
            cpPossible = and(cpPossible, ~pSmaller);
            cpPossible = and(cpPossible, ~hasPCOvr);
        end

        %%
        function nodeArray = flagPossibleLinkDirections_NodeArray(nodeArray)
            if isempty(nodeArray)
                return;
            end

            %TODO Vectorize?
            nCount = size(nodeArray, 2);
            for i = 1:nCount
                nodeStruct = nodeArray(i);
                if ~isempty(nodeStruct.linkTable)
                    ltnew = ClassFamilyAnalyzer.flagPossibleLinkDirections_SingleNode(nodeStruct);
                    nodeStruct.linkTable = ltnew;
                end
                nodeArray(i) = nodeStruct;
            end
        end

        %%
        function nodeLinks = flagPossibleLinkDirections_SingleNode(nodeStruct)
            nodeLinks = nodeStruct.linkTable;
            if isempty(nodeLinks)
                return;
            end

            nodeLinksCopy = nodeLinks;
            nodeLinksCopy{:, 'pfullLen'} = nodeStruct.len;
            [pcPossible, cpPossible] = ClassFamilyAnalyzer.flagPossibleLinkDirectionsStatic(nodeLinksCopy);
            nodeLinks{:, 'pcFlag'} = pcPossible;
            nodeLinks{:, 'cpFlag'} = cpPossible;
        end

    end
end