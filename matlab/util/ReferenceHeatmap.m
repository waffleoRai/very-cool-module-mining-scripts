%
%%
classdef ReferenceHeatmap

    properties
        addrFmtChars = 8;

        binShift = 0;
        binSize = 64;
        isSymmetric = false;
        fromRegion;
        toRegion;

        srcRefTable;
        adjSymTable;

        rawDensityMap;
        xEdges; %Local indices
        yEdges;

        %For symmetric square maps
        culledDensityMap;
        xEdgesCulled;
        yEdgesCulled;
        culledColBool; %1 if tossed, 0 if kept

        pearsonDistMap;

        cnormRadMin = 3;
        cnormRadMax = 40;
    end

    methods

        %%
        function obj = updateRawDensityMap(obj)
            addrFmtString = ['%0' num2str(obj.addrFmtChars) 'X'];
            lookupTableSymId = [obj.srcRefTable{:, 'FromSymbol'}' obj.srcRefTable{:, 'ToSymbol'}'];
            lookupTableAddr = [obj.srcRefTable{:, 'FromAddress'}' obj.srcRefTable{:, 'ToAddress'}'];

            allSyms = obj.srcRefTable{:, 'FromSymbol'}';
            allSyms = [allSyms obj.srcRefTable{:, 'ToSymbol'}'];
            allSyms = unique(allSyms);
            symCount = size(allSyms, 2);

            varNames = {'Address' 'SymbolId' 'AddressString'};
            varTypes = {'uint64' 'int32' 'string'};
            TableSize = [symCount size(varNames, 2)];
            obj.adjSymTable = table(Size=TableSize, VariableTypes=varTypes, VariableNames=varNames);

            [~, lookupIdx] = ismember(allSyms, lookupTableSymId);
            obj.adjSymTable{:, 'SymbolId'} = allSyms';
            obj.adjSymTable{:, 'Address'} = lookupTableAddr(lookupIdx)';
            obj.adjSymTable{:, 'AddressString'} = arrayfun(@(x) sprintf(addrFmtString, x), obj.adjSymTable{:, 'Address'}, 'UniformOutput', false);

            obj.adjSymTable = sortrows(obj.adjSymTable, 'Address');
            clear lookupTableSymId lookupTableAddr lookupIdx TableSize varNames varTypes allSyms symCount

            fromVec = obj.srcRefTable{:, 'FromSymbol'}';
            toVec = obj.srcRefTable{:, 'ToSymbol'}';

            %Reindex
            [~, fromVec] = ismember(fromVec, obj.adjSymTable{:, 'SymbolId'});
            [~, toVec] = ismember(toVec, obj.adjSymTable{:, 'SymbolId'});

            %Determine bins
            if obj.fromRegion.startAddress > 0
                fromMin = find(obj.adjSymTable{:, 'Address'} >= obj.fromRegion.startAddress, 1);
            else
                fromMin = min(fromVec, [], 'all', 'omitnan');
            end
            if obj.fromRegion.endAddress > 0
                fromMax = find(obj.adjSymTable{:, 'Address'} < obj.fromRegion.endAddress, 1, 'last');
            else
                fromMax = max(fromVec, [], 'all', 'omitnan');
            end
            if obj.toRegion.startAddress > 0
                toMin = find(obj.adjSymTable{:, 'Address'} >= obj.toRegion.startAddress, 1);
            else
                toMin = min(toVec, [], 'all', 'omitnan');
            end
            if obj.toRegion.endAddress > 0
                toMax = find(obj.adjSymTable{:, 'Address'} < obj.toRegion.endAddress, 1, 'last');
            else
                toMax = max(toVec, [], 'all', 'omitnan');
            end

            %Align to bins
            fromMin = alignDown(fromMin + obj.binShift, obj.binSize) + 1; 
            fromMax = alignDown(fromMax + obj.binSize + obj.binShift, obj.binSize) + 1; 
            toMin = alignDown(toMin + obj.binShift, obj.binSize) + 1; 
            toMax = alignDown(toMax + obj.binSize + obj.binShift, obj.binSize) + 1; 

            %Mirror if symmetric
            if obj.isSymmetric
                fromVecOld = fromVec;
                toVecOld = toVec;
                fromVec = [fromVecOld toVecOld];
                toVec = [toVecOld fromVecOld];

                %And clear out anything outside the trim
                fromBad = or((fromVec < fromMin), (fromVec > fromMax));
                toBad = or((toVec < toMin), (toVec > toMax));
                refBad = or(fromBad, toBad);
                fromVec = fromVec(~refBad);
                toVec = toVec(~refBad);

                clear fromBad toBad refBad
            end

            %Histogram
            obj.xEdges = fromMin:obj.binSize:fromMax;
            obj.yEdges = toMin:obj.binSize:toMax;
            [obj.rawDensityMap, obj.xEdges, obj.yEdges] = histcounts2(fromVec,toVec,obj.xEdges,obj.yEdges);
            
        end
    
        %%
        function obj = prepForClusterAnalysis(obj, maxSymRad, cleanRad, cleanThresh)
            if (nargin < 2); maxSymRad = 5000; end
            if (nargin < 3); cleanRad = 10; end
            if (nargin < 4); cleanThresh = 0.1; end

            N = size(obj.rawDensityMap, 1);
            [xx,yy] = meshgrid(1:N, 1:N);
            xydiff = xx - yy;

            %Cull low info columns
            keepvec = true(1,N);
            if ~isnan(cleanRad)
                checkCleanBool = and((xydiff < cleanRad), xydiff >= 0);
                for i = 1:N
                    checkMtx = and(checkCleanBool, (xx == i));
                    mm = mean(obj.rawDensityMap(checkMtx) > 0, 'all', 'omitnan');
                    if (mm < cleanThresh)
                        %Remove col and row
                        keepvec(i) = false;
                    end
                end

                obj.culledDensityMap = obj.rawDensityMap(keepvec, keepvec);
                keepvec = [keepvec true];
                obj.xEdgesCulled = obj.xEdges(keepvec);
                obj.yEdgesCulled = obj.yEdges(keepvec);

                N = size(obj.culledDensityMap, 1);
                [xx,yy] = meshgrid(1:N, 1:N);
                xydiff = xx - yy;
            else
                obj.culledDensityMap = obj.rawDensityMap;
                obj.xEdgesCulled = obj.xEdges;
                obj.yEdgesCulled = obj.yEdges;
            end

            obj.culledColBool = ~keepvec;

            %Pearson distance
            maxBinRad = floor(maxSymRad ./ obj.binSize);
            if obj.cnormRadMax > maxBinRad
                obj.cnormRadMax = maxBinRad - 1;
            end
            xydiff_a = abs(xydiff);
            obj.culledDensityMap(xydiff_a > maxBinRad) = NaN;

            covMtx = NaN(N,N);
            stdPMtx = NaN(N,N);
            for i = 1:N
                covMtx(i,i) = 1;
                stdPMtx(i,i) = 1;
                for j = (i+1):N
                    if isnan(obj.culledDensityMap(i,j))
                        %No overlap
                        break;
                    end
                    colA = obj.culledDensityMap(:,i);
                    colB = obj.culledDensityMap(:,j);
                    overlapReg = and(~isnan(colA), ~isnan(colB));
                    colA = colA(overlapReg);
                    colB = colB(overlapReg);

                    covRaw = cov(colA, colB, 0, 'omitrows');
                    %covRaw = cov(densityMap(:,i), densityMap(:,j), 0, 'omitrows');
                    covMtx(i,j) = covRaw(1,2);
                    covMtx(j,i) = covRaw(1,2);

                    stdPMtx(i,j) = std(colA,0,'all','omitnan') .* std(colB,0,'all','omitnan');
                    %stdPMtx(i,j) = colStdev(i) .* colStdev(j);
                    stdPMtx(j,i) = stdPMtx(i,j);
                end
            end
            clear colA colB overlapReg covRaw i j

            distMtx = 1.0 - abs(covMtx ./ stdPMtx);
            %Remove any rows with nans within CNORM_RAD_MAX of diag
            badBin = and(isnan(distMtx), (xydiff_a <= obj.cnormRadMax));
            badCol = sum(badBin, 1) > 0;
            if nnz(badCol) > 0
                keepvec = ~badCol;
                distMtx = distMtx(keepvec, keepvec);
                obj.culledDensityMap = obj.culledDensityMap(keepvec, keepvec);
                keepvec = [keepvec true];
                obj.xEdgesCulled = obj.xEdgesCulled(keepvec);
                obj.yEdgesCulled = obj.yEdgesCulled(keepvec);

                %Update overall culled flags
                prevKeptMap = find(~obj.culledColBool);
                obj.culledColBool(prevKeptMap(badCol)) = true;
            end
            obj.pearsonDistMap = distMtx;
        end

        %%
        function figHandle = renderRawDensityMap(obj, figNo)
            xbins = size(obj.xEdges, 2);
            ybins = size(obj.yEdges, 2);
            xlbl = obj.adjSymTable{obj.xEdges(1:(xbins-1)), 'AddressString'};
            ylbl = obj.adjSymTable{obj.yEdges(1:(ybins-1)), 'AddressString'};

            zmin = 0;
            zmax = max(obj.rawDensityMap, [], 'all', 'omitnan');

            figHandle = figure(figNo);
            clf;
            hm = heatmap(xlbl, ylbl, obj.rawDensityMap);
            hm.Colormap = turbo;
            hm.ColorLimits = [zmin zmax];
            hm.CellLabelColor = 'none';
            hm.GridVisible = 'off';
            title('Raw Reference Count');
        end

        %%
        function figHandle = renderPearsonDistanceMap(obj, figNo)
            figHandle = [];
            if isempty(obj.pearsonDistMap); return; end

            xbins = size(obj.xEdgesCulled, 2);
            ybins = size(obj.yEdgesCulled, 2);
            xlbl = obj.adjSymTable{obj.xEdgesCulled(1:(xbins-1)), 'AddressString'};
            ylbl = obj.adjSymTable{obj.yEdgesCulled(1:(ybins-1)), 'AddressString'};

            figHandle = figure(figNo);
            clf;
            hm = heatmap(xlbl, ylbl, obj.pearsonDistMap);
            hm.Colormap = turbo;
            hm.ColorLimits = [0 1];
            hm.CellLabelColor = 'none';
            hm.GridVisible = 'off';
            title('Pearson Distance');
        end
    
        %%
        function [startSymbol, endSymbol] = getCoverageOfXBins(obj, startBin, endBin, culledBool)
            startSymbol = 0;
            endSymbol = 0;
            if culledBool
                if (startBin > 0) & (startBin < size(obj.xEdges, 2))
                    startIdx = obj.xEdgesCulled(startBin);
                    startSymbol = obj.adjSymTable{startIdx, 'SymbolId'};
                end
                if (endBin > 0) & (endBin < size(obj.xEdges, 2))
                    endIdx = obj.xEdgesCulled(endBin);
                    endSymbol = obj.adjSymTable{endIdx, 'SymbolId'};
                end
            else
                if (startBin > 0) & (startBin < size(obj.xEdges, 2))
                    startIdx = obj.xEdges(startBin);
                    startSymbol = obj.adjSymTable{startIdx, 'SymbolId'};
                end
                if (endBin > 0) & (endBin < size(obj.xEdges, 2))
                    endIdx = obj.xEdges(endBin);
                    endSymbol = obj.adjSymTable{endIdx, 'SymbolId'};
                end
            end
        end

    end

    methods (Static)
        
        %%
        function refDensityMap = genHeatmap(refTable, binSize, binShift, asym, fromRegion, toRegion)
            if (nargin < 6); toRegion = fromRegion; end

            refDensityMap = ReferenceHeatmap;
            refDensityMap.binSize = binSize;
            refDensityMap.binShift = binShift;
            refDensityMap.isSymmetric = ~asym;

            refTable = ProgramGraph.filterFROMReferences(refTable, fromRegion);
            refTable = ProgramGraph.filterTOReferences(refTable, toRegion);
            refDensityMap.srcRefTable = refTable;

            refDensityMap.fromRegion = fromRegion;
            refDensityMap.toRegion = toRegion;
            refDensityMap = refDensityMap.updateRawDensityMap();
        end
        
    end

end