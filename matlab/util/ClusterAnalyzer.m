%
%%
classdef ClusterAnalyzer

    properties
        srcRefTable;
        regionSpec;

        minBinSizePwr = 3;
        maxBinSizePwr = 5;
        shiftsPerBin = 4;

        seedClusterSize = 4;
        maxEvalRadSyms = 5000;

        cleanRad = 10;
        cleanThresh = 0.1;

        meanGroupingFuncs = ClusterAnalyzer.genDefaultMeanGroupingFunctions();
        stdGroupingFuncs = ClusterAnalyzer.genDefaultStdGroupingFunctions();
        scoreWeights = [0.0 0.25 0.5 0.75 1.0];
        
    end

    methods

        %%
        function [symList, boundScores] = evaluate(obj, refTable, regionSpecs, verbose)

            %Filter ref table
            if ~isempty(regionSpecs)
                refTable = ProgramGraph.filterFROMReferences(refTable, regionSpecs);
                refTable = ProgramGraph.filterTOReferences(refTable, regionSpecs);
            end
            obj.srcRefTable = refTable;
            obj.regionSpec = regionSpecs;

            %Get symbol list
            symList = unique([refTable{:, 'FromSymbol'}' refTable{:, 'ToSymbol'}']);
            totalSymCount = size(symList, 2);
            boundScores = NaN(1, totalSymCount);
            runSums = zeros(1, totalSymCount);
            runTallies = zeros(1, totalSymCount);

            %Try different bin sizes and shifts
            rng('default');
            for binPwr = obj.maxBinSizePwr : -1 : obj.minBinSizePwr
                binSize = 2 ^ binPwr;
                binShift = 0;
                binShiftDone = false(1, binSize);
                binShiftDone(1) = true;

                for shroll = 1:obj.shiftsPerBin
                    if(verbose)
                        fprintf('[ClusterAnalyzer.evaluate] Now trying bin size %d with shift of %d...\n', binSize, binShift);
                    end

                    hmobj = ReferenceHeatmap.genHeatmap(obj.srcRefTable, binSize, binShift, false, obj.regionSpec);
                    hmobj = hmobj.prepForClusterAnalysis(obj.maxEvalRadSyms, obj.cleanRad, obj.cleanThresh);

                    roundScores = obj.clusterHeatmap(hmobj);
                    binCount = size(roundScores, 2);
                    for b = 1:binCount
                        [startSym, endSym] = hmobj.getCoverageOfXBins(b, b+1, true);
                        stidx = find(symList == startSym, 1);
                        edidx = find(symList == endSym, 1);
                        runSums(stidx : edidx) = runSums(stidx : edidx) + roundScores(b);
                        runTallies(stidx : edidx) = runTallies(stidx : edidx) + 1;
                    end

                    %-------------------------[DEBUG]---------------------------
                    dmcopy = hmobj.pearsonDistMap;
                    N = size(dmcopy, 1);
                    for i = 1:(N-1)
                        dmcopy(i, (i+1):N) = roundScores(i);
                    end

                    xbins = size(hmobj.xEdgesCulled, 2);
                    ybins = size(hmobj.yEdgesCulled, 2);
                    xlbl = hmobj.adjSymTable{hmobj.xEdgesCulled(1:(xbins-1)), 'AddressString'};
                    ylbl = hmobj.adjSymTable{hmobj.yEdgesCulled(1:(ybins-1)), 'AddressString'};
                    clear xbins ybins

                    figure(100);
                    clf;
                    hm = heatmap(xlbl, ylbl, dmcopy);
                    hm.Colormap = turbo;
                    hm.ColorLimits = [0 1];
                    hm.CellLabelColor = 'none';
                    hm.GridVisible = 'off';
                    title('Pearson Distance w/ Boundary Scores');

                    clear dmcopy N i xlbl ylbl hm
                    %-------------------------[DEBUG]---------------------------

                    %Reroll shift
                    shiftCheck = binShift + 1;
                    while binShiftDone(shiftCheck)
                        shiftCheck = random('unid', binSize);
                        binShift = shiftCheck - 1;
                    end
                    binShiftDone(shiftCheck) = true;
                end

                %Calculate means...
                boundScores = runSums ./ runTallies;
            end
        end

        %%
        function candSets = setTopCandAsPrincipal(obj, heatmapObj, candSets, anchorIdx)
            candSet = candSets(anchorIdx);
            clen = candSet.candTable{1, 'len'};
            cnext = anchorIdx + clen;

            %Pop from table and set as principal
            ctSize = size(candSet.candTable, 1);
            candSet.principal = table2struct(candSet.candTable(1,:));
            if(ctSize > 1)
                candSet.candTable = candSet.candTable(2:ctSize,:);
            else
                candSet.candTable = [];
            end

            %Clear out any upstream or downstream
            %candidates that split this new principal

            %Upstream
            for d = (anchorIdx - 1):-1:1
                otherCandSet = candSets(d);
                if otherCandSet.barrier > 0
                    if otherCandSet.barrier == anchorIdx
                        break;
                    end
                    if (otherCandSet.barrier < cnext) & (otherCandSet.barrier > anchorIdx)
                        otherCandSet.barrier = anchorIdx;
                        candSets(d) = otherCandSet;
                        break;
                    end
                end

                beforeMaxLen = anchorIdx - d;
                inclMinLen = beforeMaxLen + candSet.principal.len;
                addCandBool = true;
                if ~isempty(otherCandSet.principal)
                    if (otherCandSet.principal.len > beforeMaxLen) & (otherCandSet.principal.len < inclMinLen)
                        otherCandSet.principal = [];
                    end
                end
                if ~isempty(otherCandSet.candTable)
                    keepBool = otherCandSet.candTable{:, 'len'} <= beforeMaxLen;
                    keepBool = or(keepBool, (otherCandSet.candTable{:, 'len'} >= inclMinLen));
                    delBool = ~keepBool;
                    if nnz(keepBool) > 0
                        otherCandSet.candTable = otherCandSet.candTable(keepBool, :);
                        hasMe = otherCandSet.candTable{:, 'len'} == inclMinLen;
                        if nnz(hasMe) > 0
                            addCandBool = false;
                        end
                        clear hasMe
                    else
                        otherCandSet.candTable = [];
                    end
                    if nnz(delBool) < 1
                        addCandBool = false;
                    end
                    clear keepBool delBool
                else
                    addCandBool = false;
                end

                if addCandBool
                    newCand = ClusterAnalyzer.genClusterCandStruct();
                    newCand.start = d;
                    newCand.len = inclMinLen;
                    newCand.mergePoint = anchorIdx;
                    newCand = obj.initializeCandidate(newCand, heatmapObj);
                    otherCandSet.candTable = ClusterAnalyzer.addToCandTable(otherCandSet.candTable, newCand);
                end
                candSets(d) = otherCandSet;
            end
            clear d otherCandSet beforeMaxLen inclMinLen newCand addCandBool

            %Downstream
            for d = (anchorIdx + 1):(cnext - 1)
                otherCandSet = candSets(d);
                otherCandSet.principal = [];
                otherCandSet.candTable = [];
                candSets(d) = otherCandSet;
            end
            clear d otherCandSet
            clear ctSize

            candSets(anchorIdx) = candSet;
        end

        %%
        function [score, candSets, cnext] = scoreBound(obj, heatmapObj, candSets, c)
            %TODO
        end

        %%
        function boundScores = clusterHeatmap(obj, heatmapObj)
            boundScores = [];
            if isempty(heatmapObj); return; end
            if isempty(heatmapObj.pearsonDistMap); return; end

            %Create base list and seed clusters
            N = size(heatmapObj.pearsonDistMap,1);
            %boundScores = NaN(1, N);
            candSets(N) = ClusterAnalyzer.genClusterCandSetStruct();
            for i = 1:N
                candSet = candSets(N);
                candSet.anchor = i;
                cEnd = i + obj.seedClusterSize - 1;
                if (cEnd <= N)
                    seedCand = ClusterAnalyzer.genClusterCandStruct();
                    seedCand.start = i;
                    seedCand.len = obj.seedClusterSize;
                    seedCand.mergePoint = cEnd;

                    seedCand = obj.initializeCandidate(seedCand, heatmapObj);
                    if (seedCand.clusterGroup >= 4)
                        %Auto split.
                        candSets(i) = candSet;
                        candSets = addBarrier(candSets, i, cEnd);
                        candSet = candSets(i);
                    else
                        %Put as head.
                        candSet.candTable = struct2table(seedCand);
                        %candSet.sizeLookup{seedCand.len} = seedCand;
                    end
                end
                candSets(i) = candSet;
            end
            clear cEnd candSet i seedCand

            loopBreak = false;
            while ~loopBreak
                %Look for lowest scored candidate
                loopBreak = true;
                minScore = NaN;
                anchorIdx = 0;
                for c = 1:N
                    candSet = candSets(c);
                    if ~isempty(candSet.candTable)
                        if ~candSet.candTable{1, 'checkedFlag'}
                            cScore = candSet.candTable{1, 'score'};
                            if isnan(minScore) | (cScore < minScore)
                                minScore = cScore;
                                anchorIdx = c;
                            end
                            loopBreak = false;
                        end
                    end
                end
                clear cScore candSet c

                if ~isnan(minScore)
                    candSet = candSets(anchorIdx);
                    candSet.candTable{1, 'checkedFlag'} = true;
                    candSets(anchorIdx) = candSet;
                    cGroup = candSet.candTable{1, 'clusterGroup'};
                    clen = candSet.candTable{1, 'len'};
                    cnext = anchorIdx + clen;
                    if (cGroup == 0)
                        %Auto-merge
                        candSets = obj.setTopCandAsPrincipal(heatmapObj, candSets, anchorIdx);
                        candSet = candSets(anchorIdx);
                    end
                    clear cGroup

                    if (cnext <= N)
                        if (candSet.barrier < 1) | (cnext < candSet.barrier)
                            %Score candidates consuming next cluster
                            nextSet = candSets(cnext);
                            newCand = ClusterAnalyzer.genClusterCandStruct();
                            newCand.start = anchorIdx;
                            newCand.mergePoint = cnext;
                            if ~isempty(nextSet.principal)
                                newCand.len = nextSet.principal.len + clen;
                            else
                                newCand.len = clen + 1;
                            end
                            newCand = obj.initializeCandidate(newCand, heatmapObj);
                            if (newCand.clusterGroup < 4)
                                candSet.candTable = ClusterAnalyzer.addToCandTable(candSet.candTable, newCand);
                            else
                                candSets(anchorIdx) = candSet;
                                candSets = addBarrier(candSets, anchorIdx, cnext);
                                candSet = candSets(anchorIdx);
                            end
                            clear nextSet newCand
                        end
                    end
                    clear clen cnext

                    candSets(anchorIdx) = candSet;
                end
            end
            clear minScore loopBreak anchorIdx

            %Clean up sets a little bit
            candSets = ClusterAnalyzer.tidyCandSets(candSets);

            %TODO Turn results into bin bound scores --------------------
            c = 1;
            %backBarrier = 1;
            scoreSums = zeros(1, N);
            scoreTallies = zeros(1, N);
            while c < (N - 1)
                candSet = candSets(c);
                if ~isempty(candSet.principal)
                    cnext = c + candSet.principal.len; %Start of next cluster
                    scoreSums(c:(cnext - 2)) = 0.0;
                    scoreTallies(c:(cnext - 2)) = 1;
                    %backBarrier = cnext - 1;
                    c = cnext - 1;
                    clear cnext
                else
                    %Grab any cluster candidates that cover c and c+1
                    cplus = c+1;
                    for d = c:-1:1
                        candSetCheck = candSets(d);
                        if candSetCheck.barrier > 0
                            if candSetCheck.barrier == cplus
                                scoreSums(c) = scoreSums(c) + obj.scoreWeights(5);
                                scoreTallies(c) = scoreTallies(c) + 1;
                            elseif candSetCheck.barrier <= c
                                break;
                            end
                        end
                        if ~isempty(candSetCheck.candTable)
                            %TODO Only use those that have as a merge
                            %point?
                            candEnd = d + candSetCheck.candTable{:, 'len'};
                            inclBool = cplus <= candEnd;
                            inclCount = nnz(inclBool);
                            if inclCount > 0
                                subTable = candSetCheck.candTable(inclBool, :);
                                for i = 1:inclCount
                                    scoreSums(c) = scoreSums(c) + obj.scoreWeights(subTable{i, 'clusterGroup'} + 1);
                                end
                                scoreTallies(c) = scoreTallies(c) + inclCount;
                                clear subTable i
                            else
                                break;
                            end
                            clear candEnd inclBool inclCount
                        end
                    end
                    clear d candSetCheck cplus
                    c = c+1;
                end
            end
            clear c candSet

            boundScores = scoreSums ./ scoreTallies;
        end

        %%
        function candStruct = initializeCandidate(obj, candStruct, heatmapObj)
            if isempty(candStruct); return; end

            %Updates all the stats and groupings.
            endPoint = candStruct.start + candStruct.len - 1;

            submtx = heatmapObj.pearsonDistMap(candStruct.start:endPoint, candStruct.start:endPoint);
            [xx,yy] = meshgrid(1:candStruct.len, 1:candStruct.len);
            xydiff = xx - yy;
            submtx(xydiff <= 0) = NaN;

            candStruct.mean = mean(submtx, 'all', 'omitnan');
            candStruct.stdev = std(submtx, 0, 'all', 'omitnan');
            candStruct.median = median(submtx, 'all', 'omitnan');

            candStruct.score = (candStruct.mean + candStruct.stdev) ./ 2.0;

            meanGroupCount = size(obj.meanGroupingFuncs, 2);
            for i = 1:meanGroupCount
                y = ClusterAnalyzer.applyMeanGroupingFunc(candStruct.len, obj.meanGroupingFuncs(i));
                if (candStruct.mean <= y)
                    candStruct.meanGroup = i;
                    break;
                end
            end
            if candStruct.meanGroup <= 0
                candStruct.meanGroup = meanGroupCount;
            end

            stdGroupCount = size(obj.stdGroupingFuncs, 2);
            for i = 1:stdGroupCount
                y = ClusterAnalyzer.applyStdGroupingFunc(candStruct.len, obj.stdGroupingFuncs(i));
                if (candStruct.stdev <= y)
                    candStruct.stdGroup = i;
                    break;
                end
            end
            if candStruct.stdGroup <= 0
                candStruct.stdGroup = stdGroupCount;
            end

            %Overall group. 0-4. 0 is auto merge and 4 is auto split
            gSum = (candStruct.stdGroup - 1) + (candStruct.meanGroup - 1);
            if (gSum == 0)
                candStruct.clusterGroup = 0;
            elseif (gSum >= 4)
                candStruct.clusterGroup = 4;
            else
                if (candStruct.meanGroup >= 3)
                    candStruct.clusterGroup = 4;
                else
                    candStruct.clusterGroup = gSum;
                end
            end
        end

    end


    methods (Static)
        
        %%
        function candStruct = genClusterCandStruct()
            candStruct = struct();
            candStruct.start = 0;
            candStruct.len = 0;
            candStruct.mergePoint = 0;
            candStruct.checkedFlag = false;

            candStruct.mean = NaN;
            candStruct.median = NaN;
            candStruct.stdev = NaN;
            candStruct.score = NaN;

            candStruct.meanGroup = -1;
            candStruct.stdGroup = -1;
            candStruct.clusterGroup = -1;
        end

        %%
        function candSetStruct = genClusterCandSetStruct()
            candSetStruct = struct();
            candSetStruct.anchor = 0;
            candSetStruct.barrier = 0;
            candSetStruct.barrierPrimary = 0;

            candSetStruct.principal = [];
            candSetStruct.candTable = [];
            %candSetStruct.queueHead = []; %Unevaluated

            %candSetStruct.sizeLookup = cell(1,256);
        end

        %%
        function candTable = addToCandTable(candTable, candStruct)
            if isempty(candTable)
                candTable = struct2table(candStruct);
            else
                apptbl = struct2table(candStruct);
                candTable = [candTable; apptbl];
                candTable = sortrows(candTable, {'checkedFlag', 'score'});
            end
        end

        %%
        function candSets = addBarrier(candSets, c, barrier)
            %Adds barrier directrly, updates all upstream removing all
            %candidates that go through barrier
            candSet = candSets(c);
            if (candSet.barrierPrimary <= 0) | (barrier < candSet.barrierPrimary)
                candSet.barrierPrimary = barrier;
            end
            candSets(c) = candSet;

            for d = c:-1:1
                otherCandSet = candSets(d);

                %Set barrier
                if (otherCandSet.barrier <= 0) | (barrier < otherCandSet.barrier)
                    otherCandSet.barrier = barrier;
                    %Remove conflicting candidates
                    if ~isempty(otherCandSet.candTable)
                        maxLen = barrier - d;
                        keepBool = otherCandSet.candTable{:, 'len'} <= maxLen;
                        if nnz(keepBool > 0)
                            otherCandSet.candTable = otherCandSet.candTable(keepBool, :);
                        else
                            otherCandSet.candTable = [];
                        end
                    end
                    candSets(d) = otherCandSet;
                end
            end
        end

        %%
        function candSets = tidyCandSets(candSets)
            N = size(candSets, 2);
            for c = 1:N
                candSet = candSets(c);
                if ~isempty(candSet.candTable)
                    %Sort by length
                    candSet.candTable = sortrows(candSet.candTable, 'len');
                end

                if ~isempty(candSet.principal)
                    %Clear out any candidates smaller than the principal
                    if ~isempty(candSet.candTable)
                        keepBool = candSet.candTable{:, 'len'} > candSet.principal.len;
                        if nnz(keepBool) > 1
                            candSet.candTable = candSet.candTable(keepBool, :);
                        else
                            candSet.candTable = [];
                        end
                        clear keepBool
                    end
                end
                candSets(c) = candSet;
            end
        end

        %%
        function y = applyMeanGroupingFunc(x, f)
            y = x - f.x0;
            y = y .* (f.k .* -1.0);
            y = 1 + exp(y);
            y = f.L ./ y;
        end

        %%
        function y = applyStdGroupingFunc(x, f)
            y = x .* -1.0;
            y = y ./ (2.0 .* f.s .* f.s);
            y = exp(y);
            y = y .* (f.a - f.y0);
            y = y + f.y0;
        end

        %%
        function funcInfo = genMeanGroupingFunctionStruct()
            funcInfo = struct();
            funcInfo.x0 = NaN;
            funcInfo.k = NaN;
            funcInfo.L = NaN;
        end
    
        %%
        function funcInfo = genStdGroupingFunctionStruct()
            funcInfo = struct();
            funcInfo.a = NaN;
            funcInfo.s = NaN;
            funcInfo.y0 = NaN;
        end
    
        %%
        function funcParams = genDefaultMeanGroupingFunctions()
            ALL_K = [0.03 0.05];
            ALL_L = [0.9 0.95];
            ALL_X0 = [5 -15];
    
            funcCount = size(ALL_K, 2);
            funcParams(funcCount) = ClusterAnalyzer.genMeanGroupingFunctionStruct();
            for i = 1:funcCount
                funcParams(i).k = ALL_K(i);
                funcParams(i).L = ALL_L(i);
                funcParams(i).x0 = ALL_X0(i);
            end
        end

        %%
        function funcParams = genDefaultStdGroupingFunctions()
            ALL_A = [0.15, 0.25, 0.3, 0.4];
            ALL_S = [5, 5, 5.5, 5.5];
            ALL_Y0 = [0.08, 0.09, 0.11, 0.13];
    
            funcCount = size(ALL_A, 2);
            funcParams(funcCount) = ClusterAnalyzer.genStdGroupingFunctionStruct();
            for i = 1:funcCount
                funcParams(i).a = ALL_A(i);
                funcParams(i).s = ALL_S(i);
                funcParams(i).y0 = ALL_Y0(i);
            end
        end
    
    end

end