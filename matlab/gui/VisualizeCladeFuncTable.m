%
%%
function figHandle = VisualizeCladeFuncTable(cladeNodes, consensusFuncTable, emptySymId)

%Colors:
%   0 - Unused slot (Dark grey)
%   1 - NULL (Red)
%   2 - Empty Call (Magenta)
%   3 - Consensus (Cyan)
%   4 - Non-consensus unique (White)
%   5 - Non-consensus non-unique (Light grey)

COLOR_UNUSED = [0.1 0.1 0.1];
COLOR_NULL = [1 0 0];
COLOR_EMPTY = [1 0 1];
COLOR_CON = [0 1 1];
COLOR_UNQ = [1 1 1];
COLOR_SUNQ = [0.9 0.9 0.9];

COLOR_TABLE = [COLOR_UNUSED; COLOR_NULL; COLOR_EMPTY; COLOR_CON; COLOR_UNQ; COLOR_SUNQ];
COLOR_COUNT = size(COLOR_TABLE, 1);

%Create string and color matrices...
nCount = size(cladeNodes, 1);
maxFuncCount = max(cladeNodes{:, 'len'}, [], 'all', 'omitnan');
nodeNames = cladeNodes{:, 'addressString'}';
tblTypes = repmat({'string'}, [1 nCount]);

strTable = table(VariableNames=nodeNames, VariableTypes=tblTypes, Size=[maxFuncCount nCount]);
clrIdMtx = zeros(maxFuncCount, nCount);
allmtx = zeros(nCount, maxFuncCount);

conLen = size(consensusFuncTable, 2);

for i = 1:nCount
    nTable = cladeNodes{i, 'table'};
    if iscell(nTable)
        nTable = nTable{1};
    end

    tsize = size(nTable, 1);
    strTable{1:tsize, i} = nTable{:, 'AddressString'};
    allmtx(i, 1:tsize) = nTable{:, 'SymbolId'};

    if(tsize < maxFuncCount)
        strTable{(tsize+1):maxFuncCount, i} = "";
    end
end

for i = 1:nCount
    %Assign colors
    tsize = cladeNodes{i, 'len'};
    clrNull = (allmtx(i, 1:tsize) == 0);
    clrEmpty = (allmtx(i, 1:tsize) == emptySymId);

    conMax = min(conLen, tsize);
    clrConsensus = false(1, tsize);
    clrConsensus(1:conMax) = (allmtx(i, 1:conMax) == consensusFuncTable(1:conMax));
    clrConsensus = and(clrConsensus, ~clrNull);
    clrConsensus = and(clrConsensus, ~clrEmpty);

    clrUnique = and(~clrNull, ~clrEmpty);
    clrUnique = and(clrUnique, ~clrConsensus);
    clrSemiUnique = false(1, tsize);

    for j = 1:tsize
        if clrUnique(j)
            fmatch = (allmtx(:, j) == allmtx(i,j));
            if (nnz(fmatch) > 1)
                clrUnique(j) = false;
                clrSemiUnique(j) = true;
            end
        end
    end

    clrIdMtx(clrNull, i) = 1;
    clrIdMtx(clrEmpty, i) = 2;
    clrIdMtx(clrConsensus, i) = 3;
    clrIdMtx(clrUnique, i) = 4;
    clrIdMtx(clrSemiUnique, i) = 5;
end

%https://www.mathworks.com/help/matlab/ref/uitable.html
figHandle = uifigure('Name', 'Clade Function Table');
uit = uitable(figHandle, "Data", strTable, "Position", [10 10 720 480]);

for c = 1:COLOR_COUNT
    cIndex = (c-1);
    [matchRow, matchCol] = find(clrIdMtx == cIndex);

    myColor = COLOR_TABLE(c, :);
    uis = uistyle("BackgroundColor", myColor);
    addStyle(uit, uis, "cell", [matchRow, matchCol]);
end

end