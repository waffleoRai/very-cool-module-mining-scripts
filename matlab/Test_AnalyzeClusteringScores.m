%
%%
addpath('./util');

TEST_CUTOFFS = [0.75 0.7 0.65 0.6 0.55 0.5 0.45];

SAVE_DIR = 'D:\usr\bghos\code\ts3_common_re\matlab';
SAVE_PATH = [SAVE_DIR '\matgraph.mat'];
CSAVE_PATH = [SAVE_DIR '\clusterscores_10575c90_107d85c0.mat'];

load(CSAVE_PATH, 'boundScores', 'symListShort');
[symTable, refTable, funcList] = ProgramGraph.loadSavedTables(SAVE_PATH);

%Trim End
n = size(boundScores, 2);
n = n - 30;
boundScores = boundScores(1:n);

figure(300);
clf;
plot(1:n, boundScores);

boundsSm = smooth(boundScores)';
figure(301);
clf;
plot(1:n, boundsSm);

%Look for local maxima at each level
d1 = diff(boundsSm);
lastTop = 1.0;
