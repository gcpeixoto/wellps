%% CASE STUDY: cross-validation of centrality metrics from SNAP and NETWORKX
%
% This test shows that betweeness is incorrectly computed by SNAP, since 
% negative values are found. Closeness, on the other hand, seems correct.
% Degree centrality cannot be directly compared because in NETWORKX it is
% fractional.
% 
% That is why we will use NETWORKX as standard to compute betweeness.

mrstVerbose off  % turn on verbose

case_name = 'comp_SN';

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Grid processing
% See file "examples/ex0.m" for details

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');
[G,PROPS] = buildModel(f);
on = PROPS.ACTNUM == 1;
P = computeParams(G,PROPS);
S = getStats(d,P,{'DRTN_LN','DRTN_LOG10'},'n');

%% Find DRT connections

% number of significant cells
nofsc = 50;

% gets only DRT = 14
drt = S{1}(6,1);

average = 'geometric'; 
logbase = 'ln';
drtSt = findDRTConnections(d,drt, P, average, logbase, nofsc,'n', 1);

%% Computing graph metrics 
% Here, we choose parameters to compute the graph metrics over all
% the clusters previously computed with findDRTConnections 

% number of significant cells to be considered. Clusters with less 'nofsc'
% cells are ignored
opt.nofsc = nofsc;    
% linear regression slope +/- tolerance.
opt.seps = 0.05;
% minimum R2 coefficient tolerance.
opt.R2min = 0.9;

%% Compute with SNAP
[Mf,Lf] = computeDRTGraphMetrics(opt,drtSt);

% Metrics analysis 
analytics.loaddir = d.getMatDir;
analytics.drtlist = drt;
analytics.krule = average; % use the same as that for findDRTConnections! 
analytics.logbase = logbase; % use the same as that for findDRTConnections!
analytics.savedir = fullfile(d.getCsvDir,[case_name,'_snap']);
analytics.fileperf = true;
analytics.fileminmax = true;
analytics.filemetrics = true;

dataDir = metricsAnalysis(analytics);

%% Compute with NETWORKX
[MfN,LfN] = computeDRTGraphMetricsNetworkx(opt,drtSt);

% Metrics analysis 
analytics.loaddir = d.getMatDir;
analytics.drtlist = drt;
analytics.krule = average; % use the same as that for findDRTConnections! 
analytics.logbase = logbase; % use the same as that for findDRTConnections!
analytics.savedir = fullfile(d.getCsvDir,[case_name,'_networkx']);
analytics.fileperf = true;
analytics.fileminmax = true;
analytics.filemetrics = true;

dataDirN = metricsAnalysis(analytics);

%% Get files to analyze
nf = 13; % file
ds = dir(dataDir);
ds2 = importdata(fullfile(ds(nf).folder,ds(nf).name));
bet = ds2.data(:,4);
clo = ds2.data(:,5);

dsn = dir(dataDirN);
dsn2 = importdata(fullfile(dsn(nf).folder,dsn(nf).name));
betn = dsn2.data(:,4);
clon = dsn2.data(:,5);

subplot(1,2,1)
%plot(bet,betn,'o');
stem(1:length(bet),bet - betn);
aux = split(ds(nf).name,'_');
title(['SNAP vs. NETWORKX Betweeness true error :: ',aux{5},' ', aux{6}]);
xlabel('nodes');
ylabel('error');
axis equal

subplot(1,2,2)
%plot(clo,clon,'o');
stem(1:length(bet),clo - clon);
aux = split(ds(nf).name,'_');
title(['SNAP vs. NETWORKX Closeness true error:: ',aux{5},' ', aux{6}]);
xlabel('nodes');
ylabel('error');
axis equal
