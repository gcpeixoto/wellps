%% case_unisim_HP-HFU_metrics 

% Compute metrics analytics for UNISIM and 
% export .csv files to disk only for HP-HFUs 

% class instantiations 
d = DirManager(); 

% mounts standard directory tree
d.mountDir();   

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');

%% Parameters
P = computeParams(G,PROPS);

% field statistics for chosen properties
S = printStats(P,{'DRTN_LN'},'n');

%% DRT choice 

% get all DRTs 
drtlist = S{1}(2:end,1);

% number of significant cells
nofs = 30;

% compute HFUs by DRT
drtSt = findDRTConnections(drtlist, P, 'normalized','ln',nofs,'y', 1);

%% Computing graph metrics 
% Here, we choose parameters to compute the graph metrics over all
% the clusters previously computed with findDRTConnections 

% number of significant cells to be considered. Clusters with less 'nofs'
% cells are ignored
opt.nofs = nofs;    
% linear regression slope +/- tolerance.
opt.seps = 0.05;
% minimum R2 coefficient tolerance.
opt.R2min = 0.9;

% compute
[metricsSt,linregrSt] = computeDRTGraphMetrics(opt,drtSt);


%% Find HP HFUs and compute metric analytics

lrn = fieldnames(linregrSt);

for v = 1:length(lrn)    
    
    if isfield(linregrSt.(lrn{v}),'performance')
        pm = cell2mat(linregrSt.(lrn{v}).('performance'));
        c_hp = find(pm); % only high-performance clusters (pm == 1)                         
        aux = split(lrn{v},'DRT'); 
        drtValue = str2double(aux(2));
        metricsAnalyzer(drtValue,c_hp,'normalized','ln');
        
    else
        fprintf('----> No HP HFU was found for %s.\n',lrn{v});        
    end
                       
end