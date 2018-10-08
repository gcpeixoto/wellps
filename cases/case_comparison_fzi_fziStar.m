%% CASE: comparison of clustering generated from FZI and FZI*

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

% compute geometry
Gc = computeGeometry(G);

%% Parameters
% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

% field statistics for chosen properties
SN = printStats(P,{'DRTN_LN','DRTN_LOG10','DRTNStar_LN','DRTNStar_LOG10'},'n');

%% DRT list w/out zero

drt_ln = SN{1}(2:end,1);
drtStar_ln = SN{3}(2:end,1);

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

% cells with value
active = find(~isnan(Ind));

%% DRT connections

% number of significant cells (only to save info)
nofs = 50;

drt_ln = findDRTConnections(drt_ln, P, 'normalized','ln',nofs,'y', 1);
drtStar_ln = findDRTConnectionsByFZIStar(drtStar_ln, P, 'normalized','ln',nofs,'y', 1);


%% metrics
opt.nofs = nofs;    
opt.seps = 0.05;
opt.R2min = 0.9;

% compute
[M,L] = computeDRTGraphMetrics(opt,drt_ln);
[MStar,LStar] = computeDRTStarGraphMetrics(opt,drtStar_ln);


%% get clusters candidates to production

% sweep based on usual FZI
dln = fieldnames(drt_ln);
nln = numel(dln);

for i = 1:nln    
    
    % if cluster has connected components, it has the field 'performance'
    if isfield(L.(dln{i}),'performance') 
        
        % search clusters of high-performance (HP)
        p = L.(dln{i}).performance;                
        q = cell2mat(p);
        c = find(q == true); 
        
        if isempty(c)            
            continue;
        end

        % get max closeness cell, slope and R2 for HP clusters
        mcv = cell(1,numel(c));
        slopes = cell(1,numel(c));
        R2 = cell(1,numel(c));
        for cid = c                           
            mcv{cid}= M.(dln{i}).maxClosenessVoxelCoords{cid};            
            slopes{cid} = L.(dln{i}).slope{cid};            
            R2{cid} = L.(dln{i}).R2{cid};            
        end
       
        % save clusters candidates to production        
        CFZI.(dln{i}).idComp = c;             
        CFZI.(dln{i}).slope = slopes;
        CFZI.(dln{i}).R2 = R2;
        CFZI.(dln{i}).maxClosenessVoxelCoords = mcv;        
   end
        
end

% sweep based on FZI*
dlns = fieldnames(drtStar_ln);
nlns = numel(dlns);

for i = 1:nlns    
    
    % if cluster has connected components, it has the field 'performance'
    if isfield(LStar.(dlns{i}),'performance') 
        
        % search clusters of high-performance (HP)
        p = LStar.(dlns{i}).performance;                
        q = cell2mat(p);
        c = find(q == true); 
        
        if isempty(c)            
            continue;
        end

        % get max closeness cell, slope and R2 for HP clusters
        mcv = cell(1,numel(c));
        slopes = cell(1,numel(c));
        R2 = cell(1,numel(c));
        for cid = c                           
            mcv{cid}= MStar.(dlns{i}).maxClosenessVoxelCoords{cid};            
            slopes{cid} = LStar.(dlns{i}).slope{cid};            
            R2{cid} = LStar.(dlns{i}).R2{cid};            
        end
       
        % save clusters candidates to production        
        CFZIStar.(dlns{i}).idComp = c;             
        CFZIStar.(dlns{i}).slope = slopes;
        CFZIStar.(dlns{i}).R2 = R2;
        CFZIStar.(dlns{i}).maxClosenessVoxelCoords = mcv;        
   end
        
end
   

[~,fn] = csvHeader('../tmp/test','a','b','c','d');