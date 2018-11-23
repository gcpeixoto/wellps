%% CASE STUDY: FZI analysis + 6 neighbours 

% study of FZI distribution and FZI-connected clusters from 6-neighbour
% criterion

% Method:
% - Create vector of thresholds for FZI
% - Search for FZI 6-N connections and build groups
% - Find clusters  

%% Load data 
if exist('data/SCC6N.mat','file')
    load('data/SCC6N.mat');    
else
    error('You need to compute wellps:slopeConstrainedClustering6N.m before proceeding.')
end

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

% compute geometry
Gc = computeGeometry(G);

%% Parameters
% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

Log10Phiz = P.Log10PHIZ(:);
Log10RQIN = P.Log10RQIN(:);

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

% cells with value
active = find(~isnan(Ind));

%% distribution of FZI
FZIN = P.FZIN;

subplot(121)
title('FZI; LN')
FZIN = FZIN(:);
FZINA = FZIN(active);
plotCellData(G,FZINA,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(122)
histogram(FZINA)

%% w/out 0 

FZIN0 = FZIN(FZIN > 0);
subplot(121) 
histogram(FZIN0)

%% Partitions 


% partition number and thresholds
np = 20;
p = linspace(min(FZIN0),max(FZIN0),np);

% partition indices
ip = cell(np-1,1);
for k = 2:np    
    ip{k-1} = find(p(k-1) <= FZIN & FZIN < p(k));
end
last = find(abs(FZIN - p(k)) <= eps);
ip{end} = [ip{end}; last];


%% 

% create marker
for c = 1:length(SCC.clusters)
               
    members = SCC.clusters{c};
    
    marker = 0*members; 
        
    for m = 1:numel(members)
        for pp = 1:length(ip) 
            if ~isempty(find(cell2mat(ip(pp)) == members(m), 1))
                marker(m) = pp;              
            end          
        end       
    end
    
    % marker per element
    SCC.markerPerElem{c} = marker;
    
    % marker per cluster
    if all( marker/max(marker) == ones(1,length(marker)) )
        SCC.markerPerCluster{c} = max(marker);    
    end
    
end

% store clusters per partition
PC = ip; 

for i = 1:length(PC)
    
    aux = [];
    for c = 1:length(SCC.markerPerCluster)
                
        mk = SCC.markerPerCluster{c};                
        
        if ~isempty(mk) && mk == i
            aux = [aux,c];
        end
    end
    
    PC{i} = aux;
    
end

% plot clusters per partition 

R2C = PC;
MC = PC;
LOCS = PC;
for i = 1:length(PC)

    clust = PC{i};
    
    locs_pc = [];
    
    figure
    for c = 1:length(clust)        
        
        locs = SCC.clusters{clust(c)};
        
        locs_pc = [locs_pc, locs];
        
        Gfzi = extractSubgrid(G,Ind(locs));                       
        %hold on
        %plotCellData(Gfzi,FZIN(locs));    
    end
        
    [R,m,b] = regression(Log10Phiz(locs_pc),Log10RQIN(locs_pc),'one');
    R2C{i} = R*R;
    MC{i} = m;
    LOCS{i} = locs_pc;
    
end



%{

% partition voxels
vp = cell(size(ip));

% FZI connections
C = vp;

conn = [];
for k = 1:numel(ip)
    [I,J,K] = ind2sub(G.cartDims,ip{k});
    vp{k} = [I,J,K];
    C{k} = findConnectionsSimple(vp{k});
    if ~isempty(C{k})
        conn = [conn,k];
        fprintf('---> Group %d: %d components found. Maximum component has %d members.\n',k,C{k}.ncomp,C{k}.compSizes(1));
    end
end



%% Plot clusters per group

% minimum number of components to consider
nel = 10;

% loop over partitions 
for gr = conn            
    
    figure
    
    % loop over connected clusters
    for comp = 1:C{gr}.ncomp
                
        members = C{gr}.compMembers{comp};
        locs = ip{gr}(members); % get indices of cluster at reservoir    
        
        % plot all the connected clusters with nel elements for given FZI
        % range
        if numel(members) >= nel && ~isempty(locs) 
            Gfzi = extractSubgrid(G,Ind(locs));               
            hold on
            plotCellData(Gfzi,FZIN(locs));                  
        end                
    end
    axis off vis3d
    colormap(jet); colorbar;
    tit = sprintf('plot FZI; partition: %d; min. cluster: %d',gr,nel);
    title(tit);
end

%}