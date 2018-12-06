%% CASE STUDY: FZI analysis + 6 neighbours 

% study of FZI distribution and FZI-connected clusters from 6-neighbour
% criterion

% Method:
% - Create partitions from FZI thresholds 
% - Load pre-computed partitioning obtained from 
%   'slopeConstrainedClustering6N' algorithm.
% - Analyse and plot results 

%% Load data 
if exist('data/SCC6N_unisim1.mat','file')
    load('data/SCC6N_unisim1.mat'); % get structure 'SCC'    
else
    error('You need to compute wellps:slopeConstrainedClustering6N.m before proceeding.')
end

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');

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

%{
figure 
subplot(121)
title('FZI; KN-based')
FZIN = FZIN(:);
FZINA = FZIN(active);
plotCellData(G,FZINA,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(122) 
histogram(FZIN0,'Normalization','probability')
%}

%% FZI > 0 

figure
subplot(321)
locs0 = find(FZIN > 0);
FZIN0 = FZIN(locs0);
G0 = extractSubgrid(G,Ind(locs0));
plotCellData(G0,FZIN0,'EdgeColor','none')
colorbar
title('FZI > 0; KN-based')

subplot(322) 
histogram(FZIN0,'Normalization','probability')

%% FZIN25
subplot(323)
locs25 = find(FZIN > 0.25*max(FZIN(:)));
FZIN25 = FZIN(locs25);
G25 = extractSubgrid(G,Ind(locs25));
plotCellData(G25,FZIN25,'EdgeColor','none')
colorbar
title('FZI > 0.25*max(FZI); KN-based')

subplot(324) 
histogram(FZIN25,'Normalization','probability')

%% FZIN50
subplot(325)
locs50 = find(FZIN > 0.5*max(FZIN(:)));
FZIN50 = FZIN(locs50);
G50 = extractSubgrid(G,Ind(locs50));
plotCellData(G50,FZIN50,'EdgeColor','none')
colorbar
title('FZI > 0.5*max(FZI); KN-based')

subplot(326) 
histogram(FZIN50,'Normalization','probability')



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


%% Partition marker

% This will produce a marker to each individual element (point) that will
% label it with the partition number based on the FZI thresholds. It is
% possible to have unitary clusters (1 element) as well as clusters with 
% more elements. Also, it is possible to have a 'frontier' cluster in which
% there are elements with more than 1 marker, i.e. belonging to different 
% partitions. While the field 'markerPerElem' stores the marker value for
% each element, the field 'markerPerCluster' assigns a value to only those
% clusters whose elements are all in the same partition. Otherwise, the
% cluster is marked with '0'.
%
% Arrays are appended into structure 'SCC' loaded from SCC6N-algorithm.
for c = 1:length(SCC.clusters)
               
    % cluster members
    members = SCC.clusters{c};
    
    % allocate
    marker = 0*members; 
        
    % find partition where the members lie in and marks to that partition
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
    if all( marker/max(marker) == ones(1,length(marker)) ) % same value?
        SCC.markerPerCluster{c} = max(marker); 
    else
        SCC.markerPerCluster{c} = 0;    
    end
    
end

%% Clusters per partition 
%
% PC stores all the clusters per partition
%
% PC: partition cluster
% |
% |- PC_1: { C1^1, C2^1, ... Cn1^1 } % clusters of partition 1
% |- PC_2: { C1^2, C2^2, ... Cn2^2 } % clusters of partition 2
% .
% |- PC_p: { C1^p, C2^p, ... Cnp^p } % clusters of partition p
%
% Note that Ci^j are cells of arbitrary lengths

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

%% R^2, slope per partition

% To plot the clusters per partition, we need to get their original indices 

R2C = PC; % R^2 values per partition 
MC = PC; % slope values per partition
LOCS = PC; % linear indices per partition

for i = 1:length(PC)

    clust = PC{i};
    
    locs_pc = [];
    
    %figure
    for c = 1:length(clust)        
        
        locs = SCC.clusters{clust(c)};
        
        locs_pc = [locs_pc, locs];
        
        %Gfzi = extractSubgrid(G,Ind(locs));                       
        %hold on
        %plotCellData(Gfzi,FZIN(locs));    
    end
        
    [R,m,b] = regression(Log10Phiz(locs_pc),Log10RQIN(locs_pc),'one');
    R2C{i} = R*R;
    MC{i} = m;
    LOCS{i} = locs_pc;
    
end

%% Search of HFUs based only on connections, without constraints 
% Slope and R2 not guaranteed.

turn_on = true;

if turn_on

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

end