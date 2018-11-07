%% CASE STUDY: FZI analysis + 6 neighbours 

% study of FZI distribution and FZI-connected clusters from 6-neighbour
% criterion

% Method:
% - Create vector of thresholds for FZI
% - Search for FZI 6-N connections and build groups
% - Find clusters  

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

% compute geometry
Gc = computeGeometry(G);

%% Parameters
% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

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
np = 80;
p = linspace(min(FZIN0),max(FZIN0),np);

% partition indices
ip = cell(np-1,1);
for k = 2:np    
    ip{k-1} = find(p(k-1) <= FZIN & FZIN < p(k));
end
last = find(abs(FZIN - p(k)) <= eps);
ip{end} = [ip{end}; last];

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
        fprintf('---> Group %d: %d components found. Max. component: %d members.\n',k,C{k}.ncomp,C{k}.compSizes(1));
    end
end

%% Computing logs and regression of clusters

% performance parameters
seps = 1e-1;
R2min = 0.9;

for gr = conn   
    for comp = 1:C{gr}.ncomp                
        members = C{gr}.compMembers{comp};
        locs = ip{gr}(members); % get indices of cluster at reservoir    
        
        % compute logs and regression and stores into structure
        C{gr}.Log10Phiz{comp} = P.Log10PHIZ(locs);
        C{gr}.Log10RQIN{comp} = P.Log10RQIN(locs);
        [R,m,b] = regression(P.Log10PHIZ(locs),P.Log10RQIN(locs),'one');
        C{gr}.R2{comp} = R*R;
        C{gr}.slope{comp} = m;
        C{gr}.offset{comp} = b;    
        
        % performance
        if R*R >= R2min && ( 1-seps <= m && m <= 1+seps )
            C{gr}.performance{comp} = true;
        else
            C{gr}.performance{comp} = false;
        end            
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
            %hold on
            %plotCellData(Gfzi,FZIN(locs));                  
        end                
    end
    axis off vis3d
    colormap(jet); colorbar;
    tit = sprintf('plot FZI; partition: %d; min. cluster: %d',gr,nel);
    title(tit);
end

%% Plot only high-performance clusters 


% loop over partitions 
for gr = conn            
    
    figure   
    % loop over connected clusters
    for comp = 1:C{gr}.ncomp
        
        if C{gr}.performance{comp} == true
                
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
    end
    
    axis off vis3d
    colormap(jet); colorbar;
    tit = sprintf('plot FZI; gr.: %d; HP clust.: %d',gr,comp);
    title(tit);
        
end
