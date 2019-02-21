%% case_unisim_clusterFit
%
% This script will search for nonuniform 5-spot well patterns
% wrapping HFUs over the UNISIM reservoir model
%
%
% Dr. Gustavo Oliveira, @LaMEP

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
% From S, we see that the most populated DRT is 13. We then are going to
% choose it in the list to locate HP-HFUs with significant cells above 30.
% Both DRT and nofs can be freely chosen. Here, we have used them for
% testing.

% list of DRTs to consider (if array, the order is important for plotting)
drtlist = 13;

% number of significant cells
nofs = 30;

% compute HFUs by DRT
drtSt = findDRTConnections(drtlist, P, 'normalized','log10',nofs,'n', 1);

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


%% Process cluster fit for HP HFUs
% The results are saved into the cell 'clusterFit', whose lenght is the 
% same as drtlist. Each entry of 'clusterFit' is a struct containing all 
% the well patterns computed for each HP HFU. 

lrn = fieldnames(linregrSt);
clusterFit = cell(1,length(lrn)); 

for v = 1:length(lrn)    
    
    if isfield(linregrSt.(lrn{v}),'performance')
        pm = cell2mat(linregrSt.(lrn{v}).('performance'));
        c_hp = find(pm); % only high-performance clusters (pm == 1)                                                                
    
        if ~isempty(c_hp)
            % we have 4 patterns per cluster (4 angles)
            clusterFit{v} = processClusterFit(G,drtlist(v),c_hp,'normalized','log10',[0,pi/6,pi/4,pi/3]);        
        
        else 
            fprintf('----> No HP HFU was found for %s to be cluster-fittted.\n',lrn{v});      
        end
    end
                       
end


% ========================================================
%% Plot nonuniform 5-spot well pattern for a given HP HFU

% selector
tell_drtlistIndex = 1; % in this case, 1: DRT13
tell_cluster = c_hp; % cluster index in c_hp
tell_pattern = [1,2,3,4]; % pattern 1, 2, 3 or 4

for t_c = tell_cluster

    for t_p = tell_pattern

        % option struct (for future implementation as a function)
        opts.drtindex = tell_drtlistIndex; 
        opts.cluster = t_c; 
        opts.pattern = t_p;

        % inverse mapping
        Ind = nan(prod(G.cartDims),1);
        Ind(G.cells.indexMap) = 1:G.cells.num;

        % pattern checking
        if opts.pattern == 1 
            fprintf('----> Plotting well pattern for angle: 0 deg.\n')

        elseif opts.pattern == 2 
            fprintf('----> Plotting well pattern for angle: 30 deg.\n')

        elseif opts.pattern == 3 
            fprintf('----> Plotting well pattern for angle: 45 deg.\n')

        elseif opts.pattern == 4 
            fprintf('----> Plotting well pattern for angle: 60 deg.\n')

        else
            error('----> Wrong well pattern option! Must be 1, 2, 3 or 4.\n')

        end

        % structure tip prompt to pull well pattern's internal fields of interest 
        St_base = clusterFit{opts.drtindex}.(lrn{opts.drtindex}). ... 
                 (strcat('C',num2str(opts.cluster))). ...
                 (strcat('Pattern',num2str(opts.pattern)));

        % getting cell indices (injectors)
        injectors_cols = {St_base.colNeighsX1R; ... 
                          St_base.colNeighsX2R; ... 
                          St_base.colNeighsY1R; ...
                          St_base.colNeighsY2R};              

        injectors_cellInds = cell(1,4);

        for ic = 1:4    
            aux = injectors_cols{ic};    
            icol = sub2ind(G.cartDims,aux(:,1),aux(:,2),aux(:,3));
            icol = Ind(icol);    
            injectors_cellInds{ic} = icol;
        end

        % getting cell indices (producer)
        producer_col = St_base.colNeighsMaxC;
        producer_cellInds = sub2ind(G.cartDims,producer_col(:,1), ...
                                               producer_col(:,2), ...
                                               producer_col(:,3));
        producer_cellInds = Ind(producer_cellInds);                                   


        % plot cluster domain
        figure 
        plotGrid(G, Ind(drtSt.(lrn{opts.drtindex}).compVoxelInds{opts.cluster}),...
            'FaceColor',[0.6,0.6,0.6], ...
            'FaceAlpha',0.05,          ...
            'EdgeColor',[0.6,0.6,0.6], ...
            'EdgeAlpha',0.1)

        % plot 5-spot columns 

        % injectors
        for ic = 1:4
            aux = injectors_cellInds{ic};
            plotGrid(G,aux(~isnan(aux)),'FaceColor','b','EdgeColor','k')
        end

        % producer
        plotGrid(G,producer_cellInds(~isnan(producer_cellInds)), ... 
            'FaceColor','r','EdgeColor','k')

        axis off vis3d
    
    end

end

