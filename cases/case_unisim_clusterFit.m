%% case_unisim_clusterFit
%
% This script will search for nonuniform 5-spot well patterns
% wrapping HFUs over the UNISIM reservoir model
%
%
% Dr. Gustavo Oliveira, @LaMEP

%% Input parameters 

% OPTIONS:
%   -   drt_meth : 'DRT' or 'DRT*';
%   -   perm_ave : 'a', for 'arithmetic' mean 
%                  'g', for 'geometric' mean
%                  'h', for 'harmonic' mean
%                  'q', for 'quadratic' mean
%                  'n', for 'normalized'
%   -   log_base : 'ln' or 'log10';
%   -   nofs_in  : number of significant components for clustering 
%                  (arbitrary ellipsoid fit requires at least 10)
%   -   model    : 'unisim1', for UNISIM-I-D or 
%                  'spe10', for SPE 10th Project, model 2

% clustering
drt_meth = 'DRT*'; 
perm_ave = 'h'; 
log_base = 'ln'; 
nofs_in = 10; 

% analytics
slope_eps = 0.1; 
R2min = 0.9;

% number of elements
drt_ = 12; % if user inputs an invalid number, the program will suggest 
           % the possibilities
           
% model 
model = 'unisim1';

%% Standard directories 

% class instantiations 
d = DirManager(); 

% mounts standard directory tree
d.mountDir();   

%% Load grid 
switch model 
    case 'unisim1'
        [G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
        
        % compute geometry
        G = computeGeometry(G);
        
    % \TODO conversion of permeability data from m2 to mD
    case 'spe10'
        [G,PROPS] = buildModelSPE10();
        G = G;
end
%% Parameters
P = computeParams(G,PROPS);

% fieldname
bas = [];

% append method
if strcmp(drt_meth,'DRT')
    bas = 'DRT';
elseif strcmp(drt_meth,'DRT*')
    bas = 'DRTStar';
end

% append average
if strcmp(perm_ave,'a')
    bas = strcat(bas,'A');
    ave = 'arithmetic';    
elseif strcmp(perm_ave,'g')
    bas = strcat(bas,'G');
    ave = 'geometric';
elseif strcmp(perm_ave,'n')
    bas = strcat(bas,'N');
    ave = 'normalized';    
elseif strcmp(perm_ave,'q')
    bas = strcat(bas,'Q');
    ave = 'quadratic';
elseif strcmp(perm_ave,'h')
    bas = strcat(bas,'H');    
    ave = 'harmonic';
end

if strcmp(drt_meth,'DRT*')
    aux = split(split(bas,'S'),'_');
    aux2 = split(aux(2),'r');
    bas = strcat('DRT',aux2(2),'Star');
end


% append log base
if strcmp(log_base,'log10')
    bas = strcat(bas,'_LOG10');
elseif strcmp(log_base,'ln')
    bas = strcat(bas,'_LN');
end

% convert to char 
bas = char(bas);

% field statistics for chosen properties
S = printStats(P,{bas},'n');

%% DRT choice 

% list of DRTs to consider (if array, the order is important for plotting)
drtlist = S{1}(2:end,1);

if ~ismember(drt_,drtlist)    
    fprintf('!!! For your setup, only the following DRT choices are possible\n:');    
    disp(drtlist)
    error('Stopping...\n');
else
    drtlist = drt_;
end


% number of significant cells
nofs = nofs_in;

% compute HFUs by DRT
drtSt = findDRTConnections(drtlist, P, ave, log_base, nofs, 'n', 1);

if isempty(drtSt)
    fprintf('----> No components for this DRT.\n'); 
    return
end

%% Computing graph metrics 
% Here, we choose parameters to compute the graph metrics over all
% the clusters previously computed with findDRTConnections 

% number of significant cells to be considered. Clusters with less 'nofs'
% cells are ignored
opt.nofs = nofs;    
% linear regression slope +/- tolerance.
opt.seps = slope_eps;
% minimum R2 coefficient tolerance.
opt.R2min = R2min;

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
        
        %c_hp = c_hp(1:5); % uncomment this line to force the method to 
                           % get only the clusters you want. Some small
                           % clusters are returning problems with the 
                           % fitting while computing eigenvalues. \TODO
        
        % No HP HFU for this DRT
        if ~isempty(c_hp)
            % we have 4 patterns per cluster (4 angles)
            clusterFit{v} = processClusterFit(G,drtlist(v),c_hp,ave,log_base,[0,pi/6,pi/4,pi/3]);        
        
        else 
            fprintf('----> No HP HFU was found for %s to be cluster-fittted.\n',lrn{v});      
            return
        end
        
    % No HP HFU for any
    else
            fprintf('----> No HP HFU was found for any DRT to be cluster-fittted.\n');      
            return    
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

