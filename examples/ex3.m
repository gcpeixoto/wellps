%% EXAMPLE: Ellipsoid fit 
% In this tutorial, we study the shape of HFU clusters based on an 
% ellipsoid fit process. 
% See Oliveira et al. (2019); doi: 10.1016/j.petrol.2019.106680

mrstVerbose off  % turn on verbose

case_name = 'ex3'; % name for this case study

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Grid reading
f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');

[G,~] = buildModel(f);
G = computeGeometry(G);

%% Structure of options for ellipsoid fitting
% The minimum fields for input are 'grid', 'thetalist', and 'savedir',
% but the are not enough for the method to work properly.
%
% Then, in general, the following situations are expected:
%
% 1. If 'loaddir' is not given, 'metricsSt' and 'drtSt' must be provided
%    as .mat files for the metrics and DRT structures. In this case, 
%    'drtlist', 'logbase', and 'krule' are optional, since each one will be 
%    loaded in-loop. 
%    
%
% 2. If 'loaddir' is provided, 'drtlist', 'logbase', and 'krule' must also
%    be provided so that WELLPS can load the wanted files properly. 
%    If 'drtSt' or 'metricsSt' are given, they will be ignored.
%
% For both cases, the field 'comps' is optional. If provided, the cluster 
% loop only will run in the range supplied. Otherwise, all comps will be
% computed.


opt_fit.loaddir = d.getMatDir;
opt_fit.drtSt = fullfile(d.getRootDir,'examples','sample','C.mat');
opt_fit.metricsSt = fullfile(d.getRootDir,'examples','sample','CMetrics.mat');
opt_fit.grid = G; 
%opt_fit.drtlist = [12,14];
%opt_fit.logbase = 'ln';
%opt_fit.krule= 'geometric';
opt_fit.comps = 1:3;
opt_fit.thetalist = [0,pi/2,3*pi/2,2*pi];
opt_fit.savedir = fullfile(d.getMatDir,case_name);
     

%% Process cluster fit 
% This method will compute the fit parameters for given DRT,
% cluster(s) and rotation angle(s) based on a best-fit ellipsoid
% for the cluster shape according to 'opt_fit' choices. 
% These parameters are useful to build 5-spot nonuniform well patterns. 
% See the function documentation to understand what is going on behind 
% the scenes and the paper aforementioned.
clusterFitSt = processClusterFit(d,opt_fit);


%% Plot ellipsoid fit 
% Here, we plot a 3D geometric view of the ellipsoid fitting method 
% applied to clusters of the sample files. Firtsly, we set some
% parameters and then invoke 'plotEllipsoidFit', a special function
% that will draw graphical elements. 

drt = drtSt.value;
pattern = 1; % recall this number should agree with the angle
factor = 80; 
nres = 70; 
comps = [1,2,3];                    
for c = comps
[figelip,figres] = ...
    plotEllipsoidFit(G,drtSt,clusterFitSt,drt,c,pattern,factor,nres);                            
end


%% Build a nonuniform 5-spot well pattern model for a given cluster
% Here, we form the nonuniform standard 5-spot as follows:
%
% - producer well: column of all (vertical) cells neighbour to the maximum
%                  closeness centrality cell 
% - injection wells: column of all (vertical) cells neighbour to each of 
%                    the 4 cluster-fit cells +X,-X,+Y,-Y for the cluster
%
% REMARK: the z-range considered is that one of the cluster.

% parameters
drt = fieldnames(clusterFitSt); drt = drt{1};
cluster = strcat('C',num2str(comps(1)));
pat = strcat('Pattern',num2str(pattern));


% -- {+X, -X, +Y, -Y injectors} + {producer} 
injp_xy = {'colNeighsX1', ...
           'colNeighsX2', ...
           'colNeighsY1', ...
           'colNeighsY2', ...
           'colNeighsMaxC'}; 

% auxiliary function
agg = @(col) clusterFitSt.(drt).(cluster).(pat).(col); 

% fields
col = cellfun(agg,injp_xy,'UniformOutput',false);
[col_inj_1,col_inj_2,col_inj_3,col_inj_4,col_prod] = deal(col{:});

% casting to cell
cols = {col_inj_1,col_inj_2,col_inj_3,col_inj_4,col_prod}; 


% -- getting cell indices


% auxiliary function
getid = @(col) Ind( sub2ind(G.cartDims,col(:,1),col(:,2),col(:,3)) );
cols = cellfun(getid,cols,'UniformOutput',false); 
[icol_inj_1,icol_inj_2,icol_inj_3,icol_inj_4,icol_prod] = deal(cols{:});
ID = {icol_inj_1,icol_inj_2,icol_inj_3,icol_inj_4,icol_prod};


% plot 5-spot columns 
figure 
plotGrid(G, Ind(drtSt.compVoxelInds{comps(1)}),'FaceColor',[0.6,0.6,0.6], ...
    'FaceAlpha',0.05, 'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.1)

% producer
plotGrid(G,icol_prod(~isnan(icol_prod)),'FaceColor','r','EdgeColor','k')

% injectors
plt = @(id) plotGrid(G,id(~isnan(id)),'FaceColor','b','EdgeColor','k');
for k = 1:4, plt(ID{k}), axis off vis3d; end
