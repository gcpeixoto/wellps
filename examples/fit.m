%% EXAMPLE: Ellipsoid fit 

%% Load grids

mrstVerbose off  % turn on verbose

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Grid reading
f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');

[G,~] = buildModel(f);
G = computeGeometry(G);

%% Load sample cluster files 

sdir = fullfile(d.getRootDir,'examples','sample');

% Here, we load the sample .mat files 
load(fullfile(sdir,'C.mat'),'drtSt');
load(fullfile(sdir,'CLinRegr.mat'),'L');
load(fullfile(sdir,'CMetrics.mat'),'M');


%{


%% Tests

simulated = '14_23';

switch simulated
   
    % UNISIM : DRT
    case '13'
        load('../mat/DRTConnections-unisim1/DRT_harmonic_ln_13.mat');
        load('../mat/ClusterFitData/unisim1/unisim1DRT_harmonic_ln_13_clusterFitData_unisim1.mat');
        
        drt = 13;
        pattern = 1;
        factor = 60; 
        nres = 70; 
        comps = [2,4,5];            
        
        for i = 1:3                                        
            [figelip,figres] = ...
            plotEllipsoidFit(GU,drtSt,clusterFitSt,drt,comps(i),pattern,factor,nres);                            
        end
                
    case '5_95'         
        load('../mat/DRTConnections-spe10/DRT_geometric_ln_5.mat');
        load('../mat/ClusterFitData/spe10/spe10DRT_geometric_ln_5_clusterFitData_spe10.mat');
        
        drt = 5;
        comp = 95;
        pattern = 1;
        factor = 40; 
        nres = 70; 
        
        [figelip,figres] = ...
            plotEllipsoidFit(GSPE,drtSt,clusterFitSt,drt,comp,pattern,factor,nres);
        

    case '6_5'         
        load('../mat/DRTConnections-spe10/DRT_geometric_ln_6.mat');
        load('../mat/ClusterFitData/spe10/spe10DRT_geometric_ln_6_clusterFitData_spe10.mat');
        
        drt = 6;
        comp = 5;
        pattern = 1;
        factor = 40; 
        nres = 70; 
        
        [figelip,figres] = ...
            plotEllipsoidFit(GSPE,drtSt,clusterFitSt,drt,comp,pattern,factor,nres);
        
   
    case '14_23'         
        load('../mat/DRTConnections-spe10/DRT_geometric_ln_14.mat');
        load('../mat/ClusterFitData/spe10/spe10DRT_geometric_ln_14_clusterFitData_spe10.mat');
        
        drt = 14;
        comp = 23;
        pattern = 1;
        factor = 40; 
        nres = 70;         
        
        [figelip,figres] = ...
            plotEllipsoidFit(GSPE,drtSt,clusterFitSt,drt,comp,pattern,factor,nres);
        
end


%}
