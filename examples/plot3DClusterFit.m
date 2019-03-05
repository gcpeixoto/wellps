%% EXAMPLE: Ellipsoid fit 

% Cases simulated for paper 
% - Modelo UNISIM 1
%   (DRT = 13, K = harmonic, log_base = LN, cluster = 2) 
%   (DRT = 13, K = harmonic, log_base = LN, cluster = 4) 
%   (DRT = 13, K = harmonic, log_base = LN, cluster = 5) 
%
% - Modelo SPE 10/2 
%   (DRT = 5, K = geometric, log_base = LN, cluster = 95)
%   (DRT = 6, K = geometric, log_base = LN, cluster = 5) 
%   (DRT = 14, K = geometric, log_base = LN, cluster = 23)

%% Load grids

% unisim
[GU,~] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');                
GU = computeGeometry(GU);
     
% spe10
[GSPE,~] = buildModelSPE10('original');        

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
