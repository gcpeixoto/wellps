%% Case: UNISIM - Constrained Clustering
% This case applies the Peixoto's algorithms for slope constrained clustering 

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
G = computeGeometry(G);

%% Compute required parameters
P = computeParams(G,PROPS);

%% Constrained clustering

% log10(phiz) x log10(RQI); normalized averaging technique 
Log10PHIZ = P.Log10PHIZ(:);     
Log10RQIN = P.Log10RQIN(:);

% constraint slope
seps = 1e-1;

% % clustering method 
method = setSCCMethod('nl'); % 'c'; 'ct'; 'yp'; 'nl'; '6n' 

switch method

    case 'c'        
        tstart = tic;
        [Clusters,R2,M,B] = slopeConstrainedClusteringC(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCC.mat';                

    case 'ct'
        tstart = tic;
        [Clusters,R2,M,B] = slopeConstrainedClusteringCT(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCT.mat';                
        
    case 'py'
        tstart = tic;
        [Clusters,R2,M,B] = slopeConstrainedClusteringPY(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCPY.mat';                
    
    case 'nl'
        tstart = tic;
        [Clusters,R2,M,B] = slopeConstrainedClusteringNL(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCNL.mat';                
    
    case '6n'
        tstart = tic;
        [Clusters,R2,M,B] = slopeConstrainedClustering6N(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCC6N.mat';                
        
end

% saving clustering structure
SCC.method = method;
SCC.clustering = Clusters;
SCC.slope = M;
SCC.R2 = R2;
SCC.B = B;
SCC.execTime = telapsed;

save(fname,'SCC');
