%% Case: UNISIM - Constrained Clustering
% This case applies the Peixoto's algorithms for slope constrained clustering
%
%
% REMARK
% 
% All the algorithms of the family slopeConstrainedClusteringX, except 
% X = 6N were verified to be not suitable for clustering because they are
% totally dependend on a initial condition. Such 'initial condition' is due
% to the choice of the index in the global (sample) point list. Since the
% iterator can be placed on every point, this choice becomes somewhat
% arbitrary, so that the partitioning outcome is not unique. 
%
% We have decided to abort the constrained clustering approaches because 
% they were not giving reasonable results. The '6N' one is the unique that
% still might have some sense.  

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
G = computeGeometry(G);

%% Compute required parameters
P = computeParams(G,PROPS);

%% Constrained clustering

% log10(phiz) x log10(RQI); normalized averaging technique 
Log10PHIZ = P.Log10PHIZ(:);     
Log10RQIN = P.Log10RQIN(:);

%{
% TODO
% filtering
%globIndNoZero = find(Log10PHIZ~=0); 
%Log10PHIZ = Log10PHIZ(globIndNoZero);
%Log10RQIN = Log10RQIN(globIndNoZero);
%}

% constraint slope
seps = 1e-1;

% clustering method 
method = setSCCMethod('6n'); % --> READ REMARK <--

switch method

    case 'c'        
        tstart = tic;
        [partitioning,R2,M,B] = slopeConstrainedClusteringC(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCC.mat';                

    case 'ct'
        tstart = tic;
        [partitioning,R2,M,B] = slopeConstrainedClusteringCT(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCT.mat';                
        
    case 'py'
        tstart = tic;
        [partitioning,R2,M,B] = slopeConstrainedClusteringPY(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCPY.mat';                
    
    case 'nl'
        tstart = tic;
        [partitioning,R2,M,B] = slopeConstrainedClusteringNL(Log10PHIZ,Log10RQIN,1,seps);
        telapsed = toc(tstart); 
        fname = '../mat/SCCNL.mat';                
    
    case '6n'
        tstart = tic;
        [partitioning,R2,M,B] = slopeConstrainedClustering6N(Log10PHIZ,Log10RQIN,1,seps,G);
        telapsed = toc(tstart); 
        fname = '../mat/SCC6N.mat';                
        
end

% saving clustering structure
SCC.method = method;
SCC.partitioning = partitioning;
SCC.slope = M;
SCC.R2 = R2;
SCC.B = B;
SCC.execTime = telapsed;

save(fname,'SCC');
