
% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

% Compute required parameters
P = computeParams(G,PROPS);

% Constrained clustering
Log10PHIZ = P.Log10PHIZ(:);
Log10RQIN = P.Log10RQIN(:);

tstart = tic;
[Clusters,R2,m] = slopeConstrainedClustering(Log10PHIZ,Log10RQIN,1,0.05);
telapsed = toc(tstart); 
