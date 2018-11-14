%% Case: UNISIM - Constrained Clustering Analyzer
% This case analyses outcomes from the Peixoto's algorithms for slope constrained clustering 

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
G = computeGeometry(G);

%% Compute required parameters
P = computeParams(G,PROPS);

%% Constrained clustering

% log10(phiz) x log10(RQI); normalized averaging technique 
Log10PHIZ = P.Log10PHIZ(:); 
Log10RQIN = P.Log10RQIN(:);

%% Load pre-computed clusters 

[SCCC,SCCT,SCCPY,SCCNL,SCC6N] = loadClusterSCC;

% get clusters with minel or more elements
minel = 10;

getidx = @(scc) cellfun(@numel,scc.clustering);
[idx_c,idx_ct,idx_py,idx_nl,idx_6n] = ... 
deal(getidx(SCCC),getidx(SCCT),getidx(SCCPY),getidx(SCCNL),getidx(SCC6N));

getel = @(idx) idx(idx > minel);
[idx_c,idx_ct,idx_py,idx_nl,idx_6n] = ... 
deal(getel(idx_c),getel(idx_ct),getel(idx_py),getel(idx_nl),getel(idx_6n));

