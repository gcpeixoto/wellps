%% Case: UNISIM - Compute clusters
% This case compute clusters from Peixoto's algorithms 
% for slope constrained clustering. 
%
% See Remark in wellps::case_unisim_allSCC.m. Part of this methodology was 
% aborted.
%
%
% METHODOLOGY: for our purposes, SCC algorithms are for partitioning.
%              Each of them has produced partitions. Then, we go through 
%              each 'partition'and look for 'clusters' (connected components). 
%              
%              A 'partition' is also, at least theoretically, the sought 
%              HFU, since it satisfies the slope constraint
%              but it is very scattered, pulverized all over the domain
%              with disconnected elements. 
%
%              The 'clusters', however, would be like subsets of HFUs, but 
%              with connections.                 


%% INPUTS

% plot info for partitions with >= minel elements
minel = 2;

%% Reload data
fis = {'../mat/SCCC_connections.mat','../mat/SCCT_connections.mat',...
       '../mat/SCCPY_connections.mat','../mat/SCCNL_connections.mat'};
% if available. 
if unique(cellfun(@exist,fis)) == 2
    load('../mat/SCCC_connections.mat','SCCC_connections');
    load('../mat/SCCT_connections.mat','SCCT_connections');
    load('../mat/SCCPY_connections.mat','SCCPY_connections');
    load('../mat/SCCNL_connections.mat','SCCNL_connections');
    
    fprintf('Nothing to compute. Stopping...\n')    
else
    

    %% Grid reading
    [G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
    G = computeGeometry(G);

    %% Compute required parameters
    P = computeParams(G,PROPS);

    %% Log data

    % log10(phiz) x log10(RQI); normalized averaging technique 
    Log10PHIZ = P.Log10PHIZ(:); 
    Log10RQIN = P.Log10RQIN(:);

    %% Mapping
    Ind = nan(prod(G.cartDims),1);
    Ind(G.cells.indexMap) = 1:G.cells.num;

    %% Load pre-computed partitionings 

    [SCCC,SCCT,SCCPY,SCCNL] = loadClusterSCC;
    SCC = {SCCC,SCCT,SCCPY,SCCNL};
    codename = {'SCCC','SCCT','SCCPY','SCCNL'};

    getidx = @(scc) cellfun(@numel,scc.partitioning);
    [idx_c,idx_ct,idx_py,idx_nl] = ... 
    deal(getidx(SCCC),getidx(SCCT),getidx(SCCPY),getidx(SCCNL));

    getel = @(idx) idx(idx >= minel);
    [idx_c,idx_ct,idx_py,idx_nl] = ... 
    deal(getel(idx_c),getel(idx_ct),getel(idx_py),getel(idx_nl));

    % Partitions are ordered, then we can get indices using length
    getidx = @(x) 1:length(x);
    [idx_c,idx_ct,idx_py,idx_nl] = ... 
    deal(getidx(idx_c),getidx(idx_ct),getidx(idx_py),getidx(idx_nl));

    % checks if there are empty clusters 
    IDX = {idx_c,idx_ct,idx_py,idx_nl};
    aux = cellfun(@isempty,IDX);
    if any(aux)
        methodi = find(aux);
        fprintf('---> No clusters for %s clustering with min. element: %d.\n',...
                codename{methodi},minel); 
    end

    %% Clusters
    % Find connections inside each partitioning for all SCC, except 6N

    % loop clustering method
    for sc = 1:4        

        % loop partitioning
        for part = IDX{sc}        
            [I,J,K] = ind2sub(G.cartDims,SCC{sc}.partitioning{part}'); 
            cvc = [I,J,K];
            aux = findConnectionsSimple(cvc);                        

            % there might be no clusters (connected components
            if ~isempty(aux)

                aux2 = cell(1,aux.ncomp);            

                % loop connected components
                for c = 1:aux.ncomp

                    % required to get indices corresponding in global grid
                    cvi = sub2ind(G.cartDims,aux.compVoxelCoords{c}(:,1), ...
                                         aux.compVoxelCoords{c}(:,2), ...
                                         aux.compVoxelCoords{c}(:,3)); 
                    aux2{c} = cvi;            
                end            

                % append 
                aux.globalCompVoxelInds = aux2;
                SCC{sc}.connections{part} = aux;                           

            else

                SCC{sc}.connections{part} = [];   

            end
            
            SCC{sc}.minel = minel; % stores this info
        end        
    end

    % Save dataset relative to connected clusters
    SCCC_connections  = SCC{1};
    SCCT_connections  = SCC{2};
    SCCPY_connections = SCC{3};
    SCCNL_connections = SCC{4};
    save('../mat/SCCC_connections.mat','SCCC_connections');
    save('../mat/SCCT_connections.mat','SCCT_connections');
    save('../mat/SCCPY_connections.mat','SCCPY_connections');
    save('../mat/SCCNL_connections.mat','SCCNL_connections');

end