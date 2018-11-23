%% Case: UNISIM - Constrained Clustering Analyzer
% This case analyses outcomes from the Peixoto's algorithms 
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

% consider partitions with >= minel elements
minel = 10;

% limit analysis for clusters with >= minc members
minc = 5;

%% Reload data
fis = {'../mat/SCCC_connections.mat','../mat/SCCT_connections.mat',...
       '../mat/SCCPY_connections.mat','../mat/SCCNL_connections.mat'};
% if available. 
if unique(cellfun(@exist,fis)) == 2
    load('../mat/SCCC_connections.mat','SCCC_connections');
    load('../mat/SCCT_connections.mat','SCCT_connections');
    load('../mat/SCCPY_connections.mat','SCCPY_connections');
    load('../mat/SCCNL_connections.mat','SCCNL_connections');
    compute = false;
else
    compute = true;
end


%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
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

getidx = @(scc) cellfun(@numel,scc.clustering);
[idx_c,idx_ct,idx_py,idx_nl] = ... 
deal(getidx(SCCC),getidx(SCCT),getidx(SCCPY),getidx(SCCNL));

getel = @(idx) idx(idx > minel);
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

% if necessary to recompute
if compute == true
    
    % loop clustering method
    for sc = 1:4        

        % loop partitioning
        for part = IDX{sc}        
            [I,J,K] = ind2sub(G.cartDims,SCC{sc}.clustering{part}'); 
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

%% Plot 3D info

% Extract subgrid from the global grid based on the point indices of the
% connected components found by the SCC method  
% and plot property over the subgrid (e.g. FZI)
%
% EXAMPLE:
% Method: SCCC; partition: 1; connected component: 1
%{
Gi = extractSubgrid(G,Ind(SCCC_connections.connections{1}.globalCompVoxelInds{1})); % need inverse grid mapping
plotCellData(Gi,P.FZIN(SCCC_connections.connections{1}.globalCompVoxelInds{1})) % do not need grid mapping
%}


%% CLUSTER INTERSECTION ANALYSIS

% Check if there are intersections of clusters detected by each SCC method. 

% number of partitions
npc  = numel(idx_c);
npct = numel(idx_ct);
nppy = numel(idx_py);
npnl = numel(idx_nl);

% mark which method do has clusters to be analysed
goto_analysis = true(1,4);

% 'sccX_pc' are arrays storing, per partition, the indices of all the
% clusters whose number of elements is >= minc for the method 'X'.

% SCCC
sccc_pc = cell(npc,1);
for p = 1:npc
    
    if ~isempty(SCCC_connections.connections{p})
        aux = SCCC_connections.connections{p}.compSizes;        
        aux = find(aux >= minc);
        if ~isempty(aux)
            sccc_pc{p} = aux;
        end
    end                
end
% checking
if all(cellfun(@isempty,sccc_pc))
    goto_analysis(1) = false;    
else
   % store useful parts and associated clusters
   aux = find(~cellfun(@isempty,sccc_pc));
   aux2 = cell(numel(aux),2);
   for i = 1:numel(aux)       
       aux2{i,1} = aux(i);
       aux2{i,2} = sccc_pc{aux(i)};
   end
      
   sccc_pc = aux2;
end

% SCCT
scct_pc = cell(npct,1);
for p = 1:npct
    
    if ~isempty(SCCT_connections.connections{p})
        aux = SCCT_connections.connections{p}.compSizes;        
        aux = find(aux >= minc);
        if ~isempty(aux)
            scct_pc{p} = aux;
        end
    end                
end

% checking
if all(cellfun(@isempty,scct_pc))
    goto_analysis(2) = false;    
else
    % store useful parts and associated clusters
   aux = find(~cellfun(@isempty,scct_pc));
   aux2 = cell(numel(aux),2);
   for i = 1:numel(aux)
       aux2{i,1} = aux(i);
       aux2{i,2} = scct_pc{aux(i)};
   end
      
   scct_pc = aux2;
    
end

% SCCPY
sccpy_pc = cell(nppy,1);
for p = 1:nppy
    
    if ~isempty(SCCPY_connections.connections{p})
        aux = SCCPY_connections.connections{p}.compSizes;        
        aux = find(aux >= minc);
        if ~isempty(aux)
            sccpy_pc{p} = aux;
        end
    end                
end

% checking
if all(cellfun(@isempty,sccpy_pc))
    goto_analysis(3) = false;   
else
    % store useful parts and associated clusters
   aux = find(~cellfun(@isempty,sccpy_pc));
   aux2 = cell(numel(aux),2);
   for i = 1:numel(aux)
       aux2{i,1} = aux(i);
       aux2{i,2} = sccpy_pc{aux(i)};
   end
      
   sccpy_pc = aux2;
end

% SCCNL
sccnl_pc = cell(npnl,1);
for p = 1:npnl
    
    if ~isempty(SCCNL_connections.connections{p})
        aux = SCCNL_connections.connections{p}.compSizes;        
        aux = find(aux >= minc);
        if ~isempty(aux)
            sccnl_pc{p} = aux;
        end
    end                
end

% checking
if all(cellfun(@isempty,sccnl_pc))
    goto_analysis(4) = false;  
else
    % store useful parts and associated clusters
   aux = find(~cellfun(@isempty,sccnl_pc));
   aux2 = cell(numel(aux),2);
   for i = 1:numel(aux)
       aux2{i,1} = aux(i);
       aux2{i,2} = sccnl_pc{aux(i)};
   end
      
   sccnl_pc = aux2;
end

TO_ANALYSIS = {sccc_pc,scct_pc,sccpy_pc,sccnl_pc};
CONNS = {SCCC_connections,SCCT_connections,...
         SCCPY_connections,SCCNL_connections};
toa = find(goto_analysis);

% filter
stay = numel(toa);
for i = 1:stay
    fprintf('---> Keeping %s clustering to study.\n',codename{toa(i)});
    TO_STUDY.approach(i).name = codename{toa(i)};
    TO_STUDY.approach(i).partitionsID = TO_ANALYSIS{toa(i)}(:,1);
    TO_STUDY.approach(i).clustersID = TO_ANALYSIS{toa(i)}(:,2);
    TO_STUDY.approach(i).All = CONNS{toa(i)};
end

%{ 
% \TODO

TO_ANALYSIS = {sccc_pc,scct_pc,sccpy_pc,sccnl_pc};
CONNS = {SCCC_connections,SCCT_connections,...
         SCCPY_connections,SCCNL_connections};
toa = find(goto_analysis);

INTERSECTIONS = TO_ANALYSIS;

for i = toa
    for j = toa
        if i ~= j
            
            aux = [];
            
            Pi = TO_ANALYSIS{i}(:,1);
            Pj = TO_ANALYSIS{j}(:,1);
            
            np = length(Pi); nq = length(Pj); 
            
            for p = 1:np                                
                for q = 1:nq
                    
                    ip = cell2mat(Pi(p));
                    iq = cell2mat(Pj(q));
                    
                    Cp = TO_ANALYSIS{ip}(:,2);                                                                
                    Cq = TO_ANALYSIS{iq}(:,2);                                                                                                
                    
                    Clistp = cell2mat(Cp(p));
                    Clistq = cell2mat(Cq(q));
                    
                    nc = numel(Clistp); nd = numel(Clistq);
                                       
                    for c = 1:nc
                        for d = 1:nd
                            
                            cc = Clistp(c); 
                            cd = Clistq(d); 
                            
                            Ci = CONNS{i}.connections{ip}.globalCompVoxelInds{cc};
                            Cj = CONNS{j}.connections{iq}.globalCompVoxelInds{cd};
                            
                            if numel(Ci) == numel(Cj) & Ci == Cj
                                aux = [aux; [j,iq,cd]];
                                fprintf('%d %d %d\n',j,iq,cd);                                
                            end
                            
                        end
                    end
                    
                end
            end
            
        end
    end
end
%}    

%% 3D ANALYSIS FOR MAJOR CLUSTERS

napp = length(TO_STUDY.approach);

for a = 2%:napp
    
    Pa = cell2mat(TO_STUDY.approach(a).partitionsID);
       
    Ca = TO_STUDY.approach(a).clustersID;
    
    figure
    plotGrid(G, 1:G.cells.num,'FaceColor',[0.6,0.6,0.6], ...
    'FaceAlpha',0.05, 'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.)
    hold on 
    
    for p = 1:length(Pa)
        cvi = TO_STUDY.approach(a).All.connections{p}.globalCompVoxelInds{Ca{p}(1)};
        glob = Ind(cvi);
        globn = glob(~isnan(glob));
        Gi = extractSubgrid(G,globn);        
        plotCellData(Gi,P.FZIN(cvi));
        %Gi = extractSubgrid(G,Ind(cvi));        
        %plotCellData(Gi,P.FZIN(cvi));
        
    end
    

end
 





