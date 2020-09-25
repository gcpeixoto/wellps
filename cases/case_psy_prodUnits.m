%% Productivity Units 
% Script to study 3D productivity units based on proxy 
% functions of productivity potential

%
% To run this case, the MRST's function 'readGRDECL' was modified to
% include the keyword 'SO', which is a user-defined keyword for oil
% saturation.
%
% ---> MODIFICATION TO BE MADE IN 'readGRDECL' (beetween ** .. **)
%
%           case {'PORO',                         ...
%                'PERMX' , 'PERMXY', 'PERMXZ',   ...
%                'PERMYX', 'PERMY' , 'PERMYZ',   ...
%                'PERMZX', 'PERMZY', 'PERMZ' ,   ...
%                'PERMH',                        ...
%                'ACTNUM', 'SATNUM', 'ROCKTYPE', ...
%                'MULTX' , 'MULTX-',             ...
%                'MULTY' , 'MULTY-',             ...
%                'MULTZ' , 'MULTZ-',             ...
%                'NTG'   , 'VSH'   , **'SO'**    ...
%                }
%
% Moreover, for this case, we used a the file 'PSY.DATA'
% which has an INCLUDE call to the file 'SO.INC' for the distribution 
% of oil saturation obtained by the standard relative permeability curves 
% of the model. 


case_name = 'psy-prodUnits';  

%% Loading grid 

d = DirManager(); 

f = fullfile(d.getBenchMarksDir,'psy','eclipse','PSY.grdecl');

% Unable to call buildModel because SO is a nonstandard field
G = readGRDECL(f);
PROPS.PHI = G.PORO;
PROPS.KX  = G.PERMX;
PROPS.KY  = G.PERMY;
PROPS.KZ  = G.PERMZ;
PROPS.ACTNUM  = G.ACTNUM;
PROPS.SO  = G.SO;
G = processGRDECL(G);
G = computeGeometry(G);

% active cells
on = find(PROPS.ACTNUM == true);

%% Compute parameters 
P = computeParams(G,PROPS); 

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

%% RQI x SO

% This section only deals with only the local grid. 
% We use the code here to plot the fields of the proxy
% function RQI x SO for visualization.

% normalized RQI
RQI_normAll = P.RQIN./max(P.RQIN);
RQI_norm = RQI_normAll(on);
ir = find(RQI_norm > 0);
RQI_norm_pos = RQI_norm(ir);

% oil sat
SO = PROPS.SO(on);
io = find(SO > 0);
SO_pos = SO(io);

% proxy: RQInorm x oil sat
RQI_SO_norm = RQI_norm(:).*SO(:);   
in = find(RQI_SO_norm > 0);
RQI_SO_norm = RQI_SO_norm(in);

% plot proxy
%Gn = extractSubgrid(G,in);
%plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
%axis vis3d off
%colormap jet

% plot RQI_norm > 0
%figure
%Gr = extractSubgrid(G,ir);
%plotCellData(Gr,RQI_norm_pos)

% plot SO > 0
%figure
%Go = extractSubgrid(G,io);
%plotCellData(Go,SO_pos)


%% Potential ranges

% We are defining a new concept based on Productivity Units (PUs)
% that are computed on the basis of a discrete function
% that defines thresholds for the proxy function. 
% This function is called here PUC (productivity unit class).
% A priori, 4 classes are expected, but only 3 really matter.
% 
% For a given reservoir cell w and proxy function F (in this case,
% F = RQI_norm(k,phi) x SO(t)), where RQI_norm(k,phi) is the 
% normalized RQI. This way, F : w -> [0,1]
%
%          { 0, if F(w) < d*
% PUC(w) = { 1, if d1 <= F(w) <  d2
%          { 2, if d2 <= F(w) <  d3
%          { 3, if d3 <= F(w) <= 1.0
%
% where d* is the average value of F for all w and d2, d3 are
% thresholds based on fixed steps:
%
% range:  [ 0 --- d* --- d2 --- d3 --- 1 ]
% class:  |   0   |   1   |   2  |   3   |   
%
%
% Therewith, we define a PU of class X as each connected cluster 
% whose all cells have PUC(w) = X.

min = mean(RQI_SO_norm); % lim inf by mean
div = linspace(min,1,4); % divisions for 3 classes

% productivity unit class
so = reshape(PROPS.SO,G.cartDims);
RQIso = RQI_normAll.*so;

PUC = RQIso;

c0 = PUC < div(1); % class 0 (ignored)
c1 = div(1) <= PUC & PUC <  div(2); % class 1
c2 = div(2) <= PUC & PUC <  div(3); % class 2
c3 = div(3) <= PUC & PUC <= div(4); % class 3

PUC(c0) = 0;
PUC(c1) = 1;
PUC(c2) = 2;
PUC(c3) = 3;

%% Productivity Units 

% compute PUs by PUC
nofsc = 10; % only for .csv
puc = 1:3; % only 3 classes matter
pucSt = findConnectionsByPUC(d,puc,PUC,nofsc,'y',1);


% minimum number of cells to consider to get clusters
nofc = 10; 

pucx = fieldnames(pucSt);
nn = length(pucx);

%% Get info for PU clusters 
% The struct PUInfo stores information of clusters per PUC.

for f = 1:nn
    
    % id of clusters per PUC with >= nofc elements
    idsel = find ( cellfun( @numel, pucSt.(pucx{f}).compVoxelInds) >= nofc );
    
    if ~isempty(idsel)
        
        for C = idsel
        
            PUInfo.(pucx{f}).(strcat('C',num2str(C))).value = ...
            pucSt.(pucx{f}).value;
                                    
            PUInfo.(pucx{f}).(strcat('C',num2str(C))).compVoxelInds = ...
            pucSt.(pucx{f}).compVoxelInds{C};
        
            PUInfo.(pucx{f}).(strcat('C',num2str(C))).mappedCompVoxelInds = ...
            Ind( pucSt.(pucx{f}).compVoxelInds{C} ); % mapped inds
            
            PUInfo.(pucx{f}).(strcat('C',num2str(C))).compVoxelCoords = ...
            pucSt.(pucx{f}).compVoxelCoords{C};                                    
                                    
        end
    end
end

%% Computation of closeness centrality 
% In this section, we compute the closeness centrality for each 
% productivity unit (cluster) and export tailored files to setup 
% IMEX perforation tables.

for f = 1:nn % loop PUC
    
    fn = fieldnames(PUInfo.(pucx{f}));
    
    for fc = 1:length(fn) % loop clusters
                            
        % metrics 
        cvi = PUInfo.(pucx{f}).(fn{fc}).compVoxelInds;            
        cvc = PUInfo.(pucx{f}).(fn{fc}).compVoxelCoords;
                                
        indIJ = [];
        for i = 1:numel(cvi)
            for j = 1:numel(cvi)
                  if i ~= j % skipping null distance                   
                      dist = sqrt( ( cvc(i,1) - cvc(j,1) )^2 + ...
                                   ( cvc(i,2) - cvc(j,2) )^2 + ...
                                   ( cvc(i,3) - cvc(j,3) )^2 ); 

                      % detecting neighbour voxels                  
                      if dist <= 1     % connectivity criterion
                          indIJ = [ indIJ; [ i j ] ];
                          edgeList = indIJ;                      
                      end                  
                  end
             end
        end   

        % invoke graphmetrics
        aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
        indIJ = [ indIJ; aux ]; % filling       

        % !!! This 0.5 was placed here because with 1, the sparse 
        % matrix was returning entries = 2. Why is that??
        Madj = sparse( indIJ(:,1),indIJ(:,2),0.5,size(cvc,1),size(cvc,1) ); 

        edfile = saveAdjEdges(Madj);              
        outfile = fullfile(d.getTmpDir,'metrics.txt');
        exec = sprintf('graphMetrics %s %s',edfile,outfile);

        % call SNAP 
        fprintf('%s\n',repmat('=',[1,75]))
        fprintf('%sWELLPS - SNAP Interface :: running C++\n',repmat(' ',[1,20]));
        fprintf('%s\n',repmat('=',[1,75]))
        system( fullfile(d.getCppDir,exec) );

        tab = importdata(outfile);                       
        Mtab = tab.data(:,1:4);            
        nodeID = Mtab(:,1);

        % centralities            
        cln  = Mtab(:,3);            

        % delete the temporary files
        delete(outfile)                                    

        % get rid of input file
        delete(edfile);

        maxC = max(cln);        % max closeness = min farness
        iC = cln == maxC;       % network closer nodes
        iCnode = nodeID(iC);    % getting node id (there might be more than 1)
        ivC = cvc(iCnode,: );   

        % save info
        PUMetrics.(pucx{f}).(fn{fc}).maxClosenessVoxelCoords = ivC(1,:);
        PUMetrics.(pucx{f}).(fn{fc}).maxClosenessValue = maxC;
                
        % max closeness coords        
        mcx = ivC(1,1); mcy = ivC(1,2); mcz = ivC(1,3); 
                
        % maxC same-column neighbors (inclusive)
        zwide = unique(cvc(:,3)); 
        lzwide = numel(zwide);
        mcols = [ones(lzwide,1)*mcx, ones(lzwide,1)*mcy, zwide]; 
        
        % maxC + neighbors at cluster z-wide domain
        PUMetrics.(pucx{f}).(fn{fc}).maxClosenessColNeighsZ = mcols;

        aux = cvi(iCnode); aux = aux(1); % gets only first of the list
        PUMetrics.(pucx{f}).(fn{fc}).maxClosenessVoxelInd = aux;

        aux = Ind(cvi(iCnode)); aux = aux(1); % gets only first of the list
        PUMetrics.(pucx{f}).(fn{fc}).maxClosenessVoxelMappedInd = aux;                

        %save perforation table for cluster MaxC       
        ptset.savedir = fullfile(d.getTmpDir,'PSY','MaxC',pucx{f},fn{fc});
        ptset.wellname = 'W';
        ptset.geometry = 'K';
        ptset.perfs = PUMetrics.(pucx{f}).(fn{fc}).maxClosenessVoxelCoords;
        ptset = savePerfTable(ptset);
        
        %save perforation table for MaxC + column neighbors
        qtset.savedir = fullfile(d.getTmpDir,'PSY','MaxC-Column',pucx{f},fn{fc});
        qtset.wellname = 'W';
        qtset.geometry = 'K';
        qtset.perfs = PUMetrics.(pucx{f}).(fn{fc}).maxClosenessColNeighsZ;
        qtset = savePerfTable(qtset);
        
        clear ptset qtset
        
    end

end


%% Plot proxy
figure
Gn = extractSubgrid(G,in);
plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
axis off; colormap jet; colorbar
title('$RQI_n \times s_o(t)$','interpreter','latex')

%% Plot Clusters 
figure
plotGrid(G,'FaceAlpha',0.1,'FaceColor','k','EdgeColor','None'); 
colormap(jet); axis off; colorbar

% mask
mask = false(G.cells.num,1);

% example values taken from PUInfo
exi = 1; % exi-th PUC value in the list

% exic-th cluster in the PUC list
for exic = 1:length(fieldnames(PUInfo.(pucx{exi}))) 
    
    exPUC = fieldnames(PUInfo); 
    exPUC = exPUC{exi};
    exC = fieldnames(PUInfo.(exPUC)); 
    exC = exC{exic};
    excvi = PUInfo.(exPUC).(exC).compVoxelInds; % global
    excvim = PUInfo.(exPUC).(exC).mappedCompVoxelInds; % local
    exRQIso = RQIso(excvi); % proxy computed at the global index

    % background + cluster domain
    plotCellData(extractSubgrid(G,excvim),exRQIso)
    
    title(strcat('Productivity Units - PUC ',num2str(exi)))
    % max clo (cluster)
    indmax = PUMetrics.(exPUC).(exC).maxClosenessVoxelMappedInd;
    mask(indmax) = true;
    plotGrid(G,mask,'FaceColor','k','EdgeColor','w' );
end


