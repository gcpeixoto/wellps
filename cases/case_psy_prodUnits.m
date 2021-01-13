%% Productivity Units 
% Script to study 3D productivity units based on proxy 
% functions of productivity potential

%
% To run this case, the MRST's function 'readGRDECL' was modified to
% include the keywords 'SO' and 'PRESSURE', which are user-defined 
% keywords to get oil saturation and initial pressure.
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
%                'NTG'   , 'VSH'   , **'SO'**, **'PRESSURE' ...
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
PROPS.PRESS0  = G.PRESSURE;
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
RQI_normAll = P.RQIN./max(P.RQIN(:));
RQI_norm = RQI_normAll(on);
ir = find(RQI_norm > 0);
RQI_norm_pos = RQI_norm(ir);

% normalized pressure
P_normAll = PROPS.PRESS0./max(PROPS.PRESS0);
P_norm = P_normAll(on);


% oil sat
SO = PROPS.SO(on);
io = find(SO > 0);
SO_pos = SO(io);

% productivity functions
proxy = 'rqisop';
switch proxy
    case 'rqiso'
        % proxy: RQInorm x oil sat
        tit = 'RQI x SO';
        RQI_SO_norm = RQI_norm(:).*SO(:);   
        in = find(RQI_SO_norm > 0);
        RQI_SO_norm = RQI_SO_norm(in);
        flag = '';
    case 'rqisop'
        % proxy: RQInorm x oil sat x pressure_norm
        tit = 'RQI x SO x P';
        RQI_SO_norm = RQI_norm(:).*SO(:).*P_norm(:);   
        in = find(RQI_SO_norm > 0);
        RQI_SO_norm = RQI_SO_norm(in);
        flag = '-P';
end

% plot proxy
figure
Gn = extractSubgrid(G,in);
plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
axis vis3d off
colormap jet
title(tit)

% plot RQI_norm > 0
%figure
%Gr = extractSubgrid(G,ir);
%plotCellData(Gr,RQI_norm_pos)
%title(tit)

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

% productivity unit class
so = reshape(PROPS.SO,G.cartDims);
RQIso = RQI_normAll.*so;
PUC = RQIso;

nclasses = 3;
div = linspace(min,1,nclasses+1); % divisions for nclasses

switch nclasses 
    case 3       
        c0 = PUC < div(1); % class 0 (ignored)
        c1 = div(1) <= PUC & PUC <  div(2); % class 1
        c2 = div(2) <= PUC & PUC <  div(3); % class 2
        c3 = div(3) <= PUC & PUC <= div(4); % class 3
        PUC(c0) = 0;
        PUC(c1) = 1;
        PUC(c2) = 2;
        PUC(c3) = 3;
        
        puc = 1:3; % only 3 classes matter
    case 2        
        c0 = PUC < div(1); % class 0 (ignored)
        c1 = div(1) <= PUC & PUC <  div(2); % class 1
        c2 = div(2) <= PUC & PUC <  div(3); % class 2        
        PUC(c0) = 0;
        PUC(c1) = 1;
        PUC(c2) = 2;        
        puc = 1:2; % only 2 classes matter
end

%% Productivity Units 

% compute PUs by PUC
nofsc = 50; % only for .csv
pucSt = findConnectionsByPUC(d,puc,PUC,nofsc,'y',1);

% minimum number of cells to consider to get clusters
nofc = 50; 

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
    else
        % exclude empty PUC from postoperation; 
        % REMARK: improve this step, because it assumes that if I do not
        % have for instance a PUC2, I won't have a PUC3. But this is not
        % perfectly true. 
        nn = nn - 1; 
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
        ptset.savedir = fullfile(d.getTmpDir,strcat('PSY',flag),'MaxC',pucx{f},fn{fc});
        ptset.wellname = 'W';
        ptset.geometry = 'K';
        ptset.perfs = PUMetrics.(pucx{f}).(fn{fc}).maxClosenessVoxelCoords;
        ptset = savePerfTable(ptset);
        
        %save perforation table for MaxC + column neighbors
        qtset.savedir = fullfile(d.getTmpDir,strcat('PSY',flag),'MaxC-Column',pucx{f},fn{fc});
        qtset.wellname = 'W';
        qtset.geometry = 'K';
        qtset.perfs = PUMetrics.(pucx{f}).(fn{fc}).maxClosenessColNeighsZ;
        qtset = savePerfTable(qtset);
        
        clear ptset qtset
        
    end

end


%% Plot proxy
figure
subplot(211)
plotGrid(G,'FaceAlpha',0.1,'FaceColor','k','EdgeColor','None'); 
Gn = extractSubgrid(G,in);
plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
lighting phong, camproj perspective
axis tight equal off, view(101,30)
subplot(212)
plotGrid(G,'FaceAlpha',0.1,'FaceColor','k','EdgeColor','None'); 
Gn = extractSubgrid(G,in);
plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
lighting phong, camproj perspective
axis tight equal off, view(90,0)
colormap(jet(35)); 
cbar = colorbar;
cbar.Box = 'off';
cbar.Location = 'southoutside';
cbar.Label.String = '$\mathcal{I}(c,t=0)$';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 14;
cbar.FontSize = 14;
cbar.TickLabelInterpreter = 'latex';
cbar.Position = [.25 .1 .5 .05];
print('ppi-field.eps','-r0','-depsc2')

%% Plot Clusters 
figure
Gn = extractSubgrid(G,in);
plotGrid(Gn,'FaceColor','none','EdgeColor',[0.5,0.5,0.5])
axis tight equal off, view(90,90)

% mask
mask = false(G.cells.num,1);

% example values taken from PUInfo
exi = 1; % exi-th PUC value in the list

% exic-th cluster in the PUC list
cf = 1.0;
for exic = [11,13,6]%1:length(fieldnames(PUInfo.(pucx{exi}))) 
    
    cf = cf -0.06;
    exPUC = fieldnames(PUInfo); 
    exPUC = exPUC{exi};
    exC = fieldnames(PUInfo.(exPUC)); 
    exC = exC{exic};
    excvi = PUInfo.(exPUC).(exC).compVoxelInds; % global
    excvim = PUInfo.(exPUC).(exC).mappedCompVoxelInds; % local
    exRQIso = RQIso(excvi); % proxy computed at the global index

    % background + cluster domain
    %plotCellData(extractSubgrid(G,excvim),exRQIso)
    plotGrid(extractSubgrid(G,excvim),'FaceColor',[cf,cf-0.2,0.0],'EdgeColor','k')
    
    %title(strcat('Productivity Units - PUC ',num2str(exi)))
    % max clo (cluster)
    indmax = PUMetrics.(exPUC).(exC).maxClosenessVoxelMappedInd;
    mask(indmax) = true;
    %plotGrid(G,mask,'FaceColor','k','EdgeColor','w' );
end

%% Plot Histogram of RQISOP

histogram(RQI_SO_norm,'Normalization','pdf')

