%% CASE STUDY: PRODUCTIVITY UNITS
% Script to study 3D productivity units based on proxy 
% functions of productivity potential
%
% REMARK: To run this case, the MRST's function 'readGRDECL' was modified to
% include the keywords 'SO' and 'PRESSURE', which are user-defined 
% keywords to get oil saturation and initial pressure.
%
% ---> MODIFICATION TO BE MADE IN 'readGRDECL' (between ** .. **)
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
% Moreover, for this case, we used the file 'PSY.DATA'
% which has an INCLUDE call to the file 'SO.INC' for the distribution 
% of oil saturation obtained by the standard relative permeability curves 
% of the model. 

case_name = 'psy-prodUnits';  

%% Loading grid 

d = DirManager(); 

f = fullfile(d.getBenchMarksDir,'psy','eclipse','PSY.grdecl');

%% Productivity Potential Index 
% Compute PPI through one of the available methods:
m = {'rqi', 'rqip', 'kharghoria'};

% bottom hole pressure (fixed)
bhp = 200;
[J,G,PROPS,active] = computeProdProxy(f,m{2},bhp);


Ja = J(active); % J only at active cells
ip = find(Ja > 0); % positive values
Jplus = Ja(ip);
plotCellData(G,Ja)
Gn = extractSubgrid(G,ip);
plotCellData(Gn,Jplus,'FaceAlpha',1.0,'EdgeColor','none')
axis vis3d off equal
colormap(jet(35)), colorbar


%% Productivity unit classes

% We are defining a new concept based on Productivity Units (PUs)
% that are computed on the basis of a discrete function
% that defines thresholds for the proxy function. 
% This function is called here PUC (productivity unit class).

% For a given reservoir cell w and proxy function J, J : w -> [0,1]. Then
% we define
%
%          { 1, if 0.0 < F(w) <  d1
% PUC(w) = { 2, if d1 <= F(w) <  d2
%          { 3, if d2 <= F(w) <  d3
%          { ...
%          { n, if dn <= F(w) <= max(F)
%
%
% Therewith, we define a PU of class X as each connected cluster 
% whose all cells have PUC(w) = X.

% Below, there is a summary of tests performed for PSY model 
% comparing the number of PUC classes obtained per histogram 
% binning method. As seen, the Sturges method returned
% the least number of classes (15), i.e. lesser refinement.
% On the other hand, 'sqrt' returned 148 classes, which is too much.
%    SUMMARY OF TESTS 
% |----------------------|
% | method    | nclasses |
% |----------------------|
% |'auto'     | 45       |
% |'scott'    | 45       |
% |'fd'       | 45       |
% |'sturges'  | 15       | <-- Sturges returned the least refinement
% |'sqrt'     | 148      |
% |'shimazaki'| 64       |
% |----------------------|
%
%
% Now, I am reporting a sensitivity analysis summary with
% tests performed varying the number of PUCs and nofsc.
% A D-tuple entry indicates the number of clusters obtained per D value.
% with fixed NOFSC.
% 
% - D values: 1-10
% - NOFSC values: 5,10,20,35,50
%
%
% –––––––––––––––– PPI function :: J = 'rqi' –––––––––––––––– 
%
% NOFSC = 5
% –––––––––
% D  | #clusters
% 1  |  17 
% 2  | (334,358) 
% 3  | (24,571,453)  
% 4  | (2,246,332,254)
% 5  | (0,62,201,130,151)
% 6  | (0,9,102,116,46,82)
% 7  | (0,1,47,42,67,21,57)
% 8  | (0,0,13,35,44,34,8,38)
% 9  | (0,0,4,23,13,30,10,7,30)
% 10 | (0,0,0,16,14,13,16,5,2,25)
%
%
% NOFSC = 10
% –––––––––
% D  | #clusters
% 1  |  4 
% 2  | (91,167) 
% 3  | (0,216,140)  
% 4  | (0,40,67,39)
% 5  | (0,0,18,15,13)
% 6  | (0,0,8,9,3,8)
% 7  | (0,0,2,2,6,0,5)
% 8  | (0,0,0,1,1,1,0,2)
% 9  | (0,0,0,2,1,0,1,0,1)
% 10 | (0,0,0,0,0,0,1,1,0,1)
%
%
% NOFSC = 20
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (18,66) 
% 3  | (0,57,25)  
% 4  | (0,4,6,0)
% 5  | (0,0,0,1,0)
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
%
%
% NOFSC = 35
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (3,38) 
% 3  | (0,17,4)  
% 4  | (0)*4
% 5  | (0)*5
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
%
%
% NOFSC = 50
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (0,25) 
% 3  | (0,5,2)  
% 4  | (0)*4
% 5  | (0)*5
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%
% –––––––––––––––– PPI function :: J = 'rqip' –––––––––––––––– 
% 
% NOFSC = 5
% –––––––––
% D  | #clusters
% 1  |  17 
% 2  | (400,430) 
% 3  | (43,570,382)  
% 4  | (2,303,337,171)
% 5  | (0,96,240,129,87)
% 6  | (0,19,120,136,45,49)
% 7  | (0,2,72,87,55,18,26)
% 8  | (0,0,29,35,58,25,8,20)
% 9  | (0,0,9,26,21,33,9,3,11)
% 10 | (0,0,4,19,12,28,14,7,1,9)
%
%
% NOFSC = 10
% –––––––––
% D  | #clusters
% 1  |  4 
% 2  | (119,210) 
% 3  | (3,217,101)  
% 4  | (0,56,66,19)
% 5  | (0,3,32,13,6)
% 6  | (0,0,4,11,4,3)
% 7  | (0,0,3,4,3,2,2)
% 8  | (0,0,0,0,0,1,0,1)
% 9  | (0,0,0,2,1,2,0,0,1)
% 10 | (0,0,0,1,0,0,1,0,0,1)
%
%
% NOFSC = 20
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (29,92) 
% 3  | (0,73,15)  
% 4  | (0,7,9,0)
% 5  | (0)*5
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
%
%
% NOFSC = 35
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (4,52) 
% 3  | (0,20,1)  
% 4  | (0,1,1,0)
% 5  | (0)*5
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
%
%
% NOFSC = 50
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (2,33) 
% 3  | (0,12,0)  
% 4  | (0)*4
% 5  | (0)*5
% 6  | (0)*6
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%
% –––––––––––– PPI function :: J = 'kharghoria' ––––––––––––––
% NOFSC = 5
% –––––––––
% D  | #clusters
% 1  |  17 
% 2  | (93,27) 
% 3  | (252,88,1)  
% 4  | (451,88,21,0)
% 5  | (587,94,13,3,0)
% 6  | (553,115,11,3,1)
% 7  | (458,116,8,2,6,0,0)
% 8  | (367,126,11,1,0,1,0,0)
% 9  | (313,100,7,3,0,1,1,0,0)
% 10 | (254,90,15,4,0,0,0,0,0,0)
%
%
% NOFSC = 10
% –––––––––
% D  | #clusters
% 1  |  4 
% 2  | (29,1) 
% 3  | (97,6,0)  
% 4  | (210,6,0,0)
% 5  | (235,5,0,0,0)
% 6  | (185,13,0,0,0,0)
% 7  | (130,11,1,0,0,0,0)
% 8  | (101,9,0,0,0,0,0,0)
% 9  | (71,1,0,0,0,0,0,0,0)
% 10 | (49,2,0,0,0,0,0,0,0,0)
%
%
% NOFSC = 20
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (6,0) 
% 3  | (33,0,0)  
% 4  | (90,0,0)
% 5  | (87,0,0,0,0)
% 6  | (53,1,0,0,0,0)
% 7  | (28,0,0,0,0,0,0)
% 8  | (17,0,0,0,0,0,0,0)
% 9  | (8,0,0,0,0,0,0,0,0)
% 10 | (3,0,0,0,0,0,0,0,0)
%
%
% NOFSC = 35
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (1,0) 
% 3  | (18,0,0)  
% 4  | (58,0,0,0)  
% 5  | (35,0,0,0,0) 
% 6  | (18,0,0,0,0,0)
% 7  | (7,0,0,0,0,0,0)
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
%
%
% NOFSC = 50
% –––––––––
% D  | #clusters
% 1  |  1 
% 2  | (1,0) 
% 3  | (10,0,0)  
% 4  | (42,0,0,0) 
% 5  | (21,0,0,0,0) 
% 6  | (6,0,0,0,0) 
% 7  | (0)*7
% 8  | (0)*8
% 9  | (0)*9
% 10 | (0)*10
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%
% compute discrete function PUC
[PUC,nclasses,delta,div] = computePUC(J,active,'10');

% list of all PUC values except zero
puc = 1:nclasses;

% 6-neighbor connectivity criterion to find clusters from the PUC values
nofsc = 50; % only for .csv
pucSt = findConnectionsByPUC(d,1:nclasses,PUC,nofsc,'y',1);

%% minimum number of cells to consider to get clusters
nofc = nofsc; 
pucx = fieldnames(pucSt);
nn = length(pucx);

%% Get info for PU clusters 
% The struct PUInfo stores information of clusters per PUC.

% inverse mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

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
        
        %{
        % --- call SNAP 
        exec = sprintf('graphMetrics %s %s',edfile,outfile);
        fprintf('%s\n',repmat('=',[1,75]))
        fprintf('%sWELLPS - SNAP Interface :: running C++\n',repmat(' ',[1,20]));
        fprintf('%s\n',repmat('=',[1,75]))
        system( fullfile(d.getCppDir,exec) );        
        % get data
        tab = importdata(outfile);                       
        Mtab = tab.data(:,1:4);            
        nodeID = Mtab(:,1);
        % centralities            
        cln  = Mtab(:,3); % closeness      
        bet  = Mtab(:,4); % betweeness
        % delete the temporary files
        delete(outfile)  
        %}
                   
        
        % --- call NETWORKX
        fprintf('%s\n',repmat('=',[1,75]))
        fprintf('%sWELLPS - NETWORKX Interface :: running Python\n',repmat(' ',[1,20]));
        fprintf('%s\n',repmat('=',[1,75]))              
        exec = [d.getPyExec,' ',d.getPyDir,filesep, ...
                sprintf('graphMetrics.py %s %s',edfile,outfile)];
        system(exec);        
        % get data
        tab = importdata(outfile);                       
        Mtab = tab.data(:,1:4);            
        nodeID = Mtab(:,1);        
        cln = Mtab(:,3); % closeness - networkx
        bet = Mtab(:,4); % betweeness - networkx        
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

%{
%% Plot proxy
figure
subplot(211)
plotGrid(G,'FaceAlpha',0.1,'FaceColor','k','EdgeColor','None'); 
Gn = extractSubgrid(G,in);
plotCellData(Gn,PROXY,'FaceAlpha',1.0,'EdgeColor','none')
lighting phong, camproj perspective
axis tight equal off, view(101,30)
subplot(212)
plotGrid(G,'FaceAlpha',0.1,'FaceColor','k','EdgeColor','None'); 
Gn = extractSubgrid(G,in);
plotCellData(Gn,PROXY,'FaceAlpha',1.0,'EdgeColor','none')
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
%}
