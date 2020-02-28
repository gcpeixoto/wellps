%% HFU sub-clustering 
% Based on oil-saturated cells
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
% Moreover, for this case, we used a duplicata file 'UNISIM_....SO.DATA'
% which has a INCLUDE call to the file 'SO.INC' for the distribution 
% of oil saturation obtained by the standard relative permeability curves 
% of the model. 

case_name = 'unisim1-subclustering';  

%% Loading grid 

d = DirManager(); 

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE_SO.DATA');

% Unable to call buildModel because SO is a nonstandard field
G = readGRDECL(f);
PROPS.PHI = G.PORO;
PROPS.KX  = G.PERMX;
PROPS.KY  = G.PERMY;
PROPS.KZ  = G.PERMZ;
PROPS.ACTNUM  = G.ACTNUM;
PROPS.SO  = G.SO;
G = processGRDECL(G);

%% Compute parameters 
P = computeParams(G,PROPS); 

%plotCellData(G,PROPS.SO(on))

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;


%% Oil-saturated cells 
idso = find(PROPS.SO > 0);
SOP = PROPS.SO(idso);

% local indices
idsoloc = Ind(idso);

subG = extractSubgrid(G,idsoloc);
%plotCellData(subG,SOP,'EdgeColor','None')

S = getStats(d,P,{'DRTH_LN'},'n');

%% DRT connections 

% get all DRTs 
drtlist = S{1}(:,1); 
drtlist = drtlist(drtlist > 0);
drtlist = 13;

% compute HFUs by DRT
nofsc = 10;
drtSt = findDRTConnections(d,drtlist, P, 'harmonic','ln',nofsc,'n', 1);


%% Computing graph metrics 
opt.nofsc = nofsc;    
opt.seps = 0.05;
opt.R2min = 0.9;

% compute
[metricsSt,linregrSt] = computeDRTGraphMetrics(opt,drtSt);


%% Load metrics/linregr files 

mfn = fieldnames(metricsSt);

% Create structs
for f = 1:length(mfn)
    load(metricsSt.(mfn{f}),'M');
    Metrics.(mfn{f}) = M;         
    load(linregrSt.(mfn{f}),'L');
    Linregr.(mfn{f}) = L;    
    if ~isempty(find(cell2mat(L.performance), 1)) % performance
        Performance.(mfn{f}) = find(cell2mat(L.performance));
    else
        Performance.(mfn{f}) = [];
    end
    clear M L 
end

%% Get info for HP clusters 

for f = 1:length(mfn)
    p = Performance.(mfn{f});
    
    if ~isempty(p)
        
        for C = p
        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).value = ...
            drtSt.(mfn{f}).value;
            
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).averaging = ...
            drtSt.(mfn{f}).averaging;
        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).logBase = ...
            drtSt.(mfn{f}).logBase;
            
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).compVoxelInds = ...
            drtSt.(mfn{f}).compVoxelInds{C};
        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).mappedCompVoxelInds = ...
            Ind( drtSt.(mfn{f}).compVoxelInds{C} ); % mapped inds
            
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).compVoxelCoords = ...
            drtSt.(mfn{f}).compVoxelCoords{C};
            
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).compRQI = ...
            drtSt.(mfn{f}).compRQI{C};
                                    
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).closenessCentrality = ...
            Metrics.(mfn{f}).closenessCentrality{C};
        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).maxClosenessVoxelCoords = ...
            Metrics.(mfn{f}).maxClosenessVoxelCoords{C};
                
            % locate oil-saturated indices
            cvim = Ind( drtSt.(mfn{f}).compVoxelInds{C} );            
            oilSatMarker = false(length(cvim),1); 
                                    
            for i = 1:length(idsoloc)
                id = cvim == idsoloc(i);                      
                id = find(id == 1); 
                oilSatMarker(id) = true; 
            end
                                            
            % oil saturation marker
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilSatMarker = oilSatMarker;
                        
            % oil saturation value
            cvi = drtSt.(mfn{f}).compVoxelInds{C};
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilSat = PROPS.SO(cvi);
            
            if all(oilSatMarker)
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusName = 'fully saturated';
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusValue = 1;
            elseif any(oilSatMarker)
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusName = 'partially saturated';
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusValue = 0.5;                
            else
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusName = 'not saturated';
                HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilStatusValue = 0;                
            end
            
            % gammaPotential = normalizedRQI x so (local)             
            rqi = drtSt.(mfn{f}).compRQI{C};            
            so = HPInfo.(mfn{f}).(strcat('C',num2str(C))).oilSat;
            gammaPot = (rqi/max(rqi)).*so;
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).gammaPotential = gammaPot;                        
            
            % voxel coords of maximum gammaPotential
            gmax = find(gammaPot == max(gammaPot));
            cvc = drtSt.(mfn{f}).compVoxelCoords{C};
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).maxGammaPotentialCoords = cvc(gmax,:);                        
            
            % voxel indices of maximum gammaPotential                        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).maxGammaPotentialGlobalInd = cvi(gmax);                        
            HPInfo.(mfn{f}).(strcat('C',num2str(C))).maxGammaPotentialMappedInd = cvim(gmax);                                                
                                    
        end
    end
end

%% Subclustering for fully and partially oil-saturated clusters

%======
% INCORRECT!!! 
hpifn = fieldnames(HPInfo);
for f = 1:length(hpifn)
    
    hpifnc = fieldnames(HPInfo.(hpifn{f}));
    
    for fc = 1:length(hpifnc)
        
        if strcmp( HPInfo.(hpifn{f}).(hpifnc{fc}).oilStatusName, 'partially saturated' )
            
            % metrics 
            cvi = HPInfo.(hpifn{f}).(hpifnc{fc}).compVoxelInds;            
            cvc = HPInfo.(hpifn{f}).(hpifnc{fc}).compVoxelCoords;
            mark = HPInfo.(hpifn{f}).(hpifnc{fc}).oilSatMarker;   
            cvc = cvc(mark,:);            
            cvi = cvi(mark);
            
            if size(cvc,1) > 1 % bypass clusters that have only 1 cell outside oil zone
            
                indIJ = [];
                for i = 1:size(cvc,1)        
                    for j = i:size(cvc,1)                                                            
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
                aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
                indIJ = [ indIJ; aux ]; % filling                                     
                Madj = sparse( indIJ(:,1),indIJ(:,2),1,size(cvc,1),size(cvc,1) ); 

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

                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelCoords = ivC(1,:);
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessValue = maxC;
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelInd = cvi(iCnode);
            end
            
        end
        
    end

end
%=======


%% Plotting example (for visibility)

%{
exDRT = fieldnames(HPInfo); exDRT = exDRT{1};
exC = fieldnames(HPInfo.(exDRT)); exC = exC{1};
excvi = HPInfo.(exDRT).(exC).compVoxelInds;
excvim = HPInfo.(exDRT).(exC).mappedCompVoxelInds;
exmarker = HPInfo.(exDRT).(exC).oilSatMarker;
exso = HPInfo.(exDRT).(exC).oilSat;
exgamma = HPInfo.(exDRT).(exC).gammaPotential;

plotGrid(G,'FaceAlpha',0.2,'FaceColor','k','EdgeColor','None'); % background
plotCellData(subG,SOP,'FaceAlpha',0.2,'EdgeColor','None') % oil zone
plotCellData(extractSubgrid(G,excvim(exmarker)),exso(exmarker),'FaceColor','c') % cluster saturated portion
plotCellData(extractSubgrid(G,excvim(~exmarker)),exso(~exmarker),'FaceColor','g') % cluster nonsaturated portion

plotCellData(extractSubgrid(G,excvim),exgamma) % cluster nonsaturated portion

mask = false(G.cells.num,1);
mask(SubM.(exDRT).(exC).maxClosenessVoxelInd) = true;
plotGrid(G,mask)
%}

