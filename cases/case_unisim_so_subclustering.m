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
%drtlist = 13;

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
% In this step, it is necessary to run WELLPS to compute metrics for the
% region non-saturated in oil. This piece of code is included in
% computeDRTGraphMetrics and it was slightly modified for the current
% purpose.

hpifn = fieldnames(HPInfo);

for f = 1:length(hpifn) % loop DRTs
    
    hpifnc = fieldnames(HPInfo.(hpifn{f}));
    
    for fc = 1:length(hpifnc) % loop clusters
        
        % subclustering required only for partially saturated clusters
        % since the fully saturated ones are already computed.
        if strcmp( HPInfo.(hpifn{f}).(hpifnc{fc}).oilStatusName, 'partially saturated' )
            
            % metrics 
            cvi = HPInfo.(hpifn{f}).(hpifnc{fc}).compVoxelInds;            
            cvc = HPInfo.(hpifn{f}).(hpifnc{fc}).compVoxelCoords;
            mark = HPInfo.(hpifn{f}).(hpifnc{fc}).oilSatMarker;  % marker functon
            auxind = find(mark); % get indices of the oil zone cells.                        
            
            % bypass clusters that have only 1 cell outside oil zone
            if size(cvc,1) > 1 
            
                indIJ = [];
                for i = 1:numel(auxind)
                    for j = numel(auxind)
                          if i ~= j % skipping null distance                   
                              dist = sqrt( ( cvc(auxind(i),1) - cvc(auxind(j),1) )^2 + ...
                                           ( cvc(auxind(i),2) - cvc(auxind(j),2) )^2 + ...
                                           ( cvc(auxind(i),3) - cvc(auxind(j),3) )^2 ); 

                              % detecting neighbour voxels                  
                              if dist <= 1     % connectivity criterion
                                  indIJ = [ indIJ; [ i j ] ];
                                  edgeList = indIJ;                      
                              end                  
                          end
                     end
                end   
                
                % bypass non-connections (investigate here what happened!!)
                if isempty(indIJ)
                    continue
                end
                
                % invoke graphmetrics
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

                % save info
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelCoords = ivC(1,:);
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessValue = maxC;
                
                aux = cvi(iCnode); aux = aux(1); % gets only first of the list
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelInd = aux;
                
                aux = Ind(cvi(iCnode)); aux = aux(1); % gets only first of the list
                SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelMappedInd = aux;                
                
                SubM.(hpifn{f}).(hpifnc{fc}).gammaPotential = ...
                    HPInfo.(hpifn{f}).(hpifnc{fc}).gammaPotential;
                
                SubM.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialCoords = ...
                    HPInfo.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialCoords;
                
                SubM.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialGlobalInd = ...
                    HPInfo.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialGlobalInd;
                
                SubM.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialMappedInd = ...
                    HPInfo.(hpifn{f}).(hpifnc{fc}).maxGammaPotentialMappedInd;
                
                % oil saturated cells / total cells 
                SubM.(hpifn{f}).(hpifnc{fc}).oilVolumeFraction = ...
                    numel(find(HPInfo.(hpifn{f}).(hpifnc{fc}).oilSatMarker)) / ...
                    numel(HPInfo.(hpifn{f}).(hpifnc{fc}).oilSatMarker);
                
                
                %save perforation table for subcluster MaxC       
                ptset.savedir = fullfile(d.getTmpDir,'Subclustering',hpifn{f},hpifnc{fc});
                ptset.wellname = strcat('W',hpifnc{fc},'SubMaxC');
                ptset.geometry = 'K';
                ptset.perfs = SubM.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelCoords;               
                ptset = savePerfTable(ptset);
                
                %save perforation table for classical cluster MaxC                   
                qtset.savedir = fullfile(d.getTmpDir,'Subclustering',hpifn{f},hpifnc{fc});
                qtset.wellname = strcat('W',hpifnc{fc},'MaxC');
                qtset.geometry = 'K';
                qtset.perfs = HPInfo.(hpifn{f}).(hpifnc{fc}).maxClosenessVoxelCoords;               
                qtset = savePerfTable(qtset);
                
                clear ptset qtset
                
            end
            
        end
        
    end

end
%=======


%% Plotting example (for visibility)

% example values taken from HPInfo
exi = 1; % exi-th DRT value in the list
exic = 8; % exic-th cluster in the list ; 
% SubM indices: C2:1; C4:3; C6:5; C11:7; C13:8; C18:11; C22:12; C25:13; C46:18; C73:20       
exDRT = fieldnames(HPInfo); exDRT = exDRT{exi};
exC = fieldnames(HPInfo.(exDRT)); exC = exC{exic};
excvi = HPInfo.(exDRT).(exC).compVoxelInds;
excvim = HPInfo.(exDRT).(exC).mappedCompVoxelInds;
exmarker = HPInfo.(exDRT).(exC).oilSatMarker;
exso = HPInfo.(exDRT).(exC).oilSat;
exgamma = HPInfo.(exDRT).(exC).gammaPotential;
% 
% plotGrid(G,'FaceAlpha',0.2,'FaceColor','k','EdgeColor','None'); % background
% axis off
% plotCellData(subG,SOP,'FaceAlpha',0.2,'EdgeColor','None') % oil zone
 plotCellData(extractSubgrid(G,excvim(exmarker)),exso(exmarker),'FaceColor','c') % cluster saturated portion
 plotCellData(extractSubgrid(G,excvim(~exmarker)),exso(~exmarker),'FaceColor','g') % cluster nonsaturated portion
% 
% plotCellData(extractSubgrid(G,excvim),exgamma) % cluster nonsaturated portion
% 
% max clo from subcluster
mask = false(G.cells.num,1);
mask(SubM.(exDRT).(exC).maxClosenessVoxelMappedInd) = true;
plotGrid(G,mask,'FaceColor',[0.5,0.5,0.5],'EdgeColor','m' );
% 
% % max gamma from subcluster
% %mask(HPInfo.(exDRT).(exC).maxGammaPotentialMappedInd) = true;
% %plotGrid(G,mask,'FaceColor','k','EdgeColor','m' );

% max clo (cluster)
mask = false(G.cells.num,1);
cvcmax = HPInfo.(exDRT).(exC).maxClosenessVoxelCoords;
indmax = sub2ind(G.cartDims,cvcmax(1),cvcmax(2),cvcmax(3));
indmax = Ind(indmax(1));
mask(indmax) = true;
plotGrid(G,mask,'FaceColor',[0.5,0.5,0.5],'EdgeColor','r' );

% checking of zone
if exmarker(excvim == indmax)
    fprintf('==> Max closeness of cluster %s at oil zone.\n',exC)
else
    fprintf('==> Max closeness of cluster %s NOT at oil zone.\n',exC)
end

if exmarker(excvim == SubM.(exDRT).(exC).maxClosenessVoxelMappedInd(1))
    fprintf('==> Submax closeness of cluster %s at oil zone.\n',exC)
else
    warning('==> Submax closeness of cluster %s NOT at oil zone.\n',exC)
end


