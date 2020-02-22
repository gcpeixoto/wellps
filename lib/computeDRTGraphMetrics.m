function [Mf,Lf] = computeDRTGraphMetrics(opt_metrics, drtSt)
%COMPUTEDRTGRAPHMETRICS compute metrics and regression analysis for 
%                       the whole field's networks based on DRT 
%
% PARAMETERS: 
%       -   opt_metrics: structure that should contain the following 
%                        parameters:
%                           
%                        -  'nofsc' (string): number of significant cells. 
%                           This is the minimum value of cells to be 
%                           considered for a cluster (i.e., a minimum 
%                           volume threshold)
%
%                        -  'seps' (float): slope epsilon. Tolerance that
%                           allows to characterize a cluster as being of 
%                           'high-' or 'low-' performance. A cluster is 
%                           'high-performance' when slope belongs to 
%                           [1 - seps, 1 + seps].
%                           
%                        - 'R2min' (float): mininum determination
%                        coefficient. Threshold to classify cluster in the
%                        same way as 'seps'. We require that the determination 
%                        coefficient R2 belongs to [R2min,1.0].
%
%                        - 'outDir' (string): path to output directory
%                        in which the .mat files should be saved (optional).
%                        If empty, the WELLPS standard directory /mat is
%                        set.
%
%       - drtSt:         structure of DRT connections.
%
%   RETURNS:
%       -   Mf : structure having file paths to graph metrics .mat 
%
%       -   L : structure having file paths to linear regression .mat
%               
%  REMARK: see G.P. Oliveira et al. (2016), DOI: 10.1016/j.petrol.2017.06.016.

% Call 
d = DirManager;
matdir = d.getMatDir; 
cppdir = d.getCppDir;


% checking
p = {'nofsc','seps','R2min'};
if ~all(ismember(p,fieldnames(opt_metrics)))   
    error('wellps:computeDRTGraphMetrics','opt_metrics is not a struct object');
end
    
% standard output dir
if isempty(find(ismember(fieldnames(opt_metrics),'outDir'), 1))     
    warning('wellps:computeDRTGraphMetrics','Output directory was set to standard');
    outDir = matdir;
else 
    outDir = opt_metrics.outDir;
end

% fields (should result in 'DRTx')
fnames = fieldnames(drtSt);
nnames = numel(fnames);

fprintf('=> Collecting rock type clustering... \n');
tstart = tic; % timing
for n = 1:nnames
    
    S = drtSt.(fnames{n});
           
    val = S.value;     
    
    % gets rock type string to display
    [~,pos] = regexp(fnames{n},'\w\d','match');
    var = fnames{n}(1:pos);
    fprintf('==> %s = %d \n',var,val);
    
    avc = S.allVoxelCoords;
    ncomps = S.allNComps;    
    Madj = S.allAdjMatrix;    
    
    % saving options
    ave = S.averaging;
    base = S.logBase;
    M.(fnames{n}).value = val;
    M.(fnames{n}).averaging = ave;
    M.(fnames{n}).logBase = base;
    M.(fnames{n}).nofsc = opt_metrics.nofsc;
    M.(fnames{n}).slopes = [1 - opt_metrics.seps, 1 + opt_metrics.seps];
    M.(fnames{n}).R2min = opt_metrics.R2min;
    
    L.(fnames{n}).value = val;
    L.(fnames{n}).averaging = ave;
    L.(fnames{n}).logBase = base;
    L.(fnames{n}).nofsc = opt_metrics.nofsc;
    L.(fnames{n}).slopes = [1 - opt_metrics.seps, 1 + opt_metrics.seps];
    L.(fnames{n}).R2min = opt_metrics.R2min;
    
    count = 0;
    for idComp = 1:ncomps   
        
        cnn = S.compNNodes{idComp};  
        
        if cnn >= opt_metrics.nofsc % if significant components 
                          
            % performs linear regression
            logPHIZ = S.compLogPHIZ{idComp};
            logRQI  = S.compLogRQI{idComp};
            [ R, m, b ] = regression( logPHIZ, logRQI, 'one' );
                                                                                                                
            %------------------ subgraph (connected component network)
            % finding vertices in the big adjacency matrix to set up the 
            % adjacency matrix for the connected component and, then, 
            % set the subgraph of the network 
            %
            % Performance techniques introduced here: 
            % - direct dynamic allocation v(e)
            % - logical search + find for vector 'id'. A faster variant 
            %   of 'strmatch' or 'ismember'.            
            cvc = S.compVoxelCoords{idComp};            
            v = []; 
            for e = 1:size(cvc,1) 
                id = (avc(:,1)==cvc(e,1) & ...        
                      avc(:,2)==cvc(e,2) & ...
                      avc(:,3)==cvc(e,3)); 
                id = find(id ==1); 
                v(e)=id; % global indices
            end
            
            MadjComp = subgraph( Madj, v ); % component's adjacency matrix                

            %------------------ centrality metrics 
            fprintf(['==> Computing metrics for cluster ',...
                     'C%d (%d cells)...\n'], ...
                    idComp, size(v,2));

            % SNAP interface
            edfile = saveAdjEdges(MadjComp);              
            outfile = fullfile(d.getTmpDir,'metrics.txt');
            exec = sprintf('graphMetrics %s %s',edfile,outfile);
            
            % call SNAP 
            fprintf('%s\n',repmat('=',[1,75]))
            fprintf('%sWELLPS - SNAP Interface :: running C++\n',repmat(' ',[1,20]));
            fprintf('%s\n',repmat('=',[1,75]))
            system( fullfile(cppdir,exec) );
            
            % get metrics and erase file
            [nodeID,deg,clns,betw] = getMetricsData(outfile); 
            % get rid of input file
            delete(edfile);
            
            maxC = max(clns);           % max closeness = min farness
            iC = clns == maxC;          % network closer nodes
            iCnode = nodeID(iC);        % getting node id (there might be more than 1)
            ivC = avc( v(iCnode),: );   % global voxel coordinates
                        
            %disp('----> Storing structures...');
            
            count = count + 1; % component counter
            
            % store components whose volume is significant
            M.(fnames{n}).idComp{count} = idComp;                
            M.(fnames{n}).degreeCentrality{count} = deg;
            M.(fnames{n}).closenessCentrality{count} = clns;
            M.(fnames{n}).betweenessCentrality{count} = betw;                
            M.(fnames{n}).maxClosenessVoxelCoords{count} = ivC;
            M.(fnames{n}).adjMatrix{count} = MadjComp;
            
            L.(fnames{n}).idComp{count} = idComp;
            L.(fnames{n}).R2{count} = R*R;
            L.(fnames{n}).slope{count} = m;
            L.(fnames{n}).offset{count} = b;
            L.(fnames{n}).logPHIZ{count} = logPHIZ;
            L.(fnames{n}).logRQI{count} = logRQI;
            
            % linear regression criteria (performance)                 
            if ( m >= 1 - opt_metrics.seps   && ... 
                 m <= 1 + opt_metrics.seps ) && ... 
                 (R*R >= opt_metrics.R2min)
                L.(fnames{n}).performance{count} = true; % high-performance: 1 
            else
                L.(fnames{n}).performance{count} = false; % low-performance: 0
            end 
            
            
        end % nofsc loop
        
    end % components loop
        
    if count ~= 0 % saving structure to .mat, if any         
        M = M.(fnames{n});
        save( fullfile(outDir,strcat('DRT_',ave,'_',base,'_',num2str( val ),'_metrics','.mat')),'M');
        %disp('----> metrics .mat file saved.')                
                
        L = L.(fnames{n});
        save( fullfile(outDir,strcat('DRT_',ave,'_',base,'_',num2str( val ),'_linregr','.mat')),'L');
        %disp('----> regression .mat file saved.')
                
        % files
        Mf.(fnames{n}) = fullfile(outDir,strcat('DRT_',ave,'_',base,'_',num2str( val ),'_metrics','.mat'));
        Lf.(fnames{n}) = fullfile(outDir,strcat('DRT_',ave,'_',base,'_',num2str( val ),'_linregr','.mat'));
        
        
    else        
        fprintf('=> No components found. You may try to increase "nofsc" and rerun the code.\n');
    end        
    
end % DRT loop
        fprintf('=> computeDRTGraphMetrics finished after %g seconds.\n',toc(tstart));
end

%% -- HELPER

function [node,deg,cln,bet] = getMetricsData(mfile)
%GETMETRICSDATA read file with centrality data computed from SNAP and
%               retrieves in arrays.
%
% PARAMETERS:
%  mfile - SNAP output file 
%
% RETURNS:
%   node - graph node ID array
%	deg	 - degree centrality array
%	cln	 - closeness centrality array
%	bet	 - betweeness centrality array

% Assumes that SNAP has computed the centralities correctly and
% saved them to the temporary metrics.txt. 
tab = importdata(mfile);

% Getting centralities from the file computed via SNAP.
% Note that this file is a table containing all the centralities 
% computable by SNAP. However, we take only the node ID, degree, 
% closeness, and betweeness, in this order.

% interest table
Mtab = tab.data(:,1:4);

% node ID
node = Mtab(:,1);

% centralities
deg  = Mtab(:,2);
cln  = Mtab(:,3);
bet  = Mtab(:,4);

% TODO 
% This was added here to check if SNAP is producing negative betweeness. 
% I suspect that there exists an error in betweeness computation. 
% I need to check betweeness from other tools.  
aux = length(bet(bet<0));
if aux > 0
    warning('wellps:getMetricsData',['A number of %d negative values were ',...
             'found for BETWEENNESS. Metric is unreliable.\n'],aux);
end

% delete the temporary files
delete(mfile)

end
