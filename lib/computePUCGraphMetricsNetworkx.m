function Mf = computePUCGraphMetricsNetworkx(opt_metrics, pucSt)
%computePUCGraphMetricsNetworkx compute metrics for networks based 
%                       on PUC-values by using a wrap to NETWORKX package.
%
% PARAMETERS: 
%       -   opt_metrics: structure that should contain the following 
%                        parameters:
%                           
%                        - 'nofsc' (numerical): number of significant cells. 
%                          This is the minimum value of cells to be 
%                          considered for a cluster (i.e., a minimum 
%                          volume threshold)                       
%
%                        - 'outDir' (string): path to output directory
%                        in which the .mat files should be saved (optional).
%                        If empty, the WELLPS standard directory /mat is
%                        set.
%
%       - pucSt:         structure of PUC-based connections.
%
%   RETURNS:
%       -   Mf : structure having file paths to graph metrics .mat 
%               
%  REMARK: see G.P. Oliveira et al. (2016), DOI: 10.1016/j.petrol.2017.06.016.

% Call 
d = DirManager;
matdir = d.getMatDir; 
pydir = d.getPyDir;


% checking
p = {'nofsc'};
if ~all(ismember(p,fieldnames(opt_metrics)))   
    error('wellps:computePUCGraphMetricsNetworkx','opt_metrics is not a struct object');
end
    
% standard output dir
if isempty(find(ismember(fieldnames(opt_metrics),'outDir'), 1))     
    warning('wellps:computePUCGraphMetricsNetworkx','Output directory was set to standard');
    outDir = matdir;
else 
    outDir = opt_metrics.outDir;
end

% fields (should result in 'PUCx')
fnames = fieldnames(pucSt);
nnames = numel(fnames);

fprintf('=> Collecting PUC clusters... \n');
tstart = tic; % timing
for n = 1:nnames
        
    S = pucSt.(fnames{n});
           
    val = S.value;     
    
    % gets rock type string to display
    [~,pos] = regexp(fnames{n},'\w\d','match');
    var = fnames{n}(1:pos);
    fprintf('==> %s = %d \n',var,val);
    
    avc = S.allVoxelCoords;
    ncomps = S.allNComps;    
    Madj = S.allAdjMatrix;    
    
    % saving options
    M.(fnames{n}).value = val;        
    M.(fnames{n}).nofsc = opt_metrics.nofsc;        
            
    count = 0;
    for idComp = 1:ncomps   
        
        cnn = S.compNNodes{idComp};  
        
        if cnn >= opt_metrics.nofsc % if significant components 
                                                                                                                                                      
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

            % NETWORKX interface
            edfile = saveAdjEdges(MadjComp);              
            outfile = fullfile(d.getTmpDir,'metrics.txt');
            exec = [d.getPyExec,' ',pydir,filesep, ...
                    sprintf('graphMetrics.py %s %s',edfile,outfile)];
            
            % call NETWORKX 
            fprintf('%s\n',repmat('=',[1,75]))
            fprintf('%sWELLPS - NETWORKX Interface :: running Python\n',repmat(' ',[1,20]));
            fprintf('%s\n',repmat('=',[1,75]))
            system(exec);
                                                                                                                                                 
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
                                                                       
        end % nofsc loop
        
    end % components loop
        
    if count ~= 0 % saving structure to .mat, if any         
        M = M.(fnames{n});
        save( fullfile(outDir,strcat('PUC_',num2str( val ),'_metrics','.mat')),'M');
        %disp('----> metrics .mat file saved.')                
                                        
        % files
        Mf.(fnames{n}) = fullfile(outDir,strcat('PUC_',num2str( val ),'_metrics','.mat'));        
                
    else        
        fprintf('=> No components found. You may try to increase "nofsc" and rerun the code.\n');
    end        
    
end % DRT loop
        fprintf('=> computePUCGraphMetricsNetworkx finished after %g seconds.\n',toc(tstart));
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

% Getting centralities from the file computed via NETWORKX.
% Note that this file is a table containing all the centralities 
% computable by NETWORKX. However, we take only the node ID, degree, 
% closeness, and betweeness, in this order.

% interest table
Mtab = tab.data(:,1:4);

% node ID
node = Mtab(:,1);

% centralities
deg  = Mtab(:,2);
cln  = Mtab(:,3);
bet  = Mtab(:,4);

% delete the temporary files
delete(mfile)

end
