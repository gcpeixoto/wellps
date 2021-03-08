function [M,L] = computeDRTStarGraphMetrics(opt_metrics, drtSt)
%COMPUTEDRTStarGRAPHMETRICS compute metrics and regression analysis for 
%                       the whole field's networks based on DRT* found
%                       after FZI*. 
%
% PARAMETERS: 
%       -   opt_metrics:   structure that should contain the following 
%                          parameters:
%                           
%                        -  'nofs' (string): number of significant cells. 
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
%                        in which the .mat files shoul be saved (optional).
%                        If empty, the WELLPS standard directory /mat is
%                        set.
%
%   RETURNS:
%       -   M : structure having full information about the graph metrics. 
%
%       -   L : structure having linear regression information per cluster
%               
%  REMARK: see G.P. Oliveira et al. (2016), DOI: 10.1016/j.petrol.2017.06.016.


% checking
p = {'nofsc','seps','R2min'};
if ~all(ismember(p,fieldnames(opt_metrics)))   
    error('wellps:computeDRTGraphMetricsByFZIStar','opt_metrics is not a struct object');
end
    
% standard output dir
if isempty(find(ismember(fieldnames(opt_metrics),'outDir'), 1))     
    warning('wellps:computeDRTGraphMetricsByFZIStar','Output directory was set to standard');
    outDir = '../mat/';
else 
    outDir = opt_metrics.outDir;
end

% fields (should result in 'DRTx')
fnames = fieldnames(drtSt);
nnames = numel(fnames);

for n = 1:nnames
    
    S = drtSt.(fnames{n});
    
    val = S.value;     
    fprintf('----> Sweeping DRT: %d... \n',val);
    
    avc = S.allVoxelCoords;
    ncomps = S.allNComps;    
    Madj = S.allAdjMatrix;    
    
    % saving options
    ave = S.averaging;
    base = S.logBase;
    M.(fnames{n}).value = val;
    M.(fnames{n}).averaging = ave;
    M.(fnames{n}).logBase = base;
    M.(fnames{n}).nofs = opt_metrics.nofs;
    M.(fnames{n}).slopes = [1 - opt_metrics.seps, 1 + opt_metrics.seps];
    M.(fnames{n}).R2min = opt_metrics.R2min;
    
    L.(fnames{n}).value = val;
    L.(fnames{n}).averaging = ave;
    L.(fnames{n}).logBase = base;
    L.(fnames{n}).nofs = opt_metrics.nofs;
    L.(fnames{n}).slopes = [1 - opt_metrics.seps, 1 + opt_metrics.seps];
    L.(fnames{n}).R2min = opt_metrics.R2min;
    
    count = 0;
    for idComp = 1:ncomps   
        
        cnn = S.compNNodes{idComp};  
        
        if cnn >= opt_metrics.nofs % if significant components 
                          
            % performs linear regression
            logSQRTK = S.compLogSqrtK{idComp};
            logSQRTPHI  = S.compLogSqrtPHI{idComp};
            [ R, m, b ] = regression(logSQRTPHI,logSQRTK, 'one' );
                                                                                                                
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
                v(e)=id;                              % global indices
            end
            
            MadjComp = subgraph( Madj, v );           % component's adjacency matrix                

            %------------------ centrality metrics 
            fprintf('----> Computing metrics for %d nodes... \n', size(v,2) );

            % SNAP interface
            edfile = saveAdjEdges(MadjComp);  
            ! ./../cpp/graphMetrics
            [nodeID,deg,clns,betw] = getMetricsData(edfile);  
            
            maxC = max(clns);           % max closeness = min farness
            iC = clns == maxC;          % network closer nodes
            iCnode = nodeID(iC);        % getting node id (there might be more than 1)
            ivC = avc( v(iCnode),: );   % global voxel coordinates
                        
            disp('----> Storing structures...');
            
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
            L.(fnames{n}).logSQRTK{count} = logSQRTK;
            L.(fnames{n}).logSQRTPHI{count} = logSQRTPHI;
            
            % linear regression criteria (performance)                 
            if ( m >= 1 - opt_metrics.seps   && ... 
                 m <= 1 + opt_metrics.seps ) && ... 
                 (R*R >= opt_metrics.R2min)
                L.(fnames{n}).performance{count} = true; % high-performance: 1 
            else
                L.(fnames{n}).performance{count} = false; % low-performance: 0
            end 
            
            
        end % nofs loop
        
    end % components loop
        
    if count ~= 0 % saving structure to .mat, if any         
        Maux = M.(fnames{n});
        save( strcat(outDir,'DRTStar_',ave,'_',base,'_',num2str( val ),'_Metrics','.mat'),'Maux');
        disp('----> metrics .mat file saved.')                

        Laux = L.(fnames{n});
        save( strcat(outDir,'DRTStar_',ave,'_',base,'_',num2str( val ),'_LinRegr','.mat'),'Laux');
        disp('----> regression .mat file saved.')
        
        clear Maux Laux % free up memory
        
    else        
        disp('----> No components found. You may wish to change "nofn" and rerun the code.');
    end        
    
end % DRT loop




