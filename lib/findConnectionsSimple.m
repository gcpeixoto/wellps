function F = findConnectionsSimple(cvc)
%FINDCONNECTIONSSIMPLE searches for all the connected components 
%                      in a group of given cells.
%
%
% PARAMETERS:
%       cvc - list of cell logical coordinates 
%
%
% RETURNS: 
%       F - structure containing connected components for the group.
%
% 
% REMARK: see wellps::findDRTConnections for a generalization based on 
%         discrete rock typing. This method is also for clustering
%         purposes. However, it differs from the method above, which was 
%         specifically built to be used with DRT-based clustering.

indIJ = [];
for i = 1:size(cvc,1)        
    for j = i:size(cvc,1)                                                            
          if i ~= j % skipping null distance                   
              dist = sqrt( ( cvc(i,1) - cvc(j,1) )^2 + ...
                           ( cvc(i,2) - cvc(j,2) )^2 + ...
                           ( cvc(i,3) - cvc(j,3) )^2 ); 

              % detecting neighbour cells by 6-neighbor connectivity criterion              
              if dist <= 1     
                  indIJ = [ indIJ; [ i j ] ];
                  edgeList = indIJ;                      
              end                  
          end
     end
end   

if isempty(indIJ)
    fprintf('----> No connections found for DRT = %d... \n',val);    
end

aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
indIJ = [ indIJ; aux ]; % filling        

disp('----> Computing adjacency matrix...');            

% creates adjacency matrix n x n by marking 1 for connected nodes
Madj = sparse( indIJ(:,1),indIJ(:,2),1,size(cvc,1),size(cvc,1) ); 


disp('----> Finding connected components...');
[ncomp,compSizes,members] = networkComponents(Madj);

% getting cells per component
compvc = cell(1,ncomp);
adj = cell(1,ncomp);
for nm = 1:ncomp
    
    compvc{nm} = cvc(members{nm},:);
    
    
    % form component adjacency matrix 
    v = []; 
    for e = 1:size(compvc{nm},1) 
        id = (cvc(:,1) == compvc{nm}(e,1) & ...        
              cvc(:,2) == compvc{nm}(e,2) & ...
              cvc(:,3) == compvc{nm}(e,3)); 
        id = find(id == 1); 
        v(e) = id; % global indices
    end
    adj{nm} = subgraph( Madj, v ); % component's adjacency matrix
       
end

% structure
F.ncomp = ncomp;
F.compSizes = compSizes; 
F.compMembers = members;
F.allAdjMatrix = Madj;
F.compVoxelCoords = compvc;
F.compAdjMatrix = adj;
