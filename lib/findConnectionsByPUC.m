function C = findConnectionsByPUC(dobj, puc, PUC, nofsc, tocsv, dv)
%FINDCONNECTIONSBYPUC searches for all the connected components 
%                            of cells with same PUC value based on 
%                            proxy function for productivity potential.
%                            
%                            PUC = Productivity Unit Class
%
%
% PARAMETERS:
%
%       dobj    - DirManager class object.
%   
%       puc     - array listing the wanted PUC values
%
%       nofsc   - acronym for 'number of significant cells'. This is 
%                 the minimum number of cells that will be considered 
%                 when saving the .csv files.
%
%       tocsv   - logical flag to save or not files into csv/ directory. 
%
%
%       dv      - rule of Moore neighbourhood to use (1, for 6 neighbours 
%                 or sqrt(3), for 26 neighbours)
%
% RETURNS:
%
%       C       - structure containing structures of data related to each DRT.
%
%               - .mat files saved into /mat.


if ~isnumeric(puc)
    error('wellps:findConnectionsByPUC','puc must be an array containing integer values from 0 to 4.');    
end

if ~isnumeric(PUC)
    error('wellps:findConnectionsByPUC','Argument must be a matrix with proxy function values.');    
end

if ~ischar(tocsv)
    error('wellps:findConnectionsByPUC','argument tocsv must be a char: "y" [yes] or "n" [no].')    
end

if dv ~= 1 
    error('wellps:findConnectionsByPUC','Neighbourhood criterion recommended for now is 6. This means a cell distance = 1.')
end

% file header used in the loop
hdr = {'i,'; 'j,'; 'k'}'; % transposed!
hdr = sprintf('%s\t',hdr{:});
hdr(end)='';

C = []; % struct

% csv file
switch tocsv
    case 'y'
        tocsv = true;
    case 'n'
        tocsv = false;
end


% to reuse code structure, we set this
DRT = PUC;

tstart = tic; % timing
fprintf('Sweeping field...\n')

for val = puc(1):puc(end)
    
    fprintf('=> PUC = %d \n',val);
    
    % filtering grid to capture voxels with a specific PUC
    [ ~, coordsDRT, indz ] = maskVoxels(DRT,val); 
    
    
    %{
        Adjacency matrix
        ================
                    
        Strategy:        
                    
            - Find the entries whose distance is less than 1 (sqrt(3)), which
              stands for the 6-voxel (or 26-voxel) neighbourhood
    
            - Set up sparse adjacency matrix from the entries found to be
              able to create an undirected graph
    
            - Store graph edge list matrix
              
    %}
            
    indIJ = [];
    for i = 1:size(coordsDRT,1)        
        for j = i:size(coordsDRT,1)                                                            
              if i ~= j % skipping null distance                   
                  dist = sqrt( ( coordsDRT(i,1) - coordsDRT(j,1) )^2 + ...
                               ( coordsDRT(i,2) - coordsDRT(j,2) )^2 + ...
                               ( coordsDRT(i,3) - coordsDRT(j,3) )^2 ); 
                  
                  % detecting neighbour voxels                  
                  if dist <= dv     % connectivity criterion
                      indIJ = [ indIJ; [ i j ] ];
                      edgeList = indIJ;                      
                  end                  
              end
         end
    end   
    
    if isempty(indIJ)
        fprintf('==> No connections found for PUC = %d... \n',val);
        continue;
    end
    
    aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
    indIJ = [ indIJ; aux ]; % filling        
    
    %disp('----> Computing adjacency matrix...');            
    % creates adjacency matrix n x n by marking 1 for connected nodes
    MDadj = sparse( indIJ(:,1),indIJ(:,2),1,size(coordsDRT,1),size(coordsDRT,1) ); 
            
    %{ 
        Finding network components 
        ==========================     
    
        - Uses the function 'networkComponents', by Daniel Larremore 
        @MathWorks central
    
        - From adjacency matrix, computes all the components 
          (isolated + connected) for the specific DRT network
    
        - The function is based on the connected component algorithm
          for graphs (e.g., see 
         http://people.sc.fsu.edu/~jburkardt/classes/asa_2011/asa_2011_graphs.pdf)
    
       
        DRT Data structure
        ==================
    
                drtSt
                  |           
                  |- value  
                  |
                  |- all{ ... }                 }
                  .                             }
                  .                             } Global family ( all cells with such DRT )
                  .                             }
                  |                             } [ SEVERAL DATA ]
                  |
                  |- comp{ ... }                }        
                  |    |                        }
                  |    |- comp{ ... }{idcomp1}  }
                  |    |- comp{ ... }{idcomp2}  } 
                  .    .                        } Cluster Family ( K components of variable cells per DRT )
                  .    .                        }
                  .    .                        }
                  |    |                        }
                  |    |- comp{ ... }{idcompK}  } [ SEVERAL DATA ]
                                                                  
    %}      
    %disp('----> Finding connected components...');
    [ncomp,compSizes,members] = networkComponents(MDadj);
            
    % global  
    pucSt.value = val;                        % PUC value
    pucSt.allAdjMatrix = MDadj;               % graph adjacency matrix
    pucSt.allAdjEdgeList = edgeList;          % graph edge list
    pucSt.allVoxelCoords = coordsDRT;         % cell coordinates (i,j,k)     
    pucSt.allVoxelInds = indz;                % cell linear indices  
    pucSt.allNComps = ncomp;                  % number of components in the graph

    % mount matrix to export
    mat = [ coordsDRT(:,1)  ... 
            coordsDRT(:,2)  ...
            coordsDRT(:,3) ];            
            
    if tocsv == true
        
        % preparing csv file
        fname = strcat(dobj.getCsvDir,'/table-field_PUC_',num2str( val ),'.csv');    
        dlmwrite(fname,hdr,'');        
    
        % append matrix
        dlmwrite(fname,mat,'-append');
        %fprintf('----> File %s saved. \n',fname);        

    end
    
    fprintf('==> Computing %d clusters...',ncomp);
    % cluster (each component)
    cnofs = 0;
    for idcomp = 1:ncomp
                        
        aux = coordsDRT( members{idcomp},: );
        pucSt.compVoxelCoords{idcomp} = aux; 
        pucSt.compVoxelInds{idcomp} = indz( members{idcomp} );          
        pucSt.compNNodes{idcomp} = compSizes(idcomp);                  
                  
        % matrix to export
        mat = [ aux(:,1) ... 
                aux(:,2) ...
                aux(:,3) ];
                   
        if (tocsv == true) && (compSizes(idcomp) >= nofsc)              
            % preparing csv file
            fname = strcat(dobj.getCsvDir,'/table-cluster_',num2str(idcomp),'_PUC_',num2str( val ),'.csv');        
            dlmwrite(fname,hdr,'');       
            
            % append matrix 
            dlmwrite(fname,mat,'-append'); 
            
            % count clusters with >= nofsc
            cnofs = cnofs + 1;
        end                        
    end    
    fprintf(' %d clusters with >= %d cells.\n',cnofs,nofsc);            

    % saving structure to .mat     
    save( strcat(dobj.getMatDir,'/PUC_',num2str( val ),'.mat'),'pucSt'); % saving
    %fprintf('----> PUC_%s.mat file saved. \n',num2str(val));
    
    % dynamic attribution     
    f = strcat('PUC',num2str(val));
    C.(f) = pucSt;    
    
    clear pucSt;    % frees to recompute
    
end

tfinal = toc(tstart);
fprintf('findConnectionsByPUC finished after %g seconds. \n',tfinal);

end

