function C = findDRTConnections(drt, P, ave, base, nofsc, tocsv, dv)
%FINDDRTCONNECTIONS searches for all the connected components 
%                   of cells with same DRT value
%
%
% PARAMETERS:
%   
%       drt     - array listing the wanted DRT values
%
%       P       - structure of properties ('PHI','KX',etc.)
%
%       ave     - which average to use to recover values. 
%                 There are 3 options: 'arithmetic', 'geometric',
%                 'normalized'. (see ::computeParams to see definitions)
%
%       base    - which log base to use to compute properties.
%                 There are 2 options: 'log10', 'ln', where 'ln' is the
%                 Napier log.
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

if ~isnumeric(drt)
    error('wellps:findDRTConnections','drt must be an array containing DRT values to be gone through.');    
end

if ~isstruct(P)
    error('wellps:findDRTConnections','Argument must be a structure.');    
end

% copying fundamental properties 
PHI = P.PHI;    
PHIZ = P.PHIZ;  
KX = P.KX;  
KY = P.KY;  
KZ = P.KZ;

if strcmp(base,'ln')   
    
    LOGPHIZ = P.LNPHIZ;
    logphiz = 'lnphiz,'; % .csv column name
    
elseif strcmp(base,'log10')  
    
    LOGPHIZ = P.Log10PHIZ;
    logphiz = 'log10phiz,'; % .csv column name
    
end

% copying according to ave and log-base
switch ave
    
    case 'arithmetic' 
        
        KMEAN = P.KA;
        RQI = P.RQIA;
        FZI = P.FZIA;    
        
        % .csv column names
        kmean = 'ka,'; 
        rqi = 'RQIA,';
        fzi = 'FZIA,';
        
        if strcmp(base,'ln')    
            
            DRT = P.DRTA_LN;            
            LOGRQI = P.LNRQIA;
            
            logrqi = 'lnRQIA';
    
        elseif strcmp(base,'log10') 
    
            DRT = P.DRTA_LOG10;            
            LOGRQI = P.Log10RQIA;
            
            logrqi = 'log10RQIA';
    
        else
        
            error('wellps:findDRTConnections','Log base not available.');

        end

    case 'geometric' 
        
        KMEAN = P.KG;
        RQI = P.RQIG;
        FZI = P.FZIG;
        
        % .csv column names
        kmean = 'kg,'; 
        rqi = 'RQIG,';
        fzi = 'FZIG,';
        
        if strcmp(base,'ln') 
            
            DRT = P.DRTG_LN;            
            LOGRQI = P.LNRQIG;
            
            logrqi = 'lnRQIG';
    
        elseif strcmp(base,'log10') 
    
            DRT = P.DRTG_LOG10;            
            LOGRQI = P.Log10RQIG;
            
            logrqi = 'log10RQIG';
    
        else
        
            error('wellps:findDRTConnections','Log base not available.');

        end
    
    case'normalized' 
        
        KMEAN = P.KN;
        RQI = P.RQIN;
        FZI = P.FZIN;
        
        % .csv column names
        kmean = 'kn,'; 
        rqi = 'RQIN,';
        fzi = 'FZIN,';
        
        if strcmp(base,'ln')     
            
            DRT = P.DRTN_LN;            
            LOGRQI = P.LNRQIN;
            
            logrqi = 'lnRQIN';
    
        elseif strcmp(base,'log10') 
    
            DRT = P.DRTN_LOG10;            
            LOGRQI = P.Log10RQIN;
            
            logrqi = 'log10RQIN';
    
        else
        
            error('wellps:findDRTConnections','Log base not available.');

        end
        
end

if ~all(drt)
    warning('wellps:findDRTConnections','DRT = 0 was found in your list. This will result in null porosity blobs.');
end

if ~ischar(tocsv)
    error('wellps:findDRTConnections','argument tocsv must be a char: "y" [yes] or "n" [no].')    
end

if dv ~= 1 
    error('wellps:findDRTConnections','Neighbourhood criterion recommended for now is 6. This means a cell distance = 1.')
end


% file header used in the loop
hdr = {'i,'; 'j,'; 'k,'; 'phi_e,'; ...
        'kx,'; 'ky,'; 'kz,'; ...
        kmean; ...
        'phiz,'; rqi; fzi; ...
        logphiz; logrqi}'; % transposed!
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

tstart = tic; % timing
for val = drt(1):drt(end)
    
    fprintf('----> Sweeping field: DRT = %d... \n',val);
    
    % filtering grid to capture voxels with a specific DRT
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
        fprintf('----> No connections found for DRT = %d... \n',val);
        continue;
    end
    
    aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
    indIJ = [ indIJ; aux ]; % filling        
    
    disp('----> Computing adjacency matrix...');            
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
    disp('----> Finding connected components...');
    [ncomp,compSizes,members] = networkComponents(MDadj);
            
    % global  
    drtSt.value = val;                        % DRT value
    drtSt.allAdjMatrix = MDadj;               % graph adjacency matrix
    drtSt.allAdjEdgeList = edgeList;          % graph edge list
    drtSt.allVoxelCoords = coordsDRT;         % cell coordinates (i,j,k)     
    drtSt.allVoxelInds = indz;                % cell linear indices
    drtSt.allPHI = PHI( indz );               % effective porosity
    drtSt.allKX = KX( indz );                 % permeability-x
    drtSt.allKY = KY( indz );                 % permeability-y
    drtSt.allKZ = KZ( indz );                 % permeability-z
    drtSt.allKMEAN = KMEAN( indz );           % average permeability 
    drtSt.allPHIZ = PHIZ( indz );             % normalized porosity
    drtSt.allRQI = RQI( indz );               % RQI
    drtSt.allFZI = FZI( indz );               % FZI
    drtSt.allLogPHIZ = LOGPHIZ( indz );       % LogPHIZ
    drtSt.allLogRQI = LOGRQI( indz );         % LogRQI
    drtSt.allNComps = ncomp;                  % number of components in the graph

    % mount matrix to export
    mat = [ coordsDRT(:,1) ... 
            coordsDRT(:,2) ...
            coordsDRT(:,3) ...
            PHI( indz )    ...
            KX( indz )     ...
            KY( indz )     ...
            KZ( indz )     ...
            KMEAN( indz )  ...
            PHIZ( indz )   ...
            RQI( indz )    ...
            FZI( indz )    ...
        LOGPHIZ( indz )    ...
        LOGRQI( indz )   ];            
            
    if tocsv == true
        
        % preparing csv file
        fname = strcat('../csv/Table-Field_DRT_',ave,'_',base,'_',num2str( val ),'.csv');    
        dlmwrite(fname,hdr,'');        
    
        % append matrix
        dlmwrite(fname,mat,'-append');
        fprintf('----> Table-Field_DRT_%s.csv file saved. \n',...
                num2str(val));        
    end
    
    % cluster (each component)
    for idcomp = 1:ncomp
        
        fprintf('--------> Computing component = %d \n',idcomp);
        
        aux = coordsDRT( members{idcomp},: );
        drtSt.compVoxelCoords{idcomp} = aux; 
        drtSt.compVoxelInds{idcomp} = indz( members{idcomp} );          
        drtSt.compNNodes{idcomp} = compSizes(idcomp);           
        drtSt.compPHI{idcomp} = PHI( indz( members{idcomp} ) );
        drtSt.compKX{idcomp} = KX( indz( members{idcomp} ) );
        drtSt.compKY{idcomp} = KY( indz( members{idcomp} ) );
        drtSt.compKZ{idcomp} = KZ( indz( members{idcomp} ) );
        drtSt.compKMEAN{idcomp} = KMEAN( indz( members{idcomp} ) );
        drtSt.compPHIZ{idcomp} = PHIZ( indz( members{idcomp} ) );
        drtSt.compRQI{idcomp} = RQI( indz( members{idcomp} ) );
        drtSt.compFZI{idcomp} = FZI( indz( members{idcomp} ) );
        drtSt.compLogPHIZ{idcomp} = LOGPHIZ( indz( members{idcomp} ) );
        drtSt.compLogRQI{idcomp} = LOGRQI( indz( members{idcomp} ) ); 
                  
        % matrix to export
        mat = [ aux(:,1)                         ... 
                aux(:,2)                         ...
                aux(:,3)                         ...
                PHI( indz( members{idcomp} ) )   ...
                 KX( indz( members{idcomp} ) )   ...
                 KY( indz( members{idcomp} ) )   ...
                 KZ( indz( members{idcomp} ) )   ...
              KMEAN( indz( members{idcomp} ) )   ...
               PHIZ( indz( members{idcomp} ) )   ...
                RQI( indz( members{idcomp} ) )   ...
                FZI( indz( members{idcomp} ) )   ...
            LOGPHIZ( indz( members{idcomp} ) )   ...
            LOGRQI( indz( members{idcomp} ) )  ];
           
        
        if (tocsv == true) && (compSizes(idcomp) >= nofsc)              
            % preparing csv file
            fname = strcat('../csv/Table-Cluster_',num2str(idcomp),'_DRT_',ave,'_',base,'_',num2str( val ),'.csv');        
            dlmwrite(fname,hdr,'');       
            
            % append matrix 
            dlmwrite(fname,mat,'-append'); 
            fprintf('----> Table-Cluster_%s_DRT_%s.csv file saved. \n',...
                num2str(idcomp),num2str(val));            
        end                        
    end              

    % saving structure to .mat     
    save( strcat('../mat/DRT_',ave,'_',base,'_',num2str( val ),'.mat'),'drtSt'); % saving
    fprintf('----> DRT_%s.mat file saved. \n',num2str(val));
    
    % dynamic attribution     
    f = strcat('DRT',num2str(val));
    C.(f) = drtSt;    
    
    clear drtSt;    % frees to recompute
    
end

tfinal = toc(tstart);
fprintf('----> findDRTConnections finished after %g seconds. \n',tfinal);
