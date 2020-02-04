function C = findDRTConnectionsByFZIStar(dobj, drt, P, ave, base, nofsc, tocsv, dv)
%FINDDRTCONNECTIONSBYFZISTAR searches for all the connected components 
%                            of cells with same DRT value based on 
%                            FZI* by Paiaman et al (2015) 
%
%
% PARAMETERS:
%
%       dobj    - DirManager class object.
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
    error('wellps:findDRTConnectionsByFZIStar','drt must be an array containing DRT values to be gone through.');    
end

if ~isstruct(P)
    error('wellps:findDRTConnectionsByFZIStar','Argument must be a structure.');    
end

% copying fundamental properties 
PHI = P.PHI;    
PHIZ = P.PHIZ;  
KX = P.KX;  
KY = P.KY;  
KZ = P.KZ;

if strcmp(base,'ln')   
            
    LOG_SQRT_PHI = P.LNFZIStar_SQRTPHI;
    logsqrtphi = 'lnsqrtphi,'; % .csv column name
    
elseif strcmp(base,'log10')          
    
    LOG_SQRT_PHI = P.Log10FZIStar_SQRTPHI;
    logsqrtphi = 'log10sqrtphi,'; % .csv column name
    
end

% copying according to ave and log-base
switch ave
    
    case 'arithmetic' 
        
        KMEAN = P.KA;              
        FZIStar = P.FZIAStar;  
        
        % .csv column names
        kmean = 'ka,';         
        fziStar = 'FZIAStar,';
        
        if strcmp(base,'ln')    
                        
            DRTStar = P.DRTAStar_LN;                                    
            LOG_SQRT_K = P.LNFZIStarA_SQRTK;
            logsqrtk = 'lnSQRTK_A';
    
        elseif strcmp(base,'log10')                 
            
            DRTStar = P.DRTAStar_LOG10;                                    
            LOG_SQRT_K = P.Log10FZIStarA_SQRTK;
            logsqrtk = 'log10SQRTK_A';
    
        else
        
            error('wellps:findDRTConnectionsByFZIStar','Log base not available.');

        end

    case 'geometric' 
        
        KMEAN = P.KG;        
        FZIStar = P.FZIGStar;  

        
        % .csv column names
        kmean = 'kg,';         
        fziStar = 'FZIGStar,';
        
        if strcmp(base,'ln') 
                        
            DRTStar = P.DRTGStar_LN;                                    
            LOG_SQRT_K = P.LNFZIStarG_SQRTK;
            logsqrtk = 'lnSQRTK_G';
    
        elseif strcmp(base,'log10') 
                
            DRTStar = P.DRTGStar_LOG10;                                    
            LOG_SQRT_K = P.Log10FZIStarG_SQRTK;
            logsqrtk = 'log10SQRTK_G';
    
        else
        
            error('wellps:findDRTConnectionsByFZIStar','Log base not available.');

        end
    
    case 'normalized' 
        
        KMEAN = P.KN;    
        FZIStar = P.FZINStar;  
        
        % .csv column names
        kmean = 'kn,'; 
        fziStar = 'FZINStar,';

        
        if strcmp(base,'ln')     
                                    
            DRTStar = P.DRTNStar_LN;                                    
            LOG_SQRT_K = P.LNFZIStarN_SQRTK;
            logsqrtk = 'lnSQRTK_N';
    
        elseif strcmp(base,'log10') 
                
            DRTStar = P.DRTNStar_LOG10;                                    
            LOG_SQRT_K = P.Log10FZIStarN_SQRTK;
            logsqrtk = 'logSQRTK_N';
    
        else
        
            error('wellps:findDRTConnectionsByFZIStar','Log base not available.');

        end
        
    case 'quadratic'
        
        KMEAN = P.KQ;    
        FZIStar = P.FZIQStar;  
        
        % .csv column names
        kmean = 'kq,'; 
        fziStar = 'FZIQStar,';

        
        if strcmp(base,'ln')     
                                    
            DRTStar = P.DRTQStar_LN;                                    
            LOG_SQRT_K = P.LNFZIStarQ_SQRTK;
            logsqrtk = 'lnSQRTK_Q';
    
        elseif strcmp(base,'log10') 
                
            DRTStar = P.DRTQStar_LOG10;                                    
            LOG_SQRT_K = P.Log10FZIStarQ_SQRTK;
            logsqrtk = 'logSQRTK_Q';
    
        else
        
            error('wellps:findDRTConnectionsByFZIStar','Log base not available.');

        end
        
    case 'harmonic'
        
        KMEAN = P.KH;    
        FZIStar = P.FZIHStar;  
        
        % .csv column names
        kmean = 'kh,'; 
        fziStar = 'FZIHStar,';

        
        if strcmp(base,'ln')     
                                    
            DRTStar = P.DRTHStar_LN;                                    
            LOG_SQRT_K = P.LNFZIStarH_SQRTK;
            logsqrtk = 'lnSQRTK_H';
    
        elseif strcmp(base,'log10') 
                
            DRTStar = P.DRTHStar_LOG10;                                    
            LOG_SQRT_K = P.Log10FZIStarH_SQRTK;
            logsqrtk = 'logSQRTK_H';
    
        else
        
            error('wellps:findDRTConnectionsByFZIStar','Log base not available.');

        end
        
        
end

if ~all(drt)
    warning('wellps:findDRTConnectionsByFZIStar','DRT = 0 was found in your list. This will result in null porosity blobs.');
end

if ~ischar(tocsv)
    error('wellps:findDRTConnectionsByFZIStar','argument tocsv must be a char: "y" [yes] or "n" [no].')    
end

if dv ~= 1 
    error('wellps:findDRTConnectionsByFZIStar','Neighbourhood criterion recommended for now is 6. This means a cell distance = 1.')
end


% file header used in the loop
hdr = {'i,'; 'j,'; 'k,'; 'phi_e,'; ...
        'kx,'; 'ky,'; 'kz,'; ...
        kmean; ...
        fziStar; ...
        logsqrtphi; logsqrtk}'; % transposed!
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
DRT = DRTStar;

tstart = tic; % timing
fprintf('Sweeping field...\n')

for val = drt(1):drt(end)
    
    fprintf('=> DRT = %d \n',val);
    
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
        fprintf('==> No connections found for DRT = %d... \n',val);
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
    drtSt.value = val;                        % DRT value
    drtSt.averaging = ave;                    % averaging technique used
    drtSt.logBase = base;                     % log base used
    drtSt.allAdjMatrix = MDadj;               % graph adjacency matrix
    drtSt.allAdjEdgeList = edgeList;          % graph edge list
    drtSt.allVoxelCoords = coordsDRT;         % cell coordinates (i,j,k)     
    drtSt.allVoxelInds = indz;                % cell linear indices
    drtSt.allPHI = PHI( indz );               % effective porosity
    drtSt.allKX = KX( indz );                 % permeability-x
    drtSt.allKY = KY( indz );                 % permeability-y
    drtSt.allKZ = KZ( indz );                 % permeability-z
    drtSt.allKMEAN = KMEAN( indz );           % average permeability     
    drtSt.allFZIStar = FZIStar( indz );       % FZI*
    drtSt.allLogSQRTPHI = LOG_SQRT_PHI(indz); % Log(sqrt(PHI))
    drtSt.allLogSQRTK = LOG_SQRT_K( indz );   % Log(0.0314sqrt(K))
    drtSt.allNComps = ncomp;                  % number of components in the graph

    % mount matrix to export
    mat = [ coordsDRT(:,1)  ... 
            coordsDRT(:,2)  ...
            coordsDRT(:,3)  ...
            PHI( indz )     ...
            KX( indz )      ...
            KY( indz )      ...
            KZ( indz )      ...
            KMEAN( indz )   ...
            FZIStar( indz ) ...            
       LOG_SQRT_PHI( indz ) ...
       LOG_SQRT_K( indz )   ];            
            
    if tocsv == true
        
        % preparing csv file
        fname = strcat(dobj.getCsvDir,'/table-field_DRTStar_',ave,'_',base,'_',num2str( val ),'.csv');    
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
        drtSt.compVoxelCoords{idcomp} = aux; 
        drtSt.compVoxelInds{idcomp} = indz( members{idcomp} );          
        drtSt.compNNodes{idcomp} = compSizes(idcomp);           
        drtSt.compPHI{idcomp} = PHI( indz( members{idcomp} ) );
        drtSt.compKX{idcomp} = KX( indz( members{idcomp} ) );
        drtSt.compKY{idcomp} = KY( indz( members{idcomp} ) );
        drtSt.compKZ{idcomp} = KZ( indz( members{idcomp} ) );
        drtSt.compKMEAN{idcomp} = KMEAN( indz( members{idcomp} ) );        
        drtSt.compFZIStar{idcomp} = FZIStar( indz( members{idcomp} ) );
        drtSt.compLogSqrtPHI{idcomp} = LOG_SQRT_PHI( indz( members{idcomp} ) );
        drtSt.compLogSqrtK{idcomp} = LOG_SQRT_K( indz( members{idcomp} ) ); 
                  
        % matrix to export
        mat = [ aux(:,1)                         ... 
                aux(:,2)                         ...
                aux(:,3)                         ...
                PHI( indz( members{idcomp} ) )   ...
                 KX( indz( members{idcomp} ) )   ...
                 KY( indz( members{idcomp} ) )   ...
                 KZ( indz( members{idcomp} ) )   ...
              KMEAN( indz( members{idcomp} ) )   ...                               
            FZIStar( indz( members{idcomp} ) )   ...
       LOG_SQRT_PHI( indz( members{idcomp} ) )   ...
         LOG_SQRT_K( indz( members{idcomp} ) )  ];
           
        
        if (tocsv == true) && (compSizes(idcomp) >= nofsc)              
            % preparing csv file
            fname = strcat(dobj.getCsvDir,'/table-cluster_',num2str(idcomp),'_DRTStar_',ave,'_',base,'_',num2str( val ),'.csv');        
            dlmwrite(fname,hdr,'');       
            
            % append matrix 
            dlmwrite(fname,mat,'-append'); 
            
            % count clusters with >= nofsc
            cnofs = cnofs + 1;
        end                        
    end    
    fprintf(' %d clusters with >= %d cells.\n',cnofs,nofsc);            

    % saving structure to .mat     
    save( strcat(dobj.getMatDir,'/DRTStar_',ave,'_',base,'_',num2str( val ),'.mat'),'drtSt'); % saving
    %fprintf('----> DRTStar_%s.mat file saved. \n',num2str(val));
    
    % dynamic attribution     
    f = strcat('DRTStar',num2str(val));
    C.(f) = drtSt;    
    
    clear drtSt;    % frees to recompute
    
end

tfinal = toc(tstart);
fprintf('findDRTConnectionsByFZIStar finished after %g seconds. \n',tfinal);
