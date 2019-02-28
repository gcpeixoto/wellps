function clusterFitSt = processClusterFit(G,drtValue,idcomp,ave,base,theta)
% PROCESSCLUSTERFIT Determines fit parameters in order to setup nonuniform 
% (5-spot) well patterns based on a best-fit ellipsoid obtained 
% from the convex hull points of the selected cluster(s).
%
% METHODOLOGY
% ===========
%
% - Consider G a cluster with arbitrary shape and c the 
%   voxels (vertices) that belong to the convex hull H of G 
%
% - The best-fit ellipsoid based on c points returns 
%   the 3 main direction vectors, the 3 radii a,b,c of each axis 
%   and the center point 
%
%      o X1                      X1: fit-point (+x)                
%                                X2: fit-point (-x)
%        P1-c                    Y1: fit-point (+y)                
%        |  |       o Y1         Y2: fit-point (-y)
%        c  c--c-P2              Z1: fit-point (+z) (not used for vertical)
%        |    C   |              Z2: fit-point (-z) (not used for vertical)
%  o Y2  P4-c  c--c
%           |  |                  C: center
%           c-P3
%                                Pk are such that min {d(Xk,H)} in L2-norm
%                 o X2
%  RELATIONS
%  ========
%
%  X1 = C + a*E1        X2 = C - a*E1
%  Y1 = C + b*E2        Y2 = C - b*E2
%  Z1 = C + c*E3        Z2 = C - c*E3, for orthonormal basis {E1,E2,E3}
%
%  NONUNIFORM 5-SPOT PATTERNS
%  ==========================
% 
%  - Formed by the PRODUCER column placed at the max closeness 
%    vertex M = (Mx,My,Mz) z-span, i.e. (Mx,My,k), k \in klims.
%    + 4 INJECTOR columns placed at 
%    Pj = (Pj_x,Pj_y,k), k \in klims, j = 1,2,3,4. 
%  
%  ROTATED CONFIGURATIONS
%  ======================
%
%  - This script also computes this 5-spot pattern rotatated 
%    by an angle theta in relation to the reference frame 
%    {C,X1<->X2,Y1<->Y2}
%
%

%% Load cluster info

% create output dir
outDir = '../mat/ClusterFitData/';
if exist(outDir,'dir') ~= 7
    warning('wellps:processClusterFit',...
        ['Creating data directory into:',...
         '"../mat/ClusterFitData/".']);
    mkdir(outDir); 
end  

% compute grid geometry if input is not ready 
if ~any(ismember(fieldnames(G.cells),'centroids'))
    G = computeGeometry(G);
end

% inverse map 
globInd = nan(prod(G.cartDims),1);
globInd(G.cells.indexMap) = 1:G.cells.num;

% loop over DRT value
for drt = drtValue
    
    % try load .mat files required by user
    fidd = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'.mat');
    fidm = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'_Metrics.mat');
    fidr = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'_LinRegr.mat');
   
    if ~exist(fidd,'file') || ~exist(fidm,'file') || ~exist(fidr,'file')
        error('wellps:metricsAnalyzer',... 
            ['All the required files: \n %s \n %s \n were not found', ...
             'for DRT = %d. Have you run wellps:computeDRTGraphMetrics'], ...
             fidm,fidr,drt);
    else
        
        % load DRT struct
        load(fidd,'drtSt');                
        
        % load metrics
        load(fidm,'Maux');        
        
        % load linregr
        %load(fidr,'Laux');        
        
        % rename
        M = Maux; clear Maux;
        %L = Laux; clear Laux;      
        
    end

    fprintf('Sweeping DRT = %d ...\n',drt);
    
    % loop over components
    for c = idcomp
        
        fprintf('-----> Sweeping cluster = %d ...\n',c);
        
        % voxel coords and inds
        cvc = drtSt.compVoxelCoords{c};         
        cvi = drtSt.compVoxelInds{c};        
                
        % coords by column
        iG = cvc(:,1); jG = cvc(:,2); kG = cvc(:,3);
        
        % max closeness coords
        cvcmaxc = M.maxClosenessVoxelCoords{c};
        mcx = cvcmaxc(1); mcy = cvcmaxc(2); mcz = cvcmaxc(3); 
                
        % maxC same-column neighbours (inclusive)
        zz = unique(kG); 
        lz = numel(zz);
        cols = [ones(lz,1)*mcx, ones(lz,1)*mcy, zz]; 
        
        % cluster voxel's global indices
        cvgind = globInd(cvi);
        
        % centroid coordinates of cluster cells
        xc = G.cells.centroids(cvgind,1);
        yc = G.cells.centroids(cvgind,2);
        zc = G.cells.centroids(cvgind,3);
                                
        % Find convex hull in relation to centroids
        ch = convhull(xc,yc,zc);
        
        % Since 'ch' is a triangle list whose indices are 
        % over the convex hull, we get the unique list 
        % and recover their logical coords I,J,K
        ihull = unique(ch(:));
        cvchull = [iG(ihull),jG(ihull),kG(ihull)];

        % ellipsoid fit over convex hull points
        [C,R,E,P] = ellipsoid_fit(cvchull);

        % ellipsoid extrema points over principal axes
        P1 = C + R(1)*E(:,1); % +X
        P2 = C + R(2)*E(:,2); % +Y
        P3 = C + R(3)*E(:,3); % +Z

        P4 = C - R(1)*E(:,1); % -X
        P5 = C - R(2)*E(:,2); % -Y
        P6 = C - R(3)*E(:,3); % -Z

        % max closeness
        clo = M.closenessCentrality{c};          
        iMaxc = find(clo == max(clo));   
        vcMaxc = M.maxClosenessVoxelCoords{c};
        
        % store                
        clusterFitData.drtValue = drtSt.value;
        clusterFitData.ncomp = c;
        clusterFitData.iMaxc = iMaxc;
        clusterFitData.vcMaxc = vcMaxc;
        clusterFitData.colNeighsMaxC= cols;
        clusterFitData.iHull = ihull;
        clusterFitData.vcHull = cvchull;
        clusterFitData.ellipsoidFitCenter = C;
        clusterFitData.ellipsoidFitRadii = R;
        clusterFitData.ellipsoidFitEvecs = E;
        clusterFitData.ellipsoidFitParams = P;        
        
        % cluster limits
        im = min(iG);  jm = min(jG);  km = min(kG); 
        clusterFitData.clusterBoundingBoxI = [min(iG),max(iG)];
        clusterFitData.clusterBoundingBoxJ = [min(jG),max(jG)];
        clusterFitData.clusterBoundingBoxK = [min(kG),max(kG)];
        
        % vcMaxc - local coords
        aux = global2LocalCluster(mcx,mcy,mcz,im,jm,km);
        clusterFitData.vcMaxcLocalCoords = aux;
       
        % loop over angle
        pattern = 1;
        for t = theta  
            
            % rotation angle (convert to degrees)
            clusterFitData.rotationAngleDeg = 180*t/pi;
            
            %% Find 6 farthest points in the cluster from ellipsoid fit points
                       
            % ----------
            % Rotation matrix in relation to an axis
            % - matrix:
            %      R = cos(theta)I + sin(theta)UX + (1-cos(theta))*[u dyad u]
            % - angle:
            %      theta
            % - axis direction unit vector: 
            %     (ux,uy,uz)
            % ----------
            ux = E(1,3); uy = E(2,3); uz = E(3,3); % direction vector
            
            % matrices
            UX = [ 0 -uz uy; 
                   uz 0 -ux;
                  -uy ux 0 ];
              
            UU = [ ux*ux ux*uy ux*uz; 
                   ux*uy uy*uy uy*uz;
                   ux*uz uy*uz uz*uz ];
               
            % rotation matrix
            ROT = cos(t)*eye(3) + sin(t)*UX + (1-cos(t))*UU;
            
            % xy
            [RC, RP1, RP2, RP4, RP5] = deal(ROT*C, ROT*P1, ROT*P2, ROT*P4, ROT*P5);
            
            % displacing points in relation to retake center
            D = RC - C;
            RP1 = RP1 - D;
            RP2 = RP2 - D;
            RP4 = RP4 - D;
            RP5 = RP5 - D;
            
            % distances 
            dc = zeros(size(cvchull,1),1);
            dcr = zeros(size(cvchull,1),1); % rotated 

            % directions
            pos = {'+X','-X','+Y','-Y','+Z','-Z'};
            
            for p = 1:numel(pos)
                
                switch pos{p}
                    
                    case '+X'
                        
                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P1 - P' ).^2 ) );                                        
                            dcr(i) = sqrt( sum( ( RP1 - P' ).^2 ) );                                        
                        end
                        
                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);
                        
                        [dcr,auxr] = sortrows(dcr);
                        ipr = auxr(1);
                        vcr = cvchull(ipr,:);
                        
                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
                        col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 
                        
                        %---- SAVE            
                        % save colNeighsX1
                        clusterFitData.colNeighsX1 = col_neighbours;
            
                        % save local coords - colNeighsX1 
                        aux = global2LocalCluster(col_neighbours(:,1),...
                                                  col_neighbours(:,2),...
                                                  col_neighbours(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsX1LocalCoords = aux;
            
                        % save index X1
                        clusterFitData.ipX1 = ip;                        
            
                        % save voxel X1
                        clusterFitData.vcX1 = vc;
            
                        % save local coords - vcX1 
                        aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
                        clusterFitData.vcX1LocalCoords = aux;
            
                        % save colNeighsX1R
                        clusterFitData.colNeighsX1R = col_neighboursR;
            
                        % local coords - colNeighsX1R 
                        aux = global2LocalCluster(col_neighboursR(:,1),...
                                                  col_neighboursR(:,2),...
                                                  col_neighboursR(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsX1RLocalCoords = aux;
            
                        % index X1R
                        clusterFitData.ipX1R = ipr;
            
                        % voxel X1R
                        clusterFitData.vcX1R = vcr;
            
                        % local coords - vcR
                        aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
                        clusterFitData.vcX1RLocalCoords = aux;
                                    
                    case '-X'

                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P4 - P' ).^2 ) );                                        
                            dcr(i) = sqrt( sum( ( RP4 - P' ).^2 ) );                                        
                        end

                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);

                        [dcr,auxr] = sortrows(dcr);
                        ipr = auxr(1);
                        vcr = cvchull(ipr,:);

                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
                        col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 

                        %---- SAVE

                        % save colNeighsX2
                        clusterFitData.colNeighsX2 = col_neighbours;

                        % save local coords - colNeighsX2
                        aux = global2LocalCluster(col_neighbours(:,1),...
                                                  col_neighbours(:,2),...
                                                  col_neighbours(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsX2LocalCoords = aux;

                        % save index X2
                        clusterFitData.ipX2 = ip;                        

                        % save voxel X2
                        clusterFitData.vcX2 = vc;

                        % save local coords - vcX2
                        aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
                        clusterFitData.vcX2LocalCoords = aux;

                        % save colNeighsX2R
                        clusterFitData.colNeighsX2R = col_neighboursR;

                        % local coords - colNeighsX2R 
                        aux = global2LocalCluster(col_neighboursR(:,1),...
                                                  col_neighboursR(:,2),...
                                                  col_neighboursR(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsX2RLocalCoords = aux;

                        % index X2R
                        clusterFitData.ipX2R = ipr;

                        % voxel X2R
                        clusterFitData.vcX2R = vcr;

                        % local coords - vcR
                        aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
                        clusterFitData.vcX2RLocalCoords = aux;
                            
                    case '+Y'

                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P2 - P' ).^2 ) );     
                            dcr(i) = sqrt( sum( ( RP2 - P' ).^2 ) );                                        
                        end

                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);

                        [dcr,auxr] = sortrows(dcr);
                        ipr = auxr(1);
                        vcr = cvchull(ipr,:);

                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
                        col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 

                        %---- SAVE

                        % save colNeighsY1
                        clusterFitData.colNeighsY1 = col_neighbours;

                        % save local coords - colNeighsY1 
                        aux = global2LocalCluster(col_neighbours(:,1),...
                                                  col_neighbours(:,2),...
                                                  col_neighbours(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsY1LocalCoords = aux;

                        % save index Y1
                        clusterFitData.ipY1 = ip;                        

                        % save voxel Y1
                        clusterFitData.vcY1 = vc;

                        % save local coords - vcY1 
                        aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
                        clusterFitData.vcY1LocalCoords = aux;

                        % save colNeighsY1R
                        clusterFitData.colNeighsY1R = col_neighboursR;

                        % local coords - colNeighsY1R 
                        aux = global2LocalCluster(col_neighboursR(:,1),...
                                                  col_neighboursR(:,2),...
                                                  col_neighboursR(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsY1RLocalCoords = aux;

                        % index Y1R
                        clusterFitData.ipY1R = ipr;

                        % voxel Y1R
                        clusterFitData.vcY1R = vcr;

                        % local coords - vcR
                        aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
                        clusterFitData.vcY1RLocalCoords = aux;

                    case '-Y'

                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P5 - P' ).^2 ) );  
                            dcr(i) = sqrt( sum( ( RP5 - P' ).^2 ) );                                        
                        end

                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);

                        [dcr,auxr] = sortrows(dcr);
                        ipr = auxr(1);
                        vcr = cvchull(ipr,:);

                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
                        col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 

                        %---- SAVE

                        % save colNeighsY2
                        clusterFitData.colNeighsY2 = col_neighbours;

                        % save local coords - colNeighsY2
                        aux = global2LocalCluster(col_neighbours(:,1),...
                                                  col_neighbours(:,2),...
                                                  col_neighbours(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsY2LocalCoords = aux;

                        % save index Y2
                        clusterFitData.ipY2 = ip;                        

                        % save voxel Y2
                        clusterFitData.vcY2 = vc;

                        % save local coords - vcY2
                        aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
                        clusterFitData.vcY2LocalCoords = aux;

                        % save colNeighsY2R
                        clusterFitData.colNeighsY2R = col_neighboursR;

                        % local coords - colNeighsY2R 
                        aux = global2LocalCluster(col_neighboursR(:,1),...
                                                  col_neighboursR(:,2),...
                                                  col_neighboursR(:,3),...
                                                  im,jm,km);
                        clusterFitData.colNeighsY2RLocalCoords = aux;

                        % index Y2R
                        clusterFitData.ipY2R = ipr;

                        % voxel Y2R
                        clusterFitData.vcY2R = vcr;

                        % local coords - vcR
                        aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
                        clusterFitData.vcY2RLocalCoords = aux;

                    case '+Z'

                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P3 - P' ).^2 ) );                                        
                        end

                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);

                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 

                        clusterFitData.colNeighsZ1 = col_neighbours;
                        clusterFitData.ipZ1 = ip;
                        clusterFitData.vcZ1 = vc;

                    case '-Z'

                        % span convex hull
                        for i = 1:size(cvchull,1) 
                            P = cvchull(i,:);       
                            dc(i) = sqrt( sum( ( P6 - P' ).^2 ) );                                        
                        end

                        % sort distances and get minimum one
                        [dc,aux] = sortrows(dc);
                        ip = aux(1);
                        vc = cvchull(ip,:);

                        % column neighbours inside the cluster z-range
                        zz = unique(kG); 
                        lz = length(zz);
                        col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 

                        clusterFitData.colNeighsZ2 = col_neighbours;
                        clusterFitData.ipZ2 = ip;
                        clusterFitData.vcZ2 = vc;                            
                
                end
                
            end                                    
            %% Save multilevel structure
            % The multilevel structure is created as follows: 
            %
            % - 1st level: DRT{X}, where X is the DRT number            
            % - 2nd level: C{Y}, where Y is the component number 
            % - 3rd level: Pattern{Z}, where Z is the pattern number
            % according to the angle chosen. 
            % 
            % Then, we have something like:
            %
            % STRUCTURE LAYOUT
            % ================
            %
            % DRT{X},       X = 1,2,...,D
            %  |
            %  |--- C{1}
            %  |  |
            %  |  |--- Pattern{1}
            %  |  |--- Pattern{2}
            %  |    .
            %  |    .
            %  |    .
            %  |  |--- Pattern{z}
            %  |
            %  |--- C{2}
            %  |  |
            %  |  |--- Pattern{1}
            %  |  |--- Pattern{2}
            %  |    .
            %  |    .
            %  |    .            
            %  |  |--- Pattern{z}
            %  |
            %  .
            %  .
            %  .
            %  |--- C{Y}
            %  |  |
            %  |  |--- Pattern{1}
            %  |  |--- Pattern{2}
            %  |    .
            %  |    .
            %  |    .            
            %  |  |--- Pattern{z}
            %
            %
            % REMARK: since 'angle' cannot be always converted to a string
            % whose value is an integer, we are using here the field 
            % 'pattern' for which the angle is to be assigned. That is to
            % say, pattern = 1 for angle{1}; = 2 for angle{2} and so on...
            % The best way to know the angle for the pattern is to pull up
            % the field 'angle' of the minor structure
            % DRT{X}.C{Y}.Pattern{Z}.                                    
            clusterFitSt.(strcat('DRT',num2str(drt))). ... 
                         (strcat('C',num2str(c))).     ...
                         (strcat('Pattern',num2str(pattern))) = clusterFitData;
            
            % update counter for angle  
            pattern = pattern + 1;
            
        end
        
    end               
    
    % save to .mat (separated by DRT)
    save( strcat(outDir,'DRT_',ave,'_',base,'_',num2str( drt ),'_clusterFitData','.mat'),'clusterFitSt');
end

