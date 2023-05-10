function [dbmin, dbmax] = computeDistBoundary(G,varargin)
%% COMPUTEDISTBOUNDARY Compute distances of cells to the MRST 
%                       grid boundary.
%
%   This method compute the min/max Euclidean distance d(c,b) of all grid cells
%   c to all boundary cells b.
%
%   If an additional cell index c0 is passed to the function, it will
%   compute d(c0,b).
%
%   Note that d(c,b) is equivalent to a symmetric matrix of distances D 
%   with diagonal entries D(i,i) = 0. On the other hand, d(c0,b) is a 1D
%   array. However, since the distances are computer for all reservoir
%   cells, the output is a 3D (G.cartDims) array.
%
%   PARAMETERS:
%
%       G        -  reservoir grid
%       c0       -  index of reference cell from which to compute distance
%
%   RETURNS:
%
%       dbmax   - maximum distances from cell to boundary cells (3D array)
%       dbmin   - minimum distances from cell to boundary cells (3D array)



% This gets all unique boundary cells
% See function 'MRST:boundaryFaceIndices', lines 134-135 
% called by function 'MRST:pside'
bndyfac = any(G.faces.neighbors == 0, 2);
bndycells = find(accumarray(sum(G.faces.neighbors(bndyfac,:), 2), 1) > 0);

gc = G.cells.centroids;
gb = G.cells.centroids(bndycells,:);


% case for all cells
if isempty(varargin)

    % initialization
    dbmin= gc(:,1); 
    dbmax= gc(:,1); 
    
    for i = 1:size(gc,1)            
        aux = [ gc(i,1) - gb(:,1), ...
                gc(i,2) - gb(:,2), ...
                gc(i,3) - gb(:,3)  ];
        dist = sqrt(sum(aux.*aux,2));
        dbmin(i) = min(dist);                         
        dbmax(i) = max(dist);                         
    end

% case for single cells
elseif length(varargin) == 1

    c0 = varargin{1}; % reference cell index
    gc0 = G.cells.centroids(c0,:);

    aux = [ gc0(1) - gb(:,1), ...
            gc0(2) - gb(:,2), ...
            gc0(3) - gb(:,3)  ];
    dist = sqrt(sum(aux.*aux,2));
    

    dbmin = min(dist);                         
    dbmax = max(dist);                         

end

% matrix
R1 = nan(G.cartDims);
R2 = nan(G.cartDims);

% assign values for active cells
[ilog,jlog,klog] = gridLogicalIndices(G);     
R1(sub2ind(G.cartDims,ilog,jlog,klog)) = dbmin;
R2(sub2ind(G.cartDims,ilog,jlog,klog)) = dbmax;

dbmin = R1;
dbmax = R2;

end