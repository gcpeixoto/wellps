
function vcn = getVoxelNeigh(vc,neigh)
%GETVOXELNEIGH Assembles matrix of voxels belonging to the 
%               neighbourhood of a given seed voxel 
%
% input: 
%       - vc : seed voxel (1x3 array)
%       - neigh: number of neighbours (6 or 26)
%
% output: 
% 
%       - vcn : neighbour voxel array (6x3 or 26x3) 
%
% Gustavo Peixoto, @LaMEP

% check
assert(all(size(vc) == [1,3]),'seed voxel must be a 1x3 array.');

% neighbourhood
switch neigh
    
    case '26'
                                % LOCAL INDICES 
        %      xc, yc, zc ;     % 1   <-- seed voxel 
        N26 = [-1, -1, -1 ; ... % 2   ---
                0, -1, -1 ; ... % 3    |
                1, -1, -1 ; ... % 4    |
               -1,  0, -1 ; ... % 5    |
                0,  0, -1 ; ... % 6    |  upper layer (z = zc - 1)
                1,  0, -1 ; ... % 7    |
               -1,  1, -1 ; ... % 8    |
                0,  1, -1 ; ... % 9    |
                1,  1, -1 ; ... % 10  ---                
               -1, -1,  0 ; ... % 11  ---
                0, -1,  0 ; ... % 12   |
                1, -1,  0 ; ... % 13   |  
               -1,  0,  0 ; ... % 14   |  mid layer (z = zc)               
                1,  0,  0 ; ... % 15   |               
               -1,  1,  0 ; ... % 16   |
                0,  1,  0 ; ... % 17   |
                1,  1,  0 ; ... % 18  ---   
               -1, -1,  1 ; ... % 19  ---
                0, -1,  1 ; ... % 20   |
                1, -1,  1 ; ... % 21   |
               -1,  0,  1 ; ... % 22   |
                0,  0,  1 ; ... % 23   |  lower layer (z = zc + 1)
                1,  0,  1 ; ... % 24   |
               -1,  1,  1 ; ... % 25   |
                0,  1,  1 ; ... % 26   |
                1,  1,  1 ];    % 27  ---     
        
        vcn = repmat(vc,[26,1]) + N26;               
        
        
        case '6'
                               % LOCAL INDICES 
        %     xc, yc, zc ;     % 1   <-- seed voxel 
        N6 = [ 0,  0, -1 ; ... % 6   ---  upper layer (z = zc - 1)                               
               0, -1,  0 ; ... % 12  ---                
              -1,  0,  0 ; ... % 14   |  mid layer (z = zc)               
               1,  0,  0 ; ... % 15   |                              
               0,  1,  0 ; ... % 17  ---               
               0,  0,  1 ] ;   % 23   |  lower layer (z = zc + 1)                
        
        vcn = repmat(vc,[26,1]) + N6;               

end

