function vc6n = getVoxel6NeighSimple(vc)
% GETVOXEL6NEIGHSIMPLE gets the 6-neighbour coordinates of a seed voxel
%
%   input:
%          vc: seed-voxel coordinates
%
%   output:
%        vc6n: 6-neighbour voxel coordinates

assert(size(vc) ~= [1,3],'seed voxel must be a 1x3 array.');

% 6-neighbour array [+x,-x,+y,-y,+z,-z]'
vc6n = repmat(vc,[6,1]) + [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];


     
