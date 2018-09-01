function [ VB, coords, inds ] = maskVoxels( V,val )
%MASKVOXELS makes a mask for the value 'val' over the 3D array 'V'.
% 
%     input: 
%         V: 3D array
%       val: specified value 
%     
%     output:
%        VB: masked array
%    coords: coords i,j,k whose entries are equal to val
%      inds: linear indices of 'coords'.
%             

VB = V == val;
[I,J,K,inds] = logic2subind(VB);
coords = [ I J K ];

end

