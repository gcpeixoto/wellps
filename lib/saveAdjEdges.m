function [ outfile ] = saveAdjEdges( Madj )
% SAVEADJEDGES write adjacency matrix edge list to file 
%   
%   input: 
%           Madj: adjacency matrix

% temporary edge file
outfile = '../tmp/edges.txt';

if issparse(Madj)
    [i,j] = find(Madj == 1);
    edges = [i,j];
    dlmwrite(outfile,edges,'delimiter',' ','precision','%d'); % writing file
else
    error('Matrix is not sparse');
end

