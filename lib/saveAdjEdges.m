function outfile = saveAdjEdges(Madj)
% SAVEADJEDGES write adjacency matrix edge list to file 
%
% PARAMETERS: 
%       -   Madj: adjacency matrix
%
% RETURNS: 
%       - outfile: temporary edge file
%

d = DirManager; 
outfile = fullfile(d.getTmpDir,'edges.txt');
clear d;

if issparse(Madj)
    [i,j] = find(Madj == 1);    
    edges = [i,j];
    dlmwrite(outfile,edges,'delimiter',' ','precision','%d'); % writing file
else
    error('wellps:saveAdjEdges','Matrix must be sparse.');
end

end
