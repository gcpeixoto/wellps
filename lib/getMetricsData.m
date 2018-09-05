function [node,deg,cln,bet] = getMetricsData(edfile)
% GETMETRICSDATA read file with centrality data computed from SNAP.

mtfile = '../tmp/metrics.txt';
tab = importdata(mtfile);

% getting centralities from the file computed via SNAP: 
% node id, degree, closeness, betweeness
Mtab = tab.data(:,1:4);

node = Mtab(:,1);
deg  = Mtab(:,2);
cln  = Mtab(:,3);
bet  = Mtab(:,4);

% clean the temporary files
delete(edfile,mtfile)

end