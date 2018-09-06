function [node,deg,cln,bet] = getMetricsData(edfile)
%GETMETRICSDATA read file with centrality data computed from SNAP.
%
% PARAMETERS:
%	edfile	-  graph edge file (see 'saveAdjEdges')
% 
% RETURNS:
%	node	- graph node ID array
%	deg	- degree centrality array
%	cln	- closeness centrality array
%	bet	- betweeness centrality array

% Assumes that SNAP has computed the centralities correctly and
% saved them to the temporary metrics.txt. 
mtfile = '../tmp/metrics.txt';
tab = importdata(mtfile);

% Getting centralities from the file computed via SNAP.
% Note that this file is a table containing all the centralities 
% computable by SNAP. However, we take the node ID and three of them 
% to our interest, namely, degree, closeness, betweeness, 
% in this order.

% interest table
Mtab = tab.data(:,1:4);

% node ID
node = Mtab(:,1);

% centralities
deg  = Mtab(:,2);
cln  = Mtab(:,3);
bet  = Mtab(:,4);

% TODO 
% This was added here to check if SNAP is producing negative betweeness. 
% I suspect that there exists an error in betweeness computation. 
% I need to check betweeness from other tool.  
if isempty(find(bet<0))
    warning('wellps:getMetricsData','Negative values were found for BETWEENNESS. Metric not reliable.');
end

% delete the temporary files
delete(edfile,mtfile)
