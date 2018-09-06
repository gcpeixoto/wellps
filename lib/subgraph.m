function adj_sub = subgraph(adj,S)
%SUBGRAPH Outputs the adjacency matrix of a subgraph 
% 
% PARAMETERS: 
%	adj 	-   supergraph adjacency matrix 
%	S   	-   vector of subgraph node indices
% 
% RETURNS: 
%	adj_sub - adjacency matrix of the subgraph 
%
% This function is based on MIT Strategic Engineering
% toolbox
adj_sub = adj(S,S);
