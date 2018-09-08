function loc = global2LocalCluster( xglob,yglob,zglob,xm,ym,zm)
% GLOBAL2LOCALCLUSTER coordinate transform for cluster coordinates
%
% PARAMETERS:
%       - 	xglob,yglob,zglob: voxel's global coordinates (whole field)
%
%       -       xm,ym,zm: cluster's initial limits                        
%
% RETURNS:
%
%       -	loc: voxel local coordinates at cluster (1x3)
%
% REMARK:
%	This function is specific for use in CMG suite, because
%	when we make grid cut operations mainly in Builder,
%	the corner grid region is reindexed. 
%
xloc = xglob - xm + 1;
yloc = yglob - ym + 1;
zloc = zglob - zm + 1;

loc = [xloc,yloc,zloc];


