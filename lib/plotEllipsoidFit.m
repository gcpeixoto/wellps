function [figelip,figres] = plotEllipsoidFit(G,drtSt,clusterFitSt,drt,comp,pattern,factor,nres)
%PLOTELLIPSOIDFIT Special plotting for ellipsoid fitting
%
% PARAMETERS:
%
%   G     - MRST' grid structure of the reservoir model where the data 
%           will be plotted 
%   
%   drtSt - WELLPS' structure having the DRT connections for the wanted 
%           DRT
%
%   clusterFitSt - WELLPS' structure having the cluster fit information.
%                  See wellps::processClusterFit and/or file 
%                  'cases/case_clusterFit.m'.
%
%   drt   - DRT value of interest to plot 
%
%   comp  - cluster ID (connected component of the respective DRT) 
%
%   pattern - non-uniform well pattern to be considered for plotting. 
%             There are usually 4 main choices: 
%             
%               1:  0 degree
%               2: 30 degrees
%               3: 45 degrees
%               4: 60 degrees
%
%             See wellps:processClusterFit. 
%
%   factor  - number to control the stretching/shrinking of meshgrid 
%             underlying arrays to plot the 3D ellipsoid. 
%             This number should be tested from case to case to see 
%             if a smaller or larger number is good to cover all the 
%             ellipsoid' isosurface.
%
%   nres    - resolution of the surface mesh of the ellipsoid. 
%             Well tested with value = 70. 
%
%
% OUTPUT:
%   
%   figelip, figres - handles for separated figures. 
%     
% EXAMPLE:
%
% drtSt = load('../mat/DRT_harmonic_ln_13.mat');
% clusterFit = load('../mat/unisim1DRT_harmonic_ln_13_clusterFitData.mat');
% drt = 13;  comp = 2; pattern = 1; factor = 40; nres = 70; 
% [figelip,figres] = ...
% plotEllipsoidFit(G,drtSt,clusterFitSt,drt,comp,pattern,factor,nres);
%
% 
% Dr. Gustavo Oliveira, @LaMEP/UFPB
%
%

% tags
dtag = ['DRT' num2str(drt)];
ctag = ['C' num2str(comp)];
ptag = ['Pattern' num2str(pattern)];

% inverse map 
globInd = nan(prod(G.cartDims),1);
globInd(G.cells.indexMap) = 1:G.cells.num;

% voxel coords
cvc = drtSt.compVoxelCoords{comp};
cvch = clusterFitSt.(dtag).(ctag).(ptag).vcHull;

% cell indices
cvi = sub2ind(G.cartDims,cvc(:,1),cvc(:,2),cvc(:,3));
cvih = sub2ind(G.cartDims,cvch(:,1),cvch(:,2),cvch(:,3));

% max closeness cell coord
MC = clusterFitSt.(dtag).(ctag).(ptag).vcMaxc;
cvimc = sub2ind(G.cartDims,MC(1),MC(2),MC(3));

% anchor points injector coords 
P1 = clusterFitSt.(dtag).(ctag).(ptag).vcX1;
P2 = clusterFitSt.(dtag).(ctag).(ptag).vcX2;
P3 = clusterFitSt.(dtag).(ctag).(ptag).vcY1;
P4 = clusterFitSt.(dtag).(ctag).(ptag).vcY2;

cviP1 = sub2ind(G.cartDims,P1(1),P1(2),P1(3));
cviP2 = sub2ind(G.cartDims,P2(1),P2(2),P2(3));
cviP3 = sub2ind(G.cartDims,P3(1),P3(2),P3(3));
cviP4 = sub2ind(G.cartDims,P4(1),P4(2),P4(3));

% mapping
cvi = globInd(cvi);
cvih = globInd(cvih);
cvimc = globInd(cvimc);
cviP1 = globInd(cviP1);
cviP2 = globInd(cviP2);
cviP3 = globInd(cviP3);
cviP4 = globInd(cviP4);

%% Draw Ellipsoid
% REMARK: the 3D visualization of the ellipsoid and decoration 
% is built on top of the logical coordinates, but the original fitting 
% is made in relation to the centroids of the convex hull cells.

% figure handle
figelip = figure;

% ellipsoid mesh
mind = min([cvc(:,1),cvc(:,2),cvc(:,3)]);
maxd = max([cvc(:,1),cvc(:,2),cvc(:,3)]);               
step = ( maxd - mind ) / nres;

[ x, y, z ] = ...
    meshgrid(linspace( mind(1) - factor*step(1), maxd(1) + factor*step(1), nres ), ... 
             linspace( mind(2) - factor*step(2), maxd(2) + factor*step(2), nres ), ... 
             linspace( mind(3) - factor*step(3), maxd(3) + factor*step(3), nres ) );

% parameters         
P = clusterFitSt.(dtag).(ctag).(ptag).ellipsoidFitParams;
                 
Ellipsoid = P(1) *x.*x +   P(2) * y.*y + P(3) * z.*z + ...
           2*P(4) *x.*y + 2*P(5)*x.*z + 2*P(6) * y.*z + ...
           2*P(7) *x    + 2*P(8)*y    + 2*P(9) * z;
 
% isosurface 
p = patch( isosurface( x, y, z, Ellipsoid, -P(10),'verbose' ) );
p.FaceColor = [0.7,0.7,0.0];
p.EdgeColor = 'none'; 
p.FaceAlpha = 0.1;
axis off vis3d;     
camlight;      
lighting phong;
p.FaceLighting = 'gouraud';

camproj perspective

view([-85,56])


%% Decoration

hold on
% ellipsoid fit center
C = clusterFitSt.(dtag).(ctag).(ptag).ellipsoidFitCenter;

% ellipsoid fit eigenvectors and radii
E = clusterFitSt.(dtag).(ctag).(ptag).ellipsoidFitEvecs;
R = clusterFitSt.(dtag).(ctag).(ptag).ellipsoidFitRadii;

% support points
X1 = C + R(1)*E(:,1); % +X
X2 = C - R(1)*E(:,1); % -X
Y1 = C + R(2)*E(:,2); % +Y
Y2 = C - R(2)*E(:,2); % -Y
Z1 = C + R(3)*E(:,3); % +Z
Z2 = C - R(3)*E(:,3); % -Z

% plot center
scatter3(C(1),C(2),C(3),'ks');

% plot ellipsoid support points
scatter3(X1(1),X1(2),X1(3),'r^','SizeData',50);
scatter3(X2(1),X2(2),X2(3),'rv','SizeData',50);

scatter3(Y1(1),Y1(2),Y1(3),'b^','SizeData',50);
scatter3(Y2(1),Y2(2),Y2(3),'bv','SizeData',50);

% z
%scatter3(Z1(1),Z1(2),Z1(3),'k^','SizeData',50);
%scatter3(Z2(1),Z2(2),Z2(3),'kv','SizeData',50);

% max closeness point
scatter3(MC(1),MC(2),MC(3),'ko','filled','SizeData',200);

% principal axes (original ellipsoid)
lx = line([X1(1),X2(1)],[X1(2),X2(2)],[X1(3),X2(3)]); 
ly = line([Y1(1),Y2(1)],[Y1(2),Y2(2)],[Y1(3),Y2(3)]);
lz = line([Z1(1),Z2(1)],[Z1(2),Z2(2)],[Z1(3),Z2(3)]);

lx.Color = 'r';      ly.Color = 'b';      lz.Color = 'k';
lx.LineStyle = '--'; ly.LineStyle = '--'; lz.LineStyle= '--';

% plot ellipsoid anchor points
scatter3(P1(1),P1(2),P1(3),'r^','filled','SizeData',100);
scatter3(P2(1),P2(2),P2(3),'rv','filled','SizeData',100);

scatter3(P3(1),P3(2),P3(3),'b^','filled','SizeData',100);
scatter3(P4(1),P4(2),P4(3),'bv','filled','SizeData',100);

% nonuniform 5-spot lines

%{
l_mcp1 = line([MC(1),P1(1)],[MC(2),P1(2)],[MC(3),P1(3)]); 
l_mcp2 = line([MC(1),P2(1)],[MC(2),P2(2)],[MC(3),P2(3)]); 
l_mcp3 = line([MC(1),P3(1)],[MC(2),P3(2)],[MC(3),P3(3)]); 
l_mcp4 = line([MC(1),P4(1)],[MC(2),P4(2)],[MC(3),P4(3)]); 

l_mcp1.Color = 'r';      
l_mcp1.LineStyle = '-';  
l_mcp1.LineWidth = 2;

l_mcp2.Color = 'r';      
l_mcp2.LineStyle = '-';  
l_mcp2.LineWidth = 2;

l_mcp3.Color = 'b';      
l_mcp3.LineStyle = '-';  
l_mcp3.LineWidth = 2;

l_mcp4.Color = 'b';      
l_mcp4.LineStyle = '-';  
l_mcp4.LineWidth = 2;
%}

legelip = legend('$\mathcal{E}$','${\bf c}$','${\bf q}_1$','${\bf q}_2$','${\bf q}_3$','${\bf q}_4$','$M_{D,q}$','$E_1$','$E_2$','$E_3$',...
                 '$\check{\bf q}_1$','$\check{\bf q}_2$','$\check{\bf q}_3$','$\check{\bf q}_4$'); 
legelip.Location = 'bestoutside';
legelip.Interpreter = 'latex';
legelip.FontSize = 16;
legelip.EdgeColor = 'none';
legelip.Color = [0.9,0.9,0.9];


%print('elip.eps','-depsc2')

%% Cluster domain

figres = figure;
axis off;
pg1 = plotGrid(G,cvi,'FaceColor',[0.8,0.8,0.8],...
               'FaceAlpha',0.6,...
               'EdgeColor',[0.5,0.5,0.5],...
               'EdgeAlpha',0.3);
           
hold on

pg2 = plotGrid(G,cvih,'FaceColor',[0.7,0.7,0.0],...
                'FaceAlpha',0.7,...
                'EdgeColor','none');

% -------------- PLOT WELL
% Here, we're going to generate a fake rock model to plot wells 

fake_kx = (1:G.cells.num)';
rk = makeRock(G,[fake_kx,fake_kx,fake_kx],0.1);

W = addWell([],G,rk,cvimc,'Name','PROD');      
W = addWell(W,G,rk,cviP1,'Name','INJ1');  
W = addWell(W,G,rk,cviP2,'Name','INJ2');  
W = addWell(W,G,rk,cviP3,'Name','INJ3');  
W = addWell(W,G,rk,cviP4,'Name','INJ4');  

plotWell(G,W,'color','k');

pg1.FaceLighting = 'flat';  pg2.FaceLighting = 'flat';

view([-90,65])

%print('wpd.eps','-depsc2')

% ==== DISPLACED REFERENCE FRAME
%{
% displacement line 
line([C(1),MC(1)],[C(2),MC(2)],[C(3),MC(3)])

% displaced frame 
dx = C(1) - MC(1);
dy = C(2) - MC(2);
dz = C(3) - MC(3);

D = [dx,dy,dz]';
%}


end

