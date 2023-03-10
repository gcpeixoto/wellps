%% EXAMPLE: Visualization playground
%
%
%  REMARK: This example is for visualization purposes only.
%           

mrstVerbose off
mrstModule add mrst-gui

case_name = 'ex6';

%% Mounting

d = DirManager(); 

%% Reservoir model

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');
  
[G,PROPS] = buildModel(f);

G = computeGeometry(G);

%% Visualization

plotToolbar(G,G)