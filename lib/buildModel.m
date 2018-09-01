function [G,GAll,PROPS] = buildModel(f)
% Read data file and build reservoir grid
%
% SYNOPSIS: 
%   [G,PROPS] = buildModel('benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
%
% PARAMETERS:
%   f       - String holding name of grid file. 
%             Here, we are processing two source files:
%             
%               - .DAT file: a CMG IMEX file (TODO)
%
%               - .DATA file: a ECLIPSE deck file. 
%                  The workflow is the same as MRST's.
%                  No further function is implemented here. 
%
% RETURNS:
%  G     - structure already processed with processGRDECL
%         containing the GRID information (ACTNUM = 0 and 1)
%
%  GAll  - structure already processed with processGRDECL
%         containing the GRID information (only ACTNUM = 1)
%
%  PROPS - structure containing petrophysical parameters 
%        (porosity, permeability x,y,z, etc.)
%
%

%{
    Developed at LaMEP/UFPB, Brazil
    @gcpeixoto

%}

[~,~,ext] = fileparts(f);

% properties
PROPS.PHI = [];        
PROPS.KX = [];        
PROPS.KY = [];
PROPS.KZ = [];        
PROPS.ACTNUM = [];                   

fprintf('Reading file: %s\n',f);
% parsing IMEX file

switch ext

% IMEX
case '.dat'
    %TODO

% ECLIPSE
case '.DATA'

    % read and process grid with MRST
    [G,~] = readGRDECL(f);
    
    % properties
    PROPS.PHI = G.PORO;        
    PROPS.KX = G.PERMX;        
    PROPS.KY = G.PERMY;
    PROPS.KZ = G.PERMZ;        
    PROPS.ACTNUM = G.ACTNUM;                   

    % ============== ******* ==================
    % \TODO THINK ABOUT THIS. NOT SOLVED.   
    %
    % REMARK: 
    % GAll is a copy of the original G, but it is not processed
    % by 'processGRDECL' in the same way because we need to 
    % conserve all the cells activated in order to be able to 
    % plot clusters correctly. 
    %
    % Since the WELLPS is developed to accord with IMEX (I,J,K) 
    % indexing to set up simulations, we need to keep the integral 
    % grid indices. This can be done by avoiding the call of 
    % 'removeCells' inside 'processGRDECL'. 
    %
    % MRST ignores nonactive cells (ACTNUM = 0). Then we need 
    % to enforce full activation (ACTNUM = 1) and process 
    % this new grid. 
    % 
    % Once again, however, note that this step is made just for 
    % WELLPS interest. To plot other general features of the
    % reservoir, we should use the original G. Hence, GAll is 
    % like an 'abstract' grid over the original one.    
    %
    % ATTEMPTS:
    %
    % Delete kewords 'PINCH*' or set ACTNUM in brute force not
    % solved the issue. 
    GAll = G;
    GAll.ACTNUM = GAll.ACTNUM*0 + 1; 
    
    % process after to avoid overwriting of PROPS
    G = processGRDECL(G);                 
    
    % process the fully activated grid
    GAll = processGRDECL(GAll);                 
    
    
end