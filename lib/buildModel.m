function [G,PROPS] = buildModel(f)
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
%         containing the GRID information
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

    % process after to avoid overwriting of PROPS
    G = processGRDECL(G);                 
end