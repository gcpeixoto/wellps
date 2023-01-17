function [PUC,nclasses,delta,divs] = computePUC(J,active,method)
%COMPUTEPUC Compute productivity (or injecttivity) unit classes 
%           from the productivity (or injectivity) proxy function J 
%           by constructing the 3D PUC array prepared for finding 
%           connections.
%
% KNOWLEDGE BASE
%           
% # Update 2023: 
%
%   Regarding J
%   -----------
%
%   We now understand that J can be viewed as a multilinear
%   functional in space-time written in the form 
%   J(x,t) = S(x)D(x,t), where 
%
%   - S(x) holds the 'static' parameters of the reservoir
%   - D(x,t) holds the 'dynamic' parameters of the reservoir
%
%   The many faces admissible by J will lead to different 
%   physical interpretations.
%
%   See: J.W.L. Silva et al.'s conference paper
%   'ANALYSIS OF INJECTIVITY FUNCTIONALS BASED ON STATIC AND 
%   DYNAMIC PARAMETERS FOR CO2 STORAGE PROJECTS',
%   published on the XXV ENMC - Encontro Nacional de 
%   Modelagem Computacional, Brazil (Nov 2022).
%
%   Regarding 'divs'
%   ----------------
%   
%   The traditional Matlab methods determine the class partition points
%   by dividing the samples uniformly. We need to apply other class 
%   partition techniques to have nonuniformly spaced classes. For instance,
%   applying machine learning algorithms.
%
%   
% PARAMETERS:
%
%       J       - proxy function (3D double)
%       active  - active cell array indices (double)
%       method  - histogram binning method. (char)
%                 See available options below.
%
% RETURNS:
%
%       PUC     - mask with discrete values (like DRT) (3D array)
%      nclasses - number of PUC classes (double)
%       delta   - uniform width of separation (equivalent to BinWidth) (double)        
%       divs    - class partition points (equivalent to BinEdges) (1D array)


% methods
sh = 'shimazaki'; % Shimazaki
m = {'auto','scott','fd','sturges','sqrt'}; % standard Matlab 

% append
m{end+1} = sh;


% checking
if ~ismember(method,m)    
    mand = str2double(method); %convert
    assert( mand > 0,...
    sprintf(['Method ''%s'' must be one of the available binning methods.'...
    'Otherwise the input should be a positive integer number (char).\n'],method));        
    fprintf('----> PUC method set was manually set to %s.\n',method);    
end

% J only at active cells
Ja = J(active);

% positive values
Jplus = Ja(Ja > 0);

% net J: min(Jplus) <= Jplus <= max(Jplus)
Jnet = Jplus(min(Jplus) <= Jplus & Jplus <= max(Jplus));


% Shimazaki
if strcmp(method,sh)
    [nclasses, delta, divs, ~, ~] = sshist(Jnet);  
    
% Matlab methods
else
    h = histogram(Jnet,'BinMethod',method);              
    nclasses = length(h.BinCounts);
    delta = h.BinWidth;
    divs = h.BinEdges; 
        
    close all
    
end

% overwrite extrema to comply exactly with min(J), max(J) values
divs(1) = min(Jplus); divs(end) = max(Jplus); 

% finally compute 3D PUC array;
PUC = maskPUC(J,nclasses,divs);

end


% Shimazaki's method
% Extracted from: https://www.neuralengine.org/res/histogram.html#Matlab

function [optN, optD, edges, C, N] = sshist(x,N)
% [optN, optD, edges, C, N] = sshist(x,N)
%
% Function `sshist' returns the optimal number of bins in a histogram
% used for density estimation.
% Optimization principle is to minimize expected L2 loss function between 
% the histogram and an unknown underlying density function.
% An assumption made is merely that samples are drawn from the density
% independently each other.
%
% The optimal binwidth D* is obtained as a minimizer of the formula, 
% (2K-V) / D^2,
% where K and V are mean and variance of sample counts across bins with width D.
% Optimal number of bins is given as (max(x) - min(x)) / D*.
%
% For more information, visit 
% http://2000.jukuin.keio.ac.jp/shimazaki/res/histogram.html
%
% Original paper:
% Hideaki Shimazaki and Shigeru Shinomoto
% A method for selecting the bin size of a time histogram
% Neural Computation 19(6), 1503-1527, 2007
% http://dx.doi.org/10.1162/neco.2007.19.6.1503
%
% Example usage:
% optN = sshist(x); hist(x,optN);
%
% Input argument
% x:    Sample data vector.
% N (optinal):
%       A vector that specifies the number of bins to be examined. 
%       The optimal number of bins is selected from the elements of N.  
%       Default value is N = 2:500.
%       * Do not search binwidths smaller than a sampling resolution of data.
%
% Output argument
% optN: Optimal number of bins.
% N:    Bin numbers examined.
% C:    Cost function of N.
%
% See also SSKERNEL
%
%
% Copyright (c) 2009 2010, Hideaki Shimazaki All rights reserved.
% http://2000.jukuin.keio.ac.jp/shimazaki

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Setting
x = reshape(x,1,numel(x));
x_min = min(x);
x_max = max(x);

if nargin < 2
    buf = abs(diff(sort(x)));
    dx = min(buf(logical(buf ~= 0)));
    N_MIN = 2;              % Minimum number of bins (integer)
                            % N_MIN must be more than 1 (N_MIN > 1).            
    N_MAX = min(floor((x_max - x_min)/(2*dx)),500);
                            % Maximum number of bins (integer)
    N = N_MIN:N_MAX;        % # of Bins
end
    
SN = 30;                    % # of partitioning positions for shift average
D = (x_max - x_min) ./ N;   % Bin Size Vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the Cost Function
Cs = zeros(length(N),SN);
for i = 1: length(N)

       shift = linspace(0,D(i),SN);
       for p = 1 : SN
               edges = linspace(x_min+shift(p)-D(i)/2,...
                 x_max+shift(p)-D(i)/2,N(i)+1);   % Bin edges

               ki = histcounts(x,edges);               % Count # of events in bins
               ki = ki(1:end-1);

               k = mean(ki);                      % Mean of event count
               v = sum( (ki-k).^2 )/N(i);         % Variance of event count

               Cs(i,p) = ( 2*k - v ) / D(i)^2;    % The Cost Function
       end

end
C = mean(Cs,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Bin Size Selectioin
[~,idx] = min(C);
optN = N(idx);                          % Optimal number of bins
optD = D(idx);                         % *Optimal binwidth
edges = linspace(x_min,x_max,optN+1);  % Optimal segmentation

%[Cminp idxp] = min(Cs(idx,:));
%shift = linspace(0,D(idx),SN);
%edges = linspace(x_min+shift(idxp)-D(idx)/2,...
%                        x_max+shift(idxp)-D(idx)/2,N(idx)+1);
end



% --- HELPER FUNCTION 

function PUC = maskPUC(J,nclasses,divs)
%MASKPUC Mask 3D J array to PUC values according to partition points

jin = J(:);
jout = 0*jin; % assumes classes are 0.

% uses input J to collect classes and mark the positions in
% output J
for v = 1:nclasses
    jout(divs(v) < jin & jin <= divs(v+1)) = v;    
end

PUC = reshape(jout,size(J));

end
