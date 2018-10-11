function clusters = slopeConstrainedClustering(X,Y,i0,seps)
%SLOPECONSTRAINEDCLUSTERING carries out a constrained clustering procedure
%   based on unitary slope of pairwise points (Xi,Yi)
%
% PARAMETERS:
%   - X:    x-coordinate array of points
%   - Y:    y-coordinate array of points
%   - i0:   starting point used for the search
%   - seps: slope relaxation epsilon (e.g. 1e-1, 1e-2)
%
% RETURNS:
%   - disjoint clusters of variable sizes
%
% 
% TODO Extended documentation
%
%
% Dr. Gustavo Oliveira, @LaMEP/UFPB

%% Input 

% (ordered) list of indices (free indices)
I = 1:numel(X);

% number of points
n = numel(I);

% set STARTING POINT as the first of the list if user enters an 
% exceeding index; otherwise, set it to user' choice
if i0 <= 0 || i0 > n
    i = I(1); 
else
    i = i0;
end

%% Starting cluster

% cluster numbering (START WITH 1)
m = 1;

% auxiliary visited list
collected = false(1,n);

% starting cluster
Cm = i;

% mark visited
collected(i) = true;

% removes starting point
I = find(collected == false);

% pairwise checking array
paired = false(1,numel(Cm));

%% Must-link algorithm: 

% repeat until all free indices are collected 
while numel(I) > 0 % not empty list condition
                 
    % sweep free nodes
    for k = I                    
        
        % test pairwise with cluster accepted nodes
        for j = Cm
            
            m_jk = (Y(j) - Y(k))/(X(j) - X(k)); % slope
            
            % safety checking for points with same y-coordinate
            if isinf(m_jk), m_jk = 0; end;
            
            % CONSTRAINT: slope must fall into [1-seps,1+seps]
            if 1 - seps <= m_jk && m_jk <= 1 + seps
                paired(Cm == j) = true;                            
            end            
        end
        
        % transitivity: is point k must-link to all cluster points?
        % if yes, append it to cluster and checks it to collected
        if all(paired)
            Cm = union(Cm, k);        
            collected(k) = true;                               
        end
        
        % reinitialize pairwise checking 
        paired = false(1,numel(Cm));
               
    end
    
    % update list 
    I = find(collected==false);       
    
    % store cluster for output 
    clusters{m} = Cm;
            
    % new cluster 
    m = m + 1;           
    
    % update to new cluster 
    if ~isempty(I)
        i = I(1);
        Cm = i; % first point not collected    
        collected(i) = true; % remove it from the list
        I = find(collected==false); % update list
        paired = false(1,numel(Cm)); % update pairwise checking
    end
     
end

