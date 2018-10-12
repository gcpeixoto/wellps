function [clusters,R2,M] = slopeConstrainedClustering(X,Y,i0,seps)
%SLOPECONSTRAINEDCLUSTERING carries out a constrained clustering procedure
%                           based on unitary slope of pairwise points 
%                           p_i, p_j, such as p_i = (x_i,y_i) and
%                                             p_j = (x_j,y_j), with
%                           x_i = X(i); y_i = Y(i)
%                           x_j = X(j); y_j = Y(j)
%
% PARAMETERS:
%   - X:        x-coordinate array of points (double)
%   - Y:        y-coordinate array of points (double)
%   - i0:       starting point used for the algorithmic search (double)
%   - seps:     slope relaxation epsilon (e.g. 0.1, 0.01, 1e-5) (double)
%
% RETURNS:
%   - clusters:  disjoint clusters of variable sizes (cell)
%   - R2:        coefficients of determination per cluster (double)
%   - M:         slopes of best-fit lines per cluster (double)
%
% 
% DETAILS:
%
%   TODO
%
% Dr. Gustavo Oliveira, @LaMEP/UFPB

%% Input 

% (ordered) list of points (free indices)
I = 1:numel(X);

% number of points
n = numel(I);

% set STARTING POINT as the first of the list if user enters an 
% exceeding index; otherwise, set it to user's choice
if i0 <= 0 || i0 > n
    i = I(1); 
else
    i = i0;
end

%% Starting cluster

% cluster numbering (START WITH label 1)
m = 1;

% auxiliary list of collected points
collected = false(1,n);

% starting cluster begins with the sole starting point
Cm = i;

% mark as visited
collected(i) = true;

% update list by ignoring the collected starting point
I = find(collected == false);

% pairwise checking array
paired = false(1,numel(Cm));

%% Must-link algorithm

% repeat until all free indices are visited and collected 
% if they do not comply with the constraint, will be a 
% 1-element isolated group

while numel(I) > 0 % while there exists free indices
                 
    % sweep free indices
    for k = I                    
        
        % test it pairwise with all already accepted cluster points
        for j = Cm
            
            % slope 
            s_jk = (Y(j) - Y(k))/(X(j) - X(k)); 
            
            % safety checking for points with same x-coordinate.
            % This may happen and indetermine m_jk.
            if isinf(s_jk), s_jk = 0; end;
            
            % CONSTRAINT: slope must fall into [1-seps,1+seps]
            if 1 - seps <= s_jk && s_jk <= 1 + seps
                paired(Cm == j) = true; % marks as paired                            
            end            
        end
        
        % transitivity check: does free point k obeys the must-link         
        % if yes, append it to cluster and checks it to collected
        if all(paired)
            Cm = union(Cm, k); % appends the point to cluster        
            collected(k) = true; % marks as collected                               
        end
        
        % reinitialize pairwise checking array for next point
        paired = false(1,numel(Cm));
               
    end
    
    % update list by ignoring all the collected points
    I = find(collected==false);       
    
    % store cluster m for output 
    clusters{m} = Cm;
            
    % prepares new cluster 
    m = m + 1;           
    
    % update to new cluster 
    if ~isempty(I)  % required because list might be emptied at once 
        i = I(1);   % chooses first free point to continue
        Cm = i;     % starts new cluster
        collected(i) = true; % remove it from the list
        I = find(collected==false); % update list
        paired = false(1,numel(Cm)); % update pairwise checking for new cluster
    end
     
end

%% Cluster reordering 

% Reorder partitioning to output clusters in descending number of elements.
% Remark: The original cluster's label is lost, but this is irrelevant, 
% since, in the end, the points are preserved as of the clustering. 
% Only the cluster labels are changed.
[~,reordered] = sort(cellfun(@numel,clusters),2,'descend');
cReord = clusters;
for c = 1:numel(clusters)
    cReord{c} = clusters{reordered(c)};
end

% reordered cluster
clusters = cReord;

%% Linear regression over the clusters found 

% pre-allocation
R2 = zeros(1,numel(clusters));
M = R2;

% compute R2 and slopes
for c = 1:numel(clusters)  
    IC = clusters{c};   % get indices        
    [R,m,~] = regression(X(IC),Y(IC),'one');    % get coordinates
    R2(c) = R*R;    % coefficient of determination array                                 
    M(c) = m;   % slope array    
end




