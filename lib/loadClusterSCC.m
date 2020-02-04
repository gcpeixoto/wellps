function [SCCC,SCCT,SCCPY,SCCNL] = loadClusterSCC
%LOADCLUSTERSCC load data structures yielded by call of
%slopeConstrainedClustering, if they are already computed.

% C
if exist('../mat/SCCC.mat','file') == 2
    load('../mat/SCCC.mat');
    SCCC = SCC;
    clear SCC
else
    error('wellps:loadClusterSCC','Did you run "c" clustering?');
end

% CT
if exist('../mat/SCCT.mat','file') == 2
    load('../mat/SCCT.mat');
    SCCT = SCC;
    clear SCC
else
    error('wellps:loadClusterSCC','Did you run "ct" clustering?');
end

% PY
if exist('../mat/SCCPY.mat','file') == 2
    load('../mat/SCCPY.mat');
    SCCPY = SCC;
    clear SCC
else
    error('wellps:loadClusterSCC','Did you run "py" clustering?');
end

% NL
if exist('../mat/SCCNL.mat','file') == 2
    load('../mat/SCCNL.mat');
    SCCNL = SCC;
    clear SCC
else
    error('wellps:loadClusterSCC','Did you run "nl" clustering?');
end

%{
% 6N
if exist('../mat/SCC6N.mat','file') == 2
    load('../mat/SCC6N.mat');
    SCC6N = SCC;
    clear SCC
else
    error('wellps:loadClusterSCC','Did you run "6n" clustering?');
end
%}

end

