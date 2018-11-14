function sccmet = setSCCMethod(method)
%SETSCCMETHOD Simple set for Slope Constrained Clustering methods

available = {'c','ct','py','nl','6n'};
mets = {'SCCC (Planar Cone Only)',           ...
        'SCCT (Planar Cone + Transitivity)', ...
        'SCCPY (Parallel to Y axis)',        ...
        'SCCNL (Normal to 1-Slope Line)',    ...
        'SCC6N (6Neighbours + Least Square Error)'};

if ~any(ismember(available,method))
    error('wellps:setSCCMethod','Method not available.')
    disp([available',mets']);
else
    sccmet = method;
end
    
end
