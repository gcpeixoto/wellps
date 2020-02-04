function S = getStats(dobj,P,pnames,save)
%PRINTSTATS Get statistics of field property and, optionally, print to file
%   
% PARAMETERS:
%       dobj    - DirManager class object.
%       P       - structure whose fields are the property values of the reservoir.  
%
%       pnames  - property names for which to print data.
%
%       save    - flag to choose between saving or not .csv files 
%                 into the standard csv/ directory.
%
% RETURNS:
%   
%       S       - cell containing statistics for each property. Each
%                 element of S is a (mx3) matrix whose columns are  
%                 'value', 'number of occurrences' and 'percentage'.
%
%{
    Developed at LaMEP/UFPB, Brazil
    @gcpeixoto
%}

% checking
if ~isstruct(P)
    error('wellps:printStats','P must be a struct object.');
end

if ~iscell(pnames)        
    error('wellps:printStats','pnames must be a cell of strings containing valid property names.');    
end

if ~ischar(save)
    error('wellps:printStats','save must be: "y" [yes] or "n" [no].');    
end

% save to file or not
switch save 
    case 'y'
        to_file = true;
    case 'n'
        to_file = false;
end

% header
hdr = {'value,';'frequency,';'percentage'}';
hdr = sprintf('%s\t',hdr{:}); 
hdr(end)='';

np = numel(pnames); % number of chosen properties 
S = cell(1,numel(pnames)); % stats table

% printing data table to .csv file
for n = 1:np
    p = P.(pnames{n});    
    T = tabulate(p(:));    
    S{n} = T;   
    
    if to_file == true       
        aux = fullfile(dobj.getCsvDir,strcat('stats-',pnames{n},'.csv'));
        dlmwrite(aux,hdr,'');            
        dlmwrite(aux,T,'delimiter',',','-append'); 
        fprintf('----> File %s was saved.\n',aux);
    end
end

end

