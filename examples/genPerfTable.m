%% Generate automatic perforation table 

% generating 40 perforations in a column
n = 40;

ptset.wellname = 'Well1';
ptset.geometry = 'K';
ptset.perfs = [ones(n,1),2*ones(n,1),(1:n)'];               
ptset = savePerfTable(ptset);
