function ptset = savePerfTable(ptset)


if  ~isfield(ptset,'geometry') 
    ptset.geometry = 'K';
    warning('wellps:savePerfTable','Well geometry was set to "K" (vertical).');
elseif isempty(ptset.geometry)      
    % if wells are not vertical, they must be 'horizontal' ('I' or 'J')    
    assert( ~(strcmp(ptset.geometry,'I') || strcmp(ptset.geometry,'J')), 'Option not recognized.') 
end


% perforation list
if isempty(ptset.perfs) 
    error('wellps:savePerfTable','Perforation list is empty.');    
end


% well name
if isempty(ptset.wellname) 
    error('wellps:savePerfTable','Well name is empty.');    
end

% perforations
nperfs = size(ptset.perfs,1); 

% status list 
if  ~isfield(ptset,'status')
    st = cell(1,nperfs);
    st(:) = {'OPEN'};    
    ptset.status = st;
    warning('wellps:savePerfTable','All status were set to "OPEN".');
elseif isempty(ptset.status)
    error('wellps:savePerfTable','Status list is empty.');    
end

% set file name 
filen = ptset.wellname;
filen = strcat('../tmp/perfList',filen,'.txt');

% write file
if exist(filen,'file') == 2
    delete(filen);
end

fid = fopen(filen,'w');



% User Block Address
line = '** UBA  ff  Status  Connection\n';
fprintf(fid,line);

% perf table
fprintf(fid,'  PERF   GEO   %s\n',ptset.wellname);

for n = 1:nperfs
    if n == 1
        fprintf(fid,'  %d %d %d \t %.1f   %s   FLOW-TO   ''SURFACE''   REFLAYER\n', ...
                ptset.perfs(1,1), ... 
                ptset.perfs(1,2), ...
                ptset.perfs(1,3), ...
                1.0, ...
                ptset.status{1} ); 
    else
        fprintf(fid,'  %d %d %d \t %.1f   %s   FLOW-TO   %d\n', ...
                ptset.perfs(n,1), ... 
                ptset.perfs(n,2), ...
                ptset.perfs(n,3), ...
                1.0, ...
                ptset.status{1}, ...
                n-1);                         
    end
    
end

% close
fclose(fid);