function ptset = savePerfTable(ptset)
%SAVEPERFTABLE writes a .txt file containing a perforation table
%              to be pasted (MANUALLY) into .DAT CMG IMEX file.
%
% EXAMPLE:
%
%   % mounting table of 40 perforations in the column [1,2,k], k = 1:40
%   n = 40;
%   ptset.wellname = 'Well1';
%   ptset.geometry = 'K';
%   ptset.perfs = [ones(n,1),2*ones(n,1),(1:n)'];               
%   ptset = savePerfTable(ptset);
%
%
% PARAMETERS: 
%       - ptset:    settings to be used to mount the table (struct)
%
%                   'ptset' accepts the following fields:
%
%                    - 'geometry': for now, only the flag 'K' (char)
%                                  is accepted. This stands for the 
%                                  well direction and it is a parameter
%                                  understood by BUILDER/IMEX to set up a 
%                                  vertical well. See CMG manual to see 
%                                  uses and options.
%
%                    - 'perfs': logical indices of the cells to be assigned
%                               in perforation table. It is a (nx3) array,
%                               where n is the number of cells to be 
%                               perforated. For example: 
%                       
%                               ptset.perfs = [1 2 3; 1 2 4; 1 2 5];
%
%                    - 'wellname': well name (char). This is only a label.
%
%                    - 'status': the status for each perforation (cell). 
%                                It is a cell of strings. If not specified 
%                                by the user, the standard option 'OPEN' is
%                                set to all the perforations. In case of
%                                variable options, the status specification 
%                                is mandatory for each individual
%                                perforation. For example:
%
%                                ptset.status = {'OPEN','OPEN','SHUT IN'};
%
%
% RETURNS:
%
%       - file:     a .txt which is saved to standart ../tmp folder
%
%
% Dr. Gustavo Oliveira, @LaMEP/UFPB

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

fprintf('---> Perforation table saved to file: ''%s''.\n',filen);

% close
fclose(fid);