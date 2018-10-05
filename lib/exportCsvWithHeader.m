function [hdr,fn] = exportCsvWithHeader(fn,data,varargin)
%CSVHEADER  Creates a header to a nonexistent .csv file by overwriting 
%           any existent file of same name.

% PARAMETERS:
%       - fn:   physical file name (string) to which the header will be
%               written.
%
%       - data: data table to be appended into file.
%
%       - varargin: names to be used in the header row (strings). 
%                       (e.g. 'alfa','gamma','delta', and so on.)


% check
assert(isnumeric(data), ... 
    '''data'' must be an array.');

assert(all(cellfun(@ischar,varargin)), ... 
    'Not all header column names are strings.');


% add .csv extension, if none. 
[~,~,e] = fileparts(fn);
if isempty(e),  fn = strcat(fn,'.csv'); end

% write header to file
hdr = varargin;
for i = 1:length(varargin)
    aux = strcat(varargin{i},',');
    hdr{i} = aux;    
end
hdr = hdr';
dlmwrite(fn,hdr,'');

% write data
dlmwrite(fn,data,'delimiter',',');

% check if file was saved
if exist(fn,'file')
    fprintf('---> Data successfully written into ''%s''.\n',fn);    
end


