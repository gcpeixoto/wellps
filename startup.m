%% Startup
% Set environment variables and standard directories for internal use 
% of WELLPS. 
% This file must be invoked before beginning to work with WELLPS. 
% 

% WELLPS_ROOTDIR
setenv('WELLPS_ROOTDIR',fileparts(mfilename('fullpath')));
rootdir = getenv('WELLPS_ROOTDIR');

dirs = {'CSV_DIR',...
        'MAT_DIR',...
        'LOG_DIR',...
        'TMP_DIR',...
        'BENCHMARKS_DIR',...
        'CPP_DIR',...
        'PY_DIR'};

for d = dirs
    ext = split(d{1},'_'); 
    ext = lower( ext(1) );
    setenv(d{1},fullfile(rootdir,char(ext)));
end

% create standard directories
d = DirManager; 
d.mountDir;

% clean
clearvars;