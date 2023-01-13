%% Startup
%
% Set environment variables and standard directories for WELLPS to be run. 
%
% REMARKS: 
%
% - This file must be executed before any task.
% - Set the MRST and SNAP variables manually.
% - Test executions by running the codes inside the 'examples' directory.
% - Check if 'cpp/main/graphMetrics' should be recompiled for your platform.
%   The way to go is to run the companion Makefile found inside this CPP
%   directory.

%% IN-CODE VARIABLES

% WELLPS_ROOTDIR is the top-level directory
setenv('WELLPS_ROOTDIR',fileparts(mfilename('fullpath')));
rootdir = getenv('WELLPS_ROOTDIR');

% these are second-level directories
dirs = {'WELLPS_CSV_DIR',...
        'WELLPS_MAT_DIR',...
        'WELLPS_LOG_DIR',...
        'WELLPS_TMP_DIR',...
        'WELLPS_BENCHMARKS_DIR',...
        'WELLPS_CPP_DIR',...
        'WELLPS_PY_DIR',...
        'WELLPS_CASES_DIR',...
        'WELLPS_DOCS_DIR',...
        'WELLPS_THIRDPARTY_DIR',...
        'WELLPS_EXAMPLES_DIR',...
        };

% set second-level directories as env vars to be used through the code
for d = dirs
    ext = split(d{1},'_'); 
    ext = lower( ext{2} );
    setenv(d{1},fullfile(rootdir,ext));
end

% set path for lib to use DirManager
addpath(genpath(fullfile(rootdir,'lib')));

% create standard directories
d = DirManager; 
d.mountDir;

% set path for entire wellps
addpath(genpath(fullfile(rootdir)));


%% THIRD-PARTY VARIABLES

% Check if third-party environment variables are set
% This message is merely informative because the variables
% will not be found if Matlab process is running under GUI-mode call.
envs = {'MRST_DIR','SNAP_DIR'};
fenvs = cellfun(@getenv,envs,'UniformOutput',false);
if any(cellfun(@isempty,fenvs))    
    disp('Required environment variables (see ''README.md''):')
    disp(envs') 
    warning(['---> ATTENTION! Check if the environment variables above are set '...
                   'system-wide. If MATLAB is used on GUI mode, ''getenv'' is' ...
                   'unable to find global environment variables. You may wish' ...
                   'to use ''system(''env'')'' to verify if they appear on ' ...
                   'the list. Otherwise, you must set it manually to run' ...
                   'other portions of the code. Bear in mind that you cand do' ...
                   'that by using bash-like syntax with ''export'' and saving' ...
                   'them into ''~/.bash_profile'' or  ''~/.zshenv'', for instance,' ...
                   'depending on the platform you are using. If you do not ' ...
                   'understand this message, ASK FOR HELP!' ...
                   '''SNAP_DIR'' is required by ''/cpp/Makefile''.'
            ]);
else 
    disp('Required variables are set.')
end

% free up memory
clearvars;