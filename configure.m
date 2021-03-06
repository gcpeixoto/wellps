
%TODO: optimize env. variable setup for GUI mode. 
%{   
   This file is valid ONLY for execution in non-interactive mode inside
   the top of 'wellps' directory.           
   Matlab shell in interactive mode doesn't read ~/.bashrc file or alike.
   Hence, using ! or system(...) to find environmental variables 
   is not so simple. 
   
   See: https://www.mathworks.com/matlabcentral/answers/94199-can-i-use-aliases-when-using-the-and-system-commands-from-within-matlab

   To configure Matlab in non-interactive mode, 
   see: https://www.mathworks.com/help/matlab/ref/matlablinux.html  
 %}

% message
msg = ['This file is valid ONLY for execution '  ...
       'in non-interactive mode. Otherwise, use the GUI ' ...
       'command "Set Path".'];           
warning(msg);


% check if environment variables are set
envs = {'MRST_DIR','SNAP_DIR'};
fenvs = cellfun(@getenv,envs,'UniformOutput',false);
if any(cellfun(@isempty,fenvs))    
    disp('Required variables to be set (see README):')
    disp(envs') 
    error('Attention: environment variables not set.')
else
    disp('Required environment variables are set:')        
    disp([envs',fenvs']);
end

% add depend to matlab path 
cellfun(@addpath,fenvs);

% add this
addpath( genpath(fullfile(pwd,'lib')),         ...
         genpath(fullfile(pwd,'third-party')), ...
         genpath(fullfile(pwd,'benchmarks')),  ...
         genpath(fullfile(pwd,'py')),          ...
         genpath(fullfile(pwd,'examples')),    ...
         genpath(fullfile(pwd,'cases'))  );

% this performs MRST startup 
run(fullfile(fenvs{1},'startup.m'));
