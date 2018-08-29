classdef DirManager 
 %DIRMANAGER Class responsible for directory management tasks. 
 
 
%% Static methods
 % do not require class objects as arguments
 
 methods (Static)
     
     function mountDir
         %MOUNTDIR mount standard directory tree
         
         % gets working dir 
         aux = split(pwd,filesep);
         
         % check relative path of working dir to 'wellps' top 
         where = find(strcmp(aux,'wellps'));         
         level = numel(aux) - where;         
         if level == 0
             root = '.'; % this is root
         else
             root = repmat('../',1,level); % this is not
         end
                                           
         % standard dirs
         %      - csv: csv files 
         %      - tmp: temporary files
         %      - mat: MAT-file storage 
         %      - log: log files
         dirs = cellfun(@fullfile,{ root  , root  , root  , root }, ... 
                                  { 'csv' , 'tmp' , 'mat' , 'log'}, ... 
                                  'UniformOutput',false);
         
         % if any of the standard dirs doesn't exist, create it
         if ~all(cellfun(@(d) exist(d,'dir'),dirs))             
             need = find(~cellfun(@(d) exist(d,'dir'),dirs));
             for i = 1:numel(need)
                 mkdir(root,dirs{need(i)});                 
             end             
         end
         
     end
     
 end
 

end