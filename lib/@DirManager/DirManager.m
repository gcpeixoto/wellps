classdef DirManager 
 %DIRMANAGER Class responsible for directory management tasks. 
 
 
%% Static methods
 % do not require class objects as arguments
 
 methods (Static)
                         
     function mountDir
         %MOUNTDIR mount standard directory tree
                  
         root = './';
                                                    
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
               
     function tmpDir = createTempDataDir(mfilename)
        %CREATETEMPDATADIR create a temporary directory to store 
        % data output by the case in execution
        
        % Stores into '../tmp' to be ignored by .gitignore.
        % Afterwards, the user must decide if it is garbage or not
        tmpDir = strcat('../tmp/',mfilename,'-data');

        if exist(tmpDir,'dir')
            rmdir(tmpDir,'s');
            mkdir(tmpDir);    
            fprintf('----> Removing directory: "%s"\n',tmpDir);
            fprintf('----> Recreating directory: "%s"\n',tmpDir);
        else
            mkdir(tmpDir);    
            fprintf('----> Creating directory: "%s"\n',tmpDir);
        end
        
     end
   
     
     %% ---- GET methods 
     
     function d = getRootDir
         d = getenv('WELLPS_ROOTDIR');
     end
     
     function d = getCsvDir
         d = getenv('CSV_DIR');
     end
     
     function d = getMatDir
         d = getenv('MAT_DIR');
     end
     
     function d = getLogDir
         d = getenv('LOG_DIR');
     end
     
     function d = getTmpDir
         d = getenv('TMP_DIR');
     end
     
     function d = getBenchMarksDir
         d = getenv('BENCHMARKS_DIR');
     end
     
     function d = getCppDir
         d = getenv('CPP_DIR');
     end
     
     function d = getPyDir
         d = getenv('PY_DIR');
     end
     
     % This was set mainly to run beetweness centrality from Networkx
     function d = getPyExec
         d = 'python3';
     end
               
 end

end