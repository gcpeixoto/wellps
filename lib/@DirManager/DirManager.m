classdef DirManager 
 %DIRMANAGER Class responsible for directory management tasks. 
 
 
%% Static methods
 % do not require class objects as arguments
 
 methods (Static)
                         
     function mountDir
         %MOUNTDIR mount standard directory tree
                  
         root = getenv('WELLPS_ROOTDIR'); % assumes it exists from startup.m
                                                    
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
                 mkdir(dirs{need(i)});                 
             end             
         end
         
     end
     
     %{ 
     %CHECK FOR WHAT PURPOSE?
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
     %}
   
     
     %% ---- GET methods 
     
     function d = getRootDir
         d = getenv('WELLPS_ROOTDIR');
     end
     
     function d = getCsvDir
         d = getenv('WELLPS_CSV_DIR');
     end
     
     function d = getMatDir
         d = getenv('WELLPS_MAT_DIR');
     end
     
     function d = getLogDir
         d = getenv('WELLPS_LOG_DIR');
     end
     
     function d = getTmpDir
         d = getenv('WELLPS_TMP_DIR');
     end
     
     function d = getBenchMarksDir
         d = getenv('WELLPS_BENCHMARKS_DIR');
     end
     
     function d = getCppDir
         d = getenv('WELLPS_CPP_DIR');
     end
     
     function d = getPyDir
         d = getenv('WELLPS_PY_DIR');
     end

     function d = getCasesDir
         d = getenv('WELLPS_CASES_DIR');
     end

     function d = getDocsDir
         d = getenv('WELLPS_DOCS_DIR');
     end

     function d = getThirdPartyDir
         d = getenv('WELLPS_THIRDPARTY_DIR');
     end

     function d = getExamplesDir
         d = getenv('WELLPS_EXAMPLES_DIR');
     end
     
     % This was set mainly to run beetweness centrality from Networkx
     function d = getPyExec
         d = 'python3';
     end
               
 end

end