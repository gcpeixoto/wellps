function GRID = readCMGDataset(fname)
% READCMGDATASET reads permeability and porosity data 
%                from output .dat file generated from CMG IMEX
%
%
%   input: .dat file
%   output: data structure with property values by grid node
%           (only a few keywords are supported and minimal 
%            information is stored: 
%            - from 'GRID', 'FORMAT' and 'Ni' (i=I,J,K);
%            - 'PERMi' (i=I,J,K): entire data
%            - 'CMGLCustom_Poro', 'CMGLCustom_Perm': entire data. These
%            ones were included to collect the property values generated
%            through geostatistics realization in Builder. However, 
%            these keywords are not standard. They are user-defined names 
%            and need to be changed here. Moreover, the original .dat file
%            must be parsed to remove some prepended text 
%            (e.g.?'RESULTS TEMP_PROP') so that only the values can be
%            stored into output structure.%            
%
% Remark: the output GRID structure can be used as input data in MRST 
% as follows:
%
%       [I,J,K] = deal(GRID.NI,GRID.NJ,GRID.NK); 
%       rock.poro = cell2mat(GRID.POR)';
%       iactive = find(rock.poro ~= 0); 
%       rock.perm = cell2mat([GRID.PERMI;GRID.PERMJ;GRID.PERMK])';
%       
%       Cartesian
%       G = cartGrid([I,J,K]);
%       G = computeGeometry(G);
%       Kx = reshape(rock.perm(:,1),[I,J,K]); 
%       Ky = reshape(rock.perm(:,2),[I,J,K]); 
%       Kz = reshape(rock.perm(:,3),[I,J,K]); 
%
% Remark 2: if the .dat file has one or more *INCLUDE keyword, 
%           the user should provide a new file endowed with
%           the porosity and permeability tables. Otherwise, 
%           they will not be procesed by this function. 
%           This is a feature to be implementated.
%          
% Dr. Gustavo Oliveira

%%

% check file existence
[fid,msg] = fopen(fname,'r');
assert(fid ~= -1,msg);

% check file contents
finfo = dir(fname);
fsize = finfo.bytes;
assert(fsize > 0,'File is empty')
            
flag = true; % flag for fgetl

while ~feof(fid)
    
    % for first call
    if true(flag), tline = fgetl(fid); end    
    
    if strcmp(regexp(tline,'^\*INCLUDE|','match'),'*INCLUDE')   
        warning('*INCLUDE keyword detected. The output ''GRID'' variable might have missing values.');        
        disp(tline)
    end
    
    if strcmp(regexp(tline,'^GRID|','match'),'GRID')    
        sl = regexp(tline,'\s+','split');
        GRID.FORMAT = sl{2};
        GRID.NI = str2double(sl{3});
        GRID.NJ = str2double(sl{4});
        GRID.NK = str2double(sl{5});
        GRID.NCELLS = GRID.NI*GRID.NJ*GRID.NK;     
                        
    end            


    if strcmp(regexp(tline,'^KDIR|','match'),'KDIR')
        sl = regexp(tline,'\s+','split');
        GRID.KDIR = sl{2};                        
    end            
          
    
    % getting PERMEABILITY I             
    if strcmp(regexp(tline,'^PERMI|','match'),'PERMI') 
        
        PERMI = [];

        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = perm
            nps = aux(:,1);
            per = aux(:,2); 
            
            perm = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                perm = [perm;aux];
            end                
            PERMI = [PERMI;perm];                                 

            % update 'while'
            iline = fgetl(fid); 
            sil = regexp(iline,'^[\d\s]');                

        end        
        GRID.PERMI = {PERMI'}; % fill in struct             
        
        tline = iline; 
        flag = false;
        
        clear PERMI
        
    % getting PERMEABILITY J        
    elseif strcmp(regexp(tline,'^PERMJ|','match'),'PERMJ')   
        
        PERMJ = [];

        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
            % sweep cell and get numbers splitting at '*'                    
            % UNDERSTAND how IMEX deals with these format
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = perm
            nps = aux(:,1);
            per = aux(:,2); 
            
            perm = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                perm = [perm;aux];
            end                
            PERMJ = [PERMJ;perm];                                 

            % update 'while'
            iline = fgetl(fid); 
            sil = regexp(iline,'^[\d\s]');                                    

        end
        GRID.PERMJ = {PERMJ'}; % fill in struct 
        
        tline = iline; 
        flag = false;
        
        clear PERMJ
                                    
    % getting PERMEABILITY K        
    elseif strcmp(regexp(tline,'^PERMK|','match'),'PERMK')                
        
        PERMK = [];
        
        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = perm
            nps = aux(:,1);
            per = aux(:,2); 
            
            perm = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                perm = [perm;aux];
            end                
            PERMK = [PERMK;perm];                                 

            % update 'while'
            iline = fgetl(fid); 
            sil = regexp(iline,'^[\d\s]');                                    

        end
        GRID.PERMK = {PERMK'}; % fill in struct
        
        tline = iline; 
        flag = false;
    
        clear PERMK
                
    % getting POROSITY
    elseif strcmp(regexp(tline,'^POR|','match'),'POR')    
    
        POR = [];

        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
                        
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = per
            nps = aux(:,1);
            per = aux(:,2); 
            
            por = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                por = [por;aux];
            end                
            POR = [POR;por];                                 

            % update 'while'
            iline = fgetl(fid); 
            sil = regexp(iline,'^[\d\s]');                                    

        end
        GRID.POR = {POR'}; % fill in struct             
        
        tline = iline; 
        flag = false;
        
        clear POR
        
    % getting 'CMGLCustom_Poro'
    elseif strcmp(regexp(tline,'|CMGLCustom_Poro|','match'),'CMGLCustom_Poro')    
    
        GEOSTAT_POR = [];

        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
                        
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = per
            nps = aux(:,1);
            per = aux(:,2); 
            
            por = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                por = [por;aux];
            end                
            GEOSTAT_POR = [GEOSTAT_POR;por];                                 

            % update 'while'
            iline = fgetl(fid); 
            sil = regexp(iline,'^[\d\s]');                                    

        end
        GRID.GEOSTAT_POR = {GEOSTAT_POR'}; % fill in struct             
        
        tline = iline; 
        flag = false;
        
        clear GEOSTAT_POR
    
    % getting 'CMGLCustom_Perm'
    elseif strcmp(regexp(tline,'|CMGLCustom_Perm|','match'),'CMGLCustom_Perm')    
    
        GEOSTAT_PERM = [];

        % row with perm. data 
        iline = fgetl(fid);

        % allows starts with a digit or space            
        sil = regexp(iline,'^[\d\s]');  

        while ~isempty(sil)           

            data = regexp(iline,'\s+','split'); % splits the row                

            aux = [];                
                        
            for i = 1:numel(data) 
                
                % assumes there is a '*'
                a = regexp(data{i},'\*','split'); 
                
                if all(size(a) == [1,1]) && isempty(a{1}) % if a null character 
                    az = [0,0];
                    aux = [aux;az];
                elseif all(size(a) == [1,1]) && ~isempty(a{1})% if only one value (no '*')
                    az = [1,str2double(a{1})]; % multiply by itself
                    aux = [aux;az];                
                else
                    az = str2double(a);
                    aux = [aux;az];
                end                                        
            end

            % repeat data: nps blocks with value = per
            nps = aux(:,1);
            per = aux(:,2); 
            
            por = [];
            for nv = 1:numel(nps)
                aux = repmat(per(nv),[nps(nv),1]);
                por = [por;aux];
            end                
            GEOSTAT_PERM = [GEOSTAT_PERM;por];                                 

            % update 'while'
            iline = fgetl(fid); 
            if ~isa(iline,'char')
                disp('here')
            end
                
            sil = regexp(iline,'^[\d\s]');                                    

        end
        GRID.GEOSTAT_PERM = {GEOSTAT_PERM'}; % fill in struct             
        
        tline = iline; 
        flag = false;
        
        clear GEOSTAT_PERM
        
    else % go to next line
        
        flag = false;
        tline = fgetl(fid);
    
    end 

%}
end % close feof

% close file
fclose(fid); 

end
