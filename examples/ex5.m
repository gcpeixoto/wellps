%% EXAMPLE: Computation of metrics for PUCs
%
%
%  REMARK: This example has a extreme computational cost
%          in the current form. Use it for reference
%           
%           
% TODO: make available reduced structs A, B, C and D.
%
%

mrstVerbose on

case_name = 'ex5';

%% Mounting

d = DirManager(); 

%% Reservoir model

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');
  
[G,PROPS] = buildModel(f);

% compute field attributes
P = computeParams(G,PROPS); 

%% Variables

% A few static and dynamic fields

PHI = P.PHI;
KN = P.KN;
SW = P.KN*0 + 1.0; % force unitary

RQIN = P.RQIN;
RQIN_n = ( RQIN - min(RQIN(:)) ) ./ ( max(RQIN(:)) - min(RQIN(:)) );

%% Load Data

% distances
load(fullfile(d.getExamplesDir,strcat(case_name,'sample'),'Rn.mat'),'Rn');
load(fullfile(d.getExamplesDir,strcat(case_name,'sample'),'Tn.mat'),'Tn');

Rn(isnan(Rn)) = 0;
Tn(isnan(Tn)) = 0;

% pressure
load(fullfile(d.getExamplesDir,strcat(case_name,'sample'),'PRESS.mat'),'PRESS');

%% Proxy's choice

% Here, we choose the J functionals.

proxy = {'J1','J2','J3'};

for a = 1:length(proxy)

    aux = proxy{a};

    switch aux
    
        case 'J1'
            J = KN.*PHI.*SW.*PRESS;            
            A.(aux).static = {'perm','poro'};
            A.(aux).dynamic = {'sw'};
            A.(aux).functional = 'KN.*PHI.*SW';
                        
        case 'J2'
            J = RQIN_n.*SW.*PRESS.*log(Rn);
            A.(aux).static = {'perm','poro','distBoundary'};
            A.(aux).dynamic = {'sw'};
            A.(aux).functional = 'KN.*PHI.*SW.*log(Rn);';          
                        
        case 'J3'
            J = RQIN_n.*SW.*PRESS.*log(Tn);
            A.(aux).static = {'perm','poro','distTop'};
            A.(aux).dynamic = {'sw'};
            A.(aux).functional = 'KN.*PHI.*SW.*log(Tn);';   
    
    end
    
    % treatment
    J = J(:);
    J(isnan(J)) = 0;
    J(J < 0) = 0;
    J(isinf(J)) = 0;
        
    % reshape
    J = reshape(J,G.cartDims(1),G.cartDims(2),G.cartDims(3));
    
    % normalization
    J = J./max(J(:));
    
    % append to struct
    A.(aux).data = J; 
    
       
end


%% Binning

% This struct stores information for the binning methods.
% For brevity and example, we used only two methods


% load active
load(fullfile(d.getExamplesDir,strcat(case_name,'sample'),'active.mat'),'active');


binning = {'scott','sturges'}; %without integers method

for a = 1:length(proxy)
    
    auxA = proxy{a};

    % compute IUC
    for b = 1:length(binning)
        
        auxB = binning{b};
        
        [B.(auxA).(auxB).IUC,B.(auxA).(auxB).nclasses,...
        B.(auxA).(auxB).deltas, B.(auxA).(auxB).divs] ...
        = computePUC(A.(auxA).data,active,auxB);

    end

end

%% Connections

% This struct stores information on the connected components

nofsc = 30;

C = struct();

for a = 1:length(fieldnames(B))
   
    fb = fieldnames(B);
    auxA = fb{a};
    
    for b = 1:length(fieldnames(B.(auxA)))
        
        fba = fieldnames(B.(auxA));
        auxB = fba{b};      

        iuc = 1:B.(auxA).(auxB).nclasses;
        IUC = B.(auxA).(auxB).IUC;
        
        % connected components structure
        conn = findConnectionsByPUC(d, iuc, IUC, nofsc, 'n', 1);
        
        for c = 1:length(fieldnames(conn))

            fc = fieldnames(conn);
            auxC = fc{c};

            C.(auxA).(auxB).(auxC) = conn.(auxC);

            % identify components with minimum nofsc
            cvi = conn.(auxC).compVoxelInds;
            ncvi = cellfun(@length,cvi);
            ok = find(ncvi > nofsc);

            C.(auxA).(auxB).(auxC).nofscPassed = ok;

        end
    
    end

end

%% Print good clusters

% This structure gets the supposed best clusters (target of analysis)

D = struct();

fprintf('\n---> List of clusters C_{D,q} attending nofsc > %d\n',nofsc);
fprintf('(J,B,D,q)\n');

for a = 1:length(fieldnames(C))
    fb = fieldnames(C);
    auxA = fb{a};

    for b = 1:length(fieldnames(C.(auxA)))     
        fba = fieldnames(C.(auxA));
        auxB = fba{b};    
        
        for c = 1:length(fieldnames(C.(auxA).(auxB)))
            fbac = fieldnames(C.(auxA).(auxB));
            auxC = fbac{c};   

            np = C.(auxA).(auxB).(auxC).nofscPassed;
            
            
            if ~isempty(np)
                
                D.(auxA).(auxB).(auxC).nofscPassed = np;
                aux2 = auxC(4:end);
                aux3 = str2double(aux2);
                D.(auxA).(auxB).(auxC).initial_value = B.(auxA).(auxB).divs(aux3);
                D.(auxA).(auxB).(auxC).final_value = B.(auxA).(auxB).divs(aux3+1);
                
                for npi = 1:length(np)
                    fprintf('(%s, %s, %s, %d)\n',...
                        auxA,auxB,auxC,npi);                    
                end                
                                
            end
                        

        end

    end
    
end

%% Create output directories, if not existent
% Here, we have a directory tree aggreing with the (J,B,D,q) idea
%
% ex5
%  |- J1
%      |- B1
%  .   .
%  .   .
%  .   .
%      |- Bb
%          |- PUC1
%          .
%  .       .
%  .       .
%  .       |- PUCp
%  
%  |- Jj
%

% this case root dir ('ex5')
ex_dir = fullfile(d.getExamplesDir,case_name);

if ~exist(ex_dir,'dir')

    mkdir(ex_dir) % create root

    for a = 1:length(fieldnames(D))

        fb = fieldnames(D);
        auxA = fb{a};

        mkdir(ex_dir,auxA); % 1st level (functional)

        for b = 1:length(fieldnames(D.(auxA)))
        
            fba = fieldnames(D.(auxA));
            auxB = fba{b}; 

            mkdir(fullfile(ex_dir,auxA),auxB); % 2nd level (binning)

            for c = 1:length(fieldnames(D.(auxA).(auxB)))

                fbac = fieldnames(D.(auxA).(auxB));
                auxC = fbac{c};  

                mkdir(fullfile(ex_dir,auxA,auxB),auxC); % 3rd level (PUC)
        
            end

        end

    end
end


%% Compute metrics

%{
Note that we first run through the "good clusters" structure
(here it is D) to get the clusters of interest for which the metrics
need to be computed, but we get the information about the 
connections from another structure (here this is C).

The switch from D to C is made inside the 3rd-level loop,
where a dynamic reduced version of C is created to reduce
the processing cost. Otherwise, C would have to be 
swept entirely. This way, 'connSt_Reduced' is one-field only 
struct.

%}

% struct of options
opt_metrics.nofsc = nofsc;

for a = 1:length(fieldnames(D))

    fb = fieldnames(D);
    auxA = fb{a};
    
    for b = 1:length(fieldnames(D.(auxA)))
    
        fba = fieldnames(D.(auxA));
        auxB = fba{b}; 

        for c = 1:length(fieldnames(D.(auxA).(auxB)))

            fbac = fieldnames(D.(auxA).(auxB));
            auxC = fbac{c};  
            
         
            % equivalent of drtSt (note that is C)
            connSt = C.(auxA).(auxB); 

            % create volatile reduced struct           
            connSt_Reduced.(auxC) = connSt.(auxC);
            
            % output dir
            opt_metrics.outDir = fullfile(ex_dir,auxA,auxB,auxC);

            % compute            
            [~] = computePUCGraphMetrics(opt_metrics,connSt_Reduced);

            % delete to free
            clear connSt_Reduced

        end

    end

end

%% Reference list for metrics data of clusters of interest 

fprintf('\n---> Directories for metrics of clusters of interest:\n');
for a = 1:length(fieldnames(D))

    fb = fieldnames(D);
    auxA = fb{a};
    
    for b = 1:length(fieldnames(D.(auxA)))
    
        fba = fieldnames(D.(auxA));
        auxB = fba{b}; 

        for c = 1:length(fieldnames(D.(auxA).(auxB)))

            fbac = fieldnames(D.(auxA).(auxB));
            auxC = fbac{c};  
           
            % metrics directories
            fprintf('%s\n',fullfile(ex_dir,auxA,auxB,auxC));

        end
    end
end






