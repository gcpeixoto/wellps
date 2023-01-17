
%% MRST modules
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui
mrstModule list

%% Utility
d = DirManager(); 

%% Grid reading
grdecl = readGRDECL('wellps/benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
G = processGRDECL(grdecl, 'checkgrid' , true);
G = computeGeometry(G);

Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

%% Rock/fluid properties

rock.poro = grdecl.PORO;
rock.perm(:,1) = grdecl.PERMX;
rock.perm(:,2) = grdecl.PERMY;
rock.perm(:,3) = grdecl.PERMZ;

G.PORO = rock.poro;

PROPS.PHI = rock.poro;
PROPS.KX = rock.perm(:,1);
PROPS.KY = rock.perm(:,2);
PROPS.KZ = rock.perm(:,3);

p = linspace(100*barsa, 220*barsa, prod(G.cartDims));
PROPS.PRESS0  = p; % init pressure

PROPS.SW = ones(prod(G.cartDims),1);
active = find(grdecl.ACTNUM == 1);

% compute field attributes
P = computeParams(G,PROPS); 

%% proxy's underlying parameters

KN = P.KN;
PHI = reshape(PROPS.PHI,G.cartDims);
PRESS = reshape(PROPS.PRESS0,G.cartDims);
SW= reshape(PROPS.SW,G.cartDims);

%% proxy's choice

proxy = {'J1','J2'};

for a = 1:length(proxy)

    aux = proxy{a};

    switch aux
    
        case 'J1'
            J = KN.*PHI.*SW;    
            A.(aux).static = {'perm','poro'};
            A.(aux).dynamic = {'sw'};
            A.(aux).functional = 'KN.*PHI.*SW';    
        
        case 'J2'
            J = KN.*PHI.*PRESS.*SW;
            A.(aux).static = {'perm','poro'};
            A.(aux).dynamic = {'press, sw'};
            A.(aux).functional = 'KN.*PHI.*PRESS.*SW';  
    
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

% test


%% binning methods

binning = {'scott','fd','sturges','sqrt','shimazaki'};

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

%% connections

%{ 
Firstly, observe that we have an abstract pair 
(p_i,b_j), where p_i is a proxy method and b_j is a binning method.
Each pair can be accessed by using the structs A and B
we have just created above.

For example, (p_1,b_2) is the pair (J1,'fd').

%}

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

%% print good clusters

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
                for npi = 1:length(np)
                    fprintf('(%s, %s, %s, %d)\n',...
                        auxA,auxB,auxC,npi);
                end
            end

        end

    end

end


