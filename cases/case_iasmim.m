%% Localizacao de HFUs - reservatorio Iasmim/SIPGEM

model = 'w5'; % 'w4'
log_base = 'ln'; % 'log10' or 'ln'
nofs = 10; % number of significant cells (DRT connections) (only to save info)

%%
[G,PROPS] = buildModel(strcat('../benchmarks/psy-iasmim/imex/',model,'.dat'));
% forced fixing 
PROPS.KY = PROPS.KX;
PROPS.KZ = 0.1*PROPS.KX;
G = computeGeometry(G);

% plot
%rock.poro = PROPS.PHI(:);
%rock.perm = [PROPS.KX(:),PROPS.KY(:),PROPS.KZ(:)];
%plotCellData(G,rock.poro)
%plotCellData(G,rock.perm(:,1))

%% Parameters
P = computeParams(G,PROPS);

%% DRT conversion

% Required to change the adjust constant here, because DRT < 0 were found
% for all averages with the standard formula.
% CONSIDERING ONLY BASE LN FOR LOG

ave = 'n';
switch ave
    case 'a'
        FZI = P.FZIA;
        b = 8.0; % adjust constant
        r = @(FZI) round(log(FZI) + b);
        DRT = r(FZI);
        DRT(isinf(DRT)) = 0;
        P.DRTA_LN = DRT;
        
    case 'g'
        FZI = P.FZIG;
        b = 8.5;
        r = @(FZI) round(log(FZI) + b);
        DRT = r(FZI);
        DRT(isinf(DRT)) = 0;
        P.DRTG_LN = DRT;
    
    case 'n'
        FZI = P.FZIN;
        b = 7.9;
        r = @(FZI) round(log(FZI) + b);
        DRT = r(FZI);
        DRT(isinf(DRT)) = 0;
        P.DRTN_LN = DRT;
        
    case 'q'
        FZI = P.FZIQ;
        b = 8.5;
        r = @(FZI) round(log(FZI) + b);
        DRT = r(FZI);
        DRT(isinf(DRT)) = 0;
        P.DRTQ_LN = DRT;
        
    case 'h'
        FZI = P.FZIH;
        b = 9.0;
        r = @(FZI) round(log(FZI) + b);
        DRT = r(FZI);
        DRT(isinf(DRT)) = 0;
        P.DRTH_LN = DRT;
end

drtlist = unique(DRT);
drtlist = drtlist(2:end);
%% compute HFUs by DRT
drtSt = findDRTConnections(drtlist, P, 'normalized','ln',nofs,'n', 1);
