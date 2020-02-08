function dataDir = metricsAnalysis(opt_analytics)
%METRICSANALYZER compute Data Analytics of cluster metrics and export to
% file.
%
% PARAMETERS:
%       - drtValue:  array of DRT values to analyze.
%       - idcomp:   cluster components to analyze for this DRT.
%       - ave:      averaging technique option to read correspondent file.
%       - base:     log base option to read correspondent file.
%
% RETURNS:
%       - export .csv files to disk.


% checking
p = {'loaddir', ...
     'drtlist', ...
     'complist', ...
     'krule', ...
     'logbase', ...
     'savedir', ...
     'fileperf', ...
     'fileminmax', ...
     'filemetrics'};
 
if ~all(ismember(p,fieldnames(opt_analytics)))   
    error(['wellps:computeDRTGraphMetrics ', ... 
           '''opt_analytics'' must be a ''struct'' object with ', ... 
           'specified settings that are missing.', ...
           'See function help.']);
end

dataDir = opt_analytics.savedir;
D = opt_analytics.drtlist;
idcomp = opt_analytics.complist;
ld = opt_analytics.loaddir;
r = opt_analytics.krule;
b = opt_analytics.logbase;
flag_perf = opt_analytics.fileperf;
flag_mm = opt_analytics.fileminmax;
flag_metrics = opt_analytics.filemetrics;

% create savedir
if exist(dataDir,'dir') ~= 7, mkdir(dataDir); end   


% loop over DRT value
tstart = tic; % timing
for drt = D
    
    % try to load .mat files required for analysis     
    f1 = fullfile(ld,char(strcat(join({'DRT',r,b,num2str(drt)},'_'),'.mat')));
    f2 = fullfile(ld,char(strcat(join({'DRT',r,b,num2str(drt),'Metrics'},'_'),'.mat')));
    f3 = fullfile(ld,char(strcat(join({'DRT',r,b,num2str(drt),'LinRegr'},'_'),'.mat')));
           
    if ~exist(f1,'file') || ~exist(f2,'file') || ~exist(f3,'file')
        error('wellps:metricsAnalysis',... 
            ['All the required files: \n %s \n %s \n were not found', ...
             'for DRT = %d. Have you run wellps:computeDRTGraphMetrics'], ...
             f2,f3,drt);
    else
        
        % load DRT struct, metrics, and linregr
        load(f1,'drtSt');                
        load(f2,'M');                        
        load(f3,'L');                        
        
    end
    
    fprintf('=> DRT = %d ...\n',drt);
    
    % loop over components
    for c = idcomp
        
        fprintf('==> Analyzing cluster = %d ...\n',c);
                
        %-----------------------------------------
        % CLUSTER GROSS INFORMATION 
        %-----------------------------------------
        
        % global logical indices
        cvc = drtSt.compVoxelCoords{ M.idComp{c} };    
        
        % centralities        
        bet = M.betweenessCentrality{c};
        clo = M.closenessCentrality{c};
        deg = M.degreeCentrality{c};

        % matrix to save
        aux = [cvc,bet,clo,deg];

        if flag_metrics(1)
            
            % writing to file 
            hdr = {'I,';'J,';'K,';'betweeness,';'closeness,';'degree'}; 
            hdr = hdr';    
            
            fn = fullfile(dataDir,char(strcat(join({'DRT',            ...
                                                r,b,             ...
                                                num2str(drt),    ... 
                                                'Cluster',       ... 
                                                num2str(c),      ...
                                                'Metrics'},'_'), ... 
                                                '.csv')));
                        
            dlmwrite(fn,hdr,'delimiter','');
            dlmwrite(fn,aux,'-append','delimiter',',','precision','%g');   
            fprintf('==> File: ''%s'' saved.\n',fn);
        end
        
        %-----------------------------------------
        % CLUSTER POINTWISE INFORMATION 
        %-----------------------------------------                 
        
        % min bet
        minbet = min(bet);                    
        iminbet = bet == minbet;     
        cvcminbet = cvc(iminbet,:);  
        
        % max bet
        maxbet = max(bet);                    
        imaxbet = bet == maxbet;     
        cvcmaxbet = cvc(imaxbet,:);        
                                 
        % min clo
        minclo = min(clo);                    
        iminclo = clo == minclo;     
        cvcminclo = cvc(iminclo,:);        
        
        % max clo
        maxclo = max(clo);                    
        imaxclo = clo == maxclo;     
        cvcmaxclo = cvc(imaxclo,:);                 
        
        % min deg
        mindeg = min(deg);                    
        imindeg = deg == mindeg;     
        cvcmindeg = cvc(imindeg,:);        

        % max deg
        maxdeg = max(deg);                    
        imaxdeg = deg == maxdeg;     
        cvcmaxdeg = cvc(imaxdeg,:);        
              
        % COORDINATES TRANSFORM (CMG)
        %
        % REMARK: 
        % Note that we might have nonunique max/min points. 
        % When this is the case, we get the first in the list,
        % because we do not lose almost anything in terms of proximity
        % of the 'optimal' position. 
        %
        % Further, be aware that this transform is carried over the logical 
        % I,J,K indices (see wellps:global2LocalCluster).
        
        % cluster bounding box limits
        im = min(cvc(:,1)); iM = max(cvc(:,1)); ilims = [im,iM];
        jm = min(cvc(:,2)); jM = max(cvc(:,2)); jlims = [jm,jM];
        km = min(cvc(:,3)); kM = max(cvc(:,3)); klims = [km,kM];
        
        
        % ---- GET BETWEENNESS
        % transform: min point         
        cvcminbetloc = global2LocalCluster(cvcminbet(1,1),...
                                           cvcminbet(1,2),...
                                           cvcminbet(1,3),...
                                           im,jm,km);  
        % transform: max point                               
        cvcmaxbetloc = global2LocalCluster(cvcmaxbet(1,1),...
                                           cvcmaxbet(1,2),...
                                           cvcmaxbet(1,3),...
                                           im,jm,km);       

                                       
        % ---- GET DEGREE
        % transform: min point         
        cvcmindegloc = global2LocalCluster(cvcmindeg(1,1),...
                                           cvcmindeg(1,2),...
                                           cvcmindeg(1,3),...
                                           im,jm,km);  
        % transform: max point                               
        cvcmaxdegloc = global2LocalCluster(cvcmaxdeg(1,1),...
                                           cvcmaxdeg(1,2),...
                                           cvcmaxdeg(1,3),...
                                           im,jm,km);       
                                       
                                               
        % ---- GET CLOSENESS
        % transform: min point         
        cvcmincloloc = global2LocalCluster(cvcminclo(1,1),...
                                           cvcminclo(1,2),...
                                           cvcminclo(1,3),...
                                           im,jm,km);  
        % transform: max point                               
        cvcmaxcloloc = global2LocalCluster(cvcmaxclo(1,1),...
                                           cvcmaxclo(1,2),...
                                           cvcmaxclo(1,3),...
                                           im,jm,km);       

                                       
        % writing to file    
        if flag_mm(1)

            hdr = {'maxBet,';'MIglob_bet,';'MJglob_bet,';'MKglob_bet,'; ...
                              'MIloc_bet,';'MJloc_bet,';'MKloc_bet,';   ...
                   'minBet,';'mIglob_bet,';'mJglob_bet,';'mKglob_bet,'; ...
                              'mIloc_bet,';'mJloc_bet,';'mKloc_bet,';   ...
                   %%%
                   'maxClo,';'MIglob_clo,';'MJglob_clo,';'MKglob_clo,'; ...
                              'MIloc_clo,';'MJloc_clo,';'MKloc_clo,';   ...
                   'minClo,';'mIglob_clo,';'mJglob_clo,';'mKglob_clo,'; ...
                              'mIloc_clo,';'mJloc_clo,';'mKloc_clo,';   ...
                   %%%
                   'maxDeg,';'MIglob_deg,';'MJglob_deg,';'MKglob_deg,'; ...
                              'MIloc_deg,';'MJloc_deg,';'MKloc_deg,';   ...
                   'minDeg,';'mIglob_deg,';'mJglob_deg,';'mKglob_deg,'; ...
                              'mIloc_deg,';'mJloc_deg,';'mKloc_deg,';   ...
                   %%%
                   'Imin,';'Imax,';'Jmin,';'Jmax,';'Kmin,';'Kmax,';     ...
                   's,';'R2,';'performance'};            
            hdr = hdr';                                               
            aux = [maxbet(1)         ,...  % there might exist more than 1 max closeness point. We get the 1st of the list
                   cvcmaxbet(1,:)    ,...  % betweeness
                   cvcmaxbetloc(1,:) ,...  |
                   minbet(1)         ,...  |
                   cvcminbet(1,:)    ,...  |
                   cvcminbetloc(1,:) ,...  |
                   maxclo(1)         ,...  % closeness
                   cvcmaxclo(1,:)    ,...  |
                   cvcmaxcloloc(1,:) ,...  |
                   minclo(1)         ,...  |
                   cvcminclo(1,:)    ,...  |
                   cvcmincloloc(1,:) ,...  |
                   maxdeg(1)         ,...  % degree 
                   cvcmaxdeg(1,:)    ,...  |
                   cvcmaxdegloc(1,:) ,...  |
                   mindeg(1)         ,...  |
                   cvcmindeg(1,:)    ,...  |
                   cvcmindegloc(1,:) ,...  |
                   ilims             ,...  % bounds
                   jlims             ,...  |
                   klims             ,...  |
                   L.slope{c}        ,...  % curve fitting
                   L.R2{c}           ,...  |
                   L.performance{c}];      % 

            fn = fullfile(dataDir,char(strcat(join({'DRT',            ...
                                                r,b,             ...
                                                num2str(drt),    ... 
                                                'Cluster',       ... 
                                                num2str(c),      ...
                                                'MinMax'},'_'),  ... 
                                                '.csv')));

            dlmwrite(fn,hdr,'delimiter','');
            dlmwrite(fn,aux,'-append','delimiter',',','precision','%g');   
            fprintf('==> File: ''%s'' saved.\n',fn);
        end        
                
    end % close idcomp 
    
    %-----------------------------------------
    % ALL CLUSTERS STATISTICS
    %-----------------------------------------
    
    % number of elements per cluster
    nce = cell(size(M.idComp));
    for ncl = 1:numel(M.idComp)
        nce{ncl} = numel(M.degreeCentrality{ncl}); % any variable for numel
    end
        
    % performance table per DRT
    hdr = {'cluster,';'nce,';'s,';'R2,';'performance'}; 
    hdr = hdr';    
    aux = [cell2mat(L.idComp)'  ,...
           cell2mat(nce)'       ,...
           cell2mat(L.slope)'   ,...
           cell2mat(L.R2)'      ,...
           cell2mat(L.performance)'];
        
    % writing to file 
    if flag_perf(1)

        fn = fullfile(dataDir,char(strcat(join({'DRT',                     ...
                                            r,b,                      ...
                                            num2str(drt),             ... 
                                            'PerformanceTable'},'_'), ... 
                                            '.csv')));

        dlmwrite(fn,hdr,'delimiter','');
        dlmwrite(fn,aux,'-append','delimiter',',','precision','%g');   
        fprintf('==> File: ''%s'' saved.\n',fn);
    end
    
    % destroy to reallocate
    clear drtSt M L

end
   
tfinal = toc(tstart);
fprintf('metricsAnalysis finished after %g seconds. \n',tfinal);




