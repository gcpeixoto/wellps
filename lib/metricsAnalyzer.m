function metricsAnalyzer(drtValue,idcomp,ave,base)
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

% create output dir
outDir = '../csv/MetricsAnalytics/';
if exist(outDir,'dir') ~= 7
    warning('wellps:metricsAnalyzer',...
        ['Creating data analytics directory into:',...
         '"../csv/MetricsAnalytics/".']);
    mkdir(outDir); 
end   


% loop over DRT value
for drt = drtValue
    
    % try load .mat files required by user
    fidd = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'.mat');
    fidm = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'_Metrics.mat');
    fidr = strcat('../mat/DRT_',ave,'_',base,'_',num2str(drt),'_LinRegr.mat');
   
    if ~exist(fidd,'file') || ~exist(fidm,'file') || ~exist(fidr,'file')
        error('wellps:metricsAnalyzer',... 
            ['All the required files: \n %s \n %s \n were not found', ...
             'for DRT = %d. Have you run wellps:computeDRTGraphMetrics'], ...
             fidm,fidr,drt);
    else
        
        % load DRT struct
        load(fidd,'drtSt');                
        
        % load metrics
        load(fidm,'Maux');        
        
        % load linregr
        load(fidr,'Laux');        
        
        % rename
        M = Maux; clear Maux;
        L = Laux; clear Laux;      
        
    end
    
    fprintf('Sweeping DRT = %d ...\n',drt);
    
    % loop over components
    for c = idcomp
        
        fprintf('-----> Sweeping cluster = %d ...\n',c);
                
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
        aux = [cvc,deg,bet,clo];

        % writing to file 
        hdr = {'I,';'J,';'K,';'betweeness,';'closeness,';'degree'}; 
        hdr = hdr';                        
        fname = strcat('DRT_',ave,'_',base,'_',num2str(drt),'_Cluster_',num2str(c),'_MetricsAnalytics');                
        dlmwrite(strcat(outDir,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(outDir,fname,'.csv'),aux,'-append','delimiter',',','precision','%g');   
        
        fprintf('-----> File: "%s" saved. \n',fname);
        
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
        hdr = {'maxClo,';'MIglob,';'MJglob,';'MKglob,';         ...
                          'MIloc,';'MJloc,';'MKloc,';           ...
               'minClo,';'mIglob,';'mJglob,';'mKglob,';         ...
                          'mIloc,';'mJloc,';'mKloc,';           ...
               'Imin,';'Imax,';'Jmin,';'Jmax,';'Kmin,';'Kmax,'; ...
               's,';'R2,';'performance'}; 
        hdr = hdr';                                               
        aux = [maxclo(1)         ,...
               cvcmaxclo(1,:)    ,...
               cvcmaxcloloc(1,:) ,...
               minclo(1)         ,...
               cvcminclo(1,:)    ,...
               cvcmincloloc(1,:) ,...
               ilims             ,...
               jlims             ,...
               klims             ,...
               L.slope{c}        ,...
               L.R2{c}      ,...
               L.performance{c}];

        fname = strcat('DRT_',ave,'_',base,'_',num2str(drt),'_Cluster_',num2str(c),'_MetricsAnalyticsMinMax');                
        dlmwrite(strcat(outDir,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(outDir,fname,'.csv'),aux,'-append','delimiter',',','precision','%g');
        
        fprintf('-----> File: "%s" saved. \n',fname);
    end
    
    % number of elements per cluster
    nce = cell(size(M.idComp));
    for ncl = 1:numel(M.idComp)
        nce{ncl} = numel(M.degreeCentrality{ncl});
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
    fname = strcat('DRT_',ave,'_',base,'_',num2str(drt),'_PerformanceTable');                          
    dlmwrite(strcat(outDir,fname,'.csv'),hdr,'delimiter','');
    dlmwrite(strcat(outDir,fname,'.csv'),aux,'-append','delimiter',',','precision','%g'); 
    
    fprintf('File: "%s" saved. \n',fname);
            
    % destroy to reallocate
    clear D M L
end
    



