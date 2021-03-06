%% CASE STUDY: DRT distribution under different averaging

% Visualization of DRT distribution over UNISIM 1D from 
% LN and LOG10 for different permeability averaging rules

% REMARK: take the pair (KH,LN) for permeability and log base 
%         appears to be the most promising choice among the tested
%         due to its superior resolution in zoning.

warning('THIS CASE IS EXTREMELY COSTLY!!!');


d = DirManager;

%% Default
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

hnb = 6; % number of histogram bins 
FS = 16; % fontsize of labels


log_base = 'ln'; % 'log10' or 'ln'

nofs = 100; % number of significant cells (DRT connections) (only to save info)

model = 'unisim1';

%% Load grid 
switch model 
    case 'unisim1'
        [G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
        
        % compute geometry
        G = computeGeometry(G);
        
    % \TODO conversion of permeability data from m2 to mD
    case 'spe10'
        [G,PROPS] = buildModelSPE10('original');
end


%% Parameters
% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

% field statistics: triplets (DRT,perm,log_base)
SA = getStats(d,P,{'DRTA_LN','DRTA_LOG10'},'n');
SG = getStats(d,P,{'DRTG_LN','DRTG_LOG10'},'n');
SN = getStats(d,P,{'DRTN_LN','DRTN_LOG10'},'n');
SQ = getStats(d,P,{'DRTQ_LN','DRTQ_LOG10'},'n');
SH = getStats(d,P,{'DRTH_LN','DRTH_LOG10'},'n');

SAStar = getStats(d,P,{'DRTAStar_LN','DRTAStar_LOG10'},'n');
SGStar = getStats(d,P,{'DRTGStar_LN','DRTGStar_LOG10'},'n');
SNStar = getStats(d,P,{'DRTNStar_LN','DRTNStar_LOG10'},'n');
SQStar = getStats(d,P,{'DRTQStar_LN','DRTQStar_LOG10'},'n');
SHStar = getStats(d,P,{'DRTHStar_LN','DRTHStar_LOG10'},'n');


%% DRT/DRT* list w/out zero

% -- DRT list 

% arithmetic
drtlistA_ln = SA{1}(2:end,1);
drtlistA_log10 = SA{2}(2:end,1);

% geometric
drtlistG_ln = SG{1}(2:end,1);
drtlistG_log10 = SG{2}(2:end,1);

% normalised
drtlistN_ln = SN{1}(2:end,1);
drtlistN_log10 = SN{2}(2:end,1);

% quadratic
drtlistQ_ln = SQ{1}(2:end,1);
drtlistQ_log10 = SQ{2}(2:end,1);

% harmonic
drtlistH_ln = SH{1}(2:end,1);
drtlistH_log10 = SH{2}(2:end,1);

% -- DRT* list

% arithmetic
drtlistStarA_ln = SAStar{1}(2:end,1);
drtlistStarA_log10 = SAStar{2}(2:end,1);

% geometric
drtlistStarG_ln = SGStar{1}(2:end,1);
drtlistStarG_log10 = SGStar{2}(2:end,1);

% normalised
drtlistStarN_ln = SNStar{1}(2:end,1);
drtlistStarN_log10 = SNStar{2}(2:end,1);

% quadratic
drtlistStarQ_ln = SQStar{1}(2:end,1);
drtlistStarQ_log10 = SQStar{2}(2:end,1);

% harmonic
drtlistStarH_ln = SHStar{1}(2:end,1);
drtlistStarH_log10 = SHStar{2}(2:end,1);


%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

% cells with value
active = find(~isnan(Ind));


plot_all = false;

if plot_all == true

    %% plotting DRT Histograms 

    switch log_base

        case 'ln'
            % - DRT > 0 & all averages 
            b1d1a = P.DRTA_LN(:); b1d1a = b1d1a(active); b1d1a = b1d1a(b1d1a>0);
            b1d1g = P.DRTG_LN(:); b1d1g = b1d1g(active); b1d1g = b1d1g(b1d1g>0);
            b1d1n = P.DRTN_LN(:); b1d1n = b1d1n(active); b1d1n = b1d1n(b1d1n>0);
            b1d1q = P.DRTQ_LN(:); b1d1q = b1d1q(active); b1d1q = b1d1q(b1d1q>0);
            b1d1h = P.DRTH_LN(:); b1d1h = b1d1h(active); b1d1h = b1d1h(b1d1h>0);

            % fields for plotting
            flds = {'DRTA_LN','DRTG_LN','DRTN_LN','DRTQ_LN','DRTH_LN'};

        case 'log10'

            b1d1a = P.DRTA_LOG10(:); b1d1a = b1d1a(active); b1d1a = b1d1a(b1d1a>0);
            b1d1g = P.DRTG_LOG10(:); b1d1g = b1d1g(active); b1d1g = b1d1g(b1d1g>0);
            b1d1n = P.DRTN_LOG10(:); b1d1n = b1d1n(active); b1d1n = b1d1n(b1d1n>0);
            b1d1q = P.DRTQ_LOG10(:); b1d1q = b1d1q(active); b1d1q = b1d1q(b1d1q>0);
            b1d1h = P.DRTH_LOG10(:); b1d1h = b1d1h(active); b1d1h = b1d1h(b1d1h>0);

            % fields for plotting
            flds = {'DRTA_LOG10','DRTG_LOG10', ...
                    'DRTN_LOG10','DRTQ_LOG10','DRTH_LOG10'};    

    end

    f = figure;

    % -- arithmetic
    h1a = histfit(b1d1a,hnb,'lognormal');

    % patch
    h1a(1).EdgeColor = 'None';
    h1a(1).FaceColor = [0.7,0.0,0.0];
    h1a(1).FaceAlpha = 0.1;

    % line
    h1a(2).LineWidth = 2; 
    h1a(2).Color = [0.7,0.0,0.0];

    hold on 

    % -- geometric
    h1g = histfit(b1d1g,hnb,'lognormal');

    % patch
    h1g(1).EdgeColor = 'None';
    h1g(1).FaceColor = [0.0,0.7,0.0];
    h1g(1).FaceAlpha = 0.1;

    % line
    h1g(2).LineWidth = 2; 
    h1g(2).Color = [0.0,0.7,0.0];

    % -- normalised
    h1n = histfit(b1d1n,hnb,'lognormal');

    % patch
    h1n(1).EdgeColor = 'None';
    h1n(1).FaceColor = [0.0,0.0,0.7];
    h1n(1).FaceAlpha = 0.1;

    % line
    h1n(2).LineWidth = 2; 
    h1n(2).Color = [0.0,0.0,0.7];

    % -- quadratic
    h1q = histfit(b1d1q,hnb,'lognormal');

    % patch
    h1q(1).EdgeColor = 'None';
    h1q(1).FaceColor = [0.7,0.7,0.7];
    h1q(1).FaceAlpha = 0.1;

    % line
    h1q(2).LineWidth = 2; 
    h1q(2).Color = [0.7,0.7,0.7];

    % -- harmonic
    h1h = histfit(b1d1h,hnb,'lognormal');

    % patch
    h1h(1).EdgeColor = 'None';
    h1h(1).FaceColor = [0.7,0.7,0.0];
    h1h(1).FaceAlpha = 0.1;

    % line
    h1h(2).LineWidth = 2; 
    h1h(2).Color = [0.8,0.8,0.0];

    % Customising legend entries 
    set(get(get(h1a(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1g(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1n(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1q(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    leg = legend('$\mathcal{N}_A$','$\mathcal{N}_G$','$\mathcal{N}_N$',... 
           '$\mathcal{N}_Q$','$\mathcal{N}_H$');
    leg.FontSize = FS;

    % labels 
    xlabel('$DRT$','FontSize',FS);
    ylabel('$n_{e}$','FontSize',FS);

    % axis limits
    switch log_base
        case 'ln'
            ylim([0,2.2e4]);
        case 'log10'
            ylim([0,3.0e4]);
    end

    %print('../tmp/histogram-drt.eps','-depsc2');


    %% Plotting DRT distribution over reservoir


    for f = 1:5
        figure;
        prop = P.(flds{f});     
        prop = prop(:); prop = prop(active);    
        plotCellData(G,prop,'EdgeColor','none','BackFaceLighting','lit');
        axis off vis3d
        colormap(jet)
        cbar = colorbar;
        cbar.FontSize = FS;
        cbar.Box = 'off';  
        aux = split(flds{f},'_'); 
        cbar.Label.String = [aux(1),' ',aux(2)];

        fn = strcat('../tmp/3d-zoning-drt-',num2str(f),'.eps');
        %print(fn,'-depsc2');
    end      

    % - DRT* > 0 & all averages
    switch log_base

        case 'ln'

            b1d2a = P.DRTAStar_LN(:); b1d2a = b1d2a(active); b1d2a = b1d2a(b1d2a>0);
            b1d2g = P.DRTGStar_LN(:); b1d2g = b1d2g(active); b1d2g = b1d2g(b1d2g>0);
            b1d2n = P.DRTNStar_LN(:); b1d2n = b1d2n(active); b1d2n = b1d2n(b1d2n>0);
            b1d2q = P.DRTQStar_LN(:); b1d2q = b1d2q(active); b1d2q = b1d2q(b1d2q>0);
            b1d2h = P.DRTHStar_LN(:); b1d2h = b1d2h(active); b1d2h = b1d2h(b1d2h>0);

            % fields for plotting
            flds = {'DRTAStar_LN','DRTGStar_LN',... 
                    'DRTNStar_LN','DRTQStar_LN','DRTHStar_LN'};

        case 'log10'

            b1d2a = P.DRTAStar_LOG10(:); b1d2a = b1d2a(active); b1d2a = b1d2a(b1d2a>0);
            b1d2g = P.DRTGStar_LOG10(:); b1d2g = b1d2g(active); b1d2g = b1d2g(b1d2g>0);
            b1d2n = P.DRTNStar_LOG10(:); b1d2n = b1d2n(active); b1d2n = b1d2n(b1d2n>0);
            b1d2q = P.DRTQStar_LOG10(:); b1d2q = b1d2q(active); b1d2q = b1d2q(b1d2q>0);
            b1d2h = P.DRTHStar_LOG10(:); b1d2h = b1d2h(active); b1d2h = b1d2h(b1d2h>0);

              % fields for plotting
            flds = {'DRTAStar_LOG10','DRTGStar_LOG10', ...
                    'DRTNStar_LOG10','DRTQStar_LOG10','DRTHStar_LOG10'};    
    end

    f2 = figure;

    % -- arithmetic
    h1a = histfit(b1d2a,hnb,'lognormal');

    % patch
    h1a(1).EdgeColor = 'None';
    h1a(1).FaceColor = [0.7,0.0,0.0];
    h1a(1).FaceAlpha = 0.1;

    % line
    h1a(2).LineWidth = 2; 
    h1a(2).Color = [0.7,0.0,0.0];

    hold on 

    % -- geometric
    h1g = histfit(b1d2g,hnb,'lognormal');

    % patch
    h1g(1).EdgeColor = 'None';
    h1g(1).FaceColor = [0.0,0.7,0.0];
    h1g(1).FaceAlpha = 0.1;

    % line
    h1g(2).LineWidth = 2; 
    h1g(2).Color = [0.0,0.7,0.0];

    % -- normalised
    h1n = histfit(b1d2n,hnb,'lognormal');

    % patch
    h1n(1).EdgeColor = 'None';
    h1n(1).FaceColor = [0.0,0.0,0.7];
    h1n(1).FaceAlpha = 0.1;

    % line
    h1n(2).LineWidth = 2; 
    h1n(2).Color = [0.0,0.0,0.7];

    % -- quadratic
    h1q = histfit(b1d2q,hnb,'lognormal');

    % patch
    h1q(1).EdgeColor = 'None';
    h1q(1).FaceColor = [0.7,0.7,0.7];
    h1q(1).FaceAlpha = 0.1;

    % line
    h1q(2).LineWidth = 2; 
    h1q(2).Color = [0.7,0.7,0.7];

    % -- harmonic
    h1h = histfit(b1d2h,hnb,'lognormal');

    % patch
    h1h(1).EdgeColor = 'None';
    h1h(1).FaceColor = [0.7,0.7,0.0];
    h1h(1).FaceAlpha = 0.1;

    % line
    h1h(2).LineWidth = 2; 
    h1h(2).Color = [0.8,0.8,0.0];

    % Customising legend entries 
    set(get(get(h1a(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1g(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1n(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1q(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(h1h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    leg = legend('$\mathcal{N}^{*}_A$','$\mathcal{N}^{*}_G$', ... 
                 '$\mathcal{N}^{*}_N$','$\mathcal{N}^{*}_Q$', ... 
                 '$\mathcal{N}^{*}_H$');
    leg.FontSize = FS;

    % labels 
    xlabel('$DRT^{*}$','FontSize',FS);
    ylabel('$n_{e}$','FontSize',FS);

    % axis limits
    switch log_base
        case 'ln'
            ylim([0,2.2e4]);
        case 'log10'
            ylim([0,3.0e4]);
    end

    %print('../histogram-drt*.eps','-depsc2');

    %% Plotting DRT* distribution over reservoir

    for f = 1:5
        figure;
        prop = P.(flds{f});     
        prop = prop(:); prop = prop(active);    
        plotCellData(G,prop,'EdgeColor','none','BackFaceLighting','lit');
        axis off vis3d
        colormap(jet)
        cbar = colorbar;
        cbar.FontSize = FS;
        cbar.Box = 'off';      
        aux = split(flds{f},'_'); 
        cbar.Label.String = [aux(1),' ',aux(2)];

        fn = strcat('../tmp/3d-zoning-drt*-',num2str(f),'.eps');
        %print(fn,'-depsc2');
    end      

end
%% DRT connections

switch log_base
    
    case 'ln'    

        drtStA_ln = findDRTConnections(d,drtlistA_ln, P, 'arithmetic','ln',nofs,'n', 1);
        drtStG_ln = findDRTConnections(d,drtlistG_ln, P, 'geometric','ln',nofs,'n', 1);
        drtStN_ln = findDRTConnections(d,drtlistN_ln, P, 'normalized','ln',nofs,'n', 1);
        drtStQ_ln = findDRTConnections(d,drtlistQ_ln, P, 'quadratic','ln',nofs,'n', 1);
        drtStH_ln = findDRTConnections(d,drtlistH_ln, P, 'harmonic','ln',nofs,'n', 1);

    case 'log10'

        drtStA_log10 = findDRTConnections(d,drtlistA_log10, P, 'arithmetic','log10',nofs,'n', 1);
        drtStG_log10 = findDRTConnections(d,drtlistG_log10, P, 'geometric','log10',nofs,'n', 1);
        drtStN_log10 = findDRTConnections(d,drtlistN_log10, P, 'normalized','log10',nofs,'n', 1);
        drtStQ_log10 = findDRTConnections(d,drtlistQ_log10, P, 'normalized','log10',nofs,'n', 1);
        drtStH_log10 = findDRTConnections(d,drtlistH_log10, P, 'normalized','log10',nofs,'n', 1);

end


%% Computing graph metrics 

opt.nofs = nofs;    
opt.seps = 0.05;
opt.R2min = 0.9;

% compute
switch log_base
    
    case 'ln'    
        
        [metricsStA_ln,linregrStA_ln] = computeDRTGraphMetrics(opt,drtStA_ln);
        [metricsStG_ln,linregrStG_ln] = computeDRTGraphMetrics(opt,drtStG_ln);
        [metricsStN_ln,linregrStN_ln] = computeDRTGraphMetrics(opt,drtStN_ln);
        [metricsStQ_ln,linregrStQ_ln] = computeDRTGraphMetrics(opt,drtStQ_ln);
        [metricsStH_ln,linregrStH_ln] = computeDRTGraphMetrics(opt,drtStH_ln);
        
        save('../tmp/metricsStA_ln.mat','metricsStA_ln');
        save('../tmp/metricsStG_ln.mat','metricsStG_ln');
        save('../tmp/metricsStN_ln.mat','metricsStN_ln');
        save('../tmp/metricsStQ_ln.mat','metricsStQ_ln');
        save('../tmp/metricsStH_ln.mat','metricsStH_ln');
        
        save('../tmp/linregrStA_ln.mat','linregrStA_ln');
        save('../tmp/linregrStG_ln.mat','linregrStG_ln');
        save('../tmp/linregrStN_ln.mat','linregrStN_ln');
        save('../tmp/linregrStQ_ln.mat','linregrStQ_ln');
        save('../tmp/linregrStH_ln.mat','linregrStH_ln');

    case 'log10'
        
        [metricsStA_log10,linregrStA_log10] = computeDRTGraphMetrics(opt,drtStA_log10);
        [metricsStG_log10,linregrStG_log10] = computeDRTGraphMetrics(opt,drtStG_log10);
        [metricsStN_log10,linregrStN_log10] = computeDRTGraphMetrics(opt,drtStN_log10);
        [metricsStQ_log10,linregrStQ_log10] = computeDRTGraphMetrics(opt,drtStQ_log10);
        [metricsStH_log10,linregrStH_log10] = computeDRTGraphMetrics(opt,drtStH_log10);
        
        save('../tmp/metricsStA_ln.mat','metricsStA_ln');
        save('../tmp/metricsStG_ln.mat','metricsStG_ln');
        save('../tmp/metricsStN_ln.mat','metricsStN_ln');
        save('../tmp/metricsStQ_ln.mat','metricsStQ_ln');
        save('../tmp/metricsStH_ln.mat','metricsStH_ln');
        
        save('../tmp/linregrStA_log10.mat','linregrStA_log10');
        save('../tmp/linregrStG_log10.mat','linregrStG_log10');
        save('../tmp/linregrStN_log10.mat','linregrStN_log10');
        save('../tmp/linregrStQ_log10.mat','linregrStQ_log10');
        save('../tmp/linregrStH_log10.mat','linregrStH_log10');

end

fprintf('---> Execution finished for model ''%s''.\n',model);