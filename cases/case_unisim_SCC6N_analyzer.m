%% Case: UNISIM - Constrained Clustering Analyzer - 6N
% 
% Analysis of results concerning 6N algorithm. It is done separated from 
% the others

%% Input

% analyse clusters with >= minc members
minc = 2;

load('./data/SCC6N_unisim1.mat','SCC')
SCC6N = SCC; clear SCC;

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
G = computeGeometry(G);

%% Compute required parameters
P = computeParams(G,PROPS);

%% Log data

% log10(phiz) x log10(RQI); normalized averaging technique 
Log10PHIZ = P.Log10PHIZ(:); 
Log10RQIN = P.Log10RQIN(:);

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;


%% Clusters with >= minc elements 
% This is done because SCC6N delivers clusters directly, not partitions

id6n = find(cellfun(@numel,SCC6N.partitioning) >= minc); 

if ~all(cellfun(@numel,SCC6N.partitioning(id6n)) - 2*ones(size(id6n))) && ...
        minc > 2
    error(['All the partitions have 2 elements. ' ... 
            'Setting up minc > 2 is nonsense.']);    
end

% slope and R2
slopes = zeros(size(id6n));
R2 = zeros(size(id6n));
B = zeros(size(id6n));
for i = 1:length(id6n)
    [R,m,b] = regression(Log10PHIZ(SCC6N.partitioning{id6n(i)}), ...
                           Log10RQIN(SCC6N.partitioning{id6n(i)}),'one');                       
    R2(i) = R*R;
    slopes(i) = m;
    B(i) = b;
end
SCC6N.slopes = slopes;
SCC6N.R2 = R2;
SCC6N.B = B;

save('../mat/SCC6N.mat','SCC6N');

% plot of slope distribution (outlier safety)
seps = 1e-1;
points = find(1-seps <= slopes & slopes <= 1+seps);   
slopes = slopes(points);     
R2 = R2(points); 
delta = 1:numel(points);


figure
set(gca,'FontSize',14);
hold on, box on    

plot(delta,slopes,'o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])            
plot(delta,ones(1,numel(delta)),'k-')    
plot(delta,ones(1,numel(delta))*(1-seps),'r--')    
plot(delta,ones(1,numel(delta))*(1+seps),'r--')    
xlim([-1,max(delta)+1])
ylim([1-2*seps,1+2*seps])
xlabel('$\gamma$','interpreter','latex')
ylabel('$s_{\gamma}$','interpreter','latex')
xticks([1,max(delta)])
yticks([1-seps,1,1+seps])
hold off
fname = strcat('../tmp/slope_partitions6n','.eps');
print(fname,'-depsc2')

% plot of R2 distribution
figure 
set(gca,'FontSize',14);
hold on, box on         
plot(delta,R2,'o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])            
plot(delta,ones(1,numel(delta))*mean(R2),'k-')    
plot(delta,ones(1,numel(delta))*min(R2),'r--')    
plot(delta,ones(1,numel(delta))*max(R2),'r--')    
xlim([-0.5,max(delta)+0.5])
ylim([min(R2)-0.1,max(R2)+0.1])
xlabel('$\gamma$','interpreter','latex')
ylabel('$R^2_{\gamma}$','interpreter','latex')
xticks([1,max(delta)])
yticks(mean(R2))
hold off
fname = strcat('../tmp/R2_partitions6n','.eps');
print(fname,'-depsc2')


