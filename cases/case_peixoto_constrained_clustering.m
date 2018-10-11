

%% Tests

% test 1: slope = 1; 10 collinear points
% test 2: slope = 1; 30 points; 3 parallel straigth lines
% test 3: slope = 0.5,1,1.5; 3 straigth lines
% test 4: slope = -1,-0.5,0.5,1; 4 straigth lines
% test 5: random 1
% test 6: random 2
% test 7: data from unisim 

% data matrix
load ../mat/L.mat;

test = 3;

%%
switch test
    
    case 1        
        X = 1:10;
        Y = X;
        seps = 1e-2;
        i0 = 1;
        
    case 2       
        x = 1:10;
        X = [x,x,x];
        Y = [x+1,x,x-1];
        seps = 1e-2;
        i0 = 1;
        
    case 3       
        x = 1:5;
        X = [x,x,x];
        Y = [0.5*x,x,1.5*x];
        seps = 1e-2;
        i0 = 1;
    
    case 4       
        x = -5:5;
        X = [x,x,x,x];
        Y = [-x,-0.5*x,0.5*x,x];
        seps = 1e-2;
        i0 = 1;
    
    case 5       
        x = 10*rand(1,20);        
        X = [-x,x/2,2*x,x];
        Y = [log(x),sin(x),x.^2,x];
        seps = 1e-2;
        i0 = 1;
        
    case 6       
        x = 2*rand(1,50);        
        X = [x,x/2];
        Y = [log(x),log(x)];
        seps = 1e-2;
        i0 = 1;    
        
    case 7
        X = L.DRT14.logPHIZ{3};
        Y = L.DRT14.logRQI{3};        
        seps = 1e-2;
        i0 = 1;    
                
end


%% clustering
Clusters = slopeConstrainedClustering(X,Y,i0,seps);

%% Plotting best-fit lines for clusters 

for c = 1:numel(Clusters)
    
    IC = Clusters{c};        
    [R,m,b] = regression(X(IC),Y(IC),'one');
    
    x = X(IC);
    f = plot(x,m*x + b,'o-');
        
    hold on;
    labels{c} = strcat('C',num2str(c));
    
end

legend(labels,'location','bestoutside');



