%% Case - DRT comprehension 
% Discrete Rock Type (Guo, 2005) has being used as a standard tool in
% literature, although very little discussion about it is available. This
% case seeks to give a detailed comprehension of DRT for UNISIM-I. 


%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');

%% Compute FZI

P = computeParams(G,PROPS);
FZIN = P.FZIN;
FZIN = FZIN(:);
locs0 = find(FZIN > 0);
FZIN0 = FZIN(locs0);

%% DRT

% sorted DRT values
drt = @(k,C) sort(round(k*log(FZIN0) + C));

% Guo's formula 
plot(1:length(FZIN0),drt(2,10.6),'k');
axis tight
legend('Guo')
    
% Influence of amplification factor and adjust constant 0 
% with the formula DRT = round(b*log(FZI))
figure 
ns = 5;
s = cell(1,ns+1);
h = zeros(1,ns+1);
s{1} = 'Guo (2005)';
h(1) = plot(1:length(FZIN0),drt(2,10.6),'k--');
for k = 1:ns    
    hold on,
    h(k+1) = plot(1:length(FZIN0),drt(k,0));
    s{k+1} = sprintf('$C=0, d=$%d',k);
end
legend(h,s,'Location','northwest','interpreter','latex')
ylabel('$DRT$','interpreter','latex')

