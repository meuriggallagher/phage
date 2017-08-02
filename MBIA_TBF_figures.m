%MODEL-BASED IMAGE ANALYSIS OF A TETHERED BROWNIAN FIBRE FOR SHEAR STRESS
%SENSING
%
% Generates figures from the paper
%
% Requires the relevant .csv files from epapers.bham.ac.uk
%
%% Figure 2
% Load from .csv
fig2 = csvread('figure2.csv',1,0);
theta = fig2(:,1);
phi = fig2(:,2);
alpha_t = fig2(:,3);
alpha_p = fig2(:,4);

% Reshape
nTh = length(unique(theta));
nPh = length(unique(phi));
theta = reshape(theta,[nPh,nTh]);
phi = reshape(phi,[nPh,nTh]);
alpha_t = reshape(alpha_t,[nPh,nTh]);
alpha_p = reshape(alpha_p,[nPh,nTh]);

% Plot
figure(2)
set(gcf,'position',[200 900 1040 424])
subplot(1,2,1)
mesh(theta,phi,alpha_t)
shading interp
set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) ...
    diff(get(gca, 'XLim')) diff(get(gca, 'ZLim'))])
zlim([-1 1])
xlim([0 90])
ylim([0 360])
colormap(cool)
view(-43.1,10)
box on
ylabel('\phi','fontsize',14)
xlabel('\theta','fontsize',14)
zlabel('\alpha_\theta^\prime','rot',0,'fontsize',14)

subplot(1,2,2)
mesh(theta,phi,alpha_p)
shading interp
set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) ...
    diff(get(gca, 'XLim')) diff(get(gca, 'ZLim'))])
zlim([-1 1])
xlim([0 90])
ylim([0 360])
colormap(cool)
view(-43.1,10)
box on
ylabel('\phi','fontsize',14)
xlabel('\theta','fontsize',14)
zlabel('\alpha_\phi^\prime','rot',0,'fontsize',14)

suptitle(['Components of the dimensionless rotational advection'...
    ' \alpha^\prime plotted for' ...
    ' 0^\circ \leq \phi < 360^\circ, 0^\circ\leq \theta \leq 90^\circ.'])

%% Figure 3
% Load from .csv
fig3 = csvread('figure3.csv',1,0);
theta = fig3(:,1);
D_tt = fig3(:,2);
D_pp = fig3(:,3);

% Plot
figure(3)
set(gcf,'position',[225 875 1040 424])
subplot(1,2,1)
plot(theta,D_tt)
box on
ylabel('D^\prime_{\theta\theta}','rot',0)
xlabel('\theta')

subplot(1,2,2)
plot(theta,D_pp)
box on
ylabel('D^\prime_{\phi\phi}','rot',0)
xlabel('\theta')

suptitle(['Non-zero components of the dimensionless diffusion matrix ' ...
    'D^\prime plotted against 0^\circ \leq \theta \leq 90^\circ'])


%% Figure 4
% Load from .csv
fig4 = csvread('figure4.csv',1,0);
Pe = fig4(:,1);
phi = fig4(:,2);
Phi = fig4(:,3);

% Reshape
nPe = length(unique(Pe));
nPh = length(unique(phi));
Pe = reshape(Pe,[nPh,nPe]);
phi = reshape(phi,[nPh,nPe]);
Phi = reshape(Phi,[nPh,nPe]);


figure(4)
set(gcf,'position',[250 850 1040 424])
subplot(1,2,1)
surfl(Pe,phi,Phi)
view(-43.1,10)
ylim([-180,180])
zlim([0 3.5])
colormap(summer)
shading interp
ylabel('\phi','fontsize',14)
xlabel('\theta','fontsize',14)
zlabel('\Phi','rot',0,'fontsize',14)

subplot(1,2,2)
surf(Pe,phi,Phi)
view(0,90)
ylim([-180,180])
colormap(summer)
box on
shading interp
ylabel('\phi','fontsize',14,'rot',0)
xlabel('\theta','fontsize',14)

suptitle(['Marginal probability density function \Phi, the steady ' ...
    'solution to the advection diffusion equation,\newline plotted with ' ...
    '-90^\circ \leq \phi \leq 90^\circ, and 1 \leq Pe \leq 200'])

%% Figure 5
% Load from .csv
fig5 = csvread('figure5.csv',1,0);
phiSim = fig5(:,1);
phiPrep = fig5(:,2);
phiNoPrep = fig5(:,3);

% Calculate differences with and without preprocessing
dFit = phiSim - phiPrep;
dNop = phiSim - phiNoPrep;

figure(5)
subplot(1,2,1)
set(gcf,'position',[275 825 1040 424])
[nFreqFit,edges] = histcounts(dFit,'binwidth',1);
nFreqFit = nFreqFit / length(phiPrep);
hbar = bar(edges(1:end-1) + (edges(2)-edges(1))/2, nFreqFit,1,'stacked');
set(hbar, 'FaceColor', [0.3,0.75,0.83], 'EdgeColor', 'k'); 
xlim([-90 90])
xlabel('$\phi - \bar{\phi}$','fontsize',14,'interpreter','latex')
box on
title('With preprocessing')

subplot(1,2,2)
[nFreqNop,edges] = histcounts(dNop,'binwidth',1);
nFreqNop = nFreqNop / length(phiNoPrep);
hbar = bar(edges(1:end-1) + (edges(2)-edges(1))/2, nFreqNop,1,'stacked');
set(hbar, 'FaceColor', [0.3,0.75,0.83], 'EdgeColor', 'k'); 
xlim([-90 90])
xlabel('$\phi - \bar{\phi}$','fontsize',14,'interpreter','latex')
box on
title('No preprocessing')

suptitle(['Relative frequency histograms of the error between the ' ...
    'simulated M13 angles and fit angles, with 1^\circ bin widths, ' ...
    'for the analysis of the image fitting procedure.'])

%% Figure 6
% Load from .csv
fig6 = csvread('figure6.csv',1,0);
PeSim = fig6(:,1);
PePrep = fig6(:,2);
PeNoPrep = fig6(:,3);

% Mean and standard devisations
muPrep = mean(PePrep - PeSim);
muNoPrep = mean(PeNoPrep - PeSim);
sigPrep = std(PePrep - PeSim);
sigNoPrep = std(PeNoPrep - PeSim);

% Lines of best fit
bfitPrep = polyfit(PeSim,PePrep,1);
bfitNoPrep = polyfit(PeSim,PeNoPrep,1);

figure(6)
set(gcf,'position',[300 800 1040 424])
subplot(1,2,1)
plot(PeSim,PeSim,'k'),hold on
plot(PeSim,PePrep,'b.','markersize',12)
plot(PeSim,PeNoPrep,'r.','markersize',12)
plot(PeSim,bfitPrep(2) + bfitPrep(1)*PeSim,'b--','linewidth',2)
plot(PeSim,bfitNoPrep(2) + bfitNoPrep(1)*PeSim,'r--','linewidth',2)
xlim([0 200])
box on
ylabel('$\overline{\mathrm{Pe}}$','fontsize',14,'interpreter','latex')
xlabel('$\mathrm{Pe}$','fontsize',14,'interpreter','latex')

subplot(1,2,2)
plot([1,225],muPrep*[1,1],'b'),hold on
plot([1,225],muPrep*[1,1] + 1.96*sigPrep,'b--'),hold on
plot([1,225],muPrep*[1,1] - 1.96*sigPrep,'b--'),hold on
plot(0.5*(PePrep+PeSim),(PePrep-PeSim),'b.','markersize',12)
plot([1,225],muNoPrep*[1,1],'r'),hold on
plot([1,225],muNoPrep*[1,1] + 1.96*sigNoPrep,'r--'),hold on
plot([1,225],muNoPrep*[1,1] - 1.96*sigNoPrep,'r--'),hold on
plot(0.5*(PeNoPrep+PeSim),(PeNoPrep-PeSim),'r.','markersize',12)
xlim([0 225])
ylabel('$\overline{\mathrm{Pe}} - \mathrm{Pe}$','fontsize',14, ...
    'interpreter','latex')
xlabel('$\left(\mathrm{Pe} + \overline{\mathrm{Pe}}\right)/2$', ...
    'fontsize',14,'interpreter','latex')

suptitle(['(a) plots the fit Peclet number against simulated Peclet ' ...
    'for analysis of the image fitting procedure. \newlineThe blue dots shows ' ...
    'what would be perfect correspondence between the fit and ' ...
    'simulated angles, \newline with the dotted lines being the line of best ' ...
    'fit to the data. \newline (b) shows the Bland-Altman plot testing the fit ' ...
    'data against simulated data with preprocessing (blue) and without '...
    '(red).\newline The solid lines show the mean of the difference for each ' ...
    'for each case, with the dotted lines being the 95% confidence ' ...
    'interval for the difference.']);

%% Figure 7
% Load from .csv
fig7 = csvread('figure7.csv',1,0);
phiSim = fig7(:,1);
phiFull = fig7(:,2);

% Calculate differences with and without preprocessing
dFull = phiSim - phiFull;

figure(5)
set(gcf,'position',[325 775 1040 424])
[nFreqFull,edges] = histcounts(dFull,'binwidth',1);
nFreqFull = nFreqFull / length(phiFull);
hbar = bar(edges(1:end-1) + (edges(2)-edges(1))/2, nFreqFull,1,'stacked');
set(hbar, 'FaceColor', [0.3,0.75,0.83], 'EdgeColor', 'k'); 
xlim([-90 90])
xlabel('$\phi - \bar{\phi}$','fontsize',14,'interpreter','latex')
box on

suptitle(['Relative frequency histograms of the error between the ' ...
    'simulated M13 angles and fit angles, with 1^\circ bin widths, ' ...
    'for the full analysis.'])

%% Figure 8
% Load from .csv
fig8 = csvread('figure8.csv',1,0);
Pe = fig8(:,1);
PeSim = fig8(:,2);
PeFull = fig8(:,3);

% Mean and standard devisations
muFull = mean(PeFull - PeSim);
sigFull = std(PeFull - PeSim);

% Lines of best fit
bfitFull = polyfit(PeSim,PeFull,1);

figure(8)
set(gcf,'position',[350 750 1040 424])
subplot(1,2,1)
plot(Pe,Pe,'k'),hold on
plot(Pe,PeFull,'b.','markersize',12)
plot(Pe,PeSim,'r.','markersize',12)
plot(Pe,bfitFull(2) + bfitFull(1)*Pe,'b--','linewidth',2)
xlim([0 200])
box on
ylabel('$\overline{\mathrm{Pe}}$','fontsize',14,'interpreter','latex')
xlabel('$\mathrm{Pe}$','fontsize',14,'interpreter','latex')

subplot(1,2,2)
plot([1,300],muFull*[1,1],'b'),hold on
plot([1,300],muFull*[1,1] + 1.96*sigFull,'b--'),hold on
plot([1,300],muFull*[1,1] - 1.96*sigFull,'b--'),hold on
plot(0.5*(PeFull+PeSim),(PeFull-PeSim),'b.','markersize',12)
xlim([0 300])
ylabel('$\overline{\mathrm{Pe}} - \mathrm{Pe}$','fontsize',14, ...
    'interpreter','latex')
xlabel('$\left(\mathrm{Pe} + \overline{\mathrm{Pe}}\right)/2$', ...
    'fontsize',14,'interpreter','latex')

suptitle(['(a) plots the fit Peclet number against simulated Peclet ' ...
    'for the full analysis. \newlineThe blue dots shows ' ...
    'what would be perfect correspondence between the fit and ' ...
    'simulated angles, \newline while the red dots show the ' ...
    'calculated Peclet number assuming the image analysis step is ' ...
    'perfect. The black line shows what would be perfect ' ...
    'correspondence, while the\newline dotted blue line is the line of best ' ...
    'fit to the data.\newline (b) shows the Bland-Altman plot, with the '...
    'dotted lines being the 95% confidence interval for the difference.'])

%% Figure 9
% Load from .csv
fig9 = csvread('figure9_10a_10b.csv',1,0);
WSS = fig9(:,1);
t = fig9(:,2);
phi = fig9(:,3);

figure(9)
set(gcf,'position',[375 725 1040 424])
colours = parula;
nC = size(colours,1);
hold on
for ii = 1:length(WSS)
    switch WSS(ii)
        case 3.5
            nColour = 1;
        case -3.5
            nColour = 2;
        case 2.5
            nColour = 3;
        case -2.5
            nColour = 4;
        case 1.5
            nColour = 5;
        case -1.5
            nColour = 6;
        case 1
            nColour = 7;
        case -1
            nColour = 8;
        case 0.5
            nColour = 9;
        otherwise
            nColour = 0;
    end
    
    if nColour == 0
        plot(t(ii),phi(ii),'.','color',[105,105,105]/255,'markersize',8)
    else
        plot(t(ii),phi(ii),'.','color', ...
            colours(round(nC/9 * (nColour-1)) + 1,:),'markersize',8)
    end
end
hold off
axis tight
box on
ylabel('\phi','fontsize',14,'rot',0)
xlabel('time','fontsize',14)

suptitle(['Plot showing the angels obtained by fitting to the raw ' ...
    'image data of Lobo et al [11].'])


%% Figure 10
% Load from .csv
fig10 = csvread('figure10c.csv',1,0);
WSS2 = fig10(:,1);
Pe = fig10(:,2);

% Calculate distributions of angles
pdfs = cell(5,1);
for ii = 1:5
    if ii == 1
        ind = find(abs(WSS) == 1);
    elseif ii == 2
        ind = find(abs(WSS) == 3.5);
    elseif ii == 3
        ind = find(abs(WSS) == 2.5);
    elseif ii == 4
        ind = find(abs(WSS) == 1.5);
    else
        ind = find(abs(WSS) == 0.5);
    end
    pdfs{ii} = fitdist(phi(ind),'Normal');
end

figure(10)
set(gcf,'position',[400 700 1040 424])
subplot(2,1,1)
fPhi = linspace(-90,90,5e2);
for ii = 1:5
    dum = pdf(pdfs{ii},fPhi);
    plot(fPhi,dum),hold on
end
hold off
box on
xlim([-90,90])
xlabel('\phi','fontsize',14)

subplot(2,2,3)
hold on
plot(WSS,phi,'k.','markersize',0.5)
for ii = 1:10
    if ii == 1
        ind = find(WSS == 3.5);
    elseif ii == 2
        ind = find(WSS == -3.5);
    elseif ii == 3
        ind = find(WSS == 2.5);
    elseif ii == 4
        ind = find(WSS == -2.5);
    elseif ii == 5
        ind = find(WSS == 1.5);
    elseif ii == 6
        ind = find(WSS == -1.5);
    elseif ii == 7
        ind = find(WSS == 1);
    elseif ii == 8
        ind = find(WSS == -1);
    elseif ii == 9
        ind = find(WSS == 0.5);
    else
        ind = find(WSS == 0);
    end
    
    mu = mean(phi(ind));
    sig = std(phi(ind));
    
    plot(WSS(ind),mu,'bx')
    plot(WSS(ind),mu + sig, 'bx')
    plot(WSS(ind),mu - sig, 'bx')
    plot(WSS(ind(1))*[1,1],mu + [-sig,sig],'b','linewidth',2);
end
box on
ylabel('\phi','rot',0)
xlabel('Nominal WSS (dyn cm$^{-2}$)','interpreter','latex')

subplot(2,2,4)
hold on
plot(WSS2(1),Pe(1),'.','markersize',10,'color',[0.5,0.5,0.5])
for ii = 2:length(WSS2)
    if WSS2(ii) > 0
        plot(WSS2(ii),Pe(ii),'b.','markersize',10)
    else
        plot(abs(WSS2(ii)),Pe(ii),'r.','markersize',10)
    end
end
box on
ylabel('Pe','rot',0)
xlabel('Nominal WSS (dyn cm$^{-2}$)','interpreter','latex')