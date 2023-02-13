%% Flow on Mars 5 - Grain Size Distribution
% Author: Lisanne Braat (lisannebraat@gmail.com)
% Last update: 2023-02-13
% Created in Matlab version: R2022b

%% Initialize
clear variables
close all
clc

output = 'FlowOnMars_exportfig';
addpath(genpath('Checkout'));

%% Input parameters
for a = 1
%Distribution
%peak = 0.00063;
%width = 2.1;
peak = 0.00063;
width = 1.5;

%Plot
pl.width = 18;
pl.height = 8;
pl.line = 1;
pl.line_ax = 0.75;
pl.fsz = 7;
pl.fsz2 = 8.5;
pl.vertical = 1;
pl.horizantal = 3;
pl.lmarge = 0.065;
pl.rmarge = 0.025;
pl.bmarge = 0.2;
pl.tmarge = 0.045;
pl.intv = 0.065;
pl.inth = 0.065;
pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
pl.c1 = [0 0.8 0.7];
pl.c2 = [0.8 0.7 0];
pl.c1a = [0 0.8 0.7 0.3];
pl.c2a = [0.8 0.7 0 0.3];
y = 0.94;
end

%% Create distribution and bins
for a = 1
% Create distribution
scale = peak*width; %width of the ditribution
shape = -1*width;  %shape of the curve, tial and error to get symetric
data = 1000000;
X = lognormal3_rnd(peak,scale,shape,data,1);

% Create bins
% Test bins small
bins_test = logspace(-7,1,1000);
D_test = zeros(1,length(bins_test)-1);
for i = 1:length(bins_test)-1
    D_test(i) = sqrt(bins_test(i+1)/bins_test(i))*bins_test(i);
end

% Bins based on sediment classes
bins = [2e-7 6.3e-7 2e-6 6.3e-6 20e-6 63e-6 0.0002 0.00063 0.002 0.0063 0.02 0.063 0.2 0.63 2]; %ISO 14688-1:2002 (International scale)
bins_txt = {'<clay' 'clay' 'fine silt' 'silt' 'coarse silt' 'fine sand' ...
    'sand' 'course sand' 'fine gravel' 'gravel' ...
    ' course gravel' 'coble' 'boulder' '>boulder'};
D = zeros(1,length(bins)-1);
for i = 1:length(bins)-1
    D(i) = sqrt(bins(i+1)/bins(i))*bins(i);
end
end

%% Plot first part
for a = 1
close all
f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width 18],'Position',[70 -5 pl.width 18], ...
        'PaperSize',[pl.width pl.height],'visible','on');
ax(1) = subplot(3,2,1);
hold on
h_test = histogram(X,bins_test,'FaceColor',pl.c1);
title('Random lognomal sample - small bins','fontsize',14)
xlabel('Grain size bins')
ylabel('Bin Count')
set(gca,'Xscale','log');
for i = 1:length(D)
    plot([D(i) D(i)],[0 max(h_test.Values)*1.2],'color',[0.5 0.5 0.5 0.2]);
end
ylim([0 max(h_test.Values)*1.2]);

ax(2) = subplot(3,2,2);
hold on
h = histogram(X,bins,'FaceColor',pl.c2);
title('Random lognormal sample - grain size classes','fontsize',14)
xlabel('Grain size classes')
ylabel('Bin Count')
set(gca,'Xscale','log');
for i = 1:length(D)
    plot([D(i) D(i)],[0 max(h.Values)*1.2],'color',[0.5 0.5 0.5 0.2]);
end
ylim([0 max(h.Values)*1.2]);
end

%% follow up calculations
for a = 1
bins_y_test = h_test.Values;
if ~(sum(bins_y_test)==data)
    warning(['Not all data falls within bins. Outside bins = ' num2str(100-sum(bins_y_test)/data*100) '%'])
end

cum_test = zeros(size(D_test));
cum_test(1) = bins_y_test(1);
for i = 2:length(bins_y_test)
    cum_test(i) = bins_y_test(i) + cum_test(i-1);
end
cum_perc_test = cum_test/cum_test(end)*100;
bins_y_perc_test = bins_y_test/cum_test(end)*100;

D10 = InterX([D_test; ones(size(D_test))*10],[D_test; cum_perc_test]); D10 = D10(1);
D50 = InterX([D_test; ones(size(D_test))*50],[D_test; cum_perc_test]); D50 = D50(1);
D90 = InterX([D_test; ones(size(D_test))*90],[D_test; cum_perc_test]); D90 = D90(1);

bins_y = h.Values;
if ~(sum(bins_y)==data)
    warning(['Not all data falls within bins. Outside bins = ' num2str(100-sum(bins_y)/data*100) '%. ' num2str(bins_y(1)/data*100) '% in lowest class.'])
end

cum = zeros(size(D));
cum(1) = bins_y(1);
for i = 2:length(bins_y)
    cum(i) = bins_y(i) + cum(i-1);
end
cum_perc = cum/cum(end)*100;
bins_y_perc = bins_y/cum(end)*100;

D10_course = InterX([bins(2:end); ones(size(D))*10],[bins(2:end); cum_perc]); D10_course = D10_course(1);
D50_course = InterX([bins(2:end); ones(size(D))*50],[bins(2:end); cum_perc]); D50_course = D50_course(1);
D90_course = InterX([bins(2:end); ones(size(D))*90],[bins(2:end); cum_perc]); D90_course = D90_course(1);
end

%% follow up plot
for a = 1
ax(3) = subplot(3,2,3);
hold on
plot(D_test,bins_y_perc_test,'color',pl.c1);
title('Lognormal mass distribution - small bins','fontsize',14)
xlabel('Grain size bins')
ylabel('% wt')
set(gca,'Xscale','log');
for i = 1:length(D)
    plot([D(i) D(i)],[0 max(bins_y_perc_test)*1.2],'color',[0.5 0.5 0.5 0.2]);
end
ylim([0 max(bins_y_perc_test)*1.2]);
set(gca,'XTick',D,'XtickLabel',bins_txt);

ax(4) = subplot(3,2,4);
hold on
plot(D,bins_y_perc,'color',pl.c2);
title('Lognormal mass distribution - grain size classes','fontsize',14)
xlabel('Grain size classes')
ylabel('% wt')
set(gca,'Xscale','log');
for i = 1:length(D)
    plot([D(i) D(i)],[0 max(bins_y_perc)*1.2],'color',[0.5 0.5 0.5 0.2]);
end
ylim([0 max(bins_y_perc)*1.2]);
set(gca,'XTick',D,'XtickLabel',bins_txt);

ax(5) = subplot(3,2,5);
hold on
plot(bins(2:end),cum_perc,'color',pl.c2);
plot(D_test,cum_perc_test,'color',pl.c1);
plot(D_test,ones(1,length(bins_y_test))*10,'color',pl.c1a);
plot(D_test,ones(1,length(bins_y_test))*50,'color',pl.c1a);
plot(D_test,ones(1,length(bins_y_test))*90,'color',pl.c1a);
plot([D10 D10],[0 100],'color',pl.c1a);
plot([D50 D50],[0 100],'color',pl.c1a);
plot([D90 D90],[0 100],'color',pl.c1a);
plot([D10_course D10_course],[0 100],'color',pl.c2a);
plot([D50_course D50_course],[0 100],'color',pl.c2a);
plot([D90_course D90_course],[0 100],'color',pl.c2a);
title('Lognormal cumulative mass distribution','fontsize',14)
xlabel('Grain size bins/classes')
ylabel('% wt')
set(gca,'Xscale','log');
set(gca,'XTick',D,'XtickLabel',bins_txt);
ylim([0 100]);

set(ax,'xlim',[bins_test(1) bins_test(end)]);
set(ax,'box','on','Layer','top', ...
    'XMinorTick','on','YMinorTick','off', ...
    'FontSize',pl.fsz,'LineWidth',pl.line,'TickDir','in','YDir','normal', ...
    'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
    'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);
end

%% Plot paper
for a = 1
close all
xd = 10^(log10(6.3e-7) + (log10(0.63)-log10(6.3e-7))*0.02);

f2 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[70 -5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

ax2(1) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
hold on
h_test = histogram(X,bins_test,'FaceColor',pl.c1,'EdgeColor',pl.c1); %Random distribution
for i = 1:length(D) %Bin centers of grain size classes
    plot([D(i) D(i)],[0 data],'color',[0.5 0.5 0.5 0.2]);
end
set(gca,'Xscale','log');
axis([6.3e-7 0.63 0 max(h_test.Values)*1.1]);
set(gca,'XTick',[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e10]);
l(1) = xlabel('Grain size');
l(2) = ylabel('Bin Count');
title('Random lognormal sample')
t(1) = text(xd,max(h_test.Values)*1.1*y,' a');

ax2(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
hold on
h = histogram(X,bins,'FaceColor',pl.c2);
for i = 1:length(D) %Bin centers of grain size classes
    plot([D(i) D(i)],[0 data*1.1],'color',[0.5 0.5 0.5 0.2]);
end
set(gca,'Xscale','log');
axis([6.3e-7 0.63 0 max(h.Values)*1.1]);
set(gca,'YTick',(0:5:100)*10000,'YtickLabel',0:5:100);
set(gca,'XTick',D,'XtickLabel',bins_txt);
%xlabel('Grain size')
l(3) = ylabel('Weight %');
title('Grain size class bins')
t(2) = text(xd,max(h.Values)*1.1*y,' b');

ax2(3) = axes('Position',[pl.lmarge+2*pl.wt+2*pl.inth pl.bmarge pl.wt pl.ht]);
hold on
for i = 1:length(D) %Bin centers of grain size classes
    plot([D(i) D(i)],[0 data],'color',[0.5 0.5 0.5 0.2]);
end
plot(D_test,ones(1,length(bins_y_test))*10,'color','k');
plot(D_test,ones(1,length(bins_y_test))*50,'color','k');
plot(D_test,ones(1,length(bins_y_test))*90,'color','k');
plot(D_test,cum_perc_test,'color',pl.c1,'LineWidth',pl.line); %original
plot(bins(2:end),cum_perc,'color',pl.c2,'LineWidth',pl.line); %classes
set(gca,'Xscale','log');
axis([6.3e-7 0.63 0 100]);
set(gca,'XTick',D,'XtickLabel',bins_txt);
%xlabel('Grain size')
l(4) = ylabel('Weight %');
title('Cumulative distribution')
t(3) = text(xd,100*y,' c');

set(ax2,'box','on','Layer','top', ...
    'XMinorTick','off','YMinorTick','off', ...
    'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
    'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
    'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

set(t,'FontSize',pl.fsz2);
set(l,'FontSize',pl.fsz2);

xtickangle(ax2(2:3),90);
end

%%
set(gcf,'renderer','painters');
print(f2,'-dpng',[output '/FlowOnMars5_GrainSizeDistribution'],'-r400');
print(f2,'-dpdf',[output '/FlowOnMars5_GrainSizeDistribution'],'-r400');
save('FlowOnMars5_GrainSizeDistribution.mat','bins_txt','bins_y_perc','D','D10','D50','D90');

